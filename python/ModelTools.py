from __future__ import absolute_import, print_function

import os
import os.path
import re
from functools import reduce
from math import *
from sys import exit, stderr, stdout

import six
from six.moves import range

import ROOT

ROOFIT_EXPR = "expr"
ROOFIT_EXPR_PDF = "EXPR"


class SafeWorkspaceImporter:
    """Class that provides the RooWorkspace::import method, but makes sure we call the proper
    overload of it, since in ROOT 6 sometimes PyROOT calls the wrong one"""

    def __init__(self, wsp):
        self.wsp = wsp
        self.imp = getattr(wsp, "import")

    def __call__(self, *args):
        if len(args) != 1:
            self.imp(*args)
        elif (
            args[0].Class().InheritsFrom("RooAbsReal")
            or args[0].Class().InheritsFrom("RooArgSet")
            or args[0].Class().InheritsFrom("RooAbsData")
            or args[0].Class().InheritsFrom("RooCategory")
        ):
            self.imp(args[0], ROOT.RooCmdArg())  # force the proper overload to be called
        else:
            self.imp(*args)


class ModelBuilderBase:
    """This class defines the basic stuff for a model builder, and it's an interface on top of RooWorkspace::factory or HLF files"""

    def __init__(self, options):
        self.options = options
        self.out = stdout
        self.discrete_param_set = []
        if options.bin:
            if options.out == None:
                options.out = re.sub(".txt$", "", options.fileName) + ".root"
            options.baseDir = os.path.dirname(options.fileName)
            ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
            ROOT.TH1.AddDirectory(False)
            self.out = ROOT.RooWorkspace("w", "w")
            # self.out.safe_import = getattr(self.out,"import") # workaround: import is a python keyword
            self.out.safe_import = SafeWorkspaceImporter(self.out)
            self.objstore = {}
            self.out.dont_delete = []
            if options.verbose == 0:
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
            elif options.verbose < 3:
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
            if "ROOFITSYS" in os.environ:
                ROOT.gSystem.AddIncludePath(" -I%s/include " % os.environ["ROOFITSYS"])
        elif options.out != None:
            # stderr.write("Will save workspace to HLF file %s" % options.out)
            self.out = open(options.out, "w")
        if not options.bin:
            stderr.write("\nWARNING: You're not using binary mode. This is is DEPRECATED and NOT SUPPORTED anymore, and can give WRONG results.\n\n")
        if options.cexpr:
            global ROOFIT_EXPR
            ROOFIT_EXPR = "cexpr"

    def addObj(self, classtype, name, *args):
        if name not in self.objstore:
            self.objstore[name] = classtype(name, *args)
        return self.objstore[name]

    def getObj(self, name):
        if name not in self.objstore:
            raise RuntimeError("Requested object %s not found in store" % name)
        return self.objstore[name]

    def renameObj(self, currName, newName):
        if currName not in self.objstore:
            raise RuntimeError("Requested object %s not found in store" % name)
        self.objstore[currName].SetName(newName)
        self.objstore[newName] = self.objstore.pop(currName)

    def factory_(self, X):
        if self.options.verbose >= 7:
            print("RooWorkspace::factory('%s')" % X)
        if len(X) > 1000:
            print(
                "Executing factory with a string of length ",
                len(X),
                " > 1000, could trigger a bug: ",
                X,
            )
        ret = self.out.factory(X)
        if ret:
            if self.options.verbose >= 7:
                print(" ---> ", ret)
            self.out.dont_delete.append(ret)
            return ret
        else:
            print("ERROR parsing '%s'" % X)
            self.out.Print("V")
            raise RuntimeError("Error in factory statement")

    def doComment(self, X):
        if not self.options.bin:
            self.out.write("// " + X + "\n")

    def doVar(self, vardef):
        if self.options.bin:
            self.factory_(vardef)
        else:
            self.out.write(vardef + ";\n")

    def doExp(self, name, expression, vars):
        if self.options.bin:
            self.factory_('expr::%s("%s",%s)' % (name, expression, vars))
        else:
            self.out.write('%s = expr::%s("%s",%s)' % (name, name, expression, vars) + ";\n")

    def doSet(self, name, vars):
        if self.options.bin:
            self.out.defineSet(name, vars)
        else:
            self.out.write("%s = set(%s);\n" % (name, vars))

    def doObj(self, name, type, X, ignoreExisting=False):
        if self.out.obj(name) and ignoreExisting:
            return 1  # Still complain if not explicitly told to ignore the existing object
        if self.options.bin:
            return self.factory_("%s::%s(%s)" % (type, name, X))
        else:
            self.out.write("%s = %s(%s);\n" % (name, type, X))

    def addDiscrete(self, var):
        if self.options.removeMultiPdf:
            return
        self.discrete_param_set.append(var)


class ModelBuilder(ModelBuilderBase):
    """This class defines the actual methods to build a model"""

    def __init__(self, datacard, options):
        ModelBuilderBase.__init__(self, options)
        self.DC = datacard
        self.doModelBOnly = True
        self.selfNormBins = []
        self.extraNuisances = []
        self.extraGlobalObservables = []
        self.globalobs = []

    def getSafeNormName(self, n):
        # need to be careful in case user has _norm name and wants to auto-create flatPrior
        if self.options.flatParamPrior:
            if n in self.DC.pdfnorms.keys():
                return self.DC.pdfnorms[n]
        return n

    def setPhysics(self, physicsModel):
        self.physics = physicsModel
        self.physics.setModelBuilder(self)

    def doModel(self, justCheckPhysicsModel=False):
        if not justCheckPhysicsModel:
            self.doObservables()
        self.physics.doParametersOfInterest()

        # set a group attribute on POI variables
        poiIter = self.out.set("POI").createIterator()
        poi = poiIter.Next()
        while poi:
            self.out.var(poi.GetName()).setAttribute("group_POI", True)
            poi = poiIter.Next()
        self.physics.preProcessNuisances(self.DC.systs)
        self.doNuisances()
        self.doExtArgs()
        self.doRateParams()
        self.doAutoFlatNuisancePriors()
        self.doFillNuisPdfsAndSets()
        self.doExpectedEvents()
        if justCheckPhysicsModel:
            self.physics.done()
            print("Model is OK")
            exit(0)
        self.doIndividualModels()
        self.doNuisancesGroups()  # this needs to be called after both doNuisances and doIndividualModels
        self.doCombination()
        self.runPostProcesses()
        self.physics.done()
        if self.options.bin:
            self.doModelConfigs()
            if self.options.verbose > 1:
                self.out.Print("tv")
            if self.options.verbose > 2:
                self.out.pdf("model_s").graphVizTree(self.options.out + ".dot", "\\n")
                print("Wrote GraphVizTree of model_s to ", self.options.out + ".dot")

    def getRenamingParameters(self):
        toFreeze = []
        renameParamString = []
        paramString = []
        for n in self.DC.systematicsParamMap.keys():
            paramString.append(n)
            renameParamString.append(self.DC.systematicsParamMap[n])
            if n != self.DC.systematicsParamMap[n]:
                toFreeze.append(n)
        if len(renameParamString):
            renameParamString = ",".join(renameParamString)
            paramString = ",".join(paramString)
        return paramString, renameParamString, toFreeze

    def runPostProcesses(self):
        for n in self.DC.frozenNuisances:
            self.out.arg(n).setConstant(True)

    def doExtArgs(self):
        open_files = {}
        for rp in self.DC.extArgs.keys():
            if self.out.arg(rp):
                continue
            argv = self.DC.extArgs[rp][-1]
            if ":" in argv:
                split = argv.split(":")
                importargs = []
                if "RecycleConflictNodes" in split:
                    split.remove("RecycleConflictNodes")
                    importargs.append(ROOT.RooFit.RecycleConflictNodes())
                fin, wsn = split
                if (fin, wsn) in open_files:
                    wstmp = open_files[(fin, wsn)]
                    if not wstmp.arg(rp):
                        raise RuntimeError("No parameter '%s' found for extArg in workspace %s from file %s" % (rp, wsn, fin))
                    self.out.safe_import(wstmp.arg(rp), *importargs)
                else:
                    fitmp = ROOT.TFile.Open(fin)
                    if not fitmp:
                        raise RuntimeError("No File '%s' found for extArg" % fin)
                    wstmp = fitmp.Get(wsn)
                    if not wstmp:
                        raise RuntimeError("Workspace '%s' not in file %s" % (wsn, fin))
                    if not wstmp.arg(rp):
                        raise RuntimeError("No parameter '%s' found for extArg in workspace %s from file %s" % (rp, wsn, fin))
                    self.out.safe_import(wstmp.arg(rp), *importargs)
                    open_files[(fin, wsn)] = wstmp
            else:
                param_range = ""
                param_val = self.DC.extArgs[rp][-1]
                if len(self.DC.extArgs[rp]) > 3:  # range is included:
                    param_range = self.DC.extArgs[rp][-1]
                    param_val = self.DC.extArgs[rp][-2]
                    if "[" not in param_range:
                        raise RuntimeError("Expected range arguments [min,max] or [const] for extArg %s " % (rp))
                    param_range = param_range.strip("[]")

                removeRange = False
                setConst = param_range == "const"
                if param_range in ["", "const"]:
                    v = float(param_val)
                    param_range = "%g,%g" % (-2.0 * abs(v), 2.0 * abs(v))
                    removeRange = True

                self.doVar("%s[%s,%s]" % (rp, float(param_val), param_range))
                if removeRange:
                    self.out.var(rp).removeRange()
                self.out.var(rp).setConstant(False)
                if setConst:
                    self.out.var(rp).setConstant(True)
                self.out.var(rp).setAttribute("flatParam")

    def doRateParams(self):
        # First support external functions/parameters
        # keep a map of open files/workspaces
        open_files = {}

        for rp in self.DC.rateParams.keys():
            for rk in range(len(self.DC.rateParams[rp])):
                type = self.DC.rateParams[rp][rk][0][-1]
                if type != 2:
                    continue
                argu, argv = (
                    self.DC.rateParams[rp][rk][0][0],
                    self.DC.rateParams[rp][rk][0][1],
                )
                if self.out.arg(argu):
                    continue
                fin, wsn = argv.split(":")
                if (fin, wsn) in open_files:
                    wstmp = open_files[(fin, wsn)]
                    if not wstmp.arg(argu):
                        raise RuntimeError("No parameter '%s' found for rateParam in workspace %s from file %s" % (argu, wsn, fin))
                    self.out.safe_import(wstmp.arg(argu), ROOT.RooFit.RecycleConflictNodes())
                else:
                    fitmp = ROOT.TFile.Open(fin)
                    if not fitmp:
                        raise RuntimeError("No File '%s' found for rateParam" % fin)
                    wstmp = fitmp.Get(wsn)
                    if not wstmp:
                        raise RuntimeError("Workspace '%s' not in file %s" % (wsn, fin))
                    if not wstmp.arg(argu):
                        raise RuntimeError("No parameter '%s' found for rateParam in workspace %s from file %s" % (argu, wsn, fin))
                    self.out.safe_import(wstmp.arg(argu), ROOT.RooFit.RecycleConflictNodes())
                    open_files[(fin, wsn)] = wstmp
                    # fitmp.Close()

        # First do independant parameters, then expressions
        for rp in self.DC.rateParams.keys():
            for rk in range(len(self.DC.rateParams[rp])):
                type = self.DC.rateParams[rp][rk][0][-1]
                if type != 0:
                    continue
                param_range = (self.DC.rateParams[rp][rk][1]).strip("[]")
                argu, argv = (
                    self.DC.rateParams[rp][rk][0][0],
                    self.DC.rateParams[rp][rk][0][1],
                )
                if self.out.arg(argu):
                    continue

                v = float(argv)
                removeRange = len(param_range) == 0
                if param_range == "":
                    if self.options.flatParamPrior:
                        raise ValueError(
                            "Cannot create flat Prior for rateParam nuisance parameter '"
                            + argu
                            + "' without specifying a range [a,b]. Please fix in the datacard"
                        )
                    ## check range. The parameter needs to be created in range. Then we will remove it
                    param_range = "%g,%g" % (-2.0 * abs(v), 2.0 * abs(v))
                # additional check for range requested
                lo_r, hi_r = map(float, param_range.split(","))
                if v < lo_r or v > hi_r:
                    raise ValueError("Parameter: " + argu + " asked to be created out-of-range (it will lead to an error): " + argv + ":" + param_range)

                self.doVar("%s[%s,%s]" % (argu, argv, param_range))
                if removeRange:
                    self.out.var(argu).removeRange()
                self.out.var(argu).setConstant(False)
                self.out.var(argu).setAttribute("flatParam")

        # functions are tricky (functions of functions?)
        toBeCreated = []
        for rp in self.DC.rateParams.keys():
            for rk in range(len(self.DC.rateParams[rp])):
                type = self.DC.rateParams[rp][rk][0][-1]
                if type != 1:
                    continue
                argu, arge, argv = (
                    self.DC.rateParams[rp][rk][0][0],
                    self.DC.rateParams[rp][rk][0][1],
                    self.DC.rateParams[rp][rk][0][2],
                )
                if self.out.arg(argu):
                    continue
                if not reduce(
                    lambda x, y: x * y,
                    [self.out.arg(a) != None for a in argv.split(",")],
                    1,
                ):
                    toBeCreated.append([argu, arge, argv])
                else:
                    self.doExp(argu, arge, argv)

        # by now we 've probably picked up the majority of the, repeat through list until we get them all
        tbc = toBeCreated[:]
        while True:
            toBeCreated = tbc[:]
            if len(toBeCreated) == 0:
                break
            for rp in toBeCreated:
                argu, arge, argv = rp[0], rp[1], rp[2]
                if reduce(
                    lambda x, y: x * y,
                    [self.out.arg(a) != None for a in argv.split(",")],
                    1,
                ):
                    self.doExp(argu, arge, argv)
                    tbc.remove([argu, arge, argv])
            if len(tbc) == len(toBeCreated):
                print(tbc, " -> ", toBeCreated)
                raise RuntimeError("Cannot produce following rateParams (dependent parameters not found!) %s" % (",".join([t[0] for t in toBeCreated])))

    def doObservables(self):
        """create pdf_bin<X> and pdf_bin<X>_bonly for each bin"""
        raise RuntimeError("Not implemented in ModelBuilder")

    def doNuisances(self):
        for cpar in self.DC.discretes:
            self.addDiscrete(cpar)

        if len(self.DC.systs) == 0:
            return
        self.doComment(" ----- nuisances -----")
        # globalobs = []

        for n, nofloat, pdf, args, errline in self.DC.systs:
            is_func_scaled = False
            func_scaler = None
            for pn, pf in self.options.nuisanceFunctions:
                if re.match(pn, n):
                    is_func_scaled = True
                    func_scaler = pf
                    if self.options.verbose > 1:
                        print("Rescaling %s constraint as %s" % (n, pf))
            for pn, pf in self.options.nuisanceGroupFunctions:
                if pn in self.DC.groups and n in self.DC.groups[pn]:
                    is_func_scaled = True
                    func_scaler = pf
                    if self.options.verbose > 1:
                        print("Rescaling %s constraint (in group %s) as %s" % (n, pn, pf))
            if pdf == "lnN" or (pdf.startswith("shape") and pdf != "shapeU"):
                r = "-4,4" if pdf == "shape" else "-7,7"
                sig = 1.0
                for pn, pf in self.options.nuisancesToRescale:
                    if re.match(pn, n):
                        sig = float(pf)
                        sigscale = sig * (4 if pdf == "shape" else 7)
                        r = "-%g,%g" % (sigscale, sigscale)
                sig = "%g" % sig
                if is_func_scaled:
                    sig = func_scaler
                r_exp = "" if self.out.var(n) else "[%s]" % r  # Specify range to invoke factory to produce a RooRealVar only if it doesn't already exist
                if self.options.noOptimizePdf or is_func_scaled:
                    self.doObj(
                        "%s_Pdf" % n,
                        "Gaussian",
                        "%s%s, %s_In[0,%s], %s" % (n, r_exp, n, r, sig),
                        True,
                    )
                    # Use existing constraint since it could be a param
                    if is_func_scaled:
                        boundHi = self.doObj("%s_BoundHi" % n, "prod", "5, %s" % sig)
                        boundLo = self.doObj("%s_BoundLo" % n, "prod", "-5, %s" % sig)
                        self.out.var(n).setRange(boundLo, boundHi)
                else:
                    self.doObj(
                        "%s_Pdf" % n,
                        "SimpleGaussianConstraint",
                        "%s%s, %s_In[0,%s], %s" % (n, r_exp, n, r, sig),
                        True,
                    )
                    # Use existing constraint since it could be a param
                self.out.var(n).setVal(0)
                self.out.var(n).setError(1)
                self.globalobs.append("%s_In" % n)
                if self.options.bin:
                    self.out.var("%s_In" % n).setConstant(True)
                if self.options.optimizeBoundNuisances and not is_func_scaled:
                    self.out.var(n).setAttribute("optimizeBounds")
            elif pdf == "gmM":
                val = 0
                for c in errline.values():  # list channels
                    for v in c.values():  # list effects in each channel
                        if v != 0:
                            if val != 0 and v != val:
                                raise RuntimeError("Error: line %s contains two different uncertainties %g, %g, which is not supported for gmM" % (n, v, val))
                            val = v
                if val == 0:
                    raise RuntimeError("Error: line %s contains all zeroes")
                theta = val * val
                kappa = 1 / theta
                self.doObj(
                    "%s_Pdf" % n,
                    "Gamma",
                    "%s[1,%f,%f], %s_In[%g,%g,%g], %s_scaling[%g], 0"
                    % (
                        n,
                        max(0.01, 1 - 5 * val),
                        1 + 5 * val,
                        n,
                        kappa,
                        1,
                        2 * kappa + 4,
                        n,
                        theta,
                    ),
                )
                self.globalobs.append("%s_In" % n)
                if self.options.bin:
                    self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "gmN":
                if False:
                    # old version, that creates a poisson with a very large range
                    self.doObj(
                        "%s_Pdf" % n,
                        "Poisson",
                        "%s_In[%d,0,%d], %s[0,%d], 1" % (n, args[0], 2 * args[0] + 5, n, 2 * args[0] + 5),
                    )
                else:
                    # new version, that creates a poisson with a narrower range (but still +/- 7 sigmas)
                    # print "Searching for bounds for",n,"poisson with obs",args[0]
                    minExp = args[0] + 1 if args[0] > 0 else 0
                    while (ROOT.TMath.Poisson(args[0], minExp) > 1e-12) and minExp > 0:
                        # print "Poisson(%d, minExp = %f) = %g > 1e-12" % (args[0], minExp, ROOT.TMath.Poisson(args[0], minExp))
                        minExp *= 0.8
                    maxExp = args[0] + 1
                    while ROOT.TMath.Poisson(args[0], maxExp) > 1e-12:
                        # print "Poisson(%d, maxExp = %f) = %g > 1e-12" % (args[0], maxExp, ROOT.TMath.Poisson(args[0], maxExp))
                        maxExp *= 1.2
                    minObs = args[0]
                    while minObs > 0 and (ROOT.TMath.Poisson(minObs, args[0] + 1) > 1e-12):
                        # print "Poisson(minObs = %d, %f) = %g > 1e-12" % (minObs, args[0]+1, ROOT.TMath.Poisson(minObs, args[0]+1))
                        minObs -= sqrt(args[0]) if args[0] > 10 else 1
                    maxObs = args[0] + 2
                    while ROOT.TMath.Poisson(maxObs, args[0] + 1) > 1e-12:
                        # print "Poisson(maxObs = %d, %f) = %g > 1e-12" % (maxObs, args[0]+1, ROOT.TMath.Poisson(maxObs, args[0]+1))
                        maxObs += sqrt(args[0]) if args[0] > 10 else 2
                    self.doObj(
                        "%s_Pdf" % n,
                        "Poisson",
                        "%s_In[%d,%f,%f], %s[%f,%f,%f], 1" % (n, args[0], minObs, maxObs, n, args[0] + 1, minExp, maxExp),
                    )
                self.globalobs.append("%s_In" % n)
                if self.options.bin:
                    self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "trG":
                trG_min = -7
                trG_max = +7
                for b in errline.keys():
                    for v in errline[b].values():
                        if v > 0 and 1.0 + trG_min * v < 0:
                            trG_min = -1.0 / v
                        if v < 0 and 1.0 + trG_max * v < 0:
                            trG_max = -1.0 / v
                r = "%f,%f" % (trG_min, trG_max)
                self.doObj("%s_Pdf" % n, "Gaussian", "%s[0,%s], %s_In[0,%s], 1" % (n, r, n, r))
                self.globalobs.append("%s_In" % n)
                if self.options.bin:
                    self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "lnU" or pdf == "shapeU":
                self.doObj("%s_Pdf" % n, "Uniform", "%s[-1,1]" % n)
            elif pdf == "unif":
                self.doObj("%s_Pdf" % n, "Uniform", "%s[%f,%f]" % (n, args[0], args[1]))
            elif pdf == "flatParam" and self.options.flatParamPrior:
                c_param_name = self.getSafeNormName(n)
                if self.out.var(c_param_name):
                    v, x1, x2 = self.out.var(c_param_name).getVal(), self.out.var(c_param_name).getMin(), self.out.var(c_param_name).getMax()
                    self.DC.toCreateFlatParam[c_param_name] = [v, x1, x2]
                else:
                    self.DC.toCreateFlatParam[c_param_name] = []

            elif pdf == "dFD" or pdf == "dFD2":
                dFD_min = -(1 + 8 / args[0])
                dFD_max = +(1 + 8 / args[0])
                for b in errline.keys():
                    for v in errline[b].values():
                        if v > 0 and 1.0 + dFD_min * v < 0:
                            dFD_min = -1.0 / v
                        if v < 0 and 1.0 + dFD_max * v < 0:
                            dFD_max = -1.0 / v
                r = "%f,%f" % (dFD_min, dFD_max)
                # r = "%f,%f" % (-(1+8/args[0]), +(1+8/args[0]));
                # r = "-1,1"
                if pdf == "dFD":
                    self.doObj(
                        "%s_Pdf" % n,
                        ROOFIT_EXPR_PDF,
                        "'1/(2*(1+exp(%f*((@0-@1)-1)))*(1+exp(-%f*((@0-@1)+1))))', %s[0,%s], %s_In[0,%s]" % (args[0], args[0], n, r, n, r),
                    )
                else:
                    self.doObj(
                        "%s_Pdf" % n,
                        ROOFIT_EXPR_PDF,
                        "'1/(2*(1+exp(%f*(@0-1)))*(1+exp(-%f*(@0+1))))', %s[0,%s], %s_In[0,%s]" % (args[0], args[0], n, r, n, r),
                    )
                self.globalobs.append("%s_In" % n)
                if self.options.bin:
                    self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "constr":
                print("-------------- WARNING, constraint found --> make sure that you know what you are doing!")
                ## I want to construct this line
                ## constr1_In[0.],RooFormulaVar::fconstr1("r_Bin0+r_Bin2-2*r_Bin1",{r_Bin0,r_Bin1,r_Bin2}),constr1_S[0.001000]
                ##  Assuming args=
                ##   r_Bin0+r_Bin2-2*r_Bin1 0.001
                ## or r_Bin0+r_Bin2-2*r_Bin1 {r_Bin0,r_Bin1,r_Bin2} 0.001000
                ## the parameter can be a number or a variable
                d = {
                    "pdf": "%s_Pdf" % n,
                    "name": n,
                    "function": "%s_Func" % n,
                    "in": "%s_In" % n,
                    "sigma": "%s_S" % n,
                    "formula": args[0],
                    "param": args[-1],
                }
                if len(args) > 2:
                    d["depend"] = args[1] if args[1][0] == "{" else "{" + args[1] + "}"
                else:
                    remove = set(["TMath", "Exp", "::", ""])
                    l = list(set(re.split("\\+|-|\\*|/|\\(|\\)", d["formula"])) - remove)
                    l2 = []  ## remove all the non-float expressions
                    for x in l:
                        try:
                            float(x)
                        except ValueError:
                            l2.append(x)
                    d["depend"] = "{" + ",".join(l2) + "}"

                ## derve the constrain strength
                try:
                    float(d["param"])  # raise exception
                    if self.options.verbose > 2:
                        print("DEBUG constr", "param is a number")
                    d["fullsigma"] = "%(sigma)s[%(param)s]" % d
                except ValueError:
                    if self.options.verbose > 2:
                        print("DEBUG constr", "param is a variable")
                    d["fullsigma"] = d["param"]

                if self.options.verbose > 2:
                    print("DEBUG constr", "args", args)
                if self.options.verbose > 2:
                    print("DEBUG constr", "Dictionary", d)

                full = '%(in)s[0.],RooFormulaVar::%(function)s("%(formula)s",%(depend)s),%(fullsigma)s' % d
                # self.doObj("%s_Pdf"%n, "Gaussian"," ".join(args),True)
                v = self.options.verbose
                self.options.verbose = 10  # force debug this line
                self.doObj("%s_Pdf" % n, "Gaussian", full, True)
                self.options.verbose = v
            elif pdf == "param":
                mean = float(args[0])
                if "/" in args[1]:
                    sigmaL, sigmaR = args[1].split("/")
                    if sigmaL[0] != "-" or sigmaR[0] != "+":
                        raise RuntimeError("Asymmetric parameter uncertainties should be entered as -x/+y")
                    sigmaL = sigmaL[1:]
                    sigmaR = sigmaR[1:]
                    if len(args) == 3:  # mean, sigma, range
                        if self.out.var(n):
                            bounds = [float(x) for x in args[2][1:-1].split(",")]
                            self.out.var(n).setConstant(False)
                            if self.out.var(n).getMin() != bounds[0] or self.out.var(n).getMax() != bounds[1]:
                                print(
                                    "Resetting range for %s to be [%s,%s] from param statement (was [%s,%s])"
                                    % (
                                        n,
                                        bounds[0],
                                        bounds[1],
                                        self.out.var(n).getMin(),
                                        self.out.var(n).getMax(),
                                    )
                                )
                            self.out.var(n).setRange(bounds[0], bounds[1])
                        else:
                            self.doVar("%s%s" % (n, args[2]))
                    else:
                        if self.out.var(n):
                            self.out.var(n).setConstant(False)
                            self.out.var(n).setRange(mean - 4 * float(sigmaL), mean + 4 * float(sigmaR))
                        else:
                            self.doVar(
                                "%s[%g,%g]"
                                % (
                                    n,
                                    mean - 4 * float(sigmaL),
                                    mean + 4 * float(sigmaR),
                                )
                            )
                    self.out.var(n).setVal(mean)
                    self.out.var(n).setError(0.5 * (float(sigmaL) + float(sigmaR)))

                    sigmaStrL = sigmaL
                    sigmaStrR = sigmaR
                    if is_func_scaled:
                        sigmaStrL = "%s_WidthScaledL" % n
                        sigmaStrR = "%s_WidthScaledR" % n
                        self.doObj(sigmaStrL, "prod", "%g, %s" % (float(sigmaL), func_scaler))
                        self.doObj(sigmaStrR, "prod", "%g, %s" % (float(sigmaR), func_scaler))
                    self.doObj(
                        "%s_Pdf" % n,
                        "BifurGauss",
                        "%s, %s_In[%s,%g,%g], %s, %s"
                        % (
                            n,
                            n,
                            args[0],
                            self.out.var(n).getMin(),
                            self.out.var(n).getMax(),
                            sigmaStrL,
                            sigmaStrR,
                        ),
                        True,
                    )
                    self.out.var("%s_In" % n).setConstant(True)
                    if is_func_scaled:
                        self.doExp(
                            "%s_BoundHi" % n,
                            "%g+%g*@0" % (mean, self.out.var(n).getMax() - mean),
                            "%s" % (func_scaler),
                        )
                        self.doExp(
                            "%s_BoundLo" % n,
                            "%g-%g*@0" % (mean, mean - self.out.var(n).getMin()),
                            "%s" % (func_scaler),
                        )
                        self.out.var(n).setRange(
                            self.out.function("%s_BoundLo" % n),
                            self.out.function("%s_BoundHi" % n),
                        )
                else:
                    if len(args) == 3:  # mean, sigma, range
                        sigma = float(args[1])
                        if self.out.var(n):
                            bounds = [float(x) for x in args[2][1:-1].split(",")]
                            self.out.var(n).setConstant(False)
                            if self.out.var(n).getMin() != bounds[0] or self.out.var(n).getMax() != bounds[1]:
                                print(
                                    "Resetting range for %s to be [%s,%s] from param statement (was [%s,%s])"
                                    % (
                                        n,
                                        bounds[0],
                                        bounds[1],
                                        self.out.var(n).getMin(),
                                        self.out.var(n).getMax(),
                                    )
                                )
                                self.out.var(n).setRange(bounds[0], bounds[1])
                        else:
                            self.doVar("%s%s" % (n, args[2]))
                    else:
                        sigma = float(args[1])
                        if self.out.var(n):
                            self.out.var(n).setConstant(False)
                            self.out.var(n).setRange(mean - 4 * sigma, mean + 4 * sigma)
                        else:
                            self.doVar("%s[%g,%g]" % (n, mean - 4 * sigma, mean + 4 * sigma))
                    self.out.var(n).setVal(mean)
                    self.out.var(n).setError(sigma)
                    sigmaStr = args[1]
                    if is_func_scaled:
                        sigmaStr = "%s_WidthScaled" % n
                        self.doObj(sigmaStr, "prod", "%g, %s" % (float(args[1]), func_scaler))
                    if self.options.noOptimizePdf or is_func_scaled:
                        self.doObj(
                            "%s_Pdf" % n,
                            "Gaussian",
                            "%s, %s_In[%s,%g,%g], %s"
                            % (
                                n,
                                n,
                                args[0],
                                self.out.var(n).getMin(),
                                self.out.var(n).getMax(),
                                sigmaStr,
                            ),
                            True,
                        )
                    else:
                        self.doObj(
                            "%s_Pdf" % n,
                            "SimpleGaussianConstraint",
                            "%s, %s_In[%s,%g,%g], %s"
                            % (
                                n,
                                n,
                                args[0],
                                self.out.var(n).getMin(),
                                self.out.var(n).getMax(),
                                sigmaStr,
                            ),
                            True,
                        )
                    self.out.var("%s_In" % n).setConstant(True)
                    if is_func_scaled:
                        boundHi = self.doExp(
                            "%s_BoundHi" % n,
                            "%g+%g*@0" % (mean, self.out.var(n).getMax() - mean),
                            "%s" % (func_scaler),
                        )
                        boundLo = self.doExp(
                            "%s_BoundLo" % n,
                            "%g-%g*@0" % (mean, mean - self.out.var(n).getMin()),
                            "%s" % (func_scaler),
                        )
                        self.out.var(n).setRange(
                            self.out.function("%s_BoundLo" % n),
                            self.out.function("%s_BoundHi" % n),
                        )
                self.globalobs.append("%s_In" % n)
                # if self.options.optimizeBoundNuisances: self.out.var(n).setAttribute("optimizeBounds")
            elif pdf == "extArg":
                continue

            else:
                raise RuntimeError("Unsupported pdf %s" % pdf)
            if nofloat:
                self.out.var(n).setAttribute("globalConstrained", True)
            # self.out.var(n).Print('V')
            if n in self.DC.frozenNuisances:
                self.out.var(n).setConstant(True)

    def doFillNuisPdfsAndSets(self):
        if self.options.bin:
            # avoid duplicating  _Pdf in list
            setNuisPdf = []
            nuisPdfs = ROOT.RooArgList()
            nuisVars = ROOT.RooArgSet()
            for n, nf, p, a, e in self.DC.systs:
                c_param_name = self.getSafeNormName(n)
                if p != "constr":
                    nuisVars.add(self.out.var(c_param_name))
                setNuisPdf.append(c_param_name)
            setNuisPdf = list(dict.fromkeys((setNuisPdf)))
            for n in setNuisPdf:
                nuisPdfs.add(self.out.pdf(n + "_Pdf"))
            self.out.defineSet("nuisances", nuisVars)
            self.out.nuisPdf = ROOT.RooProdPdf("nuisancePdf", "nuisancePdf", nuisPdfs)
            self.out.safe_import(self.out.nuisPdf)
            self.out.nuisPdfs = nuisPdfs
            gobsVars = ROOT.RooArgSet()
            for g in self.globalobs:
                gobsVars.add(self.out.var(g))
            self.out.defineSet("globalObservables", gobsVars)
        else:  # doesn't work for too many nuisances :-(
            # avoid duplicating  _Pdf in list
            setNuisPdf = list(dict.fromkeys(keywords([self.getSafeNormName(n) for (n, nf, p, a, e) in self.DC.systs])))
            self.doSet("nuisances", ",".join(["%s" % self.getSafeNormName(n) for (n, nf, p, a, e) in self.DC.systs]))
            self.doObj("nuisancePdf", "PROD", ",".join(["%s_Pdf" % n for n in setNuisPdf]))
            self.doSet("globalObservables", ",".join(self.globalobs))

    def doAutoFlatNuisancePriors(self):
        if len(self.DC.toCreateFlatParam.keys()) > 0:
            for flatNP in self.DC.toCreateFlatParam.items():
                c_param_name = flatNP[0]
                c_param_details = flatNP[1]
                if len(c_param_details):
                    v, x1, x2 = c_param_details
                else:
                    v, x1, x2 = self.out.var(c_param_name).getVal(), self.out.var(c_param_name).getMin(), self.out.var(c_param_name).getMax()
                if self.options.verbose > 2:
                    print("Will create flat prior for parameter ", c_param_name, " with range [", x1, x2, "]")
                self.doExp(
                    "%s_diff_expr" % c_param_name, "%s-%s_In" % (c_param_name, c_param_name), "%s,%s_In[%g,%g,%g]" % (c_param_name, c_param_name, v, x1, x2)
                )
                self.doObj("%s_Pdf" % c_param_name, "Uniform", "%s_diff_expr" % c_param_name)
                self.out.var("%s_In" % c_param_name).setConstant(True)
                self.globalobs.append("%s_In" % c_param_name)

    def doNuisancesGroups(self):
        # Prepare a dictionary of which group a certain nuisance belongs to
        groupsFor = {}
        # existingNuisanceNames = tuple(set([syst[0] for syst in self.DC.systs]+self.DC.flatParamNuisances.keys()+self.DC.rateParams.keys()+self.DC.extArgs.keys()+self.DC.discretes))
        existingNuisanceNames = self.DC.getAllVariables()
        for groupName, nuisanceNames in six.iteritems(self.DC.groups):
            for nuisanceName in nuisanceNames:
                if nuisanceName not in existingNuisanceNames:
                    raise RuntimeError(
                        'Nuisance group "%(groupName)s" refers to nuisance "%(nuisanceName)s" but it does not exist. Perhaps you misspelled it.' % locals()
                    )
                if nuisanceName in groupsFor:
                    groupsFor[nuisanceName].append(groupName)
                else:
                    groupsFor[nuisanceName] = [groupName]

        # print self.DC.groups
        # print groupsFor
        for n in existingNuisanceNames:
            # set an attribute related to the group(s) this nuisance belongs to
            if n in groupsFor:
                groupNames = groupsFor[n]
                if self.options.verbose > 1:
                    print('Nuisance "%(n)s" is assigned to the following nuisance groups: %(groupNames)s' % locals())
                for groupName in groupNames:
                    var = self.out.var(n)
                    if not var:
                        var = self.out.cat(n)
                    if not var:
                        raise RuntimeError('Nuisance group "%(groupName)s" refers to nuisance but it is not an independent parameter.' % locals())
                    var.setAttribute("group_" + groupName, True)

        for groupName, nuisanceNames in six.iteritems(self.DC.groups):
            nuisanceargset = ROOT.RooArgSet()
            for nuisanceName in nuisanceNames:
                nuisanceargset.add(self.out.var(nuisanceName))
            self.out.defineSet("group_%s" % groupName, nuisanceargset)

    def doExpectedEvents(self):
        self.doComment(" --- Expected events in each bin, for each process ----")
        for b in self.DC.bins:
            for p in self.DC.exp[b].keys():  # so that we get only self.DC.processes contributing to this bin
                # if it's a zero background, write a zero and move on
                if self.DC.exp[b][p] == 0:
                    self.doVar("n_exp_bin%s_proc_%s[%g]" % (b, p, self.DC.exp[b][p]))
                    continue
                # get model-dependent scale factor
                scale = self.physics.getYieldScale(b, p)
                if scale == 0:
                    self.doVar("n_exp_bin%s_proc_%s[%g]" % (b, p, 0))
                    continue
                # collect multiplicative corrections
                nominal = self.DC.exp[b][p]
                gamma = None
                #  gamma normalization (if present, DC.exp[b][p] is ignored)
                factors = []  # RooAbsReal multiplicative factors (including gmN)
                logNorms = []  # (kappa, RooAbsReal) lnN (or lnN)
                alogNorms = []  # (kappaLo, kappaHi, RooAbsReal) asymm lnN
                if scale == 1:
                    pass
                elif type(scale) == str:
                    factors.append(scale)
                else:
                    raise RuntimeError("Physics model returned something that is neither a name, nor 0, nor 1.")

                # look for rate param for this bin
                if "%sAND%s" % (b, p) in list(self.DC.rateParams.keys()):
                    for rk in range(len(self.DC.rateParams["%sAND%s" % (b, p)])):
                        argu = self.DC.rateParams["%sAND%s" % (b, p)][rk][0][0]
                        if self.out.arg(argu):
                            factors.append(argu)
                        else:
                            raise RuntimeError("No rate parameter found %s, are you sure you defined it correctly in the datacard?" % (argu))
                selfNormRate = 1.0
                for n, nofloat, pdf, args, errline in self.DC.systs:
                    if pdf == "param":
                        continue
                    if pdf == "constr":
                        continue
                    if pdf == "rateParam" or pdf == "flatParam":
                        continue
                    if p not in errline[b]:
                        continue
                    if errline[b][p] == 0.0:
                        continue
                    if pdf.startswith("shape") and pdf.endswith("?"):  # might be a lnN in disguise
                        if not self.isShapeSystematic(b, p, n):
                            pdf = "lnN"
                    if pdf.startswith("shape"):
                        continue
                    if pdf == "lnN" and errline[b][p] == 1.0:
                        continue
                    if pdf == "lnN" or pdf == "lnU":
                        if type(errline[b][p]) == list:
                            elow, ehigh = errline[b][p]
                            alogNorms.append((elow, ehigh, n))
                        else:
                            logNorms.append((errline[b][p], n))
                    elif pdf == "gmM":
                        factors.append(n)
                    # elif pdf == "trG" or pdf == "unif" or pdf == "flatParam" or pdf == "dFD" or pdf == "dFD2":
                    elif pdf == "trG" or pdf == "unif" or pdf == "dFD" or pdf == "dFD2":
                        myname = "n_exp_shift_bin%s_proc_%s_%s" % (b, p, n)
                        self.doObj(myname, ROOFIT_EXPR, "'1+%f*@0', %s" % (errline[b][p], n))
                        factors.append(myname)
                    elif pdf == "gmN":
                        factors.append(n)
                        if abs(errline[b][p] * args[0] - self.DC.exp[b][p]) > max(0.05 * max(self.DC.exp[b][p], 1), errline[b][p]):
                            raise RuntimeError(
                                "Values of N = %d, alpha = %g don't match with expected rate %g for systematics %s "
                                % (args[0], errline[b][p], self.DC.exp[b][p], n)
                            )
                        if gamma != None:
                            raise RuntimeError("More than one gmN uncertainty for the same bin and process (second one is %s)" % n)
                        gamma = n
                        nominal = errline[b][p]
                        # The case with N=0 isn't relevant if the process provides its own normalisation,
                        # so we don't need to do anything special to handle it here.
                        if args[0] > 0:
                            selfNormRate = selfNormRate / args[0]
                    else:
                        raise RuntimeError("Unsupported pdf %s" % pdf)
                # optimize constants
                if len(factors) + len(logNorms) + len(alogNorms) == 0:
                    norm = selfNormRate if b in self.selfNormBins else self.DC.exp[b][p]
                    self.doVar("n_exp_bin%s_proc_%s[%g]" % (b, p, norm))
                else:
                    norm = selfNormRate if b in self.selfNormBins else nominal
                    # print "Process %s of bin %s depends on:\n\tlog-normals: %s\n\tasymm log-normals: %s\n\tother factors: %s\n" % (p,b,logNorms, alogNorms, factors)
                    procNorm = ROOT.ProcessNormalization("n_exp_bin%s_proc_%s" % (b, p), "", norm)
                    for kappa, thetaName in logNorms:
                        procNorm.addLogNormal(kappa, self.out.function(thetaName))
                    for kappaLo, kappaHi, thetaName in alogNorms:
                        procNorm.addAsymmLogNormal(kappaLo, kappaHi, self.out.function(thetaName))
                    for factorName in factors:
                        if self.out.function(factorName):
                            procNorm.addOtherFactor(self.out.function(factorName))
                        elif self.out.var(factorName):
                            procNorm.addOtherFactor(self.out.var(factorName))
                        elif self.out.arg(factorName):
                            raise RuntimeError(
                                "Factor %s for process %s, bin %s is a %s (not supported)"
                                % (
                                    factorName,
                                    p,
                                    b,
                                    self.out.arg(factorName).ClassName(),
                                )
                            )
                        else:
                            raise RuntimeError("Cannot add nonexistent factor %s for process %s, bin %s" % (factorName, p, b))

                    # take care of any variables which were renamed (eg for "param")
                    (
                        paramString,
                        renameParamString,
                        toFreeze,
                    ) = self.getRenamingParameters()
                    if len(renameParamString):
                        self.out.safe_import(
                            procNorm,
                            ROOT.RooFit.RecycleConflictNodes(),
                            ROOT.RooFit.RenameVariable(paramString, renameParamString),
                        )
                    else:
                        self.out.safe_import(procNorm)

    def doIndividualModels(self):
        """create pdf_bin<X> and pdf_bin<X>_bonly for each bin"""
        raise RuntimeError("Not implemented in ModelBuilder")

    def doCombination(self):
        """create model_s and model_b pdfs"""
        raise RuntimeError("Not implemented in ModelBuilder")

    def doModelConfigs(self):
        if not self.options.bin:
            raise RuntimeException
        if self.options.out == None:
            raise RuntimeException
        for nuis, warn in six.iteritems(self.DC.flatParamNuisances):
            if self.out.var(nuis):
                self.out.var(nuis).setAttribute("flatParam")
            elif warn:
                stderr.write("Missing variable %s declared as flatParam, will create one!\n" % nuis)
        mc_s = ROOT.RooStats.ModelConfig("ModelConfig", self.out)
        mc_b = ROOT.RooStats.ModelConfig("ModelConfig_bonly", self.out)
        for l, mc in [("s", mc_s), ("b", mc_b)]:
            if self.doModelBOnly:
                mc.SetPdf(self.out.pdf("model_" + l))
            else:
                mc.SetPdf(self.out.pdf("model_s"))
            # if l == 's' or mc.GetPdf().dependsOnValue(self.out.set("POI")):
            mc.SetParametersOfInterest(self.out.set("POI"))
            mc.SetObservables(self.out.set("observables"))
            nuisancesSet = ROOT.RooArgSet()
            if self.out.set("nuisances"):
                nuisancesSet = self.out.set("nuisances")
            for nuis in self.extraNuisances:
                nuisancesSet.add(nuis)
            if nuisancesSet.getSize():
                mc.SetNuisanceParameters(nuisancesSet)
            gObsSet = ROOT.RooArgSet()
            if self.out.set("globalObservables"):
                gObsSet = self.out.set("globalObservables")
            for gobs in self.extraGlobalObservables:
                gObsSet.add(gobs)
            if gObsSet.getSize():
                mc.SetGlobalObservables(gObsSet)
            if self.options.verbose > 2:
                mc.Print("V")
            self.out.safe_import(mc, mc.GetName())
            if self.options.noBOnly:
                break
        discparams = ROOT.RooArgSet("discreteParams")
        for cpar in self.discrete_param_set:
            discparams.add(self.out.cat(cpar))
        self.out.safe_import(discparams, discparams.GetName())
        self.out.writeToFile(self.options.out)

    def isShapeSystematic(self, channel, process, syst):
        return False


class CountingModelBuilder(ModelBuilder):
    """ModelBuilder to make a counting experiment"""

    def __init__(self, datacard, options):
        ModelBuilder.__init__(self, datacard, options)
        if datacard.hasShapes:
            raise RuntimeError("You're using a CountingModelBuilder for a model that has shapes")

    def doObservables(self):
        if len(self.DC.obs):
            self.doComment(" ----- observables (already set to observed values) -----")
            for b in self.DC.bins:
                self.doVar("n_obs_bin%s[%f]" % (b, self.DC.obs[b]))
                self.out.var("n_obs_bin%s" % b).setMin(0)
        else:
            self.doComment(" ----- observables -----")
            for b in self.DC.bins:
                self.doVar("n_obs_bin%s[1]" % b)
                self.out.var("n_obs_bin%s" % b).setMin(0)

        self.doSet("observables", ",".join(["n_obs_bin%s" % b for b in self.DC.bins]))
        if len(self.DC.obs):
            if self.options.bin:
                self.out.data_obs = ROOT.RooDataSet(self.options.dataname, "observed data", self.out.set("observables"))
                self.out.data_obs.add(self.out.set("observables"))
                self.out.safe_import(self.out.data_obs)

    def doIndividualModels(self):
        self.doComment(" --- Expected events in each bin, total (S+B and B) ----")
        for b in self.DC.bins:
            self.doObj(
                "n_exp_bin%s_bonly" % b,
                "sum",
                ", ".join(["n_exp_bin%s_proc_%s" % (b, p) for p in self.DC.exp[b].keys() if not self.DC.isSignal[p]]),
            )
            self.doObj(
                "n_exp_bin%s" % b,
                "sum",
                ", ".join(["n_exp_bin%s_proc_%s" % (b, p) for p in self.DC.exp[b].keys()]),
            )
            self.doObj("pdf_bin%s" % b, "Poisson", "n_obs_bin%s, n_exp_bin%s, 1" % (b, b))
            self.doObj(
                "pdf_bin%s_bonly" % b,
                "Poisson",
                "n_obs_bin%s, n_exp_bin%s_bonly, 1" % (b, b),
            )

    def doCombination(self):
        prefix = "modelObs" if len(self.DC.systs) else "model"  # if no systematics, we build directly the model
        nbins = len(self.DC.bins)
        if nbins > 50:
            from math import ceil

            nblocks = int(ceil(nbins / 10.0))
            for i in range(nblocks):
                self.doObj(
                    "%s_s_%d" % (prefix, i),
                    "PROD",
                    ",".join(["pdf_bin%s" % self.DC.bins[j] for j in range(10 * i, min(nbins, 10 * i + 10))]),
                )
                self.doObj(
                    "%s_b_%d" % (prefix, i),
                    "PROD",
                    ",".join(["pdf_bin%s_bonly" % self.DC.bins[j] for j in range(10 * i, min(nbins, 10 * i + 10))]),
                )
            self.doObj(
                "%s_s" % prefix,
                "PROD",
                ",".join([prefix + "_s_%d" % i for i in range(nblocks)]),
            )
            self.doObj(
                "%s_b" % prefix,
                "PROD",
                ",".join([prefix + "_b_%d" % i for i in range(nblocks)]),
            )
        else:
            self.doObj(
                "%s_s" % prefix,
                "PROD",
                ",".join(["pdf_bin%s" % b for b in self.DC.bins]),
            )
            self.doObj(
                "%s_b" % prefix,
                "PROD",
                ",".join(["pdf_bin%s_bonly" % b for b in self.DC.bins]),
            )
        if len(self.DC.systs):  # multiply by nuisances if needed
            self.doObj("model_s", "PROD", "modelObs_s, nuisancePdf")
            self.doObj("model_b", "PROD", "modelObs_b, nuisancePdf")
