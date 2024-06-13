import collections, itertools, math, re
import numpy as np
from HiggsAnalysis.CombinedLimit.PhysicsModel import  CanTurnOffBkgModel, MultiSignalModel, PhysicsModelBase_NiceSubclasses

### This is the base python class to study the SpinZero structure

class SpinZeroHiggsBaseAC(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        super(SpinZeroHiggsBaseAC, self).__init__()

        self.faiphiaistatus = collections.defaultdict(lambda: "fix")
        self.fairelative = collections.defaultdict(lambda: False)
        self.allowPMF = collections.defaultdict(lambda: False)
        self.allowPMF[0] = False

        self.HWWcombination = False

        self.offshell = False

    def setModelBuilder(self, modelBuilder):
        super(SpinZeroHiggsBaseAC, self).setModelBuilder(modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        result = super(SpinZeroHiggsBaseAC, self).getYieldScale(bin, process)
        if result not in (0, 1): print "Process {0} will scale by {1}".format(process,result)
        return result

    def processPhysicsOptions(self,physOptions):
        processed = super(SpinZeroHiggsBaseAC, self).processPhysicsOptions(physOptions)

        if {self.faidefinitionorder(i) for i in xrange(self.numberoffais)} != set(xrange(self.numberoffais)):
            raise ValueError("faidefinitionorder is not defined right, should go from 0 to {0} for inputs from 0 to {0}\n{1}\n{2}".format(self.numberoffais-1, {self.faidefinitionorder(i) for i in xrange(self.numberoffais)}, set(xrange(self.numberoffais))))

        for po in physOptions:
            newpo = po.lower().replace("a1", "ai0")
            match = re.match("(f|phi)ai([0-9]+)(fixed|notpoi|floating|aspoi)((?:relative)?)$", newpo)
            if match:
                parametertype = match.group(1)
                i = int(match.group(2))
                whattodo = match.group(3)
                relative = bool(match.group(4))

                if not 0 <= i < self.numberoffais:
                    raise ValueError("There are only {} fais available, so can't do anything with {}".format(self.numberoffais-1, po))

                if relative:
                    if whattodo == "fixed": raise ValueError("relative doesn't make sense for fixed fais")
                    if parametertype == "phi": raise ValueError("relative doesn't make sense for phi")

                key = parametertype, i

                if key in self.faiphiaistatus:
                    raise ValueError("Specified multiple physics options for {}ai{}".format(parametertype, i).replace("ai0", "a1"))
                if self.faidefinitionorder(i) == self.numberoffais-1:
                    if parametertype == "phi":
                        raise ValueError("fai{} is the last parameter to be defined, so it doesn't have a phi.  Can't do anything with {}".format(i, po).replace("ai0", "a1"))

                if whattodo == "fixed":
                    self.faiphiaistatus[key] = "fix"
                    if parametertype == "f" and i == 1:
                        print "Will fix {} to 0".format(self.parametername(parametertype, i, relative))
                elif whattodo == "floating" or whattodo == "notpoi":
                    self.faiphiaistatus[key] = "float"
                    if parametertype == "f" and i == 1:
                        print "{} is NOT A POI".format(self.parametername(parametertype, i, relative))
                    else:
                        print "Will float {}".format(self.parametername(parametertype, i, relative))
                elif whattodo == "aspoi":
                    self.faiphiaistatus[key] = "POI"
                    print "Will consider {} as a parameter of interest".format(self.parametername(parametertype, i, relative))
                else:
                    assert False, whattodo

                if relative:
                    self.fairelative[key] = True

                processed.append(po)

            if po.lower() == 'allowpmf':
                self.allowPMF.default_factory = lambda: True
                processed.append(po)

            match = re.match("(allow|forbid)pmfai([0-9]+)", newpo)
            if match:
                i = int(match.group(2))
                toset = {"allow": True, "forbid": False}[match.group(1)]
                if i == 0:
                    if toset == True:
                        raise ValueError("fa1 has to be fixed to positive")
                else:
                    if i in self.allowPMF: raise ValueError("Two different options set for allow/forbidPMfai{}".format(i).replace("ai0", "a1"))
                    self.allowPMF[i] = toset
                processed.append(po)

            if po.lower() == 'hwwcombination':
                self.HWWcombination = True
                processed.append(po)

            if po.lower() == "offshell":
                self.offshell = True
                processed.append(po)

        return processed

    @property
    def numberoffais(self):
        return 3 #including fa1

    def faidefinitionorder(self, i):
        if i == 0: return self.numberoffais-1
        return i-1
    def faidefinitionorderinverse(self, i):
        for _ in xrange(self.numberoffais):
            if self.faidefinitionorder(_) == i:
                return _
        raise ValueError("faidefinitionorder doesn't have anything that gives {}".format(i))
    def parametername(self, parametertype, i, relative):
      if parametertype == "phi" and self.faidefinitionorder(i) == self.numberoffais-1: raise ValueError("Can't have phiai for the last fai")
      if i >= self.numberoffais: raise ValueError("Only have {} fais".format(self.numberoffais-1))
      if i < 0: raise ValueError("Only have positive integer fais, plus 0 for fa1")
      return "CMS_zz4l_{}ai{:d}{}".format(parametertype, i, "_relative" if relative else "").replace("ai0", "a1")

    def getPOIList(self):
        poi = []
        poi += super(SpinZeroHiggsBaseAC, self).getPOIList()

        for i in xrange(self.numberoffais):
            varname = self.parametername("f", i, False)
            if self.faidefinitionorder(i) == self.numberoffais-1:
                if self.faiphiaistatus["f", i] == "POI":
                    if not self.allowPMF[i]:
                        raise ValueError("fai{} is the last parameter to be defined, so it's defined as a function of the others.  It's sign is fixed.  Can't set it as a POI")
                    self.faiphiaistatus["f", i] = "lastPOI"
                else:
                    self.faiphiaistatus["f", i] = "last"
        if ("f", 1) not in self.faiphiaistatus: self.faiphiaistatus["f", 1] = "POI"

        for parametertype in "f", "phi":
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1: continue

                varname = self.parametername(parametertype, i, self.fairelative[parametertype, i])

                if not self.modelBuilder.out.var(varname):
                    self.modelBuilder.doVar(varname+"[0,0,1]") #will set the range later

                self.modelBuilder.out.var(varname).setVal(1 if i==0 else 0) #set fa1 to 1, anomalous couplings to 0

        parametertype = "f"
        done = {i: False for i in xrange(self.numberoffais)}
        while not all(done.values()):
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1: done[i] = True
                if done[i]: continue
                if self.fairelative[parametertype, i]:
                    otheris = [j for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i)]
                    if not all(done[j] for j in otheris): continue
                    expr = "-".join(["1"] + ["abs(@{})".format(k) for k, j in enumerate(otheris)])
                    expr = "(" + expr + ")" + " * @{}".format(len(otheris))
                    fais = ", ".join(
                      [self.parametername("f", j, False) for j in otheris]
                      + [self.parametername("f", i, True)]
                    )
                    self.modelBuilder.doVar('expr::{}("{}", {})'.format(self.parametername("f", i, False), expr, fais))
                done[i] = True


        for parametertype in "f", "phi":
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1 and parametertype == "phi": continue
                relative = self.fairelative[parametertype, i]
                varname = self.parametername(parametertype, i, relative)
                status = self.faiphiaistatus[parametertype, i]
                if status in ("float", "POI"):
                    if parametertype == "f":
                        if self.faidefinitionorder(i) == 0 or relative:
                            parameterrange = (-1 if self.allowPMF[i] else 0), 1
                        else:
                            expr = "-".join(["1"] + ["abs(@{})".format(k) for k, j in enumerate(j for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i))])
                            fais = ", ".join(self.parametername("f", j, False) for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i))
                            self.modelBuilder.doVar('expr::max_'+varname+'("{}", {})'.format(expr, fais))
                            if self.allowPMF[i]:
                                self.modelBuilder.doVar('expr::min_{0}("-@0", max_{0})'.format(varname))
                            else:
                                self.modelBuilder.doVar('expr::min_{0}("0")'.format(varname))
                            parameterrange = (
                                self.modelBuilder.out.obj("min_"+varname),
                                self.modelBuilder.out.obj("max_"+varname),
                            )
                    elif parametertype == "phi":
                        parameterrange = -math.pi, math.pi
                    else:
                        assert False
                    self.modelBuilder.out.var(varname).setRange(*parameterrange)
                    self.modelBuilder.out.var(varname).setConstant(False)
                    if status == "POI":
                        print "Treating "+varname+" as a POI"
                        poi.append(varname)
                    else:
                        print "Floating "+varname
                        self.modelBuilder.out.var(varname).setAttribute("flatParam")
                    if parametertype == "f" and self.allowPMF[i]: print "Allowing negative "+varname
                elif status == "fix":
                    self.modelBuilder.out.var(varname).setConstant()
                    print "Fixing "+varname
                elif status in ("last", "lastPOI"):
                    expr = "-".join(["1"] + ["abs(@{})".format(j) for j in xrange(self.numberoffais-1)])
                    fais = ", ".join(self.parametername(parametertype, j, relative) for j in xrange(self.numberoffais) if j != i)

                    if self.allowPMF[i]:
                        self.modelBuilder.doVar('sgn{}[negative,positive]'.format(varname))
                        expr = "({}) * 2 * (@{}-0.5)".format(expr, self.numberoffais-1)
                        fais += ", sgn"+varname
                        self.modelBuilder.out.cat("sgn"+varname).setAttribute("flatParam")
                        self.modelBuilder.addDiscrete("sgn"+varname)
                        print "Floating the sign of "+varname

                    self.modelBuilder.doVar('expr::{}("{}", {})'.format(varname, expr, fais))
                    print "Setting "+varname+" to 1 - the sum of the other fais"
                else:
                    assert False, status

        if self.HWWcombination:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                print "Found CMS_zz4l_alpha; setting to 0.5"
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0.5)
            else:
                print "Creating CMS_zz4l_alpha; setting to 0.5"
                self.modelBuilder.doVar("CMS_zz4l_alpha[0.5,-1,1]")
            print "Treating alpha (Rai1) as a POI"
            self.modelBuilder.out.var("CMS_zz4l_alpha").setConstant(False)
            poi.append("CMS_zz4l_alpha")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                print "Found CMS_zz4l_alpha; setting to constant 0"
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_alpha").setConstant()

        # GGsm offshell variable
        if self.modelBuilder.out.var("CMS_zz4l_GGsm"):
            print "CMS_zz4l_GGsm is being renamed to GGsm"
            self.modelBuilder.out.var("CMS_zz4l_GGsm").SetName("GGsm")
        if self.offshell:
            if self.modelBuilder.out.var("GGsm"):
                print "Found GGsm, setting range to [1., 0., 50.]"
                self.modelBuilder.out.var("GGsm").setRange(0, 50)
                self.modelBuilder.out.var("GGsm").setVal(1)
            else:
                print "Creating GGsm; setting range to [1., 0., 50.]"
                self.modelBuilder.doVar("GGsm[1,0,50]")
            print "GGsm is a POI."
            self.modelBuilder.out.var("GGsm").setConstant(False)
            poi.append("GGsm")
        else:
            if self.modelBuilder.out.var("GGsm"):
                print "Found GGsm, fixing to 1"
                self.modelBuilder.out.var("GGsm").setVal(1)
                self.modelBuilder.out.var("GGsm").setConstant()

        return poi

class SpinZeroHiggsAC(SpinZeroHiggsBaseAC):
    def __init__(self):
        super(SpinZeroHiggsAC, self).__init__()
        self.muFloating = True
        self.muAsPOI = False

    def processPhysicsOptions(self,physOptions):
        processed = super(SpinZeroHiggsAC, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.lower() == 'mufixed':
                print "Will consider the signal strength as a fixed parameter"
                self.muFloating = False
                processed.append(po)
            elif po.lower() == 'muaspoi':
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                processed.append(po)

        if self.muAsPOI and not self.muFloating:
            raise ValueError("Specified both muFixed and muAsPOI!")

        if not self.muAsPOI and not self.fai1POI and not self.fai2POI and not self.phiai1POI and not self.phiai2POI:
            print "No POIs detected: Switching to default configuration: Floating nuisance mu, floating POI fai1, eveything else fixed"
            self.muFloating = True
            self.muAsPOI = False
            self.fai1Floating = True
            self.fai1POI = True
            self.phiai1Floating = False
            self.phiai1POI = False
            self.fai2Floating = False
            self.fai2POI = False
            self.phiai2Floating = False
            self.phiai2POI = False
            self.allowPMF = defaultdict(lambda: False)

        return processed

    def getPOIList(self):
        poi = super(SpinZeroHiggsAC, self).getPOIList()
        if self.muFloating:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setRange(0.,400.)
                self.modelBuilder.out.var("r").setVal(1)
            else:
                self.modelBuilder.doVar("r[1,0,400]")
            if self.HWWcombination:
                self.modelBuilder.out.var("r").removeMax()
                print "Removed maximum of r"

            self.modelBuilder.out.var("r").setConstant(False)

            if self.muAsPOI:
                print "Treating r as a POI"
                poi.append("r")
            else:
                self.modelBuilder.out.var("r").setAttribute("flatParam")
        else:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setVal(1)
                self.modelBuilder.out.var("r").setConstant()
            else:
                self.modelBuilder.doVar("r[1]")
        return poi

class SpinZeroHiggsAiBase(SpinZeroHiggsBaseAC):
    def __init__(self):
        super(SpinZeroHiggsAiBase, self).__init__()
        self.aistatus = collections.defaultdict(lambda: "fix")
        self.floatgamma = True

    def processPhysicsOptions(self,physOptions):
        processed = super(SpinZeroHiggsAiBase, self).processPhysicsOptions(physOptions)

        for po in physOptions:
            newpo = po.lower().replace("a1", "ai0")
            match = re.match("ai([0-9]+)(fixed|notpoi|floating|aspoi)$", newpo)
            if match:
                i = int(match.group(1))
                whattodo = match.group(2)

                if not 0 <= i < self.numberoffais:
                    raise ValueError("There are only {} fais available, so can't do anything with {}".format(self.numberoffais-1, po))

                if i in self.aistatus:
                    raise ValueError("Specified multiple physics options for ai{}".format(i).replace("ai0", "a1"))

                if whattodo == "fixed":
                    self.aistatus[i] = "fix"
                    print "Will fix ai{} to 0".format(i)
                elif whattodo == "floating" or whattodo == "notpoi":
                    self.aistatus[i] = "float"
                    print "Will float ai{}".format(i)
                elif whattodo == "aspoi":
                    self.aistatus[i] = "POI"
                    print "Will consider ai{} as a parameter of interest".format(i)
                else:
                    assert False, whattodo

                processed.append(po)

            if po.lower() == "fixgamma":
                self.fixgamma = True
                processed.append(po)


        if 1 not in self.aistatus: self.aistatus[1] = "POI"

        return processed

    def getPOIList(self):
        poi = []
        poi += super(SpinZeroHiggsAiBase, self).getPOIList()

        if 1 not in self.aistatus: self.aistatus[1] = "POI"

        for i in xrange(self.numberoffais):
            varname = self.getcouplingname(self.sortedcouplings[i])

            if not self.modelBuilder.out.var(varname):
                self.modelBuilder.doVar(varname+"[0,0,1]") #will set the range later                                                                 

            self.modelBuilder.out.var(varname).setVal(1 if varname == "g1" else 0) #set a1 to 1, anomalous couplings to 0                            

        for i in xrange(self.numberoffais):
            varname = self.getcouplingname(self.sortedcouplings[i])
            status = self.aistatus[i]
            if status in ("float", "POI"):
                self.modelBuilder.out.var(varname).removeMax()
                if self.allowPMF[i]: self.modelBuilder.out.var(varname).removeMin()
                self.modelBuilder.out.var(varname).setConstant(False)
                if status == "POI":
                    print "Treating "+varname+" as a POI"
                    poi.append(varname)
                else:
                    print "Floating "+varname
                    self.modelBuilder.out.var(varname).setAttribute("flatParam")
                if self.allowPMF[i]: print "Allowing negative "+varname
            elif status == "fix":
                print "Fixing "+varname
                self.modelBuilder.out.var(varname).setConstant()
            else:
                assert False, status

        return poi





class MultiSignalSpinZeroHiggsAC(SpinZeroHiggsBaseAC,CanTurnOffBkgModel,MultiSignalModel):
    def __init__(self):
        super(MultiSignalSpinZeroHiggsAC, self).__init__()

        self.scalemuvfseparately = True
        self.scaledifferentsqrtsseparately = False
        self.uservoverrf = False
        self.sqrts = None
        self.fixed = []
        self.floated = []
        self.noRV = self.noRF = False
        #not doing muAsPOI or fixMu, there are too many permutations.
        #should just be set when running combine.

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("map=") for po in physOptions):
            #no po started with map --> no manual overriding --> use the defaults
            #can still override with e.g. turnoff=ZH,WH
            physOptions = ["map=.*/(gg|qq|Z|W|V|tt|bb)H$:1"] + physOptions
        super(MultiSignalSpinZeroHiggsAC, self).setPhysicsOptions(physOptions)

    def processPhysicsOptions(self, physOptions):
        processed = super(MultiSignalSpinZeroHiggsAC, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.lower() == "scalemuvmuftogether":
                self.scalemuvfseparately = False
                processed.append(po)
            if po.lower() == "scaledifferentsqrtsseparately":
                self.scaledifferentsqrtsseparately = True
                processed.append(po)
            if po.lower().startswith("sqrts="):
                if self.sqrts is not None: raise ValueError("Duplicate physicsoption sqrts=?? provided")
                self.sqrts = [int(_) for _ in po.replace("sqrts=", "").split(",")]
                processed.append(po)
            if po.lower() == "uservoverrf":
                self.uservoverrf = True
                processed.append(po)
            if po.lower() == "norv":
                self.noRV = True
                processed.append(po)
            if po.lower() == "norf":
                self.noRF = True
                processed.append(po)

        if self.uservoverrf and not self.scalemuvfseparately:
            raise ValueError("can't specify both uservoverrf and scalemuvmuftogether")

        if self.sqrts is None:
            raise ValueError("PhysicsOption sqrts=?? is mandatory.  example: sqrts=7,8,13")

        if self.scaledifferentsqrtsseparately and self.scalemuvfseparately:
            if self.uservoverrf:
                self.fixed = ["RV", "RF", "R"] + ["RF_{}TeV".format(_) for _ in self.sqrts]
                self.floated = ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "") for _2 in self.sqrts]
            else:
                self.fixed = ["RV", "RF", "R"] + ["R_{}TeV".format(_) for _ in self.sqrts]
                self.floated = ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
        elif self.scaledifferentsqrtsseparately and not self.scalemuvfseparately:
            self.fixed = ["RV", "RF", "R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
            self.floated = ["R_{}TeV".format(_) for _ in self.sqrts]
        elif not self.scaledifferentsqrtsseparately and self.scalemuvfseparately:
            if self.uservoverrf:
                self.fixed = ["RF"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.floated = ["RV", "R"]
            else:
                self.fixed = ["R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.floated = ["RV", "RF"]
        elif not self.scaledifferentsqrtsseparately and not self.scalemuvfseparately:
            self.fixed = ["RV", "RF"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
            self.floated = ["R"]
        else:
            assert False, "?????"

        if self.noRV:
            self.fixed = [_ for _ in self.fixed if "RV" not in _]
            self.floated = [_ for _ in self.floated if "RV" not in _]
        if self.noRF:
            self.fixed = [_ for _ in self.fixed if "RF" not in _]
            self.floated = [_ for _ in self.floated if "RF" not in _]

        return processed

    def getPOIList(self):
        result = super(MultiSignalSpinZeroHiggsAC, self).getPOIList()

        fixedorfloated = self.fixed+self.floated
        for variable in fixedorfloated:
            if not self.modelBuilder.out.var(variable):
                if variable in self.fixed: continue
                raise RuntimeError("{} does not exist in the workspace!  Check:\n - your datacard maker\n - your sqrts option".format(variable))
            else:
                if 'r' in variable.lower():
                    print "Setting {} range to [1.,0.,400.]".format(variable)
                    self.modelBuilder.out.var(variable).setRange(0.,400.)
                    self.modelBuilder.out.var(variable).setVal(1)
                elif "ggsm" in variable.lower():
                    print "Setting {} range to [1.,0.,50.]".format(variable)
                    self.modelBuilder.out.var(variable).setRange(0.,50.)
                    self.modelBuilder.out.var(variable).setVal(1)
                else:
                    print "Setting {} value to 0".format(variable)
                    self.modelBuilder.out.var(variable).setVal(0)
        for variable in self.fixed:
            if self.modelBuilder.out.var(variable):
                print "Fixing {}".format(variable)
                self.modelBuilder.out.var(variable).setConstant()
        for variable in self.floated:
            if self.modelBuilder.out.var(variable):
                print "Floating {} and assigning attribute flatParam".format(variable)
                self.modelBuilder.out.var(variable).setConstant(False)
                self.modelBuilder.out.var(variable).setAttribute("flatParam")

        return result

class HZZAnomalousCouplingsFromHistogramsBase(SpinZeroHiggsBaseAC):
    """                                                                                                                                                                                                                              
    This class expects histograms (which could be TH1 or RooDataHist) for each component of the PDF.                                                                                                                                 
    The histograms should be normalized to ai = 1.                                                                                                                                                                                   
                                                                                                                                                                                                                                     
    The pure components should be named ggH_0PM, qqH_0L1Zg, etc.                                                                                                                                                                     
                                                                                                                                                                                                                                     
    The interference components should be split in two, one for positive bins and one for negative bins flipped.                                                                                                                     
    This way all bins are positive.                                                                                                                                                                                                  
    They should be named ZH_g13g1prime21_positive or ttH_g11g41_negative, for example.                                                                                                                                               
    And they should also be normalized to ai = 1, for all ais involved in that term.                                                                                                                                                 
                                                                                                                                                                                                                                     
    The only exception is the L1 and L1Zg terms, which should be scaled to g1prime2=10000 for HZZ and HZgamma                                                                                                                        
                                                                                                                                                                                                                                     
    For anomalous fermion couplings, call the histogram, for example, ttH_0PMff_0L1 or ggH_0Mff_g41g21_negative                                                                                                                      
    Interference between fermion couplings is not implemented yet                                                                                                                                                                    
    In that case the histograms should be normalized to kappa=1 or kappa_tilde=1 for ttH, or g2=1 or g4=1 for ggH                                                                                                                    
                                                                                                                                                                                                                                     
    By default, ttH and ggH anomalous couplings will be related.  The physicsmodel will take care of the scaling from                                                                                                                
    kappa_tilde to g4                                                                                                                                                                                                                
    """

    kappa_tilde_ttH = 1.6

    def __init__(self):
        self.anomalouscouplings = []
        self.turnoff = []
        self.scalegL1by10000 = True
        self.useHffanomalous = False
        self.separateggHttH = False
        self.fixgamma = False
        super(HZZAnomalousCouplingsFromHistogramsBase, self).__init__()

    def processPhysicsOptions(self,physOptions):
        processed = []
        for po in physOptions:
            if po in ("fa3", "fa2", "fL1", "fL1Zg", "fa1"):
                if po in self.anomalouscouplings: raise ValueError("Provided physOption "+po+" twice")
                self.anomalouscouplings.append(po)
                processed.append(po)

            if po.lower() == "separategghtth":
                self.separateggHttH = True
                self.noRF = True
                processed.append(po)

        if "fa1" not in self.anomalouscouplings: self.anomalouscouplings.append("fa1")

        for po in physOptions[:]:
            for i, fai in enumerate(self.anomalouscouplings):
                if fai == "fa1": continue #handled in the base class                                                                                                                                                                 
                ai = fai[1:]
                if re.search("{}(fixed|notpoi|floating|aspoi)".format(ai).lower(), po.lower()):
                    physOptions.append(po.replace(ai, "ai{}".format(self.faidefinitionorderinverse(i))))
                    processed.append(po)

        processed += super(HZZAnomalousCouplingsFromHistogramsBase, self).processPhysicsOptions(physOptions)

        if not self.anomalouscouplings: raise ValueError("Have to provide an anomalous coupling as a physOption (fa3, fa2, fL1, fL1Zg)")
        return processed

    @property
    def numberoffais(self):
        return len(self.anomalouscouplings)  #including fa1                                                                                                                                                                          

    @property
    def sortedcouplings(self):
        return sorted(self.anomalouscouplings, key=["fa1", "fa3", "fa2", "fL1", "fL1Zg"].index)

    def faidefinitionorder(self, i):
        #CMS_zz4l_fai1, CMS_zz4l_fai2, etc. correspond to fa3, fa2, fL1, fL1Zg in that order                                                                                                                                         
        #However they might not be defined in that order, e.g. CMS_zz4l_fai1 might be restricted to (0, 1-CMS_zz4l_fai2)                                                                                                             
        return self.anomalouscouplings.index(self.sortedcouplings[i])

    signalprocessregex = (
        "(?P<production>gg|tt|bb|qq|Z|W|V)H_"
        "(?:(?P<Hffpure>0(?:PM|M)ff)_)?"
        "(?:"
          "(?P<HVVpure>0(?:PM|M|PH|L1|L1Zg))|"
          "(?P<HVVint>(?:g(?:1|2|4|1prime2|hzgs1prime2)[1234])*)_(?P<HVVintsign>positive|negative)"
        ")$"
    )

    @staticmethod
    def getcouplingname(processorfai, production=None):
        if processorfai == "0PMff": return {"gg": "ghg2", "tt": "kappa", "bb": "kappa"}[production]
        if processorfai == "0Mff": return {"gg": "ghg4", "tt": "kappatilde", "bb": "kappatilde"}[production]

        return {
            "0PM": "g1",
            "0PH": "g2",
            "0M": "g4",
            "0L1": "g1prime2",
            "0L1Zg": "ghzgs1prime2",
            "fa1": "g1",
            "fa2": "g2",
            "fa3": "g4",
            "fL1": "g1prime2",
            "fL1Zg": "ghzgs1prime2",
        }[processorfai]

    def tellAboutProcess(self, bin, process):
        match = re.match(self.signalprocessregex, process)
        if match and match.group("Hffpure"):
            self.useHffanomalous = True

    def getYieldScale(self,bin,process):
        match = re.match(self.signalprocessregex, process)

        if not match:
            if any(process.startswith(_) for _ in ("ggH", "ttH", "bbH", "qqH", "ZH", "WH", "VH")):
                raise ValueError("Your signal process "+process+" doesn't match the pattern")
            return super(HZZAnomalousCouplingsFromHistogramsBase, self).getYieldScale(bin, process)

        if match.group("production")+"H" in self.turnoff: return 0

        if match.group("production") == "gg": maxpower = 2; production = "ggHVV"
        elif match.group("production") in ("tt", "bb"): maxpower = 2; production = "ttHVV"
        elif match.group("production") in ("qq", "Z", "W", "V"): maxpower = 4; production = "VVHVV"

        result = production

        if match.group("Hffpure") is not None:
            if match.group("production") not in ("gg", "tt", "bb"): raise ValueError("Don't put fermion couplings for {}H: {}".format(match.group("production"), process))
            Hffpowerdict = {self.getcouplingname(match.group("Hffpure"), match.group("production")): 2}
            result += "_" + "".join("{}{}".format(k, v) for k, v in Hffpowerdict.iteritems())

        if match.group("HVVpure") is not None:
            powerdict = {self.getcouplingname(match.group("HVVpure")): maxpower}
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems())
        elif match.group("HVVint") is not None:
            powerdict = {coupling: int(power) for coupling, power in re.findall("(g(?:1|2|4|1prime2|hzgs1prime2))([1234])", match.group("HVVint"))}

            if sum(powerdict.values()) != maxpower:
                raise ValueError("power dict doesn't add up properly!  Sum should be {}\n{}\n{}".format(maxpower, process, powerdict))

            powerdict = collections.OrderedDict(
                sorted(powerdict.iteritems(), key = lambda x: "g1 g4 g2 g1prime2 ghzgs1prime2".index(x[0]))
            )

            sign = match.group("HVVintsign")
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems()) + "_" + sign
        else:
            assert False

        if self.verbose:
            print "Process {0} will scale by {1}".format(process,result)

        return result


class HZZAnomalousCouplingsFromHistogramsNonSMEFT(MultiSignalSpinZeroHiggsAC):
    """
    This class expects histograms (which could be TH1 or RooDataHist) for each component of the PDF.
    The histograms should be normalized to ai = 1.

    The pure components should be named ggH_0PM, qqH_0L1Zg, etc.

    The interference components should be split in two, one for positive bins and one for negative bins flipped.
    This way all bins are positive.
    They should be named ZH_g13g1prime21_positive or ttH_g11g41_negative, for example.
    And they should also be normalized to ai = 1, for all ais involved in that term.

    The only exception is the L1 and L1Zg terms, which should be scaled to g1prime2=10000 for HZZ and HZgamma

    For anomalous fermion couplings, call the histogram, for example, ttH_0PMff_0L1 or ggH_0Mff_g41g21_negative
    Interference between fermion couplings is not implemented yet
    In that case the histograms should be normalized to kappa=1 or kappa_tilde=1 for ttH, or g2=1 or g4=1 for ggH

    By default, ttH and ggH anomalous couplings will be related.  The physicsmodel will take care of the scaling from
    kappa_tilde to g4
    """

    '''

    HERE FOR EFT XSECTIONS
    eft xsections
    aidecay = {
      "g2": 0.394465808268,
      "g4": 2.55052,
      "g1prime2": -4363.84210717,
      "ghzgs1prime2": -7613.351302119843,
    }

    '''
    #non eft crossections
    aidecay = {
      "g2": 1.65684,
      "g4": 2.55052,
      "g1prime2": -12100.42,
      "ghzgs1prime2": -7613.351302119843,
    }


    
    kappa_tilde_ttH = 1.6

    def __init__(self):
        self.anomalouscouplings = []
        self.turnoff = []
        self.scalegL1by10000 = True
        self.useHffanomalous = False
        self.separateggHttH = False
        super(HZZAnomalousCouplingsFromHistogramsNonSMEFT, self).__init__()

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("sqrts=") for po in physOptions):
            physOptions = physOptions + ["sqrts=13"]
        for po in physOptions:
            if po.startswith("turnoff="):
                self.turnoff += po.replace("turnoff=", "").split(",")
                #po gets removed in super
        super(MultiSignalSpinZeroHiggsAC, self).setPhysicsOptions(physOptions)
        if self.sqrts != [13]:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is set up for 13 TeV only")
        if self.scaledifferentsqrtsseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scaledifferentsqrtsseparately")
        if not self.scalemuvfseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scalemuvmuftogether")

    def processPhysicsOptions(self,physOptions):
        processed = []
        for po in physOptions:
            if po in ("fa3", "fa2", "fL1", "fL1Zg", "fa1"):
                if po in self.anomalouscouplings: raise ValueError("Provided physOption "+po+" twice")
                self.anomalouscouplings.append(po)
                processed.append(po)

            if po.lower() == "separategghtth":
                self.separateggHttH = True
                self.noRF = True
                processed.append(po)

        if "fa1" not in self.anomalouscouplings: self.anomalouscouplings.append("fa1")

        for po in physOptions[:]:
            for i, fai in enumerate(self.anomalouscouplings):
                if fai == "fa1": continue #handled in the base class
                ai = fai[1:]
                if re.match("(f|phi){}(fixed|notpoi|floating|aspoi)(?:relative)?$".format(ai).lower(), po.lower()):
                    physOptions.append(po.replace(ai, "ai{}".format(self.faidefinitionorderinverse(i))))
                    processed.append(po)

        processed += super(HZZAnomalousCouplingsFromHistogramsNonSMEFT, self).processPhysicsOptions(physOptions)

        if not self.anomalouscouplings: raise ValueError("Have to provide an anomalous coupling as a physOption (fa3, fa2, fL1, fL1Zg)")
        return processed

    @property
    def numberoffais(self):
        return len(self.anomalouscouplings)  #including fa1

    @property
    def sortedcouplings(self):
        return sorted(self.anomalouscouplings, key=["fa1", "fa3", "fa2", "fL1", "fL1Zg"].index)

    def faidefinitionorder(self, i):
        #CMS_zz4l_fai1, CMS_zz4l_fai2, etc. correspond to fa3, fa2, fL1, fL1Zg in that order
        #However they might not be defined in that order, e.g. CMS_zz4l_fai1 might be restricted to (0, 1-CMS_zz4l_fai2)
        return self.anomalouscouplings.index(self.sortedcouplings[i])

    def getPOIList(self):
        if self.useHffanomalous:
            self.modelBuilder.doVar("fCP_Htt[0.,-1,1]")
            self.modelBuilder.out.var("fCP_Htt").setConstant(False)
            self.modelBuilder.out.var("fCP_Htt").setAttribute("flatParam")

        if self.separateggHttH:
            self.modelBuilder.doVar("Rg[1.0,0,400]")
            self.modelBuilder.doVar("Rt[1.0,0,400]")
            self.modelBuilder.out.var("Rg").setConstant(False)
            self.modelBuilder.out.var("Rg").setAttribute("flatParam")
            self.modelBuilder.out.var("Rt").setConstant(False)
            self.modelBuilder.out.var("Rt").setAttribute("flatParam")

            if self.useHffanomalous:
                self.modelBuilder.doVar("fa3_ggH[0.,-1,1]")
                self.modelBuilder.out.var("fa3_ggH").setConstant(False)
                self.modelBuilder.out.var("fa3_ggH").setAttribute("flatParam")
        else:
            self.modelBuilder.doVar("RF[1.0,0,10]")
            if self.useHffanomalous:
                self.modelBuilder.doVar('expr::fa3_ggH("@0 == 0 ? 0 : (@0>0 ? 1 : -1) * 1. / (1. + 4./9. * (1.0 / abs(@0) - 1.))", fCP_Htt)')

        if self.useHffanomalous:
            self.modelBuilder.doVar('expr::kappa("sqrt(1-abs(@0))", fCP_Htt)')
            self.modelBuilder.doVar('expr::kappa_tilde("(@0>0 ? 1 : -1) * sqrt(abs(@0)) * {kappatildettH}", fCP_Htt)'.format(kappatildettH=self.kappa_tilde_ttH))
            self.modelBuilder.doVar('expr::ghg2("sqrt(1-abs(@0))", fa3_ggH)')
            self.modelBuilder.doVar('expr::ghg4("(@0>0 ? 1 : -1) * sqrt(abs(@0))", fa3_ggH)')


        self.modelBuilder.doVar("RV[1.0,0,10]")
        self.modelBuilder.doVar("R[1.0,0,10]")

        pois = super(HZZAnomalousCouplingsFromHistogramsNonSMEFT, self).getPOIList()

        if not self.modelBuilder.out.var("g1"):
            self.modelBuilder.doVar('expr::g1("sqrt(@0)", CMS_zz4l_fa1)')

        couplings = ["g1"]
        i = 0
        for fai in self.sortedcouplings:
            if fai == "fa1": continue
            ai = self.getcouplingname(fai)
            i += 1

            if self.scalegL1by10000:
                divideby = {
                    "g4": 1,
                    "g2": 1,
                    "g1prime2": 10000,
                    "ghzgs1prime2": 10000,
                }[ai]
            else:
                divideby = 1

            kwargs = {
              "i": i,
              "ai": ai,
              "aidecay": self.aidecay[ai] / divideby,
            }
            self.modelBuilder.doVar('expr::{ai}("(@0>0 ? 1 : -1) * sqrt(abs(@0))*{aidecay}", CMS_zz4l_fai{i})'.format(**kwargs))
            couplings.append(ai)

        if self.scaledifferentsqrtsseparately: raise ValueError("HZZAnomalousCouplingsFromHistograms is not compatible with scaledifferentsqrtsseparately")

        for g in couplings:
            if self.separateggHttH:
                self.modelBuilder.doVar('expr::ggHVV_{g}2("@0*@1*@2*@2", R, Rg, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_{g}2("@0*@1*@2*@2", R, Rt, {g})'.format(g=g))
                if self.useHffanomalous:
                    self.modelBuilder.doVar('expr::ggHVV_ghg22_{g}2("@0*@1*@2*@2*@3*@3", R, Rg, ghg2, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ggHVV_ghg42_{g}2("@0*@1*@2*@2*@3*@3", R, Rg, ghg4, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappa2_{g}2("@0*@1*@2*@2*@3*@3", R, Rt, kappa, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g}2("@0*@1*@2*@2*@3*@3", R, Rt, kappa_tilde, {g})'.format(g=g))
            else:
                self.modelBuilder.doVar('expr::ffHVV_{g}2("@0*@1*@2*@2", R, RF, {g})'.format(g=g))
                if self.useHffanomalous:
                    self.modelBuilder.doVar('expr::ggHVV_ghg22_{g}2("@0*@1*@2*@2*@3*@3", R, RF, ghg2, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ggHVV_ghg42_{g}2("@0*@1*@2*@2*@3*@3", R, RF, ghg4, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappa2_{g}2("@0*@1*@2*@2*@3*@3", R, RF, kappa, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g}2("@0*@1*@2*@2*@3*@3", R, RF, kappa_tilde, {g})'.format(g=g))
            self.modelBuilder.doVar('expr::VVHVV_{g}4("@0*@1*@2*@2*@2*@2", R, RV, {g})'.format(g=g))

        kwargs = {}
        for kwargs["signname"], kwargs["sign"] in ("positive", ""), ("negative", "-"):
            for kwargs["g1"], kwargs["g2"] in itertools.combinations(couplings, 2):
                if self.separateggHttH:
                    self.modelBuilder.doVar('expr::ggHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, Rg, {g1}, {g2})'.format(**kwargs))
                    self.modelBuilder.doVar('expr::ttHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, Rt, {g1}, {g2})'.format(**kwargs))
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rg, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rg, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rt, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rt, kappa_tilde, {g1}, {g2})'.format(**kwargs))
                else:
                    self.modelBuilder.doVar('expr::ffHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, RF, {g1}, {g2})'.format(**kwargs))
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, kappa_tilde, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}3_{signname}("{sign}@0*@1*@2*@3*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}2_{signname}("{sign}@0*@1*@2*@2*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}3{g2}1_{signname}("{sign}@0*@1*@2*@2*@2*@3", R, RV, {g1}, {g2})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"] in itertools.combinations(couplings, 3):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}2_{signname}("{sign}@0*@1*@2*@3*@4*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}2{g3}1_{signname}("{sign}@0*@1*@2*@3*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}1{g3}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"], kwargs["g4"] in itertools.combinations(couplings, 4):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}1{g4}1_{signname}("{sign}@0*@1*@2*@3*@4*@5", R, RV, {g1}, {g2}, {g3}, {g4})'.format(**kwargs))

        return pois

    signalprocessregex = (
        "(?P<production>gg|tt|bb|qq|Z|W|V)H_"
        "(?:(?P<Hffpure>0(?:PM|M)ff)_)?"
        "(?:"
          "(?P<HVVpure>0(?:PM|M|PH|L1|L1Zg))|"
          "(?P<HVVint>(?:g(?:1|2|4|1prime2|hzgs1prime2)[1234])*)_(?P<HVVintsign>positive|negative)"
        ")$"
    )

    @staticmethod
    def getcouplingname(processorfai, production=None):
        if processorfai == "0PMff": return {"gg": "ghg2", "tt": "kappa", "bb": "kappa"}[production]
        if processorfai == "0Mff": return {"gg": "ghg4", "tt": "kappatilde", "bb": "kappatilde"}[production]

        return {
            "0PM": "g1",
            "0PH": "g2",
            "0M": "g4",
            "0L1": "g1prime2",
            "0L1Zg": "ghzgs1prime2",
            "fa1": "g1",
            "fa2": "g2",
            "fa3": "g4",
            "fL1": "g1prime2",
            "fL1Zg": "ghzgs1prime2",
        }[processorfai]

    def tellAboutProcess(self, bin, process):
        match = re.match(self.signalprocessregex, process)
        if match and match.group("Hffpure"):
            self.useHffanomalous = True

    def getYieldScale(self,bin,process):
        match = re.match(self.signalprocessregex, process)

        if not match:
            if any(process.startswith(_) for _ in ("ggH", "ttH", "bbH", "qqH", "ZH", "WH", "VH")):
                raise ValueError("Your signal process "+process+" doesn't match the pattern")
            return super(HZZAnomalousCouplingsFromHistogramsNonSMEFT, self).getYieldScale(bin, process)

        if match.group("production")+"H" in self.turnoff: return 0

        if (self.separateggHttH or match.group("Hffpure") is not None) and match.group("production") == "gg": maxpower = 2; production = "ggHVV"
        elif (self.separateggHttH or match.group("Hffpure") is not None) and match.group("production") in ("tt", "bb"): maxpower = 2; production = "ttHVV"
        elif match.group("production") in ("gg", "tt", "bb"): maxpower = 2; production = "ffHVV"
        elif match.group("production") in ("qq", "Z", "W", "V"): maxpower = 4; production = "VVHVV"

        result = production

        if match.group("Hffpure") is not None:
            if match.group("production") not in ("gg", "tt", "bb"): raise ValueError("Don't put fermion couplings for {}H: {}".format(match.group("production"), process))
            Hffpowerdict = {self.getcouplingname(match.group("Hffpure"), match.group("production")): 2}
            result += "_" + "".join("{}{}".format(k, v) for k, v in Hffpowerdict.iteritems())

        if match.group("HVVpure") is not None:
            powerdict = {self.getcouplingname(match.group("HVVpure")): maxpower}
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems())
        elif match.group("HVVint") is not None:
            powerdict = {coupling: int(power) for coupling, power in re.findall("(g(?:1|2|4|1prime2|hzgs1prime2))([1234])", match.group("HVVint"))}

            if sum(powerdict.values()) != maxpower:
                raise ValueError("power dict doesn't add up properly!  Sum should be {}\n{}\n{}".format(maxpower, process, powerdict))

            powerdict = collections.OrderedDict(
                sorted(powerdict.iteritems(), key = lambda x: "g1 g4 g2 g1prime2 ghzgs1prime2".index(x[0]))
            )

            sign = match.group("HVVintsign")
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems()) + "_" + sign
        else:
            assert False

        if self.verbose:
            print "Process {0} will scale by {1}".format(process,result)

        return result

class HZZAnomalousCouplingsFromHistograms(MultiSignalSpinZeroHiggsAC):
    """
    This class expects histograms (which could be TH1 or RooDataHist) for each component of the PDF.
    The histograms should be normalized to ai = 1.

    The pure components should be named ggH_0PM, qqH_0L1Zg, etc.

    The interference components should be split in two, one for positive bins and one for negative bins flipped.
    This way all bins are positive.
    They should be named ZH_g13g1prime21_positive or ttH_g11g41_negative, for example.
    And they should also be normalized to ai = 1, for all ais involved in that term.

    The only exception is the L1 and L1Zg terms, which should be scaled to g1prime2=10000 for HZZ and HZgamma

    For anomalous fermion couplings, call the histogram, for example, ttH_0PMff_0L1 or ggH_0Mff_g41g21_negative
    Interference between fermion couplings is not implemented yet
    In that case the histograms should be normalized to kappa=1 or kappa_tilde=1 for ttH, or g2=1 or g4=1 for ggH

    By default, ttH and ggH anomalous couplings will be related.  The physicsmodel will take care of the scaling from
    kappa_tilde to g4
    """

    '''

    HERE FOR EFT XSECTIONS

    
    
    aidecay = {
      "g2": 1.65684,
      "g4": 2.55052,
      "g1prime2": -12100.42,
      "ghzgs1prime2": -7613.351302119843,
    }

    '''

    aidecay = {
      "g2": 0.394465808268,
      "g4": 2.55052,
      "g1prime2": -4363.84210717,
      "ghzgs1prime2": -7613.351302119843,
    }

    


    
    kappa_tilde_ttH = 1.6

    def __init__(self):
        self.anomalouscouplings = []
        self.turnoff = []
        self.scalegL1by10000 = True
        self.useHffanomalous = False
        self.separateggHttH = False
        super(HZZAnomalousCouplingsFromHistograms, self).__init__()

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("sqrts=") for po in physOptions):
            physOptions = physOptions + ["sqrts=13"]
        for po in physOptions:
            if po.startswith("turnoff="):
                self.turnoff += po.replace("turnoff=", "").split(",")
                #po gets removed in super
        super(MultiSignalSpinZeroHiggsAC, self).setPhysicsOptions(physOptions)
        if self.sqrts != [13]:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is set up for 13 TeV only")
        if self.scaledifferentsqrtsseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scaledifferentsqrtsseparately")
        if not self.scalemuvfseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scalemuvmuftogether")

    def processPhysicsOptions(self,physOptions):
        processed = []
        for po in physOptions:
            if po in ("fa3", "fa2", "fL1", "fL1Zg", "fa1"):
                if po in self.anomalouscouplings: raise ValueError("Provided physOption "+po+" twice")
                self.anomalouscouplings.append(po)
                processed.append(po)

            if po.lower() == "separategghtth":
                self.separateggHttH = True
                self.noRF = True
                processed.append(po)

        if "fa1" not in self.anomalouscouplings: self.anomalouscouplings.append("fa1")

        for po in physOptions[:]:
            for i, fai in enumerate(self.anomalouscouplings):
                if fai == "fa1": continue #handled in the base class
                ai = fai[1:]
                if re.match("(f|phi){}(fixed|notpoi|floating|aspoi)(?:relative)?$".format(ai).lower(), po.lower()):
                    physOptions.append(po.replace(ai, "ai{}".format(self.faidefinitionorderinverse(i))))
                    processed.append(po)

        processed += super(HZZAnomalousCouplingsFromHistograms, self).processPhysicsOptions(physOptions)

        if not self.anomalouscouplings: raise ValueError("Have to provide an anomalous coupling as a physOption (fa3, fa2, fL1, fL1Zg)")
        return processed

    @property
    def numberoffais(self):
        return len(self.anomalouscouplings)  #including fa1

    @property
    def sortedcouplings(self):
        return sorted(self.anomalouscouplings, key=["fa1", "fa3", "fa2", "fL1", "fL1Zg"].index)

    def faidefinitionorder(self, i):
        #CMS_zz4l_fai1, CMS_zz4l_fai2, etc. correspond to fa3, fa2, fL1, fL1Zg in that order
        #However they might not be defined in that order, e.g. CMS_zz4l_fai1 might be restricted to (0, 1-CMS_zz4l_fai2)
        return self.anomalouscouplings.index(self.sortedcouplings[i])

    def getPOIList(self):
        if self.useHffanomalous:
            self.modelBuilder.doVar("fCP_Htt[0.,-1,1]")
            self.modelBuilder.out.var("fCP_Htt").setConstant(False)
            self.modelBuilder.out.var("fCP_Htt").setAttribute("flatParam")

        if self.separateggHttH:
            self.modelBuilder.doVar("Rg[1.0,0,400]")
            self.modelBuilder.doVar("Rt[1.0,0,400]")
            self.modelBuilder.out.var("Rg").setConstant(False)
            self.modelBuilder.out.var("Rg").setAttribute("flatParam")
            self.modelBuilder.out.var("Rt").setConstant(False)
            self.modelBuilder.out.var("Rt").setAttribute("flatParam")

            if self.useHffanomalous:
                self.modelBuilder.doVar("fa3_ggH[0.,-1,1]")
                self.modelBuilder.out.var("fa3_ggH").setConstant(False)
                self.modelBuilder.out.var("fa3_ggH").setAttribute("flatParam")
        else:
            self.modelBuilder.doVar("RF[1.0,0,10]")
            if self.useHffanomalous:
                self.modelBuilder.doVar('expr::fa3_ggH("@0 == 0 ? 0 : (@0>0 ? 1 : -1) * 1. / (1. + 4./9. * (1.0 / abs(@0) - 1.))", fCP_Htt)')

        if self.useHffanomalous:
            self.modelBuilder.doVar('expr::kappa("sqrt(1-abs(@0))", fCP_Htt)')
            self.modelBuilder.doVar('expr::kappa_tilde("(@0>0 ? 1 : -1) * sqrt(abs(@0)) * {kappatildettH}", fCP_Htt)'.format(kappatildettH=self.kappa_tilde_ttH))
            self.modelBuilder.doVar('expr::ghg2("sqrt(1-abs(@0))", fa3_ggH)')
            self.modelBuilder.doVar('expr::ghg4("(@0>0 ? 1 : -1) * sqrt(abs(@0))", fa3_ggH)')


        self.modelBuilder.doVar("RV[1.0,0,10]")
        self.modelBuilder.doVar("R[1.0,0,10]")

        pois = super(HZZAnomalousCouplingsFromHistograms, self).getPOIList()

        if not self.modelBuilder.out.var("g1"):
            self.modelBuilder.doVar('expr::g1("sqrt(@0)", CMS_zz4l_fa1)')

        couplings = ["g1"]
        i = 0
        for fai in self.sortedcouplings:
            if fai == "fa1": continue
            ai = self.getcouplingname(fai)
            i += 1

            if self.scalegL1by10000:
                divideby = {
                    "g4": 1,
                    "g2": 1,
                    "g1prime2": 10000,
                    "ghzgs1prime2": 10000,
                }[ai]
            else:
                divideby = 1

            kwargs = {
              "i": i,
              "ai": ai,
              "aidecay": self.aidecay[ai] / divideby,
            }
            self.modelBuilder.doVar('expr::{ai}("(@0>0 ? 1 : -1) * sqrt(abs(@0))*{aidecay}", CMS_zz4l_fai{i})'.format(**kwargs))
            couplings.append(ai)

        if self.scaledifferentsqrtsseparately: raise ValueError("HZZAnomalousCouplingsFromHistograms is not compatible with scaledifferentsqrtsseparately")

        for g in couplings:
            if self.separateggHttH:
                self.modelBuilder.doVar('expr::ggHVV_{g}2("@0*@1*@2*@2", R, Rg, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_{g}2("@0*@1*@2*@2", R, Rt, {g})'.format(g=g))
                if self.useHffanomalous:
                    self.modelBuilder.doVar('expr::ggHVV_ghg22_{g}2("@0*@1*@2*@2*@3*@3", R, Rg, ghg2, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ggHVV_ghg42_{g}2("@0*@1*@2*@2*@3*@3", R, Rg, ghg4, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappa2_{g}2("@0*@1*@2*@2*@3*@3", R, Rt, kappa, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g}2("@0*@1*@2*@2*@3*@3", R, Rt, kappa_tilde, {g})'.format(g=g))
            else:
                self.modelBuilder.doVar('expr::ffHVV_{g}2("@0*@1*@2*@2", R, RF, {g})'.format(g=g))
                if self.useHffanomalous:
                    self.modelBuilder.doVar('expr::ggHVV_ghg22_{g}2("@0*@1*@2*@2*@3*@3", R, RF, ghg2, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ggHVV_ghg42_{g}2("@0*@1*@2*@2*@3*@3", R, RF, ghg4, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappa2_{g}2("@0*@1*@2*@2*@3*@3", R, RF, kappa, {g})'.format(g=g))
                    self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g}2("@0*@1*@2*@2*@3*@3", R, RF, kappa_tilde, {g})'.format(g=g))
            self.modelBuilder.doVar('expr::VVHVV_{g}4("@0*@1*@2*@2*@2*@2", R, RV, {g})'.format(g=g))

        kwargs = {}
        for kwargs["signname"], kwargs["sign"] in ("positive", ""), ("negative", "-"):
            for kwargs["g1"], kwargs["g2"] in itertools.combinations(couplings, 2):
                if self.separateggHttH:
                    self.modelBuilder.doVar('expr::ggHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, Rg, {g1}, {g2})'.format(**kwargs))
                    self.modelBuilder.doVar('expr::ttHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, Rt, {g1}, {g2})'.format(**kwargs))
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rg, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rg, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rt, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, Rt, kappa_tilde, {g1}, {g2})'.format(**kwargs))
                else:
                    self.modelBuilder.doVar('expr::ffHVV_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, RF, {g1}, {g2})'.format(**kwargs))
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RF, kappa_tilde, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}3_{signname}("{sign}@0*@1*@2*@3*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}2_{signname}("{sign}@0*@1*@2*@2*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}3{g2}1_{signname}("{sign}@0*@1*@2*@2*@2*@3", R, RV, {g1}, {g2})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"] in itertools.combinations(couplings, 3):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}2_{signname}("{sign}@0*@1*@2*@3*@4*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}2{g3}1_{signname}("{sign}@0*@1*@2*@3*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}1{g3}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"], kwargs["g4"] in itertools.combinations(couplings, 4):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}1{g4}1_{signname}("{sign}@0*@1*@2*@3*@4*@5", R, RV, {g1}, {g2}, {g3}, {g4})'.format(**kwargs))

        return pois

    signalprocessregex = (
        "(?P<production>gg|tt|bb|qq|Z|W|V)H_"
        "(?:(?P<Hffpure>0(?:PM|M)ff)_)?"
        "(?:"
          "(?P<HVVpure>0(?:PM|M|PH|L1|L1Zg))|"
          "(?P<HVVint>(?:g(?:1|2|4|1prime2|hzgs1prime2)[1234])*)_(?P<HVVintsign>positive|negative)"
        ")$"
    )

    @staticmethod
    def getcouplingname(processorfai, production=None):
        if processorfai == "0PMff": return {"gg": "ghg2", "tt": "kappa", "bb": "kappa"}[production]
        if processorfai == "0Mff": return {"gg": "ghg4", "tt": "kappatilde", "bb": "kappatilde"}[production]

        return {
            "0PM": "g1",
            "0PH": "g2",
            "0M": "g4",
            "0L1": "g1prime2",
            "0L1Zg": "ghzgs1prime2",
            "fa1": "g1",
            "fa2": "g2",
            "fa3": "g4",
            "fL1": "g1prime2",
            "fL1Zg": "ghzgs1prime2",
        }[processorfai]

    def tellAboutProcess(self, bin, process):
        match = re.match(self.signalprocessregex, process)
        if match and match.group("Hffpure"):
            self.useHffanomalous = True

    def getYieldScale(self,bin,process):
        match = re.match(self.signalprocessregex, process)

        if not match:
            if any(process.startswith(_) for _ in ("ggH", "ttH", "bbH", "qqH", "ZH", "WH", "VH")):
                raise ValueError("Your signal process "+process+" doesn't match the pattern")
            return super(HZZAnomalousCouplingsFromHistograms, self).getYieldScale(bin, process)

        if match.group("production")+"H" in self.turnoff: return 0

        if (self.separateggHttH or match.group("Hffpure") is not None) and match.group("production") == "gg": maxpower = 2; production = "ggHVV"
        elif (self.separateggHttH or match.group("Hffpure") is not None) and match.group("production") in ("tt", "bb"): maxpower = 2; production = "ttHVV"
        elif match.group("production") in ("gg", "tt", "bb"): maxpower = 2; production = "ffHVV"
        elif match.group("production") in ("qq", "Z", "W", "V"): maxpower = 4; production = "VVHVV"

        result = production

        if match.group("Hffpure") is not None:
            if match.group("production") not in ("gg", "tt", "bb"): raise ValueError("Don't put fermion couplings for {}H: {}".format(match.group("production"), process))
            Hffpowerdict = {self.getcouplingname(match.group("Hffpure"), match.group("production")): 2}
            result += "_" + "".join("{}{}".format(k, v) for k, v in Hffpowerdict.iteritems())

        if match.group("HVVpure") is not None:
            powerdict = {self.getcouplingname(match.group("HVVpure")): maxpower}
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems())
        elif match.group("HVVint") is not None:
            powerdict = {coupling: int(power) for coupling, power in re.findall("(g(?:1|2|4|1prime2|hzgs1prime2))([1234])", match.group("HVVint"))}

            if sum(powerdict.values()) != maxpower:
                raise ValueError("power dict doesn't add up properly!  Sum should be {}\n{}\n{}".format(maxpower, process, powerdict))

            powerdict = collections.OrderedDict(
                sorted(powerdict.iteritems(), key = lambda x: "g1 g4 g2 g1prime2 ghzgs1prime2".index(x[0]))
            )

            sign = match.group("HVVintsign")
            result += "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems()) + "_" + sign
        else:
            assert False

        if self.verbose:
            print "Process {0} will scale by {1}".format(process,result)

        return result



class HZZAnomalousCouplingsFromHistogramsAim125p38(HZZAnomalousCouplingsFromHistogramsBase, SpinZeroHiggsAiBase, CanTurnOffBkgModel, MultiSignalModel):
    def __init__(self):
        super(HZZAnomalousCouplingsFromHistogramsAim125p38, self).__init__()
        print self

    def getPOIList(self):
        self.modelBuilder.doVar("ghg2[1.0,0,10]")
        self.modelBuilder.out.var("ghg2").removeRange()
        self.modelBuilder.out.var("ghg2").setAttribute("flatParam")
        if self.separateggHttH:
            self.modelBuilder.doVar("kappa[1.0,0,10]")
            self.modelBuilder.out.var("kappa").removeRange()
            self.modelBuilder.out.var("kappa").setAttribute("flatParam")
        else:
            self.modelBuilder.doVar('expr::kappa("@0", ghg2)')

        if self.useHffanomalous:
            self.modelBuilder.doVar("ghg4[0.0,0,10]")
            self.modelBuilder.out.var("ghg4").removeRange()
            self.modelBuilder.out.var("ghg4").setAttribute("flatParam")

            if self.separateggHttH:
                self.modelBuilder.doVar("kappa_tilde[1.0,0,10]")
                self.modelBuilder.out.var("kappa_tilde").removeRange()
                self.modelBuilder.out.var("kappa_tilde").setAttribute("flatParam")
            else:
                self.modelBuilder.doVar('expr::kappa_tilde("0.6482*@0", ghg4)') # 0.6482 comes from sqrt( sigma_(kf =1 )/ sigma_(kftilde =1))        

        self.modelBuilder.doVar("RV[1.0,0,10]")
        self.modelBuilder.doVar("R[1.0,0,10]")


        pois = super(HZZAnomalousCouplingsFromHistogramsAim125p38, self).getPOIList()

        couplings = ["g1"]
        i = 0
        for fai in self.sortedcouplings:
            if fai == "fa1": continue
            ai = self.getcouplingname(fai)
            i += 1
            couplings.append(ai)



        #define EFT constants                                                                                                                        
        self.modelBuilder.doVar('expr::cosW("0.87681811112",)')
        self.modelBuilder.doVar('expr::sinW("0.48082221247",)')
        self.modelBuilder.doVar('expr::mZ("91.2",)')
        self.modelBuilder.doVar('expr::Lambda1("100.0",)')

        self.modelBuilder.doVar('expr::EFT_g1WW("@0",g1)')
        self.modelBuilder.doVar('expr::EFT_g2WW("@0*@0*@1",cosW,g2)')
        self.modelBuilder.doVar('expr::EFT_g4WW("@0*@0*@1",cosW,g4)')

        self.modelBuilder.doVar('expr::EFT_L1WW("(@2 / (@0*@0 - @1*@1) - 2*@1*@1*@3*@4*@4 /(@5*@5*(@0*@0 - @1*@1)))",cosW,sinW,g1prime2,g2,Lambda1,m\
Z)')
        self.modelBuilder.doVar('expr::EFT_L1Zg_L1("2*@0*@1*@2/(@0*@0 - @1*@1)",cosW,sinW,g1prime2)')
        self.modelBuilder.doVar('expr::EFT_L1Zg_g2("-2*@0*@1*@3*@4*@4/((@2*@2)*(@0*@0 - @1*@1))",cosW,sinW,mZ,g2,Lambda1)')
        self.modelBuilder.doVar('expr::EFT_L1Zg("@0 + @1",EFT_L1Zg_L1,EFT_L1Zg_g2)')




        #define expressions                                                                                                                          



        zz_expr      = '"4*(@0*@0/4. + 0.1695*@3*@3 + 0.09076*@1*@1 + 0.03809*@2*@2 + 0.8095*@0*@3/2. + 0.5046*@0*@1/2. + 0.2092*@1*@3 + 0.1023*@4*@\
4 + 0.1901*@0*@4/2. + 0.07429*@3*@4 + 0.04710*@1*@4) ",g1,g2,g4,g1prime2,EFT_L1Zg'
        ww_expr      = '"4*(@0*@0/4. + 0.1320*@3*@3 + 0.1944*@1*@1 + 0.08075*@2*@2 + 0.7204*@0*@3/2. + 0.7437*@0*@1/2. + 0.2774*@3*@1) ",EFT_g1WW,EF\
T_g2WW,EFT_g4WW,EFT_L1WW'

        zgamma_expr  = '"4*(1.118600*@0*@0/4. +0.0035*@1*@1 -  0.125010*@0*@1/2. + 0.000003*@1*@1 - 0.00018*@1*@1 + 0.003100*@0*@1/2. +0.00126*@2*@2\
 + 0.000005*@2*@2 -0.00047*@2*@2)",EFT_g1WW,kappa,kappa_tilde'
        gg_expr      = '"(1.1068*@0*@0 + 0.0082*@0*@0 - 0.1150*@0*@0 + 2.5717*@1*@1 + 0.0091*@1*@1 - 0.1982*@1*@1)",kappa,kappa_tilde'

        bb_expr      = '"(@0*@0 + @1*@1)",kappa,kappa_tilde'
        cc_expr      = '"(@0*@0 + @1*@1)",kappa,kappa_tilde'
        tautau_expr  = '"(@0*@0 + @1*@1)",kappa,kappa_tilde'
        mumu_expr    = '"(@0*@0 + @1*@1)",kappa,kappa_tilde'
        gmgm_expr    = '"4*(1.6054*@0*@0/4. + 0.07312*@1*@1 - 0.6854*@0*@1/2. + 0.00002*@1*@1 - 0.0018*@1*@1 + 0.0085*@0*@1/2. + 0.1699*@2*@2 + 0.00\
002*@2*@2 - 0.0031*@2*@2)",EFT_g1WW,kappa,kappa_tilde'
        
        self.modelBuilder.doVar('expr::R_WW('+str(ww_expr)+')')
        self.modelBuilder.doVar('expr::R_ZZ('+str(zz_expr)+')')
        self.modelBuilder.doVar('expr::R_Zgamma('+str(zgamma_expr)+')')
        self.modelBuilder.doVar('expr::R_gg('+str(gg_expr)+')')
        self.modelBuilder.doVar('expr::R_bb('+str(bb_expr)+')')
        self.modelBuilder.doVar('expr::R_cc('+str(cc_expr)+')')
        self.modelBuilder.doVar('expr::R_tautau('+str(tautau_expr)+')')
        self.modelBuilder.doVar('expr::R_mumu('+str(mumu_expr)+')')
        self.modelBuilder.doVar('expr:R_gammagamma('+str(gmgm_expr)+')')

        if self.fixgamma :
            self.modelBuilder.doVar('expr::gammaH("1",)')
        else :
            #m125.00                                                                                                                                 
            #self.modelBuilder.doVar('expr::gammaH("(0.5824*@0 + 0.2137*@1 + 0.08187*@2 + 0.06272*@3 + 0.02891*@4 + 0.02619*@5 + 0.002270*@6 +  0.001533*@7 + 0.0002176*@8 )/0.9998",R_bb,R_WW,R_gg,R_tautau,R_cc,R_ZZ,R_gammagamma,R_Zgamma,R_mumu)')                                                   
            #m 125.38                                                                                                                                
            self.modelBuilder.doVar('expr::gammaH("(0.5760*@0 + 0.2203*@1 + 0.08154*@2 + 0.06208*@3 + 0.02860*@4 + 0.02716*@5 + 0.002270*@6 +  0.001567*@7 + 0.0002153*@8 )/0.99973230",R_bb,R_WW,R_gg,R_tautau,R_cc,R_ZZ,R_gammagamma,R_Zgamma,R_mumu)')
            

            


        self.modelBuilder.doVar('expr::alpha(0.39,)')


        for g in couplings:
            if self.useHffanomalous:
                self.modelBuilder.doVar('expr::ggHVV_{g}2("(@1*@1+@2*@2) * @3*@3 / @0", gammaH, ghg2, ghg4, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_{g}2("(@1*@1+@4*@2*@2) * @3*@3 / @0", gammaH, kappa, kappa_tilde, {g}, alpha)'.format(g=g))
                self.modelBuilder.doVar('expr::ggHVV_ghg22_{g}2("@1*@1*@2*@2 / @0", gammaH, ghg2, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ggHVV_ghg42_{g}2("@1*@1*@2*@2 / @0", gammaH, ghg4, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_kappa2_{g}2("@1*@1*@2*@2 / @0", gammaH, kappa, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g}2("@1*@1*@2*@2 / @0", gammaH, kappa_tilde, {g})'.format(g=g))
            else:
                self.modelBuilder.doVar('expr::ggHVV_{g}2("@1*@1*@2*@2 / @0", gammaH, ghg2, {g})'.format(g=g))
                self.modelBuilder.doVar('expr::ttHVV_{g}2("@1*@1*@2*@2 / @0", gammaH, kappa, {g})'.format(g=g))
            self.modelBuilder.doVar('expr::VVHVV_{g}4("@1*@1*@1*@1 / @0", gammaH, {g})'.format(g=g))

        kwargs = {}
        for kwargs["signname"], kwargs["sign"] in ("positive", ""), ("negative", "-"):


            for kwargs["g1"], kwargs["g2"] in itertools.combinations(couplings, 2):



                if self.separateggHttH:
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_{g1}1{g2}1_{signname}("{sign}(@1*@1+@2*@2) * @3*@4 / @0", gammaH, ghg2, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_{g1}1{g2}1_{signname}("{sign}(@1*@1+@5*@2*@2) * @3*@4 / @0", gammaH, kappa, kappa_tilde, {g1}, {g2}, alpha)'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, kappa_tilde, {g1},{g2})'.format(**kwargs))
                    else:
                        self.modelBuilder.doVar('expr::ggHVV_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, kappa, {g1}, {g2})'.format(**kwargs))
                else:
                    self.modelBuilder.doVar('expr::ggHVV_{g1}1{g2}1_{signname}("{sign}@1*@2 / @0", gammaH, {g1}, {g2})'.format(**kwargs))
                    self.modelBuilder.doVar('expr::ttHVV_{g1}1{g2}1_{signname}("{sign}@1*@2 / @0", gammaH, {g1}, {g2})'.format(**kwargs))
                    if self.useHffanomalous:
                        self.modelBuilder.doVar('expr::ggHVV_ghg22_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, ghg2, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ggHVV_ghg42_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, ghg4, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappa2_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, kappa, {g1}, {g2})'.format(**kwargs))
                        self.modelBuilder.doVar('expr::ttHVV_kappatilde2_{g1}1{g2}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, kappa_tilde, {g1},{g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}3_{signname}("{sign}@1*@2*@2*@2 / @0", gammaH, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}2_{signname}("{sign}@1*@1*@2*@2 / @0", gammaH, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}3{g2}1_{signname}("{sign}@1*@1*@1*@2 / @0", gammaH, {g1}, {g2})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"] in itertools.combinations(couplings, 3):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}2_{signname}("{sign}@1*@2*@3*@3 / @0", gammaH, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}2{g3}1_{signname}("{sign}@1*@2*@2*@3 / @0", gammaH, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVHVV_{g1}2{g2}1{g3}1_{signname}("{sign}@1*@1*@2*@3 / @0", gammaH, {g1}, {g2}, {g3})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"], kwargs["g4"] in itertools.combinations(couplings, 4):
                self.modelBuilder.doVar('expr::VVHVV_{g1}1{g2}1{g3}1{g4}1_{signname}("{sign}@1*@2*@3*@4 / @0", gammaH, {g1}, {g2}, {g3}, {g4})'.format(**kwargs))

        return pois







spinZeroHiggsAC = SpinZeroHiggsAC()
multiSignalSpinZeroHiggs = MultiSignalSpinZeroHiggsAC()
hzzAnomalousCouplingsFromHistograms = HZZAnomalousCouplingsFromHistograms() 
hzzAnomalousCouplingsFromHistogramsNonSMEFT = HZZAnomalousCouplingsFromHistograms() 
hzzAnomalousCouplingsFromHistogramsAim125p38 = HZZAnomalousCouplingsFromHistogramsAim125p38()
