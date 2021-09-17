from sys import stdout, stderr, exit
import os.path
import ROOT
from collections import defaultdict
from math import *

RooArgSet_add_original = ROOT.RooArgSet.add
def RooArgSet_add_patched(self, obj, *args, **kwargs):
    if isinstance(obj, ROOT.RooAbsCollection):
        return ROOT.RooAbsCollection.add(self, obj, *args, **kwargs)
    else:
        return RooArgSet_add_original(self, obj, *args, **kwargs)
ROOT.RooArgSet.add = RooArgSet_add_patched

from HiggsAnalysis.CombinedLimit.ModelTools import ModelBuilder

class FileCache:
    def __init__(self, basedir, maxsize=250):
        self._basedir = basedir
        self._maxsize = maxsize
        self._files = {}
        self._hits  = defaultdict(int)
        self._total = 0
    def __getitem__(self, fname):
        self._total += 1
        if fname not in self._files:
            if len(self._files) >= self._maxsize:
                print "Flushing file cache of size %d" % len(self._files)
                keys = self._files.keys()
                keys.sort(key = lambda k : self._files[k][1] + 10*self._hits[k])
                for k in keys[:self._maxsize/2]:
                    self._files[k][0].Close()
                    del self._files[k]
            trueFName = fname 
            if not os.path.exists(trueFName) and not os.path.isabs(trueFName) and os.path.exists(self._basedir+"/"+trueFName):
                trueFName = self._basedir+"/"+trueFName;
            self._files[fname] = [ ROOT.TFile.Open(trueFName), self._total ]
        else:
            self._files[fname][1] = self._total
        self._hits[fname] += 1
        return self._files[fname][0]

class ShapeBuilder(ModelBuilder):
    def __init__(self,datacard,options):
        ModelBuilder.__init__(self,datacard,options) 
        if not datacard.hasShapes: 
            raise RuntimeError, "You're using a ShapeBuilder for a model that has no shapes"
        if options.libs:
            for lib in options.libs:
                ROOT.gSystem.Load(lib)
    	self.wspnames = {}
    	self.wsp = None
    	self.extraImports = []
	self.norm_rename_map = {}
        self._fileCache = FileCache(self.options.baseDir)
    ## ------------------------------------------
    ## -------- ModelBuilder interface ----------
    ## ------------------------------------------
    def doObservables(self):
        if (self.options.verbose > 2): stderr.write("Using shapes: qui si parra' la tua nobilitate\n")
        self.prepareAllShapes();
        if len(self.DC.bins) > 1 or not self.options.forceNonSimPdf:
            ## start with just a few channels
            strexpr="CMS_channel[" + ",".join(["%s=%d" % (l,i) for i,l in enumerate(self.DC.bins[:5])]) + "]";
            self.doVar(strexpr);
            self.out.binCat = self.out.cat("CMS_channel");
            ## then add all the others, to avoid a too long factory string
            for i,l in enumerate(self.DC.bins[5:]): self.out.binCat.defineType(l,i+5)   
            if self.options.verbose > 1: stderr.write("Will use category 'CMS_channel' to identify the %d channels\n" % self.out.binCat.numTypes())
            self.out.obs = ROOT.RooArgSet()
            self.out.obs.add(self.out.binVars)
            self.out.obs.add(self.out.binCat)
        else:
            self.out.obs = self.out.binVars
        self.doSet("observables",self.out.obs)
        if len(self.DC.obs) != 0 and not self.options.noData:
            self.doCombinedDataset()
    def doIndividualModels(self):
        if self.options.verbose:
            stderr.write("Creating pdfs for individual modes (%d): " % len(self.DC.bins));
            stderr.flush()
        bbb_names = []
        for i,b in enumerate(self.DC.bins):
            #print "  + Getting model for bin %s" % (b)
            pdfs   = ROOT.RooArgList(); bgpdfs   = ROOT.RooArgList()
            coeffs = ROOT.RooArgList(); bgcoeffs = ROOT.RooArgList()
            sigcoeffs = []
            binconstraints = ROOT.RooArgList()
            bbb_args = None
            channelBinParFlag = b in self.DC.binParFlags.keys()
            if channelBinParFlag:
                print 'Channel %s will use autoMCStats with settings: event-threshold=%g, include-signal=%i, hist-mode=%i' % ((b,)+self.DC.binParFlags[b])
            for p in self.DC.exp[b].keys(): # so that we get only self.DC.processes contributing to this bin
                if self.DC.exp[b][p] == 0: continue
                if self.physics.getYieldScale(b,p) == 0: continue # exclude really the pdf
                #print "  +--- Getting pdf for %s in bin %s" % (p,b)
                (pdf,coeff) = (self.getPdf(b,p), self.out.function("n_exp_bin%s_proc_%s" % (b,p)))
                if self.options.optimizeExistingTemplates:
                    pdf1 = self.optimizeExistingTemplates(pdf)
                    if (pdf1 != pdf):
                        self.out.dont_delete.append(pdf1)
                        pdf = pdf1
                extranorm = self.getExtraNorm(b,p)
                if extranorm:
                    if self.options.packAsymPows:
                        if coeff.ClassName() == "ProcessNormalization":
                            pass # nothing to do
                        elif coeff.ClassName() == "RooRealVar":
                            coeff = self.addObj(ROOT.ProcessNormalization, "n_exp_final_bin%s_proc_%s" % (b,p), "", coeff.getVal())
                        else:
                            raise RuntimeError("packAsymPows: can't work with a coefficient of kind %s for %s %s" % (coeff.ClassName(), b, p))
                        for X in extranorm:
                            if type(X) == tuple:
                                (klo, khi, syst) = X
                                coeff.addAsymmLogNormal(klo,khi, self.out.var(syst))
                            else:
                                if self.out.function(X):
                                    coeff.addOtherFactor(self.out.function(X))
                                else:
                                    coeff.addOtherFactor(self.getObj(X))
                    else:   
                        prodset = ROOT.RooArgList(self.out.function("n_exp_bin%s_proc_%s" % (b,p)))
                        for X in extranorm:
                            # X might already be in the workspace (e.g. _norm term)...
                            if self.out.function(X):
                                    prodset.add(self.out.function(X))
                            # ... but usually it's only in our object store (e.g. AsymPow for shape systs)
                            else:
                    		prodset.add(self.getObj(X))
                        coeff = self.addObj(ROOT.RooProduct, "n_exp_final_bin%s_proc_%s" % (b,p), "", prodset)
                pdf.setStringAttribute("combine.process", p)
                pdf.setStringAttribute("combine.channel", b)
                pdf.setAttribute("combine.signal", self.DC.isSignal[p])
                if channelBinParFlag and self.DC.isSignal[p] and not self.DC.binParFlags[b][1]:
                    pdf.setAttribute('skipForErrorSum')
                coeff.setStringAttribute("combine.process", p)
                coeff.setStringAttribute("combine.channel", b)
                coeff.setAttribute("combine.signal", self.DC.isSignal[p])
                pdfs.add(pdf); coeffs.add(coeff)
                if not self.DC.isSignal[p]:
                    bgpdfs.add(pdf); bgcoeffs.add(coeff)
                else:
                    sigcoeffs.append(coeff)
            if self.options.verbose > 1: print "Creating RooAddPdf %s with %s elements" % ("pdf_bin"+b, coeffs.getSize())
            if channelBinParFlag:
                if self.options.useCMSHistSum:
                    prop = self.addObj(ROOT.CMSHistSum, "prop_bin%s" % b, "", pdfs.at(0).getXVar(), pdfs, coeffs)
                    prop.setAttribute('CachingPdf_NoClone', True)
                else:
                    prop = self.addObj(ROOT.CMSHistErrorPropagator, "prop_bin%s" % b, "", pdfs.at(0).getXVar(), pdfs, coeffs)
                prop.setAttribute('CachingPdf_Direct', True)
                if self.DC.binParFlags[b][0] >= 0.:
                    bbb_args = prop.setupBinPars(self.DC.binParFlags[b][0])
                    for bidx in range(bbb_args.getSize()):
                        arg = bbb_args.at(bidx)
                        n = arg.GetName()
                        bbb_names.append(n)
                        parname = n
                        self.out._import(arg)
                        if arg.getAttribute("createGaussianConstraint"):
                            if self.options.noOptimizePdf:
                                self.doObj("%s_Pdf" % n, "Gaussian", "%s, %s_In[0,%s], %s" % (n, n, '-7,7', '1.0'), True)
                            else:
                                self.doObj("%s_Pdf" % n, "SimpleGaussianConstraint", "%s, %s_In[0,%s], %s" % (n, n, '-7,7', '1.0'), True)
                            self.out.var(n).setVal(0)
                            self.out.var(n).setError(1)
                            if self.options.optimizeBoundNuisances: self.out.var(n).setAttribute("optimizeBounds")
                        elif arg.getAttribute("createPoissonConstraint"):
                            nom = arg.getVal()
                            pval = ROOT.Math.normal_cdf_c(7)
                            minObs = nom
                            while minObs > 0 and (ROOT.TMath.Poisson(minObs, nom + 1) > pval):
                                minObs -= (sqrt(nom) if nom > 10 else 1)
                            maxObs = nom + 2
                            while (ROOT.TMath.Poisson(maxObs, nom + 1) > pval):
                                #print "Poisson(maxObs = %d, %f) = %g > 1e-12" % (maxObs, args[0]+1, ROOT.TMath.Poisson(maxObs, args[0]+1))
                                maxObs += (sqrt(nom) if nom > 10 else 2)
                            self.doObj("%s_Pdf" % n, "Poisson", "%s_In[%d,%f,%f], %s, 1" % (n, nom, minObs, maxObs, n))
                            if n.endswith('_prod'):
                                parname = n[:-5]
                        binconstraints.add(self.out.pdf('%s_Pdf' % n))
                        self.out.var("%s_In" % n).setConstant(True)
                        self.extraNuisances.append(self.out.var("%s" % parname))
                        self.extraGlobalObservables.append(self.out.var("%s_In" % n))
                if not self.out.var('ONE'):
                    self.doVar('ONE[1.0]')
                sum_s = self.addObj(ROOT.RooRealSumPdf, "pdf_bin%s"       % b,  "", ROOT.RooArgList(prop),   ROOT.RooArgList(self.out.var('ONE')), True)
                if not self.options.noBOnly:
                    if not self.out.var('ZERO'):
                        self.doVar('ZERO[0.0]')
                    customizer = ROOT.RooCustomizer(prop, "")
                    for arg in sigcoeffs:
                        customizer.replaceArg(arg, self.out.var('ZERO'))
                    prop_b = customizer.build(True)
                    if len(sigcoeffs):
                        prop_b.SetName("prop_bin%s_bonly" % b)
                    self.objstore[prop_b.GetName()] = prop_b
                    sum_b = self.addObj(ROOT.RooRealSumPdf, "pdf_bin%s_bonly"       % b,  "", ROOT.RooArgList(prop_b),   ROOT.RooArgList(self.out.var('ONE')), True)
            else:
                sum_s = self.addObj(ROOT.RooAddPdf,"pdf_bin%s"       % b,  "",  pdfs,   coeffs)
                if not self.options.noBOnly: sum_b = self.addObj(ROOT.RooAddPdf, "pdf_bin%s_bonly" % b, "", bgpdfs, bgcoeffs)
            sum_s.setAttribute("MAIN_MEASUREMENT") # useful for plain ROOFIT optimization on ATLAS side
            if b in self.pdfModes: 
                sum_s.setAttribute('forceGen'+self.pdfModes[b].title())
                if not self.options.noBOnly: sum_b.setAttribute('forceGen'+self.pdfModes[b].title())
            addSyst = False
            if    self.options.moreOptimizeSimPdf == "none":   addSyst = True
            # New behaviour for "lhchcg" mode - since autoMCStats constraints are only added to their respective channels,
            # we can't get away with only adding a constraint production to the first channel. Instead we will always enter
            # the code block below (addSyst=True), and only add the normal nuisPdfs on i==0, while the binconstraints will
            # always be added.
            elif  self.options.moreOptimizeSimPdf == "lhchcg": addSyst = True
            elif  self.options.moreOptimizeSimPdf == "cms":
                if self.options.noOptimizePdf: raise RuntimeError, "--optimize-simpdf-constraints=cms is incompatible with --no-optimize-pdfs"
                addSyst = False
            if (len(self.DC.systs) or binconstraints.getSize()) and addSyst:
                ## rename the pdfs
                self.renameObj("pdf_bin%s" % b, "pdf_bin%s_nuis" % b)
                if not self.options.noBOnly:
                	self.renameObj("pdf_bin%s_bonly" % b, "pdf_bin%s_bonly_nuis" % b)
                # now we multiply by all the nuisances, but avoiding nested products
                # so we first make a list of all nuisances plus the RooAddPdf
                if len(self.DC.systs) and not (self.options.moreOptimizeSimPdf == "lhchcg" and i > 0):
                    sumPlusNuis_s = ROOT.RooArgList(self.out.nuisPdfs)
                else:
                    sumPlusNuis_s = ROOT.RooArgList()
                sumPlusNuis_s.add(sum_s)
                pdf_bins = self.addObj(ROOT.RooProdPdf, 'pdfbins_bin%s' % b, "", binconstraints)
                sumPlusNuis_s.add(pdf_bins)
                # then make RooProdPdf and import it
                pdf_s = self.addObj(ROOT.RooProdPdf, "pdf_bin%s"       % b, "", sumPlusNuis_s)
                if not self.options.noBOnly:
                    if len(self.DC.systs):
                        sumPlusNuis_b = ROOT.RooArgList(self.out.nuisPdfs)
                    else:
                        sumPlusNuis_b = ROOT.RooArgList()
                    sumPlusNuis_b.add(sum_b)
                    sumPlusNuis_b.add(pdf_bins)
                    pdf_b = self.addObj(ROOT.RooProdPdf, "pdf_bin%s_bonly" % b, "", sumPlusNuis_b)
                if b in self.pdfModes: 
                    pdf_s.setAttribute('forceGen'+self.pdfModes[b].title())
                    if not self.options.noBOnly: pdf_b.setAttribute('forceGen'+self.pdfModes[b].title())
                if self.options.verbose:
                    if i > 0: stderr.write("\b\b\b\b\b");
                    stderr.write(". %4d" % (i+1))
                    stderr.flush()
            else:
                if self.options.verbose:
                    if i > 0: stderr.write("\b\b\b\b\b");
                    stderr.write(". %4d" % (i+1))
                    stderr.flush()
            if channelBinParFlag and not self.options.noHistFuncWrappers:
                for idx in xrange(pdfs.getSize()):
                    wrapper = ROOT.CMSHistFuncWrapper(pdfs[idx].GetName() + '_wrapper', '', pdfs.at(idx).getXVar(), pdfs.at(idx), prop, idx)
                    wrapper.setStringAttribute("combine.process", pdfs.at(idx).getStringAttribute("combine.process"))
                    wrapper.setStringAttribute("combine.channel", pdfs.at(idx).getStringAttribute("combine.channel"))
                    self.extraImports.append(wrapper)

        if len(bbb_names) > 0 :
            bbb_nuisanceargset = ROOT.RooArgSet()
            for nuisanceName in bbb_names:
               bbb_nuisanceargset.add(self.out.var(nuisanceName))
               # FIXME - should restructure to have autoMCStats pdfs added
               # to nuisPdfs *before* creating the final channel pdfs
               if self.options.moreOptimizeSimPdf == "cms":
                   self.out.nuisPdfs.add(self.out.pdf(nuisanceName + '_Pdf'))
            self.out.defineSet("group_autoMCStats",bbb_nuisanceargset)
        if self.options.verbose:
            stderr.write("\b\b\b\bdone.\n"); stderr.flush()
    def doCombination(self):
        ## Contrary to Number-counting models, here each channel PDF already contains the nuisances
        ## So we just have to build the combined pdf
        dupObjs = set()
        dupNames = set()
        if len(self.DC.bins) > 1 or not self.options.forceNonSimPdf:
            if self.options.doMasks:
                maskList = ROOT.RooArgList()
                for b in self.DC.bins:
                    maskList.add(self.out.arg(self.physics.getChannelMask(b)))
            for (postfixIn,postfixOut) in [ ("","_s"), ("_bonly","_b") ]:
                simPdf = ROOT.RooSimultaneous("model"+postfixOut, "model"+postfixOut, self.out.binCat) if self.options.noOptimizePdf else ROOT.RooSimultaneousOpt("model"+postfixOut, "model"+postfixOut, self.out.binCat)
                for b in self.DC.bins:
                    pdfi = self.getObj("pdf_bin%s%s" % (b,postfixIn))
                    self.RenameDupObjs(dupObjs, dupNames, pdfi, b)
                    simPdf.addPdf(pdfi, b)
                if (not self.options.noOptimizePdf) and self.options.doMasks:
                    simPdf.addChannelMasks(maskList)
                if len(self.DC.systs) and (not self.options.noOptimizePdf) and self.options.moreOptimizeSimPdf == "cms":
                    simPdf.addExtraConstraints(self.out.nuisPdfs)
                if self.options.verbose:
                    stderr.write("Importing combined pdf %s\n" % simPdf.GetName()); stderr.flush()

		# take care of any variables which were renamed (eg for "param")
		paramString,renameParamString,toFreeze = self.getRenamingParameters()
		if len(renameParamString): 
                  self.out._import(simPdf, ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.RenameVariable(paramString,renameParamString))
                else: self.out._import(simPdf, ROOT.RooFit.RecycleConflictNodes())
		for pfreeze in toFreeze:
		  if self.out.var(pfreeze) : self.out.var(pfreeze).setConstant(True)
                if self.options.noBOnly: break
        else:
            self.out._import(self.getObj("pdf_bin%s"       % self.DC.bins[0]).clone("model_s"), ROOT.RooFit.Silence())
            if not self.options.noBOnly: 
                self.out._import(self.getObj("pdf_bin%s_bonly" % self.DC.bins[0]).clone("model_b"), ROOT.RooFit.Silence())
        for arg in self.extraImports:
            #print 'Importing extra arg: %s' % arg.GetName()
            self.out._import(arg, ROOT.RooFit.RecycleConflictNodes())
        if self.options.fixpars:
            pars = self.out.pdf("model_s").getParameters(self.out.obs)
            iter = pars.createIterator()
            while True:
                arg = iter.Next()
                if arg == None: break;
                if arg.InheritsFrom("RooRealVar") and arg.GetName() != "r": 
                    arg.setConstant(True);

    def RenameDupObjs(self, dupObjs, dupNames, newObj, postFix):
        #print 'Checking for duplicates in %s' % newObj.GetName()
        branchNodes = ROOT.RooArgList()
        newObj.branchNodeServerList(branchNodes)
        # branchNodes.Print('v')
        for i in xrange(1, branchNodes.getSize()):
            arg = branchNodes.at(i)
            if arg.GetName() in dupNames and arg not in dupObjs:
                if self.options.verbose > 1 : stderr.write('Object %s is duplicated, will rename to %s_%s\n'%(arg.GetName(),arg.GetName(),postFix))
                arg.SetName(arg.GetName() + '_%s' % postFix)
            # if arg.GetName() in dupNames and arg in dupObjs:
                # print 'Objected %s is repeated' % arg.GetName()
            dupObjs.add(arg)
            dupNames.add(arg.GetName())
    ## --------------------------------------
    ## -------- High level helpers ----------
    ## --------------------------------------
    def prepareAllShapes(self):
        shapeTypes = []; shapeBins = {}; shapeObs = {}
        self.pdfModes = {}
        for ib,b in enumerate(self.DC.bins):
            databins = {}; bgbins = {}
            channelBinParFlag = b in self.DC.binParFlags.keys()
            for p in [self.options.dataname]+self.DC.exp[b].keys():
                if len(self.DC.obs) == 0 and p == self.options.dataname: continue
                if p != self.options.dataname and self.DC.exp[b][p] == 0: continue
                shape = self.getShape(b,p); norm = 0;
                if shape == None: # counting experiment
                    if not self.out.var("CMS_fakeObs"): 
                        self.doVar("CMS_fakeObs[0,1]");
                        self.out.var("CMS_fakeObs").setBins(1);
                        self.doSet("CMS_fakeObsSet","CMS_fakeObs");
			self.doVar("CMS_fakeWeight[0,1]");self.out.var("CMS_fakeWeight").removeRange()
                        shapeObs["CMS_fakeObsSet"] = self.out.set("CMS_fakeObsSet")
                    if p == self.options.dataname:
                        self.pdfModes[b] = 'binned'
                        shapeTypes.append("RooDataHist")
                    else:
                        shapeTypes.append("RooAbsPdf");
                elif shape.ClassName().startswith("TH1"):
                    shapeTypes.append("TH1"); shapeBins[b] = shape.GetNbinsX()
                    if channelBinParFlag:
                        self.selfNormBins.append(b)
                    norm = shape.Integral()
                    if p == self.options.dataname: 
                        if self.options.poisson > 0 and norm > self.options.poisson:
                            self.pdfModes[b] = 'poisson'
                        else:
                            self.pdfModes[b] = 'binned'
                        for i in xrange(1, shape.GetNbinsX()+1):
                            if shape.GetBinContent(i) > 0: databins[i] = True
                    elif not self.DC.isSignal[p]:
                        for i in xrange(1, shape.GetNbinsX()+1):
                            if shape.GetBinContent(i) > 0: bgbins[i] = True
                elif shape.InheritsFrom("RooDataHist"):
                    shapeTypes.append("RooDataHist"); 
                    #if doPadding: shapeBins[b] = shape.numEntries() --> Not clear this is needed at all for RooDataHists so just ignore
                    shapeObs[self.argSetToString(shape.get())] = shape.get()
                    norm = shape.sumEntries()
                    if p == self.options.dataname: 
                        if self.options.poisson > 0 and norm > self.options.poisson:
                            self.pdfModes[b] = 'poisson'
                        else:
                            self.pdfModes[b] = 'binned'
                elif shape.InheritsFrom("RooDataSet"):
                    shapeTypes.append("RooDataSet"); 
                    shapeObs[self.argSetToString(shape.get())] = shape.get()
                    norm = shape.sumEntries()
                    if p == self.options.dataname: self.pdfModes[b] = 'unbinned'
                elif shape.InheritsFrom("TTree"):
                    shapeTypes.append("TTree"); 
                    if p == self.options.dataname: self.pdfModes[b] = 'unbinned'
                elif shape.InheritsFrom("RooAbsPdf"):
                    shapeTypes.append("RooAbsPdf");
                elif shape.InheritsFrom("CMSHistFunc"):
                    shapeTypes.append("CMSHistFunc");
                else: raise RuntimeError, "Currently supporting only TH1s, RooDataHist and RooAbsPdfs"
                if norm != 0:
                    if p == self.options.dataname:
                        if len(self.DC.obs):
                            if self.DC.obs[b] == -1: self.DC.obs[b] = norm
                            elif self.DC.obs[b] == 0 and norm > 0.01:
                                if not self.options.noCheckNorm: raise RuntimeError, "Mismatch in normalizations for observed data in bin %s: text %f, shape %f" % (b,self.DC.obs[b],norm)
                            elif self.DC.obs[b] >0 and abs(norm/self.DC.obs[b]-1) > 0.005:
                                if not self.options.noCheckNorm: raise RuntimeError, "Mismatch in normalizations for observed data in bin %s: text %f, shape %f" % (b,self.DC.obs[b],norm)
                    else:
                        if self.DC.exp[b][p] == -1: self.DC.exp[b][p] = norm
                        elif self.DC.exp[b][p] > 0 and abs(norm-self.DC.exp[b][p]) > 0.01*max(1,self.DC.exp[b][p]): 
                            if not self.options.noCheckNorm: raise RuntimeError, "Mismatch in normalizations for bin %s, process %s: rate %f, shape %f" % (b,p,self.DC.exp[b][p],norm)
            if len(databins) > 0:
                for i in databins.iterkeys():
                    if i not in bgbins: stderr.write("Channel %s has bin %d fill in data but empty in all backgrounds\n" % (b,i))
        if shapeTypes.count("TH1"):
	    self.TH1Observables = {}
	    self.out.binVars = ROOT.RooArgSet()
            self.out.maxbins = max([shapeBins[k] for k in shapeBins.keys()])
	    if self.options.optimizeTemplateBins:
              if self.options.verbose > 1: stderr.write("Will use binning variable CMS_th1x with %d bins\n" % self.out.maxbins)
	      self.doVar("CMS_th1x[0,%d]" % self.out.maxbins); self.out.var("CMS_th1x").setBins(self.out.maxbins)
              self.out.binVars.add(self.out.var("CMS_th1x"))
              shapeObs['CMS_th1x'] = self.out.var("CMS_th1x")
	      for b in shapeBins: self.TH1Observables[b] = "CMS_th1x"
	    else:
	      for b in shapeBins:
		binVar = "CMS_th1x_%s"%b 
                if self.options.verbose > 1: stderr.write("Will use binning variable %s with %d bins\n" %(binVar,shapeBins[b]))
                self.doVar("%s[0,%d]" %(binVar,shapeBins[b])); self.out.var(binVar).setBins(shapeBins[b])
                self.out.binVars.add(self.out.var(binVar))
                shapeObs['CMS_th1x_%s'%b] = self.out.var(binVar)
	        self.TH1Observables[b] = binVar
        if shapeTypes.count("TH1") == len(shapeTypes):
            self.out.mode    = "binned"
        elif shapeTypes.count("RooDataSet") > 0 or shapeTypes.count("TTree") > 0 or len(shapeObs.keys()) > 1: # remake RooArgSet for binVars with all Variables inside 
            self.out.mode = "unbinned"
            if self.options.verbose > 1: stderr.write("Will work with unbinned datasets\n")
            if self.options.verbose > 1: stderr.write("Observables: %s\n" % str(shapeObs.keys()))
            if len(shapeObs.keys()) != 1:
                self.out.binVars = ROOT.RooArgSet()
                for obs in shapeObs.values():
                     self.out.binVars.add(obs, True)
            else:
                self.out.binVars = shapeObs.values()[0]
            self.out._import(self.out.binVars)
        else:
            self.out.mode = "binned"
            if self.options.verbose > 1: stderr.write("Will make a binned dataset\n")
            if self.options.verbose > 1: stderr.write("Observables: %s\n" % str(shapeObs.keys()))
            if len(shapeObs.keys()) != 1:
                raise RuntimeError, "There's more than once choice of observables: %s\n" % str(shapeObs.keys())
            self.out.binVars = shapeObs.values()[0]
            self.out._import(self.out.binVars)
    def doCombinedDataset(self):
        if len(self.DC.bins) == 1 and self.options.forceNonSimPdf:
            data = self.getData(self.DC.bins[0],self.options.dataname).Clone(self.options.dataname)
            self.out._import(data)
            return

        """ Combine is able to handle the binned/vs unbinned properly so no need for separate commands
	, commenting this switch helps with avoiding creating an n-dim dataset (i.e padding rows for the actual data)
	, not clear why the "binned" version was ever needed, but should be checked. 

	if self.out.mode == "binned":
            combiner = ROOT.CombDataSetFactory(self.out.obs, self.out.binCat)
            for b in self.DC.bins: 
	    	combiner.addSetBin(b, self.getData(b,self.options.dataname))
            self.out.data_obs = combiner.done(self.options.dataname,self.options.dataname)
            self.out._import(self.out.data_obs)
        elif self.out.mode == "unbinned":
            combiner = ROOT.CombDataSetFactory(self.out.obs, self.out.binCat)
            for b in self.DC.bins: combiner.addSetAny(b, self.getData(b,self.options.dataname))
            self.out.data_obs = combiner.doneUnbinned(self.options.dataname,self.options.dataname)
            self.out._import(self.out.data_obs)
        else: raise RuntimeException, "Only combined datasets are supported"
	"""
        combiner = ROOT.CombDataSetFactory(self.out.obs, self.out.binCat)
        for b in self.DC.bins: 
		combiner.addSetAny(b, self.getData(b,self.options.dataname))
        self.out.data_obs = combiner.doneUnbinned(self.options.dataname,self.options.dataname)
        self.out._import(self.out.data_obs)
	if self.options.verbose>2:
          print "Created combined dataset with ",self.out.data_obs.numEntries()," entries, out of:"
          for b in self.DC.bins: print "  bin", b, ": entries = ", self.getData(b,self.options.dataname).numEntries()
    ## -------------------------------------
    ## -------- Low level helpers ----------
    ## -------------------------------------
    def getShape(self,channel,process,syst="",_cache={},allowNoSyst=False):
        if _cache.has_key((channel,process,syst)): 
            if self.options.verbose > 2: print "recyling (%s,%s,%s) -> %s\n" % (channel,process,syst,_cache[(channel,process,syst)].GetName())
            return _cache[(channel,process,syst)];
        postFix="Sig" if (process in self.DC.isSignal and self.DC.isSignal[process]) else "Bkg"
        bentry = None
        if self.DC.shapeMap.has_key(channel): bentry = self.DC.shapeMap[channel]
        elif self.DC.shapeMap.has_key("*"):   bentry = self.DC.shapeMap["*"]
        else: raise KeyError, "Shape map has no entry for channel '%s'" % (channel)
        names = []
        if bentry.has_key(process): names = bentry[process]
        elif bentry.has_key("*"):   names = bentry["*"]
        elif self.DC.shapeMap["*"].has_key(process): names = self.DC.shapeMap["*"][process]
        elif self.DC.shapeMap["*"].has_key("*"):     names = self.DC.shapeMap["*"]["*"]
        else: raise KeyError, "Shape map has no entry for process '%s', channel '%s'" % (process,channel)
        if len(names) == 1 and names[0] == "FAKE": return None
        if syst != "": 
            if len(names) == 2:
                if allowNoSyst: return None
                raise RuntimeError, "Can't find systematic "+syst+" for process '%s', channel '%s'" % (process,channel)
            names = [names[0], names[2]]
        else:   
            names = [names[0], names[1]]
        strmass = "%d" % self.options.mass if self.options.mass % 1 == 0 else str(self.options.mass)
        finalNames = [ x.replace("$PROCESS",process).replace("$CHANNEL",channel).replace("$SYSTEMATIC",syst).replace("$MASS",strmass) for x in names ]
	for mp in self.options.modelparams:
	   if len(mp.split('='))!=2 : raise RuntimeError, "No value found for keyword in %s (use --keyword-value WORD=VALUE)"%mp 
	   mpname, mpv = mp.split('=')
	   protected_kwords =  ["PROCESS","CHANNEL","SYSTEMATIC","MASS"]
	   if mpname in protected_kwords: raise RuntimeError, "Cannot use the following keywords (already assigned in combine): $"+" $".join(protected_kwords) 
           finalNames = [ fn.replace("$%s"%mpname,mpv) for fn in finalNames ]
        file = self._fileCache[finalNames[0]]; objname = finalNames[1]
        if not file: raise RuntimeError, "Cannot open file %s (from pattern %s)" % (finalNames[0],names[0])
        if ":" in objname: # workspace:obj or ttree:xvar or th1::xvar
            (wname, oname) = objname.split(":")
            if (file,wname) not in self.wspnames : 
		self.wspnames[(file,wname)] = file.Get(wname)
	    self.wsp = self.wspnames[(file,wname)]
            if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
            if self.wsp.ClassName() == "RooWorkspace":
                ret = self.wsp.data(oname)
                if not ret: ret = self.wsp.pdf(oname)
                if not ret: ret = self.wsp.function(oname)
                if not ret:
                    if allowNoSyst: return None
                    raise RuntimeError, "Object %s in workspace %s in file %s does not exist or it's neither a data nor a pdf" % (oname, wname, finalNames[0])
                # Fix the fact that more than one entry can refer to the same object
                ret = ret.Clone("shape%s_%s_%s%s" % (postFix,process,channel, "_"+syst if syst else ""))
                if self.options.removeMultiPdf and ret.InheritsFrom("RooMultiPdf"):
                    print ("removeMultiPdf",process,channel,oname,"current index",ret.getCurrentIndex())
                    ret=ret.getCurrentPdf().Clone(ret.GetName())
                if self.options.optimizeMHDependency and ret.InheritsFrom("RooAbsReal"):
                    ret = self.optimizeMHDependency(ret,self.wsp)
                _cache[(channel,process,syst)] = ret
                if not syst:
                  normname = "%s_norm" % (oname)
                  norm = self.wsp.arg(normname)
		  if norm==None: 
			if normname in self.norm_rename_map.keys(): norm = self.wsp.arg(self.norm_rename_map[normname])
                  if norm: 
                    if normname in self.DC.flatParamNuisances: 
                        self.DC.flatParamNuisances[normname] = False # don't warn if not found
                        norm.setAttribute("flatParam")
                    elif self.options.optimizeMHDependency:
                        norm = self.optimizeMHDependency(norm,self.wsp)
                    norm.SetName("shape%s_%s_%s%s_norm" % (postFix,process,channel, "_"))
		    self.norm_rename_map[normname]=norm.GetName()

		    # take care of any variables which were renamed (eg for "param")
		    paramString,renameParamString,toFreeze = self.getRenamingParameters()
		    if len(renameParamString):   self.out._import(norm, ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.RenameVariable(paramString,renameParamString))
                    else : self.out._import(norm, ROOT.RooFit.RecycleConflictNodes()) 
                if self.options.verbose > 2: print "import (%s,%s) -> %s\n" % (finalNames[0],objname,ret.GetName())
                return ret;
            elif self.wsp.ClassName() == "TTree":
                ##If it is a tree we will convert it in RooDataSet . Then we can decide if we want to build a
                ##RooKeysPdf or if we want to use it as an unbinned dataset 
                if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
                self.doVar("%s[%f,%f]" % (oname,self.wsp.GetMinimum(oname),self.wsp.GetMaximum(oname)))
                #Check if it is weighted
                self.doVar("__WEIGHT__[0.,1000.]")
                rds = ROOT.RooDataSet("shape%s_%s_%s%s" % (postFix,process,channel, "_"+syst if syst else ""), "shape%s_%s_%s%s" % (postFix,process,channel, "_"+syst if syst else ""),self.wsp,ROOT.RooArgSet(self.out.var(oname)),"","__WEIGHT__")
                rds.var = oname
                _cache[(channel,process,syst)] = rds
                if self.options.verbose > 2: print "import (%s,%s) -> %s\n" % (finalNames[0],wname,rds.GetName())
                return rds
            elif self.wsp.InheritsFrom("TH1"):
                ##If it is a Histogram we will convert it in RooDataSet preserving the bins 
                if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
                name = "shape%s_%s_%s%s" % (postFix,process,channel, "_"+syst if syst else "")
                # don't make it twice
                for X in _neverDelete:
                    if X.InheritsFrom("TNamed") and X.GetName() == name: return X
                self.doVar("%s[%f,%f]" % (oname,self.wsp.GetXaxis().GetXmin(),self.wsp.GetXaxis().GetXmax()))
                rds = ROOT.RooDataHist(name, name, ROOT.RooArgList(self.out.var(oname)), self.wsp)
                rds.var = oname
                if self.options.verbose > 2: stderr.write("import (%s,%s) -> %s\n" % (finalNames[0],wname,rds.GetName()))
                _neverDelete.append(rds)
                return rds
            else:
                raise RuntimeError, "Object %s in file %s has unrecognized type %s" (wname, finalNames[0], self.wsp.ClassName())
        else: # histogram
            ret = file.Get(objname);
            if not ret: 
                if allowNoSyst: return None
                raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (objname,finalNames[0],names[1],names[0])
            ret.SetName("shape%s_%s_%s%s" % (postFix,process,channel, "_"+syst if syst else ""))
            if self.options.verbose > 2: print "import (%s,%s) -> %s\n" % (finalNames[0],objname,ret.GetName())
            _cache[(channel,process,syst)] = ret
            return ret
    def getData(self,channel,process,syst="",_cache={}):
        return self.shape2Data(self.getShape(channel,process,syst),channel,process)
    def getPdf(self,channel,process,_cache={}):
        postFix="Sig" if (process in self.DC.isSignal and self.DC.isSignal[process]) else "Bkg"
        if _cache.has_key((channel,process)): return _cache[(channel,process)]
        shapeNominal = self.getShape(channel,process)
        nominalPdf = self.shape2Pdf(shapeNominal,channel,process) if (self.options.useHistPdf == "always" or shapeNominal == None) else shapeNominal
        if shapeNominal == None: return nominalPdf # no point morphing a fake shape
        morphs = []; shapeAlgo = None
	channelBinParFlag = channel in self.DC.binParFlags.keys()
        for (syst,nofloat,pdf,args,errline) in self.DC.systs:
            if not "shape" in pdf: continue
            if errline[channel][process] == 0: continue
            allowNoSyst = (pdf[-1] == "?")
            pdf = pdf.replace("?","")
            if pdf[-1] == "U": pdf = pdf[:-1]
            if shapeAlgo == None:  shapeAlgo = pdf
            elif pdf != shapeAlgo: 
                errmsg =  "ERROR for channel %s, process %s. " % (channel,process)
                errmsg += "Requesting morphing %s  for systematic %s after having requested %s. " % (pdf, syst, shapeAlgo)
                raise RuntimeError, errmsg+" One can use only one morphing algorithm for a given shape";
            if errline[channel][process] != 0:
                if allowNoSyst and not self.isShapeSystematic(channel,process,syst): continue
		systShapeName = syst
		if (syst,channel,process) in self.DC.systematicsShapeMap.keys(): systShapeName = self.DC.systematicsShapeMap[(syst,channel,process)]
                shapeUp   = self.getShape(channel,process,systShapeName+"Up")
                shapeDown = self.getShape(channel,process,systShapeName+"Down")
                if shapeUp.ClassName()   != shapeNominal.ClassName() and nominalPdf.ClassName()!="RooParametricHist": raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst %s" % (channel,process,syst)
                if shapeDown.ClassName() != shapeNominal.ClassName() and nominalPdf.ClassName()!="RooParametricHist": raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst %s" % (channel,process,syst)
                if self.options.useHistPdf == "always":
                    morphs.append((syst,errline[channel][process],self.shape2Pdf(shapeUp,channel,process),self.shape2Pdf(shapeDown,channel,process)))
                else:
                    morphs.append((syst,errline[channel][process],shapeUp,shapeDown))
        if len(morphs) == 0:
            if self.options.useHistPdf == "always":
                return nominalPdf
            else:
                return self.shape2Pdf(shapeNominal,channel,process)
        if shapeAlgo == "shapeN": stderr.write("Warning: the shapeN implementation in RooStats and L&S are different\n")
        pdfs = ROOT.RooArgList(nominalPdf) if self.options.useHistPdf == "always" else ROOT.TList()
        if self.options.useHistPdf != "always": pdfs.Add(nominalPdf)
        coeffs = ROOT.RooArgList()
        minscale = 1
        for (syst,scale,pdfUp,pdfDown) in morphs:
            if self.options.useHistPdf == "always": 
                pdfs.add(pdfUp); pdfs.add(pdfDown);
            else:
                pdfs.Add(pdfUp); pdfs.Add(pdfDown);
            if scale == 1:
                coeffs.add(self.out.var(syst))
            else: # must scale it :-/
                coeffs.add(self.doObj("%s_scaled_%s_%s" % (syst,channel,process), "prod","%s, %s" % (scale,syst)))
                if scale < minscale: minscale = scale
        qrange = minscale; qalgo = 0;
        if shapeAlgo[-1] == "*": 
            qalgo = 100
            shapeAlgo = shapeAlgo[:-1]
        if shapeAlgo == "shape": shapeAlgo = self.options.defMorph
        if "shapeL"   in shapeAlgo:
                raise RuntimeError, "No algorithm shapeL - this mode is depricated" 
		#qrange = 0;
        elif "shapeN" in shapeAlgo: qalgo = -1;
        if self.options.useHistPdf != "always":
            if nominalPdf.InheritsFrom("TH1"):
                rebins = ROOT.TList()
                maxbins = 0
                for i in xrange(pdfs.GetSize()):
                    rebinned = self.rebinH1(pdfs.At(i))
                    rebins.Add(rebinned)
                    maxbins = max(maxbins, rebinned._original_bins)
                if channelBinParFlag:
                    rhp = ROOT.CMSHistFunc("shape%s_%s_%s_morph" % (postFix,channel,process), "", self.out.var(self.TH1Observables[channel]), rebins[0])
                    rhp.setVerticalMorphs(coeffs)
                    rhp.setVerticalType(ROOT.CMSHistFunc.QuadLinear if qalgo >= 0 else ROOT.CMSHistFunc.LogQuadLinear)
                    rhp.setVerticalSmoothRegion(qrange)
                    rhp.prepareStorage()
                    rhp.setShape(0, 0, 0, 0, rebins[0])
                    for i in xrange(len(coeffs)):
                        if self.DC.binParFlags[channel][2] in [2]:
                            rhp.setShape(0, 0, i+1, 0, rebins[2 + i*2])
                            rhp.setShape(0, 0, i+1, 1, rebins[1 + i*2])
                        elif self.DC.binParFlags[channel][2] in [1]:
                            renormLo = rebins[2 + i*2].Clone()
                            if renormLo.Integral() > 0.:
                                renormLo.Scale(rebins[0].Integral() / renormLo.Integral())
                            renormHi = rebins[1 + i*2].Clone()
                            if renormHi.Integral() > 0.:
                                renormHi.Scale(rebins[0].Integral() / renormHi.Integral())
                            rhp.setShape(0, 0, i+1, 0, renormLo)
                            rhp.setShape(0, 0, i+1, 1, renormHi)
                    if self.options.optimizeTemplateBins and maxbins < self.out.maxbins:
                        #print "Optimizing binning: %d -> %d for %s " % (self.out.maxbins, maxbins, rhp.GetName())
                        rhp.setActiveBins(maxbins)

                else:
                    rhp = ROOT.FastVerticalInterpHistPdf2("shape%s_%s_%s_morph" % (postFix,channel,process), "", self.out.var(self.TH1Observables[channel]), rebins, coeffs, qrange, qalgo)
                    if self.options.optimizeTemplateBins and maxbins < self.out.maxbins:
                        #print "Optimizing binning: %d -> %d for %s " % (self.out.maxbins, maxbins, rhp.GetName())
                        rhp.setActiveBins(maxbins) 
                _cache[(channel,process)] = rhp
                return rhp
            elif nominalPdf.InheritsFrom("RooHistPdf") or nominalPdf.InheritsFrom("RooDataHist"):
                nominalPdf = self.shape2Pdf(shapeNominal,channel,process)
                pdfs.Clear(); ## now it contains "RooDataHist", it should contain "RooHistPdf"
                pdfs.Add(nominalPdf)
                for (syst,scale,shapeUp,shapeDown) in morphs:
                    pdfs.Add(self.shape2Pdf(shapeUp,channel,process))
                    pdfs.Add(self.shape2Pdf(shapeDown,channel,process))
                histpdf =  nominalPdf if nominalPdf.InheritsFrom("RooDataHist") else nominalPdf.dataHist()
                xvar = histpdf.get().first()
                rhp = ROOT.FastVerticalInterpHistPdf2("shape%s_%s_%s_morph" % (postFix,channel,process), "", xvar, pdfs, coeffs, qrange, qalgo)
                _cache[(channel,process)] = rhp
                return rhp
	    elif nominalPdf.InheritsFrom("RooParametricHist") : 
	        # Add the shape morphs to it. Cannot pass a collection of DataHists so we have to convert to PDFs first?! 
                for (syst,scale,shapeUp,shapeDown) in morphs:
		  nominalPdf.addMorphs(shapeUp,shapeDown,coeffs.find(syst),qrange)
                _cache[(channel,process)] = nominalPdf
		return nominalPdf
	
            else:
                pdflist = ROOT.RooArgList()
                nominalPdf = self.shape2Pdf(shapeNominal,channel,process)
                pdflist.add(nominalPdf)
                for (syst,scale,shapeUp,shapeDown) in morphs:
                    pdflist.add(self.shape2Pdf(shapeUp,channel,process))
                    pdflist.add(self.shape2Pdf(shapeDown,channel,process))
                pdfs = pdflist
		
        if "2a" in shapeAlgo: # old shape2
            if not nominalPdf.InheritsFrom("RooHistPdf"):  raise RuntimeError, "Algorithms 'shape2',  'shapeN2' only work with histogram templates"
            if nominalPdf.dataHist().get().getSize() != 1: raise RuntimeError, "Algorithms 'shape2',  'shapeN2' only work in one dimension"
            xvar = nominalPdf.dataHist().get().first()
            _cache[(channel,process)] = ROOT.VerticalInterpHistPdf("shape%s_%s_%s_morph" % (postFix,channel,process), "", xvar, pdfs, coeffs, qrange, qalgo)
        elif "2" in shapeAlgo:  # new faster shape2
            if not nominalPdf.InheritsFrom("RooHistPdf"):  raise RuntimeError, "Algorithms 'shape2',  'shapeN2' only work with histogram templates"
            if nominalPdf.dataHist().get().getSize() != 1: raise RuntimeError, "Algorithms 'shape2',  'shapeN2' only work in one dimension"
            xvar = nominalPdf.dataHist().get().first()
            _cache[(channel,process)] = ROOT.FastVerticalInterpHistPdf("shape%s_%s_%s_morph" % (postFix,channel,process), "", xvar, pdfs, coeffs, qrange, qalgo)
        else:
            _cache[(channel,process)] = ROOT.VerticalInterpPdf("shape%s_%s_%s_morph" % (postFix,channel,process), "", pdfs, coeffs, qrange, qalgo)
        return _cache[(channel,process)]
    def isShapeSystematic(self,channel,process,syst):
    	systShapeName = syst
	if (syst,channel,process) in self.DC.systematicsShapeMap.keys(): systShapeName = self.DC.systematicsShapeMap[(syst,channel,process)]
        shapeUp = self.getShape(channel,process,systShapeName+"Up",allowNoSyst=True)    
        return shapeUp != None
    def getExtraNorm(self,channel,process):
        if channel in self.selfNormBins and self.DC.binParFlags[channel][2] in [2]:
            if self.options.verbose > 1:
                print 'Skipping getExtraNorm for (%s,%s)' % (channel, process)
            return None
        postFix="Sig" if (process in self.DC.isSignal and self.DC.isSignal[process]) else "Bkg"
        terms = []
        shapeNominal = self.getShape(channel,process)
        if shapeNominal == None: 
            # FIXME no extra norm for dummy pdfs (could be changed)
            return None
        if shapeNominal.InheritsFrom("RooAbsPdf") and not shapeNominal.InheritsFrom("RooParametricHist") or shapeNominal.InheritsFrom("CMSHistFunc"): 
            # return nominal multiplicative normalization constant
            normname = "shape%s_%s_%s%s_norm" % (postFix,process,channel, "_")
            if self.out.arg(normname): return [ normname ]
            else: return None
        normNominal = 0
        if shapeNominal.InheritsFrom("TH1"): normNominal = shapeNominal.Integral()
        elif shapeNominal.InheritsFrom("RooDataHist"): normNominal = shapeNominal.sumEntries()
	elif shapeNominal.InheritsFrom("RooParametricHist"): 
	   normNominal = shapeNominal.quickSum()
           normname = "shape%s_%s_%s%s_norm" % (postFix,process,channel, "_")
           if self.out.arg(normname): terms.append(normname) 
        else: return None    
        if normNominal == 0: raise RuntimeError, "Null norm for channel %s, process %s" % (channel,process)
        for (syst,nofloat,pdf,args,errline) in self.DC.systs:
            if "shape" not in pdf: continue
            if errline[channel][process] != 0:
                if pdf[-1] == "?" and not self.isShapeSystematic(channel,process,syst): continue
		systShapeName = syst
	        if (syst,channel,process) in self.DC.systematicsShapeMap.keys(): systShapeName = self.DC.systematicsShapeMap[(syst,channel,process)]
                shapeUp   = self.getShape(channel,process,systShapeName+"Up")
                shapeDown = self.getShape(channel,process,systShapeName+"Down")
                if shapeUp.ClassName()   != shapeNominal.ClassName() and shapeNominal.ClassName()!="RooParametricHist": raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst %s" % (channel,process,syst)
                if shapeDown.ClassName() != shapeNominal.ClassName() and shapeNominal.ClassName()!="RooParametricHist": raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst %s" % (channel,process,syst)
                kappaUp,kappaDown = 1,1
                if shapeNominal.InheritsFrom("TH1"):
                    kappaUp,kappaDown = shapeUp.Integral(),shapeDown.Integral()
                elif shapeNominal.InheritsFrom("RooDataHist") or shapeNominal.InheritsFrom("RooParametricHist"):
                    kappaUp,kappaDown = shapeUp.sumEntries(),shapeDown.sumEntries()
                if not kappaUp > 0: raise RuntimeError, "Bogus norm %r for channel %s, process %s, systematic %s Up" % (kappaUp, channel,process,syst)
                if not kappaDown > 0: raise RuntimeError, "Bogus norm %r for channel %s, process %s, systematic %s Down" % (kappaDown, channel,process,syst)
                kappaUp /=normNominal; kappaDown /= normNominal
                if abs(kappaUp-1) < 1e-3 and abs(kappaDown-1) < 1e-3: continue
                # if errline[channel][process] == <x> it means the gaussian should be scaled by <x> before doing pow
                # for convenience, we scale the kappas
                kappasScaled = [ pow(x, errline[channel][process]) for x in kappaDown,kappaUp ]
                if self.options.packAsymPows:
                    terms.append( (kappasScaled[0], kappasScaled[1], syst) )
                else:
                    obj_kappaDown = self.addObj(ROOT.RooConstVar, '%f' %  kappasScaled[0], "", float('%f' %  kappasScaled[0]))
                    obj_kappaUp = self.addObj(ROOT.RooConstVar, '%f' %  kappasScaled[1], "", float('%f' %  kappasScaled[1]))
                    obj_var = self.out.var(syst)
                    self.addObj(ROOT.AsymPow, "systeff_%s_%s_%s" % (channel,process,syst), "", obj_kappaDown, obj_kappaUp, obj_var)
                    terms.append( "systeff_%s_%s_%s" % (channel,process,syst) )
        return terms if terms else None;

    def rebinH1(self,shape):
    	
	if self.options.optimizeTemplateBins:
          rebinh1 = ROOT.TH1F(shape.GetName()+"_rebin", "", self.out.maxbins, 0.0, float(self.out.maxbins))
          for i in range(1,min(shape.GetNbinsX(),self.out.maxbins)+1): 
            rebinh1.SetBinContent(i, shape.GetBinContent(i))
            rebinh1.SetBinError(i, shape.GetBinError(i))
          rebinh1._original_bins = shape.GetNbinsX()
	else :
	  shapeNbins = shape.GetNbinsX()
          rebinh1 = ROOT.TH1F(shape.GetName()+"_rebin", "", shapeNbins, 0.0, float(shapeNbins))
          for i in range(1,shapeNbins+1): 
            rebinh1.SetBinContent(i, shape.GetBinContent(i))
            rebinh1.SetBinError(i, shape.GetBinError(i))
          rebinh1._original_bins = shapeNbins
        return rebinh1;
	   
    def shape2Data(self,shape,channel,process,_cache={}):
        postFix="Sig" if (process in self.DC.isSignal and self.DC.isSignal[process]) else "Bkg"
        if shape == None:
            name = "shape%s_%s_%s" % (postFix,channel,process)
            if not _cache.has_key(name):
                obs = ROOT.RooArgSet(self.out.var("CMS_fakeObs"))
                obs.setRealValue("CMS_fakeObs",0.5);
                if self.out.mode == "binned":
                    self.out.var("CMS_fakeObs").setBins(1)
                    rdh = ROOT.RooDataHist(name, name, obs)
                    rdh.set(obs, self.DC.obs[channel])
                    _cache[name] = rdh
                else:
		    obs.add(self.out.var("CMS_fakeWeight"))
                    rds = ROOT.RooDataSet(name, name, obs,"CMS_fakeWeight")
		    obs.setRealValue("CMS_fakeWeight", self.DC.obs[channel])
                    rds.add(obs, self.DC.obs[channel])
                    _cache[name] = rds
            return _cache[name]
        if not _cache.has_key(shape.GetName()):
            if shape.ClassName().startswith("TH1"):
                rebinh1 = self.rebinH1(shape)
                rdh = ROOT.RooDataHist(shape.GetName(), shape.GetName(), ROOT.RooArgList(self.out.var(self.TH1Observables[channel])), rebinh1)
                #self.out._import(rdh)
                _cache[shape.GetName()] = rdh
            elif shape.ClassName() in ["RooDataHist", "RooDataSet"]:
                return shape
            else: raise RuntimeError, "shape2Data not implemented for %s" % shape.ClassName()
        return _cache[shape.GetName()]
    def shape2Pdf(self,shape,channel,process,_cache={}):
        postFix="Sig" if (process in self.DC.isSignal and self.DC.isSignal[process]) else "Bkg"
	channelBinParFlag = channel in self.DC.binParFlags.keys()
        if shape == None:
            name = "shape%s_%s_%s" % (postFix,channel,process)
            if not _cache.has_key(name):
                _cache[name] = ROOT.RooUniform(name, name, ROOT.RooArgSet(self.out.var("CMS_fakeObs")))
            return _cache[name]
        if not _cache.has_key(shape.GetName()+"Pdf"):
            if shape.ClassName().startswith("TH1"):
                if self.options.useHistPdf == "never":
                    shape = self.rebinH1(shape)
                    list = ROOT.TList(); list.Add(shape);
                    if channelBinParFlag:
                        rhp = ROOT.CMSHistFunc("%sPdf" % shape.GetName(), "", self.out.var(self.TH1Observables[channel]), shape)
                        rhp.prepareStorage()
                        rhp.setShape(0, 0, 0, 0, shape)
                        if self.options.optimizeTemplateBins:
                            rhp.setActiveBins(shape._original_bins)
                    else:
                        rhp = ROOT.FastVerticalInterpHistPdf2("%sPdf" % shape.GetName(), "", self.out.var(self.TH1Observables[channel]), list, ROOT.RooArgList())
                    _cache[shape.GetName()+"Pdf"] = rhp
                else:
                    rdh = self.shape2Data(shape,channel,process)
                    rhp = ROOT.RooHistPdf("%sPdf" % shape.GetName(), "", ROOT.RooArgSet(self.out.var(self.TH1Observables[channel])), rdh)
                    rhp.rdh = rdh # so it doesn't get deleted
                    _cache[shape.GetName()+"Pdf"] = rhp
            elif shape.InheritsFrom("RooAbsPdf"):
                if shape.ClassName() == "RooExtendPdf": 
                    raise RuntimeError, "Error in channel %s, process %s: pdf %s is a RooExtendPdf, this is not supported" % (channel,process,shape.GetName())
                elif shape.ClassName() == "RooAddPdf":
                    self.checkRooAddPdf(channel,process,shape)
                _cache[shape.GetName()+"Pdf"] = shape
            elif shape.InheritsFrom("RooDataHist"):
                rhp = ROOT.RooHistPdf("%sPdf" % shape.GetName(), "", shape.get(), shape) 
                _cache[shape.GetName()+"Pdf"] = rhp
            elif shape.InheritsFrom("RooDataSet"):
                rkp = ROOT.RooKeysPdf("%sPdf" % shape.GetName(), "", self.out.var(shape.var), shape,3,1.5); 
                _cache[shape.GetName()+"Pdf"] = rkp
            elif shape.InheritsFrom("CMSHistFunc"):
                _cache[shape.GetName()+"Pdf"] = shape
            else:
                raise RuntimeError, "shape2Pdf not implemented for %s" % shape.ClassName()
        return _cache[shape.GetName()+"Pdf"]
    def checkRooAddPdf(self,channel,process,pdf):
        coeflist = pdf.coefList()
        if (coeflist.getSize() == pdf.pdfList().getSize()-1):
            return True
        sum = 0.0
        for i in xrange(coeflist.getSize()):
            sum += coeflist.at(i).getVal()
        if abs(sum-1.0) > 1e-4:
            raise RuntimeError, "Error in channel %s, process %s: RooAddPdf %s has coefficients that sum up to %g, and not to unity. This is not supported (but it could be supported on request).\n" % (channel,process,pdf.GetName(),sum)
    def argSetToString(self,argset):
        names = []
        it = argset.createIterator()
        while True:
            arg = it.Next()
            if not arg: break
            names.append(arg.GetName())
        return ",".join(names)
    def optimizeExistingTemplates(self,pdf):
        if pdf.ClassName() == "FastVerticalInterpHistPdf2D":
            return ROOT.FastVerticalInterpHistPdf2D2(pdf, "%s_opt" % pdf.GetName())
        elif pdf.ClassName() == "RooProdPdf" and pdf.pdfList().getSize() == 2:
            f1 = pdf.pdfList().at(0)
            f2 = pdf.pdfList().at(1)
            #pdf.Print("")
            if f2.ClassName() == "FastVerticalInterpHistPdf2D" and f2.conditional():
                f1 = self.optimizeExistingTemplates(f1)
                f2 = ROOT.FastVerticalInterpHistPdf2D2(f2, "%s_opt" % f2.GetName()) 
                ret = ROOT.RooProdPdf("%s_opt" % pdf.GetName(), "", ROOT.RooArgSet(f1), ROOT.RooFit.Conditional(ROOT.RooArgSet(f2),ROOT.RooArgSet(f2.y())))
                ret.optf2 = f2
                ret.optf1 = f1
                #print "Optimize %s in \t" % (pdf.GetName()),; ret.Print("")
                return ret
        return pdf
    def optimizeMHDependency(self,arg,wsp,MH=None,indent=""):
        if arg.isFundamental(): return arg
        if not MH:
            MH = wsp.var("MH")
            if not MH: return arg
        if not arg.dependsOn(MH):
            #print "%s%s does not depend on MH" % (indent,arg.GetName())
            return arg
        depvars = arg.getVariables()
        if depvars.getSize() == 1:
            if not depvars.find("MH"):
                depvars.Print("")
                raise RuntimeError("???")
            depvars.setRealValue("MH", self.options.mass) # be safe
            if self.options.optimizeMHDependency == "fixed":
                print "%s%s depends only on MH, will freeze to its value at MH=%g, %g" % (indent, arg.GetName(), MH.getVal(), arg.getVal())
                ret = ROOT.RooConstVar("%s__frozenMH" % arg.GetName(), "", arg.getVal())
                self.out.dont_delete.append(ret)
            elif self.options.optimizeMHDependency in ("pol0", "pol1", "pol2", "pol3", "pol4"):
                order = int(self.options.optimizeMHDependency[3:])
                ret = ROOT.SimpleTaylorExpansion1D("%s__MHpol%d" % (arg.GetName(),order), "", arg, MH, 0.1, order)
                self.out.dont_delete.append(ret)
            else:
                raise RuntimeError("Unknown option value %r for optimizeMHDependency" % self.options.optimizeMHDependency)
            return ret
        else:
            print "%s%s depends on MH and other %d variables." % (indent, arg.GetName(), depvars.getSize())
            #depvars.Print("")
            srviter = arg.serverMIterator()
            servers = []
            while True:
                a = srviter.next()
                if not a: break
                servers.append(a)
            #print "%sFound %d servers: %s" % (indent, len(servers), ", ".join(a.GetName() for a in servers))
            newservers = []
            for a in servers:
                aopt = self.optimizeMHDependency(a,wsp,MH,indent=indent+"   ")
                if aopt != a:
                    newservers.append((a,aopt))
            if newservers:
                print "%sCan do replacements of %d servers: %s" % (indent, len(newservers), ", ".join(a.GetName() for (a,aopt) in newservers))
                #print "-- Before --"
                #arg.Print("t")
                cust = ROOT.RooCustomizer(arg,"__optMH")
                for (a,aopt) in newservers:
                    cust.replaceArg(a,aopt)
                ret = cust.build(True)
                self.out.dont_delete.append(ret)
                #print "-- After --"
                #ret.Print("t")
                #print "-- End --"
                arg = ret
            return arg

