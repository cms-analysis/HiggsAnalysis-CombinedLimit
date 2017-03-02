from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import ROOT, os

class TrilinearHiggs(SMLikeHiggsModel):
    "Float independently cross sections and branching ratios"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.mHRange = []
        self.poiNames = []
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
	
	# trilinear Higgs couplings modified 
	self.modelBuilder.doVar("k_lambda[1,-20,20]")
	self.poiNames="k_lambda"

        # --- Higgs Mass as other parameter ----
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                self.poiNames += ",MH"
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                self.poiNames += ",MH" 
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

	# now do the scaling - taken from Tab. 1 of https://arxiv.org/pdf/1607.04251v1.pdf
	cGammap = {"hgg":0.49e-2,"hzz":0.83e-2,"hww":0.73e-2,"hgluglu":0.66e-2,"htt":0,"hbb":0,"hcc":0,"hmm":0}
	cGTot   =  2.5e-3

	# taken from Tab. 2 of https://arxiv.org/pdf/1607.04251v1.pdf - need to recalculate these for specific analyses 
	cXSmap_7   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.06e-2,"ZH":1.23e-2,"ttH":3.87e-2}
	cXSmap_8   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.05e-2,"ZH":1.22e-2,"ttH":3.78e-2}
	cXSmap_13  = {"ggH":0.66e-2,"qqH":0.64e-2,"WH":1.03e-2,"ZH":1.19e-2,"ttH":3.51e-2}
	cXSmaps = {"7TeV":cXSmap_7, "8TeV":cXSmap_8, "13TeV":cXSmap_13}
	dZH = -1.536e-3
	
        self.modelBuilder.factory_("expr::C2(\"%g/(1-@0*@0*%g)\",k_lambda)" % (dZH,dZH))

	for sqrts in [7,8,13]: 
	  for proc in ["ggH","qqH","WH","ZH","ttH"]: 
	    valC1 = cXSmaps["%dTeV"%sqrts][proc]
            self.modelBuilder.factory_("expr::XSscal_%s_%dTeV(\"1+(@0-1)*%g + (@0*@0-1)*@1\",k_lambda,C2)" % (proc,sqrts,valC1))

	for dec in ["hgg","hzz","hww","hgluglu","htt","hbb","hcc","hmm"]: 
	    valC1 = cGammap[dec]
            self.modelBuilder.factory_("expr::BRscal_%s(\"1+((@0-1)*(%g-%g)/(1+(@0-1)*%g))\",k_lambda)" % (dec,valC1,cGTot,cGTot))
	    
	# make a dummy scaler of one ?
        self.modelBuilder.doVar("ONE[1,1,1]")
	self.modelBuilder.out.var("ONE").setConstant(True)

        print self.poiNames
        self.modelBuilder.doSet("POI",self.poiNames)

    def getHiggsSignalYieldScale(self,production,decay, energy):

        name = "XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name): return name

	prodname  = "XSscal_%s_%s"%(production,energy)
	decayname = "BRscal_%s"%decay

	if not self.modelBuilder.out.function("XSscal_%s_%s"%(production,energy)): prodname = "ONE"
	if not self.modelBuilder.out.function("BRscal_%s"%(decay)): decayname = "ONE"

	print " Building ---  XSBRscal_%s_%s_%s" %(production,decay,energy), " Using ", prodname,decayname
        self.modelBuilder.factory_("expr::XSBRscal_%s_%s_%s(\"@0*@1\",%s,%s)" % (production,decay,energy,prodname,decayname))
	   
	return name

trilinearHiggs = TrilinearHiggs()
