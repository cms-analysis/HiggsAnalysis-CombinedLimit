from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the Higgs width

class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.RVRFfixed = False
        self.GGsmRVRFfixed = False
        self.is2l2nu = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        if process == "ggH_s": return "ggH_s_func"
        elif process == "ggH_b": return "ggH_b_func"
        elif process == "ggH_sbi": return "ggH_sbi_func"
        if process == "qqH_s": return "qqH_s_func"
        elif process == "qqH_b": return "qqH_b_func"
        elif process == "qqH_sbi": return "qqH_sbi_func"
        elif process in ["ggH","ttH"]:
            if self.RVRFfixed or self.GGsmRVRFfixed:
                return "R"
            else:
                return "RF"
                
        elif process in ["qqH","WH","ZH","VH"]:
            if self.RVRFfixed or self.GGsmRVRFfixed:
                return "R"
            else:
                return "RV"
        else:
            return 1
            
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "GGsmfixed":
                print "Will fix CMS_zz4l_GGsm to 1 and float muV, muF"
                self.GGsmfixed = True
            if po == "RVRFfixed":
                print "Will fix muV, muF to 1 and float mu"
                self.RVRFfixed = True
            if po == "GGsmRVRFfixed":
                print "Will fix muV, muF, CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmRVRFfixed = True
            if po == "is2l2nu":
                print "Will consider cards in 2l2nu style (separated S, B, S+B+I)"
                self.is2l2nu = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("R[1.,0.,4]")
        self.modelBuilder.doVar("RV[1.,0.,4]")
        self.modelBuilder.doVar("RF[1.,0.,4]")
        self.modelBuilder.doVar("CMS_zz4l_GGsm[1.,0.,50.]")
        self.modelBuilder.doVar("CMS_widthH_kbkg[1.,0.,2.]")
	self.modelBuilder.out.var("R").setVal(1)
	self.modelBuilder.out.var("RV").setVal(1)
	self.modelBuilder.out.var("RF").setVal(1)
	self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(1)
	self.modelBuilder.out.var("CMS_widthH_kbkg").setVal(1)

	if self.GGsmfixed:
	    self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(True)
            print "Fixing CMS_zz4l_GGsm and R"
            self.modelBuilder.out.var("R").setConstant(True)
            poi = "RV,RF"
	elif self.GGsmRVRFfixed:
	    self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(True)
            print "Fixing CMS_zz4l_GGsm and RV, RF"
            self.modelBuilder.out.var("RV").setConstant(True)
            self.modelBuilder.out.var("RF").setConstant(True)
            poi = "R"
        else:
            if self.RVRFfixed:
                self.modelBuilder.out.var("RV").setConstant(True)
                self.modelBuilder.out.var("RF").setConstant(True)
            else:
                self.modelBuilder.out.var("R").setConstant(True)
            poi = "CMS_zz4l_GGsm"

	self.modelBuilder.factory_("expr::ggH_s_func(\"@0*@1*@3-sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")
        self.modelBuilder.factory_("expr::ggH_b_func(\"@2-sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")
        self.modelBuilder.factory_("expr::ggH_sbi_func(\"sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")

        self.modelBuilder.factory_("expr::qqH_s_func(\"@0*@1*@2-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        self.modelBuilder.factory_("expr::qqH_b_func(\"1-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        self.modelBuilder.factory_("expr::qqH_sbi_func(\"sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        
        
        self.modelBuilder.doSet("POI",poi)
        
higgswidth = Higgswidth()
