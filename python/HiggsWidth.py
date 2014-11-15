from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the Higgs width

class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.RVRFfixed = False
        self.GGsmRVRFfixed = False
        self.hasACfai1 = False
        self.forbidPMF = False
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
            if 'GGsmfixed' in po:
                print "Will fix CMS_zz4l_GGsm to 1 and float muV, muF"
                self.GGsmfixed = True
            if 'RVRFfixed' in po:
                print "Will fix muV, muF to 1 and float mu"
                self.RVRFfixed = True
            if 'GGsmRVRFfixed' in po:
                print "Will fix muV, muF, CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmRVRFfixed = True
            if 'is2l2nu' in po:
                print "Will consider cards in 2l2nu style (separated S, B, S+B+I)"
                self.is2l2nu = True
            if 'ACfai1' in po:
                print "Model will consider fai1 for anomalous couplings onshell and offshell. Notice that it is not going to be the single POI"
                self.hasACfai1 = True
                if 'forbidPMF' in po:
                    print "Negative real phase will be forbidden for anomalous couplings"
                    self.forbidPMF = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if not self.modelBuilder.out.var("R"):
            self.modelBuilder.doVar("R[1.,0.,100.]")
        if not self.modelBuilder.out.var("RV"):
            self.modelBuilder.doVar("RV[1.,0.,100.]")
        if not self.modelBuilder.out.var("RF"):
            self.modelBuilder.doVar("RF[1.,0.,100.]")
        if not self.modelBuilder.out.var("CMS_zz4l_GGsm"):
            self.modelBuilder.doVar("CMS_zz4l_GGsm[1.,0.,50.]")
        if not self.modelBuilder.out.var("CMS_widthH_kbkg"):
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

        if self.hasACfai1:
            if not self.modelBuilder.out.var("CMS_zz4l_fai1"):
                print "Could not detect fai1, building a new one"
                self.modelBuilder.doVar("CMS_zz4l_fai1[0,-1.,1.]")
            self.modelBuilder.out.var("CMS_zz4l_fai1").setVal(0)
            if self.forbidPMF:
                print "fai1 cannot fall below 0"
                self.modelBuilder.out.var("CMS_zz4l_fai1").setRange(0,1)
            poi += ",CMS_zz4l_fai1"
        else:
            if self.modelBuilder.out.var("CMS_zz4l_fai1"):
                print "Found fai1 but will fix it to 0"
                self.modelBuilder.out.var("CMS_zz4l_fai1").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_fai1").setConstant()

        self.modelBuilder.factory_("expr::ggH_s_func(\"@0*@1*@3-sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")
        self.modelBuilder.factory_("expr::ggH_b_func(\"@2-sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")
        self.modelBuilder.factory_("expr::ggH_sbi_func(\"sqrt(@0*@1*@2*@3)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg,RF)")

        self.modelBuilder.factory_("expr::qqH_s_func(\"@0*@1*@2-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        self.modelBuilder.factory_("expr::qqH_b_func(\"1-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        self.modelBuilder.factory_("expr::qqH_sbi_func(\"sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,RV)")
        
        
        self.modelBuilder.doSet("POI",poi)
        
higgswidth = Higgswidth()
