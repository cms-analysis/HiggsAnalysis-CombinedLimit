from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class FA2_Interference_JHU_rw_HZZ(PhysicsModelBase_NiceSubclasses):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        xsecs = {
            "sigma1_HZZ": 290.58626,
            #          "sigma3_HZZ": 581.17253,
            "sigma3_HZZ": 44.670158,
            "sigma1_VBF": 968.674,
            "sigma3_VBF": 10909.54,
            "sigma1_ZH":  9022.36,
            "sigma3_ZH":  434763.7,
            "sigma1_WH":  30998.54,
            "sigma3_WH":  2028656,            
            "sigma2_HZZ": 105.85594,
            "sigmaL1_HZZ": 1.9846071e-06,
            "sigma2_VBF": 13102.71,
            "sigmaL1_VBF": 2.08309E-4,
            "sigma2_ZH": 713123,
            "sigmaL1_ZH": 33652.46e-6,
            "sigma2_WH": 3106339,
            "sigmaL1_WH": 11234.91e-5,
            "sigmaa1a3int_VBF": 0,
            "sigmaa1a3int_ZH": 0,
            "sigmaa1a3int_WH": 0,
            "sigmaa1a2int_VBF": 2207.73 ,
            "sigmaa1a2int_ZH": 4258.966,
            "sigmaa1a2int_WH": 16486.68,
            "sigmaa1L1int_VBF": 923.86500,
            "sigmaa1L1int_ZH": -11192.413,
            "sigmaa1L1int_WH": -36694.71,
            "yield_Powheg_ggH": 5.127128e+00,
            "yield_SM_ggH": 7.937751e+01,
            "yield_BSM_ggH": 8.952848e+01,
            "sigma_Powheg_ggH": 48.58,
            "sigma_SM_ggH": 15980,
            "sigma_BSM_ggH": 15981,
            "BR_H_tt": 0.0627


        }

        if not self.modelBuilder.out.var("CMS_zz4l_fai1"):
            self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");

#        self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");
        self.modelBuilder.doVar("RF[1.0,0,10]");
        self.modelBuilder.doVar("muTT[1.0,0,10]");
        self.modelBuilder.doVar("RV[1.0,0.0,10.0]");
#        self.modelBuilder.doSet("POI","CMS_zz4l_fai1,RV,RF,muTT")
        
        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::ai("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma2_HZZ})", CMS_zz4l_fai1)'.format(**xsecs))


        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigma2_VBF}/{sigma1_VBF})", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigma2_ZH}/{sigma1_ZH})", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigma2_WH}/{sigma1_WH})", RV,a1,ai,muTT)'.format(**xsecs))


        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@3*@1**2*{sigma2_VBF}/{sigma1_VBF} - @0*@1*@2*@3*sqrt({sigma2_VBF}/{sigma1_VBF})", RV,ai,a1,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@3*@1**2*{sigma2_ZH}/{sigma1_ZH} - @0*@1*@2*@3*sqrt({sigma2_ZH}/{sigma1_ZH})", RV,ai,a1,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@3*@1**2*{sigma2_WH}/{sigma1_WH} - @0*@1*@2*@3*sqrt({sigma2_WH}/{sigma1_WH})", RV,ai,a1,muTT)'.format(**xsecs))



        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*@3*sqrt({sigma2_VBF}/{sigma1_VBF})*{sigmaa1a2int_VBF}/{sigma1_VBF}", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*@3*sqrt({sigma2_ZH}/{sigma1_ZH})*{sigmaa1a2int_ZH}/{sigma1_ZH}", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*@3*sqrt({sigma2_WH}/{sigma1_WH})*{sigmaa1a2int_WH}/{sigma1_WH}", RV,a1,ai,muTT)'.format(**xsecs))

        return ["CMS_zz4l_fai1","RV","RF","muTT"] + super(FA2_Interference_JHU_rw_HZZ, self).getPOIList()


    def getYieldScale(self,bin,process):
        if process in ["GGH2Jets_sm_M",]:
            return 'RF*muTT'
#        if process in ["GGH2Jets_pseudoscalar_M",]:
#            return 'bsmCoupling_ggH'
        if process in ["reweighted_qqH_htt_0PM",]:
            return 'smCoupling_VBF'
        if process in ["reweighted_WH_htt_0PM",]:
            return 'smCoupling_WH'
        if process in ["reweighted_ZH_htt_0PM",]:
            return 'smCoupling_ZH'
        if process in ["reweighted_qqH_htt_0PH",]:
            return 'bsmCoupling_VBF'
        if process in ["reweighted_WH_htt_0PH",]:
            return 'bsmCoupling_WH'
        if process in ["reweighted_ZH_htt_0PH",]:
            return 'bsmCoupling_ZH'
        if process in ["reweighted_qqH_htt_0PHf05ph0"]:
            return 'intCoupling_VBF'
        if process in ["reweighted_ZH_htt_0PHf05ph0"]:
            return 'intCoupling_ZH'
        if process in ["reweighted_WH_htt_0PHf05ph0"]:
            return 'intCoupling_WH'
        #return 1
        return super(FA2_Interference_JHU_rw_HZZ, self).getYieldScale(bin, process)


FA2_Interference_JHU_rw_HZZ = FA2_Interference_JHU_rw_HZZ()
