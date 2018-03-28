from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class FA3_Interference_JHU_ggHSyst_rw(PhysicsModel):
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
            "sigmaa1a3int_VBF": 1937.15,
            "sigmaa1a3int_ZH": 18044.72,
            "sigmaa1a3int_WH": 61997.07,
            # for ggH what is listed are yields in histos (all normliazed to powheg xsection):
            "yield_Powheg_ggH": 5.127128e+00,
            "yield_SM_ggH": 7.937751e+01,
            "yield_BSM_ggH": 8.952848e+01,
            "sigma_Powheg_ggH": 48.58,
            "sigma_SM_ggH": 15980,
            "sigma_BSM_ggH": 15981,
            "BR_H_tt": 0.0627

        }

        self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");
        self.modelBuilder.doVar("fa3_ggH[0.0,0.0,1.0]");
        self.modelBuilder.doVar("muf[1.0,0,10]");
        self.modelBuilder.doVar("muV[1.0,0.0,10.0]");
        #self.modelBuilder.doSet("POI","CMS_zz4l_fai1,muV,muf")
        self.modelBuilder.doSet("POI","CMS_zz4l_fai1,fa3_ggH,muV,muf")
#        self.modelBuilder.out.var("muf").setAttribute("flatParam")
        
        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::a3("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma3_HZZ})", CMS_zz4l_fai1)'.format(**xsecs))

        self.modelBuilder.doVar('expr::a1_ggH("sqrt(1-abs(@0))", fa3_ggH)')
        self.modelBuilder.doVar('expr::a3_ggH("sqrt(abs(@0))", fa3_ggH)'.format(**xsecs))


        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a1,a3)'.format(**xsecs))

# use yield from Powheg sample to fix the normalization
#        self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@1**2*{yield_Powheg_ggH}/{yield_SM_ggH}", muf,a1_ggH)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@1**2*1./13.45", muf,a1_ggH)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@1**2", muf,a1_ggH)'.format(**xsecs))

        
#        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a3,a1)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a3,a1)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a3,a1)'.format(**xsecs))

        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a3,a1)'.format(**xsecs))


#        self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@1**2*{yield_Powheg_ggH}/{yield_SM_ggH}", muf,a3_ggH)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@1**2*1./13.45", muf,a3_ggH)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@1**2", muf,a3_ggH)'.format(**xsecs))


#        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*2", muV,a1,a3)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*2", muV,a1,a3)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*2", muV,a1,a3)'.format(**xsecs))

        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*{sigmaa1a3int_VBF}/{sigma1_VBF}", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*{sigmaa1a3int_ZH}/{sigma1_ZH}", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*{sigmaa1a3int_WH}/{sigma1_WH}", muV,a1,a3)'.format(**xsecs))


    def getYieldScale(self,bin,process):
#        if process in ["ggH_htt",]:
#            return 'muf'
        if process in ["GGH2Jets_sm_M",]:
            return 'smCoupling_ggH'
        if process in ["GGH2Jets_pseudoscalar_M",]:
            return 'bsmCoupling_ggH'
        if process in ["reweighted_qqH_htt_0PM",]:
            return 'smCoupling_VBF'
        if process in ["reweighted_WH_htt_0PM",]:
            return 'smCoupling_WH'
        if process in ["reweighted_ZH_htt_0PM",]:
            return 'smCoupling_ZH'
        if process in ["reweighted_qqH_htt_0M",]:
            return 'bsmCoupling_VBF'
        if process in ["reweighted_WH_htt_0M",]:
            return 'bsmCoupling_WH'
        if process in ["reweighted_ZH_htt_0M",]:
            return 'bsmCoupling_ZH'
        if process in ["reweighted_qqH_htt_0Mf05ph0"]:
            return 'intCoupling_VBF'
        if process in ["reweighted_ZH_htt_0Mf05ph0"]:
            return 'intCoupling_ZH'
        if process in ["reweighted_WH_htt_0Mf05ph0"]:
            return 'intCoupling_WH'
        return 1


FA3_Interference_JHU_ggHSyst_rw = FA3_Interference_JHU_ggHSyst_rw()
