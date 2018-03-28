from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class FA3_Interference_wrongS3HZZ(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        xsecs = {
            "sigma1_HZZ": 290.58626,
            "sigma3_HZZ": 581.17253,
#            "sigma3_HZZ": 44.670158,
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
            "sigmaL1_WH": 11234.91e-5
        }

        self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");
        self.modelBuilder.doVar("muf[1.0,0,10]");
        self.modelBuilder.doVar("muV[1.0,0.0,10.0]");
        #self.modelBuilder.doSet("POI","CMS_zz4l_fai1,muV,muf")
        self.modelBuilder.doSet("POI","CMS_zz4l_fai1,muV,muf")
#        self.modelBuilder.out.var("muf").setAttribute("flatParam")
        
        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::a3("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma3_HZZ})", CMS_zz4l_fai1)'.format(**xsecs))

        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a1,a3)'.format(**xsecs))
        
        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a3,a1)'.format(**xsecs))

        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*2", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*2", muV,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*2", muV,a1,a3)'.format(**xsecs))


    def getYieldScale(self,bin,process):
        if process in ["ggH_htt",]:
            return 'muf'
        if process in ["qqH_htt",]:
            return 'smCoupling_VBF'
        if process in ["WH_htt",]:
            return 'smCoupling_WH'
        if process in ["ZH_htt",]:
            return 'smCoupling_ZH'
        if process in ["qqH_htt_0M",]:
            return 'bsmCoupling_VBF'
        if process in ["WH_htt_0M",]:
            return 'bsmCoupling_WH'
        if process in ["ZH_htt_0M",]:
            return 'bsmCoupling_ZH'
        if process in ["qqH_htt_0Mf05ph0"]:
            return 'intCoupling_VBF'
        if process in ["ZH_htt_0Mf05ph0"]:
            return 'intCoupling_ZH'
        if process in ["WH_htt_0Mf05ph0"]:
            return 'intCoupling_WH'
        return 1


FA3_Interference_wrongS3HZZ = FA3_Interference_wrongS3HZZ()
