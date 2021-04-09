from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from HiggsAnalysis.CombinedLimit.LHCHCGModels import * 

import ROOT

class KappaVKappaT(LHCHCGBaseModel):
    """
    Copy of Kappas model with a combined kappa_V (for kappa_W and kappa_Z),
    and where hcc is independent of kappa_t.
   
    For tHq multilepton analysis (HIG-17-005)

    NOTE - Do not use this model for a generic analysis, 
    instead use the LHCHCGModels:K3 or K7 models   and freeze POIs accordingly
    """
    def __init__(self,resolved=True,BRU=True,addInvisible=False,coupleTopTau=False):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU
        self.resolved = resolved
        self.addInvisible = addInvisible
        self.coupleTopTau = coupleTopTau
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = (po.replace("BRU=","") in [ "yes", "1", "Yes", "True", "true" ])
        print "BR uncertainties in partial widths: %s " % self.doBRU
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("r[1,0.0,10.0]")
        self.modelBuilder.doVar("kappa_V[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_t[1,-10.0,10.0]")
        self.modelBuilder.doVar("kappa_mu[1,0.0,5.0]")
        if not self.coupleTopTau:
            self.modelBuilder.doVar("kappa_tau[1,0.0,3.0]")
            self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")
        else:
            self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_t)")
        self.modelBuilder.doVar("kappa_b[1,0.0,3.0]")
        self.modelBuilder.doVar("kappa_c[1,0.0,3.0]") # treat hcc independently from kappa_t
        if not self.resolved:
            self.modelBuilder.doVar("kappa_g[1,0.0,2.0]")
            self.modelBuilder.doVar("kappa_gam[1,0.0,2.5]")
        self.modelBuilder.doVar("BRinv[0,0,1]")
        if not self.addInvisible: self.modelBuilder.out.var("BRinv").setConstant(True)
        pois = 'kappa_V,kappa_t,kappa_b,kappa_c'
        if not self.coupleTopTau:
            pois += ',kappa_tau'
        if not self.resolved:
            pois += ',kappa_g,kappa_gam'
        if self.addInvisible: pois+=",BRinv"
        self.doMH()
        self.modelBuilder.doSet("POI",pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        # SM BR
        for d in SM_HIGG_DECAYS + [ "hss" ]:
            self.SMH.makeBR(d)
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()
        else:
            for d in SM_HIGG_DECAYS:
                self.modelBuilder.factory_('HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % d)
        # get VBF, tHq, tHW, ggZH cross section
        self.SMH.makeScaling('qqH', CW='kappa_V', CZ='kappa_V')
        self.SMH.makeScaling("tHq", CW='kappa_V', Ctop="kappa_t")
        self.SMH.makeScaling("tHW", CW='kappa_V', Ctop="kappa_t")
        self.SMH.makeScaling("ggZH", CZ='kappa_V', Ctop="kappa_t",Cb="kappa_b")
        # resolve loops
        if self.resolved:
            self.SMH.makeScaling('ggH', Cb='kappa_b', Ctop='kappa_t', Cc="kappa_t")
            self.SMH.makeScaling('hgluglu', Cb='kappa_b', Ctop='kappa_t')
            if not self.coupleTopTau:
                self.SMH.makeScaling('hgg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_V', Ctau='kappa_tau')
                self.SMH.makeScaling('hzg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_V', Ctau='kappa_tau')
            else:
                self.SMH.makeScaling('hgg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_V', Ctau='kappa_t')
                self.SMH.makeScaling('hzg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_V', Ctau='kappa_t')
        else:
            self.modelBuilder.factory_('expr::Scaling_hgluglu("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_hgg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_hzg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_ggH_7TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_8TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_13TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_14TeV("@0*@0", kappa_g)')

        ## partial witdhs, normalized to the SM one
        self.modelBuilder.factory_('expr::c7_Gscal_Z("@0*@0*@1*@2", kappa_V, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)')
        self.modelBuilder.factory_('expr::c7_Gscal_W("@0*@0*@1*@2", kappa_V, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)')
        if not self.coupleTopTau:
            self.modelBuilder.factory_('expr::c7_Gscal_tau("@0*@0*@1*@4+@2*@2*@3*@5", kappa_tau, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)')
        else:
            self.modelBuilder.factory_('expr::c7_Gscal_tau("@0*@0*@1*@4+@2*@2*@3*@5", kappa_t, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)')
        self.modelBuilder.factory_('expr::c7_Gscal_top("@0*@0 * @1*@2", kappa_c, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_Gscal_bottom("@0*@0 * (@1*@3+@2)", kappa_b, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)')
        self.modelBuilder.factory_('expr::c7_Gscal_gluon("  @0  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)')
        self.modelBuilder.factory_('expr::c7_Gscal_gamma("@0*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)')
        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::c7_SMBRs(%s)" %  (",".join("SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("c7_SMBRs").Print("")

        ## total witdh, normalized to the SM one
        self.modelBuilder.factory_('expr::c7_Gscal_tot("(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0)", BRinv, c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs)')

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM)
        self.modelBuilder.factory_('expr::c7_BRscal_hww("@0*@0*@2/@1", kappa_V, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzz("@0*@0*@2/@1", kappa_V, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)')
        if not self.coupleTopTau:
            self.modelBuilder.factory_('expr::c7_BRscal_htt("@0*@0*@2/@1", kappa_tau, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)')
        else:
            self.modelBuilder.factory_('expr::c7_BRscal_htt("@0*@0*@2/@1", kappa_t, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)')
        self.modelBuilder.factory_('expr::c7_BRscal_hmm("@0*@0*@2/@1", kappa_mu_expr, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)')
        self.modelBuilder.factory_('expr::c7_BRscal_hbb("@0*@0*@2/@1", kappa_b, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)')
        self.modelBuilder.factory_('expr::c7_BRscal_hcc("@0*@0*@2/@1", kappa_c, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgg("@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzg("@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgluglu("@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)')

        self.modelBuilder.factory_('expr::c7_BRscal_hinv("@0", BRinv)')

    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c7_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name) == None:
            if production in [ "ggH", "qqH", "ggZH", "tHq", "tHW"]:
                XSscal = ("@0", "Scaling_%s_%s" % (production,energy) )
            elif production == "WH":  XSscal = ("@0*@0", "kappa_V")
            elif production == "ZH":  XSscal = ("@0*@0", "kappa_V")
            elif production == "ttH": XSscal = ("@0*@0", "kappa_t")
            elif production == "bbH": XSscal = ("@0*@0", "kappa_b")
            else: raise RuntimeError, "Production %s not supported" % production
            BRscal = decay
            if not self.modelBuilder.out.function("c7_BRscal_"+BRscal):
                raise RuntimeError, "Decay mode %s not supported" % decay
            if decay == "hss": BRscal = "hbb"
            if production in ['tHq', 'tHW', 'ttH']:
                self.modelBuilder.factory_('expr::%s("%s*@1*@2", %s, c7_BRscal_%s, r)' % (name, XSscal[0], XSscal[1], BRscal))
            elif production == "ggH" and (decay in self.add_bbH) and energy in ["7TeV","8TeV","13TeV","14TeV"]:
                b2g = "CMS_R_bbH_ggH_%s_%s[%g]" % (decay, energy, 0.01)
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_('expr::%s("(%s + @1*@1*@2*@3)*@4", %s, kappa_b, %s, %s, c7_BRscal_%s)' % (name, XSscal[0], XSscal[1], b2g, b2gs, BRscal))
            else:
                self.modelBuilder.factory_('expr::%s("%s*@1*@2", %s, c7_BRscal_%s,r)' % (name, XSscal[0], XSscal[1], BRscal))
            print '[LHC-HCG Kappas]', name, production, decay, energy,": ",
            self.modelBuilder.out.function(name).Print("")
        return name


K4 = KappaVKappaT(resolved=True)
K5 = KappaVKappaT(resolved=False)
K6 = KappaVKappaT(resolved=False, coupleTopTau=True)
K7 = KappaVKappaT(resolved=True, coupleTopTau=True)
