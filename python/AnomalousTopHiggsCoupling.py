# This file contains to classes needed for anomalous top-Higgs coupling fit
#
# AnomalousTopHiggsBuilder replaces the makeScaling function of SMHiggsBuilder
# with one that includes the pseudoscalar coupling modifier Ctildetop
#
# KappasAnomalousTopHiggs is a rewrite of the Kappas model in LHCHCGModels
# which includes the pseudoscalar coupling modifier kappa_tilde_top

from SMHiggsBuilder import *
from LHCHCGModels import *

class AnomalousTopHiggsBuilder(SMHiggsBuilder):

    def makeScaling(self,what, Cb='Cb', Ctop='Ctop', Ctildetop='Ctildetop', CW='CW', CZ='CZ', Ctau='Ctau', suffix=''):
        prefix = 'SM_%(what)s_' % locals()
        if suffix:
            suffix += '_'
            prefix += suffix
        if what.startswith('qqH'): # no change w.r.t. SM: qqH not affected by top-Higgs coupling
            for sqrts in ('7TeV', '8TeV'):
                rooName = prefix+'RVBF_'+sqrts
                self.textToSpline(rooName, os.path.join(self.coupPath, 'R_VBF_%(sqrts)s.txt'%locals()), ycol=1 )
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s(\
                "(@0*@0 + @1*@1 * @2 )/(1+@2)",\
                %(CW)s, %(CZ)s,\
                %(rooName)s\
)'%locals()
                self.modelBuilder.factory_(rooExpr)

        ## gluon fusion scaling
        # pseudoscalar cross section fraction defined in HiggsAnalysis/CombinedLimit/data/lhc-hxswg/couplings/ggA_*TeV.txt    
        elif what.startswith('ggH'): #
            structure = {'sigma_tt':2, 'sigma_bb':3, 'sigma_tb':4}
            structure_CPVtth = {'sigma_tildettildet':1}
            for sqrts in ('7TeV', '8TeV'):
                for qty, column in structure.iteritems():
                    rooName = prefix+qty+'_'+sqrts
                    self.textToSpline(rooName, os.path.join(self.coupPath, 'ggH_%(sqrts)s.txt'%locals()), ycol=column )
                # CP violating pseudoscalar admixture
                for qty, column in structure_CPVtth.iteritems():
                    rooName = prefix+qty+'_'+sqrts
                    self.textToSpline(rooName, os.path.join(self.coupPath, 'ggA_%(sqrts)s.txt'%locals()), ycol=column )
                scalingName = 'Scaling_'+what+'_'+sqrts
                # added pseudoscalar cross section to scaling
                rooExpr = 'expr::%(scalingName)s(\
                "(@0*@0)*@3  + (@1*@1)*@4 + (@0*@1)*@5 + (@2*@2)*@6",\
                %(Ctop)s, %(Cb)s, %(Ctildetop)s,\
                %(prefix)ssigma_tt_%(sqrts)s, %(prefix)ssigma_bb_%(sqrts)s, %(prefix)ssigma_tb_%(sqrts)s, %(prefix)ssigma_tildettildet_%(sqrts)s\
                )'%locals()
                self.modelBuilder.factory_(rooExpr)

        ## scaling for ttH
        # pseudoscalar cross section fraction defined in HiggsAnalysis/CombinedLimit/data/lhc-hxswg/couplings/ttA_*TeV.txt
        # does not exist in standard Higgs builder, as scaling is trivial there
        elif what.startswith('ttH'): #
            structure_CPVtth = {'sigma_tildettildet':1}
            for sqrts in ('7TeV', '8TeV'):
                for qty, column in structure_CPVtth.iteritems():
                    rooName = prefix+qty+'_'+sqrts
                    self.textToSpline(rooName, os.path.join(self.coupPath, 'ttA_%(sqrts)s.txt'%locals()), ycol=column )
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s(\
                "(@0*@0)  + (@1*@1)*@2",\
                %(Ctop)s, %(Ctildetop)s,\
                %(prefix)ssigma_tildettildet_%(sqrts)s\
                )'%locals()       
                self.modelBuilder.factory_(rooExpr)

        ## scaling for H -> gluon gluon
        # pseudoscalar cross section fraction defined in HiggsAnalysis/CombinedLimit/data/lhc-hxswg/couplings/Gamma_Agluongluon.txt
        elif what.startswith('hgluglu'):
            structure = {'Gamma_tt':2, 'Gamma_bb':3, 'Gamma_tb':4}
            structure_CPVtth = {'Gamma_tildettildet':1}
            for qty, column in structure.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, 'Gamma_Hgluongluon.txt'), ycol=column )
            # CP violating pseudoscalar admixture
            for qty, column in structure_CPVtth.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, 'Gamma_Agluongluon.txt'), ycol=column )
            scalingName = 'Scaling_'+what
            # modified scaling
            rooExpr = 'expr::%(scalingName)s(\
            "(@0*@0)*@3  + (@1*@1)*@4 + (@0*@1)*@5 + (@2*@2)*@6",\
            %(Ctop)s, %(Cb)s,%(Ctildetop)s,\
            %(prefix)sGamma_tt, %(prefix)sGamma_bb, %(prefix)sGamma_tb, %(prefix)sGamma_tildettildet\
            )'%locals()
            print  rooExpr
            self.modelBuilder.factory_(rooExpr)

        ## scaling for H -> gamma gamma
        # pseudoscalar cross section fraction defined in HiggsAnalysis/CombinedLimit/data/lhc-hxswg/couplings/Gamma_Agammagamma.txt
        elif what.startswith('hgg') or what.startswith('hzg'): #in ['hgg', 'hzg']:
            fileFor = {'hgg':'Gamma_Hgammagamma.txt',
                       'hzg':'Gamma_HZgamma.txt'}
            structure = {'Gamma_tt':2, 'Gamma_bb':3, 'Gamma_WW':4,
                         'Gamma_tb':5, 'Gamma_tW':6, 'Gamma_bW':7,
                         'Gamma_ll':8,
                         'Gamma_tl':9, 'Gamma_bl':10, 'Gamma_lW':11}
            fileForCPVtth = {'hgg':'Gamma_Agammagamma.txt',
                             'hzg':'Gamma_AZgamma.txt'}
            structureCPVtth = {'Gamma_tildettildet':1}
            for qty, column in structure.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, fileFor['hgg' if what.startswith('hgg') else 'hzg']), ycol=column )
            for qty, column in structureCPVtth.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, fileForCPVtth['hgg' if what.startswith('hgg') else 'hzg']), ycol=column )
            scalingName = 'Scaling_'+what
            # modified gamma gamma scaling
            rooExpr = 'expr::%(scalingName)s(\
            "( (@0*@0)*@5 + (@1*@1)*@6 + (@2*@2)*@7 + (@0*@1)*@8 + (@0*@2)*@9 + (@1*@2)*@10 + (@3*@3)*@11 + (@0*@3)*@12 + (@1*@3)*@13 + (@2*@3)*@14 + (@4*@4)*@15) / (@5+@6+@7+@8+@9+@10+@11+@12+@13+@14)",\
            %(Ctop)s, %(Cb)s, %(CW)s, %(Ctau)s, %(Ctildetop)s,\
            %(prefix)sGamma_tt, %(prefix)sGamma_bb, %(prefix)sGamma_WW,\
            %(prefix)sGamma_tb, %(prefix)sGamma_tW, %(prefix)sGamma_bW,\
            %(prefix)sGamma_ll,\
            %(prefix)sGamma_tl, %(prefix)sGamma_bl, %(prefix)sGamma_lW,\
            %(prefix)sGamma_tildettildet)'%locals()
            self.modelBuilder.factory_(rooExpr)

        elif what.startswith('ggZH'):
            for sqrts in ('7TeV', '8TeV','13TeV'):
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s( "(@0*@0)*2.27  + (@1*@1)*0.37 - (@0*@1)*1.64", %(CZ)s, %(Ctop)s)'%locals()
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('tHq'):
            for sqrts in ('7TeV', '8TeV','13TeV'):
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s( "(@0*@0)*3.4  + (@1*@1)*3.56 - (@0*@1)*5.96", %(Ctop)s, %(CW)s)'%locals()
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('tHW'):
            for sqrts in ('7TeV', '8TeV','13 TeV'):
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s( "(@0*@0)*1.84  + (@1*@1)*1.57 - (@0*@1)*2.41", %(Ctop)s, %(CW)s)'%locals()
                self.modelBuilder.factory_(rooExpr)
        else:
            raise RuntimeError, "There is no scaling defined for %(what)s" % locals()


# copy of the Kappas model from KHCHCGModels.py, modifications commented inline
class KappasAnomalousTopHiggs(LHCHCGBaseModel):
    "kappa model with additional pseudoscalar top-Higgs coupling"
    def __init__(self,resolved=True,BRU=True):
        LHCHCGBaseModel.__init__(self)
        self.doBRU = BRU
        self.resolved = resolved
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = (po.replace("BRU=","") in [ "yes", "1", "Yes", "True", "true" ])
        print "BR uncertainties in partial widths: %s " % self.doBRU
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kappa_W[1,0.0,2.0]") 
        self.modelBuilder.doVar("kappa_Z[1,0.0,2.0]") 
        self.modelBuilder.doVar("kappa_tau[1,0.0,3.0]")
        self.modelBuilder.doVar("kappa_mu[1,0.0,5.0]") 
        self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")
        self.modelBuilder.doVar("kappa_t[1,0.0,4.0]")
        # additional kappa for the anomalous coupling
        self.modelBuilder.doVar("kappa_tilde_t[0.0,0.0,4.0]")
        self.modelBuilder.doVar("kappa_b[1,0.0,3.0]")
        if not self.resolved:
            self.modelBuilder.doVar("kappa_g[1,0.0,2.0]")
            self.modelBuilder.doVar("kappa_gam[1,0.0,2.5]")
	self.modelBuilder.doVar("BRinv[0,0,1]")
        self.modelBuilder.out.var("BRinv").setConstant(True)
        # adding additional kappa to list of parameters of interest
        pois = 'kappa_W,kappa_Z,kappa_tau,kappa_t,kappa_tilde_t,kappa_b'
        if not self.resolved:
            pois += ',kappa_g,kappa_gam'
        self.doMH()
        self.modelBuilder.doSet("POI",pois)
        # use modified Higgs Builder
        self.SMH = AnomalousTopHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        self.dobbH()
        for d in SM_HIGG_DECAYS + [ "hss" ]: 
            self.SMH.makeBR(d)
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()
        else:
            for d in SM_HIGG_DECAYS: 
                self.modelBuilder.factory_('HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % d)
        self.SMH.makeScaling('qqH', CW='kappa_W', CZ='kappa_Z')
        # cross section depending on kappa_t and kappa_tilde_t
        self.SMH.makeScaling("tHq", CW='kappa_W', Ctop="kappa_t", Ctildetop="kappa_tilde_t")
        self.SMH.makeScaling("tHW", CW='kappa_W', Ctop="kappa_t", Ctildetop="kappa_tilde_t")
        self.SMH.makeScaling("ggZH", CZ='kappa_Z', Ctop="kappa_t", Ctildetop="kappa_tilde_t")
        # resolved loops: all cross section depending on kappa_t and kappa_tilde_t
        if self.resolved:
            self.SMH.makeScaling('ggH', Cb='kappa_b', Ctop='kappa_t', Ctildetop="kappa_tilde_t")
            self.SMH.makeScaling('hgluglu', Cb='kappa_b', Ctop='kappa_t', Ctildetop="kappa_tilde_t")
            self.SMH.makeScaling('hgg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_W', Ctau='kappa_tau', Ctildetop="kappa_tilde_t")
            self.SMH.makeScaling('hzg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_W', Ctau='kappa_tau', Ctildetop="kappa_tilde_t")
            self.SMH.makeScaling('ttH', Ctop='kappa_t', Ctildetop="kappa_tilde_t")
        # using unresolved mode possible, too
        else:
            self.modelBuilder.factory_('expr::Scaling_hgluglu("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_hgg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_hzg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_ggH_7TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_8TeV("@0*@0", kappa_g)')

        self.modelBuilder.factory_('expr::c7_Gscal_Z("@0*@0*@1*@2", kappa_Z, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)')
        self.modelBuilder.factory_('expr::c7_Gscal_W("@0*@0*@1*@2", kappa_W, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)')
        self.modelBuilder.factory_('expr::c7_Gscal_tau("@0*@0*@1*@4+@2*@2*@3*@5", kappa_tau, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)')
        # h--> cc BR uncertainty is now independent of kappa_t
        self.modelBuilder.factory_('expr::c7_Gscal_top("@1*@2", kappa_t, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_Gscal_bottom("@0*@0 * (@1*@3+@2)", kappa_b, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)')
        self.modelBuilder.factory_('expr::c7_Gscal_gluon("  @0  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)')
        self.modelBuilder.factory_('expr::c7_Gscal_gamma("@0*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_("sum::c7_SMBRs(%s)" %  (",".join("SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("c7_SMBRs").Print("")        

        ## total witdh, normalized to the SM one
        self.modelBuilder.factory_('expr::c7_Gscal_tot("(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0)", BRinv, c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs)')

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
        self.modelBuilder.factory_('expr::c7_BRscal_hww("@0*@0*@2/@1", kappa_W, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzz("@0*@0*@2/@1", kappa_Z, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)')
        self.modelBuilder.factory_('expr::c7_BRscal_htt("@0*@0*@2/@1", kappa_tau, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)')
        self.modelBuilder.factory_('expr::c7_BRscal_hmm("@0*@0*@2/@1", kappa_mu_expr, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)')
        self.modelBuilder.factory_('expr::c7_BRscal_hbb("@0*@0*@2/@1", kappa_b, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)')
        # h--> cc BR fixed to SM
        self.modelBuilder.factory_('expr::c7_BRscal_hcc("@2/@1", kappa_t, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgg("@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzg("@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgluglu("@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)')

    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c7_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name) == None:
            # implemented more complicated scaling for ttH, too
            if production in [ "ggH", "qqH", "ggZH", "tHq", "tHW", "ttH"]: 
                XSscal = ("@0", "Scaling_%s_%s" % (production,energy) )
            elif production == "WH":  XSscal = ("@0*@0", "kappa_W")
            elif production == "ZH":  XSscal = ("@0*@0", "kappa_Z")
            # ttH no longer scales with kappa_t only
#            elif production == "ttH": XSscal = ("@0*@0", "kappa_t")
            elif production == "bbH": XSscal = ("@0*@0", "kappa_b")
            else: raise RuntimeError, "Production %s not supported" % production
            BRscal = decay
            if not self.modelBuilder.out.function("c7_BRscal_"+BRscal):
                raise RuntimeError, "Decay mode %s not supported" % decay
            if decay == "hss": BRscal = "hbb"
            if production == "ggH" and (decay in self.add_bbH) and energy in ["7TeV","8TeV"]:
                b2g = "CMS_R_bbH_ggH_%s_%s[%g]" % (decay, energy, 0.01) 
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_('expr::%s("(%s + @1*@1*@2*@3)*@4", %s, kappa_b, %s, %s, c7_BRscal_%s)' % (name, XSscal[0], XSscal[1], b2g, b2gs, BRscal))
            else:
                self.modelBuilder.factory_('expr::%s("%s*@1", %s, c7_BRscal_%s)' % (name, XSscal[0], XSscal[1], BRscal))
            print '[LHC-HCG Kappas]', name, production, decay, energy,": ",
            self.modelBuilder.out.function(name).Print("")
        return name

anomTTH=KappasAnomalousTopHiggs()
