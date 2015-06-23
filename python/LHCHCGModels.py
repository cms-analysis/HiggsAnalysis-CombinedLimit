from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import ROOT

## Naming conventions
CMS_to_LHCHCG_Dec = { 
    'hww': 'WW',
    'hzz': 'ZZ',
    'hgg': 'gamgam',
    'hbb': 'bb',
    'hcc': 'cc',
    'htt': 'tautau',
    'hmm': 'mumu',
    'hzg': 'Zgam',
    'hgluglu': 'gluglu',
}
CMS_to_LHCHCG_DecSimple = { 
    'hww': 'WW',
    'hzz': 'ZZ',
    'hgg': 'gamgam',
    'hbb': 'bb',
    'hcc': 'bb',
    'htt': 'tautau',
    'hmm': 'tautau',
    'hzg': 'gamgam',
    'hgluglu': 'bb',
}
CMS_to_LHCHCG_Prod = { 
    'ggH': 'ggF',
    'qqH': 'VBF',
    'WH': 'WH',
    'ZH': 'qqZH',
    'ggZH': 'ggZH',
    'ttH': 'ttH',
    'tHq': 'tHjb',
    'tHW': 'WtH',
    'bbH': 'bbH',
 } 

class LHCHCGBaseModel(SMLikeHiggsModel):
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
        self.add_bbH = [ "hzz", "hgg" ]
        self.bbH_pdf = "pdf_gg" 
    def preProcessNuisances(self,nuisances):
        if self.add_bbH and not any(row for row in nuisances if row[0] == "QCDscale_bbH"):
            nuisances.append(("QCDscale_bbH",False, "param", [ "0", "1"], [] ) )
        if self.add_bbH and not any(row for row in nuisances if row[0] == self.bbH_pdf):
            nuisances.append((pdf_bbH,False, "param", [ "0", "1"], [] ) )
    def setPhysicsOptionsBase(self,physOptions):
        for po in physOptions:
            if po.startswith("bbh="):
                self.add_bbH = [d.strip() for d in po.replace("bbh=","").split(",")]
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
        print "Will add bbH to signals in the following Higgs boson decay modes: %s" % (", ".join(self.add_bbH))
    def dobbH(self):
        self.modelBuilder.doVar("QCDscale_bbH[-7,7]")
        self.modelBuilder.doVar(self.bbH_pdf+"[-7,7]")
        scaler7 = ROOT.ProcessNormalization("CMS_bbH_scaler_7TeV","",1.0)
        scaler8 = ROOT.ProcessNormalization("CMS_bbH_scaler_8TeV","",1.0)
        self.modelBuilder.out._import(scaler7)
        self.modelBuilder.out._import(scaler8)
        self.modelBuilder.out.function("CMS_bbH_scaler_7TeV").addAsymmLogNormal(1.0/114.5, 1.106, self.modelBuilder.out.var("QCDscale_bbH"))
        self.modelBuilder.out.function("CMS_bbH_scaler_8TeV").addAsymmLogNormal(1.0/114.8, 1.103, self.modelBuilder.out.var("QCDscale_bbH"))
        self.modelBuilder.out.function("CMS_bbH_scaler_7TeV").addLogNormal(1.061, self.modelBuilder.out.var(self.bbH_pdf))
        self.modelBuilder.out.function("CMS_bbH_scaler_8TeV").addLogNormal(1.062, self.modelBuilder.out.var(self.bbH_pdf))
    def doMH(self):
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

class SignalStrengths(LHCHCGBaseModel):
    "Allow different signal strength fits"
    def __init__(self):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.POIs = "mu"
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=","")
    def doVar(self,x,constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*","",x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print "SignalStrengths:: declaring %s as %s, and set to constant" % (vname,x)
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.doMH()
        self.doVar("mu[1,0,5]")
        for X in CMS_to_LHCHCG_Dec.values():
            self.doVar("mu_BR_%s[1,0,5]" % X)
        for X in CMS_to_LHCHCG_Prod.values() + [ "ZH", "tH", "ggFbbH", "ttHtH", "VH" ]:
            self.doVar("mu_XS_%s[1,0,5]" % X)
            self.doVar("mu_XS7_%s[1,0,5]" % X)
            self.doVar("mu_XS8_%s[1,0,5]" % X)
        for X in CMS_to_LHCHCG_DecSimple.values():
            self.doVar("mu_V_%s[1,0,5]" % X)
            self.doVar("mu_F_%s[1,0,5]" % X)
        print "Default parameters of interest: ", self.POIs
        self.modelBuilder.doSet("POI",self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()
    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH": continue # skip aggregates 
            for D in SM_HIGG_DECAYS:
                for E in 7, 8:
                    terms = [ "mu", "mu_BR_"+CMS_to_LHCHCG_DecSimple[D] ]
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01) 
                        ggs = ",".join([ "mu_XS_ggF", "mu_XS%d_ggF"%E ])
                        bbs = ",".join([ "mu_XS_bbH", "mu_XS%d_bbH"%E, "CMS_bbH_scaler_%dTeV"%E ])
                        ## FIXME should include the here also logNormal for QCDscale_bbH
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV(\"@1*@2+@0*@3*@4*@5\",%s,%s,%s)' % (D,E,b2g,ggs,bbs))
                        terms += [ 'ggH_bbH_sum_%s_%dTeV' % (D,E),  "mu_XS_ggFbbH", "mu_XS%d_ggFbbH"%E ]
                    else:
                        if P in [ "ggH", "bbH" ]:
                            terms += [ "mu_XS_ggFbbH", "mu_XS%d_ggFbbH"%E ]
                        terms += [ "mu_XS_"+CMS_to_LHCHCG_Prod[P],  "mu_XS%d_%s"%(E,CMS_to_LHCHCG_Prod[P])  ]
                    # Summary modes
                    if P in [ "tHW", "tHq" ]:
                        terms += [ "mu_XS_tH", "mu_XS%d_tH"%E ]
                    if P in [ "tHW", "tHq", "ttH" ]:
                        terms += [ "mu_XS_ttHtH", "mu_XS%d_ttHtH"%E ]
                    if P in [ "ggZH", "ZH" ]:
                        terms += [ "mu_XS_ZH", "mu_XS%d_ZH"%E ]
                    if P in [ "WH", "ZH", "ggZH" ]:
                        terms += [ "mu_XS_VH", "mu_XS%d_VH"%E ]
                    # for 2D scans
                    if P in [ "ggH", "ttH", "bbH", "tHq", "tHW" ]:
                        terms += [ "mu_F_"+CMS_to_LHCHCG_DecSimple[D] ]
                    else:
                        terms += [ "mu_V_"+CMS_to_LHCHCG_DecSimple[D] ]
                    self.modelBuilder.factory_('prod::scaling_%s_%s_%dTeV(%s)' % (P,D,E,",".join(terms)))
                    self.modelBuilder.out.function('scaling_%s_%s_%dTeV' % (P,D,E)).Print("")

    def getHiggsSignalYieldScale(self,production,decay,energy):
        return "scaling_%s_%s_%s" % (production,decay,energy)
class SignalStrengthRatios(LHCHCGBaseModel):
    "Allow for fits of ratios of signal strengths"
    def __init__(self):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=","")
    def doVar(self,x,constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*","",x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print "SignalStrengthRatios:: declaring %s as %s, and set to constant" % (vname,x)
    def doParametersOfInterest(self):
        """Create POI out of signal strength ratios and MH"""
        self.doMH()
        self.doVar("mu_V_r_F[1,0,5]")
        for X in CMS_to_LHCHCG_DecSimple.values():
            self.doVar("mu_F_%s[1,0,5]" % X)
            self.doVar("mu_V_r_F_%s[1,0,5]" % X)
        self.POIs = "mu_V_r_F"
        print "Default parameters of interest: ", self.POIs
        self.modelBuilder.doSet("POI",self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()
    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH": continue # skip aggregates 
            for D in SM_HIGG_DECAYS:
                for E in 7, 8:
                    terms = [ "mu_F_"+CMS_to_LHCHCG_DecSimple[D] ]
                    if P in [ "qqH","VBF", "VH", "WH", "ZH","qqZH", "ggZH" ]:
                        terms += [ "mu_V_r_F", "mu_V_r_F_"+CMS_to_LHCHCG_DecSimple[D] ]

                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        bbs = "CMS_bbH_scaler_%dTeV"%E
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV(\"1+@0*@1\",%s,%s)' % (D,E,b2g,bbs))
                        terms += [ 'ggH_bbH_sum_%s_%dTeV' % (D,E) ]

                    self.modelBuilder.factory_('prod::scaling_%s_%s_%dTeV(%s)' % (P,D,E,",".join(terms)))
                    self.modelBuilder.out.function('scaling_%s_%s_%dTeV' % (P,D,E)).Print("")

    def getHiggsSignalYieldScale(self,production,decay,energy):
        return "scaling_%s_%s_%s" % (production,decay,energy)

class XSBRratios(LHCHCGBaseModel):
    "Model with the ratios of cross sections and branching ratios "
    def __init__(self):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("poi="):
                self.POIs = po.replace("poi=","")
    def doVar(self,x,constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*","",x)
        self.modelBuilder.out.var(vname).setConstant(constant)
        print "XSBRratios:: declaring %s as %s, and set to constant" % (vname,x)

    def doParametersOfInterest(self):
        """Create POI out of signal strength ratios and MH"""
        self.doMH()
        self.doVar("mu_XS7_r_XS8_ggF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_VBF[1,0,5]")
        self.doVar("mu_XS7_r_XS8_WH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ZH[1,0,5]")
        self.doVar("mu_XS7_r_XS8_ttH[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ggF_x_BR_WW[1,0,5]")
        self.modelBuilder.doVar("mu_XS_VBF_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_WH_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ZH_r_XS_ggF[1,0,5]")
        self.modelBuilder.doVar("mu_XS_ttH_r_XS_ggF[1,0,5]")
        self.POIs = "mu_XS_ggF_x_BR_WW,mu_XS_VBF_r_XS_ggF,mu_XS_ttH_r_XS_ggF,mu_XS_WH_r_XS_ggF,mu_XS_ZH_r_XS_ggF"
        for X in ["ZZ","tautau","bb","gamgam"]:
            self.modelBuilder.doVar("mu_BR_%s_r_BR_WW[1,0,5]" % X)
            self.POIs += ",mu_BR_%s_r_BR_WW"%X
        print "Default parameters of interest: ", self.POIs
        self.modelBuilder.doSet("POI",self.POIs)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()
    def setup(self):
        self.dobbH()
        for P in ALL_HIGGS_PROD:
            if P == "VH": continue # skip aggregates 
            for D in SM_HIGG_DECAYS:
                for E in 7, 8:
                    terms = [ "mu_XS_ggF_x_BR_WW" ]
                    if CMS_to_LHCHCG_DecSimple[D] != "WW":
                        terms += [ "mu_BR_%s_r_BR_WW"%CMS_to_LHCHCG_DecSimple[D] ]
                    if P in ["ttH","tHq","tHW"]:
                        terms += [ "mu_XS_ttH_r_XS_ggF" ]
                    if P in ["ZH","ggZH"]:
                        terms += [ "mu_XS_ZH_r_XS_ggF" ]
                    if P in ["WH"]:
                        terms += [ "mu_XS_WH_r_XS_ggF" ]
                    if P in ["qqH"]:
                        terms += [ "mu_XS_VBF_r_XS_ggF" ]
                    if E == 7:
                        if P in ["ggH","bbH"]:
                            terms += [ "mu_XS7_r_XS8_ggF" ]
                        if P in ["ttH","tHq","tHW"]:
                            terms += [ "mu_XS7_r_XS8_ttH" ]
                        if P in ["ZH","ggZH"]:
                            terms += [ "mu_XS7_r_XS8_ZH" ]
                        if P in ["WH"]:
                            terms += [ "mu_XS7_r_XS8_WH" ]
                        if P in ["qqH"]:
                            terms += [ "mu_XS7_r_XS8_VBF" ]
                    # Hack for ggH
                    if D in self.add_bbH and P == "ggH":
                        b2g = "CMS_R_bbH_ggH_%s_%dTeV[%g]" % (D, E, 0.01)
                        bbs = "CMS_bbH_scaler_%dTeV"%E
                        self.modelBuilder.factory_('expr::ggH_bbH_sum_%s_%dTeV(\"1+@0*@1\",%s,%s)' % (D,E,b2g,bbs))
                        terms += [ 'ggH_bbH_sum_%s_%dTeV' % (D,E) ]

                    self.modelBuilder.factory_('prod::scaling_%s_%s_%dTeV(%s)' % (P,D,E,",".join(terms)))
                    self.modelBuilder.out.function('scaling_%s_%s_%dTeV' % (P,D,E)).Print("")

    def getHiggsSignalYieldScale(self,production,decay,energy):
        return "scaling_%s_%s_%s" % (production,decay,energy)


class Kappas(LHCHCGBaseModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self,resolved=True,BRU=True):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
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
        self.modelBuilder.doVar("kappa_b[1,0.0,3.0]")
        if not self.resolved:
            self.modelBuilder.doVar("kappa_g[1,0.0,2.0]")
            self.modelBuilder.doVar("kappa_gam[1,0.0,2.5]")
	self.modelBuilder.doVar("BRinv[0,0,1]")
        self.modelBuilder.out.var("BRinv").setConstant(True)
        pois = 'kappa_W,kappa_Z,kappa_tau,kappa_t,kappa_b'
        if not self.resolved:
            pois += ',kappa_g,kappa_gam'
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
        self.SMH.makeScaling('qqH', CW='kappa_W', CZ='kappa_Z')
        self.SMH.makeScaling("tHq", CW='kappa_W', Ctop="kappa_t")
        self.SMH.makeScaling("tHW", CW='kappa_W', Ctop="kappa_t")
        self.SMH.makeScaling("ggZH", CZ='kappa_Z', Ctop="kappa_t")
        # resolve loops
        if self.resolved:
            self.SMH.makeScaling('ggH', Cb='kappa_b', Ctop='kappa_t')
            self.SMH.makeScaling('hgluglu', Cb='kappa_b', Ctop='kappa_t')
            self.SMH.makeScaling('hgg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_W', Ctau='kappa_tau')
            self.SMH.makeScaling('hzg', Cb='kappa_b', Ctop='kappa_t', CW='kappa_W', Ctau='kappa_tau')
        else:
            self.modelBuilder.factory_('expr::Scaling_hgluglu("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_hgg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_hzg("@0*@0", kappa_gam)')
            self.modelBuilder.factory_('expr::Scaling_ggH_7TeV("@0*@0", kappa_g)')
            self.modelBuilder.factory_('expr::Scaling_ggH_8TeV("@0*@0", kappa_g)')

        ## partial witdhs, normalized to the SM one
        self.modelBuilder.factory_('expr::c7_Gscal_Z("@0*@0*@1*@2", kappa_Z, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)')
        self.modelBuilder.factory_('expr::c7_Gscal_W("@0*@0*@1*@2", kappa_W, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)')
        self.modelBuilder.factory_('expr::c7_Gscal_tau("@0*@0*@1*@4+@2*@2*@3*@5", kappa_tau, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)')
        self.modelBuilder.factory_('expr::c7_Gscal_top("@0*@0 * @1*@2", kappa_t, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_Gscal_bottom("@0*@0 * (@1*@3+@2)", kappa_b, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)')
        self.modelBuilder.factory_('expr::c7_Gscal_gluon("  @0  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)')
        self.modelBuilder.factory_('expr::c7_Gscal_gamma("@0*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)')
        # fix to have all BRs add up to unity
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
        self.modelBuilder.factory_('expr::c7_BRscal_hcc("@0*@0*@2/@1", kappa_t, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgg("@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hzg("@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgluglu("@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)')

    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c7_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name) == None:
            if production in [ "ggH", "qqH", "ggZH", "tHq", "tHW"]: 
                XSscal = ("@0", "Scaling_%s_%s" % (production,energy) )
            elif production == "WH":  XSscal = ("@0*@0", "kappa_W")
            elif production == "ZH":  XSscal = ("@0*@0", "kappa_Z")
            elif production == "ttH": XSscal = ("@0*@0", "kappa_t")
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


class Lambdas(LHCHCGBaseModel):
    def __init__(self,BRU=True):
        LHCHCGBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.doBRU = BRU
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = (po.replace("BRU=","") in [ "yes", "1", "Yes", "True", "true" ])
        print "BR uncertainties in partial widths: %s " % self.doBRU
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.doMH()
        self.modelBuilder.doVar("lambda_WZ[1,0.0,2.0]") 
        self.modelBuilder.doVar("lambda_Zg[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_bZ[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_gamZ[1,0.0,2.0]")
        self.modelBuilder.doVar("lambda_tauZ[1,0.0,4.0]")
        self.modelBuilder.doVar("lambda_tg[1,0.0,4.0]")
        self.modelBuilder.doVar("kappa_gZ[1,0.0,3.0]")
        self.modelBuilder.doSet("POI",'lambda_WZ,lambda_Zg,lambda_bZ,lambda_gamZ,lambda_tauZ,lambda_tg,kappa_gZ')
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()
    def setup(self):
        self.dobbH()
        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()

        self.modelBuilder.doVar("lambda_one[1]")

        self.modelBuilder.factory_("expr::lambda_tZ(\"@0/@1\",lambda_tg,lambda_Zg)");
        self.SMH.makeScaling("qqH", CW='lambda_WZ', CZ="lambda_one")
        self.SMH.makeScaling("tHq", CW='lambda_WZ', Ctop="lambda_tZ")
        self.SMH.makeScaling("tHW", CW='lambda_WZ', Ctop="lambda_tZ")
        self.SMH.makeScaling("ggZH", CZ='lambda_one', Ctop="lambda_tZ")
        self.SMH.makeScaling('hzg', Cb='lambda_bZ', Ctop='lambda_tZ', CW='lambda_WZ', Ctau='lambda_tauZ')
        self.modelBuilder.factory_("expr::lambda_gZ(\"1/@0\",lambda_Zg)");
        self.modelBuilder.factory_("expr::sqrt_zgamma(\"sqrt(@0)\",Scaling_hzg)");

        for E in "7TeV", "8TeV":
            for P in "qqH", "tHq", "tHW", "ggZH":
                self.modelBuilder.factory_("expr::PW_XSscal_%s_%s(\"@0*@1*@1\",Scaling_%s_%s,lambda_Zg)"%(P,E,P,E))
            self.modelBuilder.factory_("expr::PW_XSscal_WH_%s(\"@0*@0*@1*@1\",lambda_Zg,lambda_WZ)" % E)
            self.modelBuilder.factory_("expr::PW_XSscal_ZH_%s(\"@0*@0\",lambda_Zg)" % E)
            self.modelBuilder.factory_("expr::PW_XSscal_ttH_%s(\"@0*@0\",lambda_tg)" % E)
        self.decayMap_ = {
            'hww' : 'lambda_WZ',
            'hzz' : 'lambda_one',
            'hgg' : 'lambda_gamZ',
            'hbb' : 'lambda_bZ',
            'htt' : 'lambda_tauZ',
            'hmm' : 'lambda_tauZ',
            'hcc' : 'lambda_tZ',     # charm scales as top
            'hgluglu' : 'lambda_gZ', # glu scales as 1/Zgky
            'hzg'     : 'lambda_gamZ',   # fancier option: 'sqrt_zgamma',
            #'hss' : 'lambda_bZ', # strange scales as bottom # not used
        } 
    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c7_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name):
            return name
        dscale = self.decayMap_[decay]
        if self.doBRU:
            name += "_noBRU"
        if production == "ggH": 
            if decay in self.add_bbH:
                b2g = "CMS_R_bbH_ggH_%s_%s[%g]" % (decay, energy, 0.01) 
                b2gs = "CMS_bbH_scaler_%s" % energy
                self.modelBuilder.factory_("expr::%s(\"@0*@0*@1*@1*(1+@2*@3*@4*@4*@5*@5)\",kappa_gZ,%s,%s,%s,lambda_bZ,lambda_Zg)" % (name, dscale, b2g, b2gs))
            else:
                self.modelBuilder.factory_("expr::%s(\"@0*@0*   @1*@1\",kappa_gZ,%s)" % (name, dscale))
        else:
            self.modelBuilder.factory_("expr::%s(\"@0*@0*@1*@2*@2\",kappa_gZ,PW_XSscal_%s_%s,%s)" % (name, production, energy, dscale))
        if self.doBRU:
            name = name.replace("_noBRU","")
            if decay == "hzz":
                self.modelBuilder.factory_("prod::%s(%s_noBRU, HiggsDecayWidth_UncertaintyScaling_%s)" % (name, name, "hzz"))
            else:
                self.modelBuilder.factory_("expr::%s(\"@0*(@1/@2)\", %s_noBRU, HiggsDecayWidth_UncertaintyScaling_%s, HiggsDecayWidth_UncertaintyScaling_%s)" % (name, name, decay, "hzz"))
        print '[LHC-HCG Lambdas]', name, production, decay, energy,": ",
        self.modelBuilder.out.function(name).Print("")
        return name


A1 = SignalStrengths()
A2 = SignalStrengthRatios()
B1 = XSBRratios()
K1 = Kappas(resolved=True)
K2 = Kappas(resolved=False)
L1 = Lambdas()
