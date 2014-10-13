from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import ROOT, os
##  think about BRUndet  and BRscal_hcc/gg/ss  
class C8(SMLikeHiggsModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
        self.doHZg = False
        self.doHInv = False
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po == 'doHZg':
                self.doHZg = True
            if po == 'doHInv':
                self.doHInv = True
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kV[1,0.0,1.0]") # bounded to 1
        self.modelBuilder.doVar("ktau[1,0.0,2.0]")
        self.modelBuilder.doVar("ktop[1,0.0,4.0]")
        self.modelBuilder.doVar("kbottom[1,0.0,3.0]")
        self.modelBuilder.doVar("kgluon[1,0.0,2.0]")
        self.modelBuilder.doVar("kgamma[1,0.0,2.5]")
	self.modelBuilder.doVar("BRInv[0,0,1]")
	self.modelBuilder.doVar("BRUndet[0,0,1]")
        pois = 'kV,ktau,ktop,kbottom,kgluon,kgamma,BRInv,BRUndet' 
        if self.doHZg:
            self.modelBuilder.doVar("kZgamma[1,0.0,30.0]")
            pois += ",kZgamma"
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
            self.modelBuilder.doSet("POI",pois+',MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
            self.modelBuilder.doSet("POI",pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):

        # SM BR
        for d in [ "htt", "hbb", "hcc", "hww", "hzz", "hgluglu", "htoptop", "hgg", "hzg", "hmm", "hss" ]: self.SMH.makeBR(d)

        ## total witdhs, normalized to the SM one
        self.modelBuilder.factory_('expr::c8_Gscal_Vectors("@0*@0 * (@1+@2)", kV, SM_BR_hzz, SM_BR_hww)')
        self.modelBuilder.factory_('expr::c8_Gscal_tau("@0*@0 * (@1+@2)", ktau, SM_BR_htt, SM_BR_hmm)')
        self.modelBuilder.factory_('expr::c8_Gscal_top("@0*@0 * (@1+@2)", ktop, SM_BR_htoptop, SM_BR_hcc)')
        self.modelBuilder.factory_('expr::c8_Gscal_bottom("@0*@0 * (@1+@2)", kbottom, SM_BR_hbb, SM_BR_hss)')
        self.modelBuilder.factory_('expr::c8_Gscal_gluon("@0*@0 * @1", kgluon, SM_BR_hgluglu)')
        if not self.doHZg:
            self.modelBuilder.factory_('expr::c8_Gscal_gamma("@0*@0 * (@1+@2)", kgamma, SM_BR_hgg, SM_BR_hzg)')
        else:
            self.modelBuilder.factory_('expr::c8_Gscal_gamma("@0*@0 *@1+@2*@2*@3", kgamma, SM_BR_hgg, kZgamma, SM_BR_hzg)')


        self.modelBuilder.factory_('expr::c8_Gscal_tot("(@1+@2+@3+@4+@5+@6)/(1-@0-@7)",BRInv,c8_Gscal_Vectors, c8_Gscal_tau, c8_Gscal_top, c8_Gscal_bottom, c8_Gscal_gluon, c8_Gscal_gamma,BRUndet)')

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM)^2 / (total/total_SM)^2 
        self.modelBuilder.factory_('expr::c8_BRscal_hvv("@0*@0/@1", kV, c8_Gscal_tot)')
        self.modelBuilder.factory_('expr::c8_BRscal_htt("@0*@0/@1", ktau, c8_Gscal_tot)')
        self.modelBuilder.factory_('expr::c8_BRscal_hbb("@0*@0/@1", kbottom, c8_Gscal_tot)')
        self.modelBuilder.factory_('expr::c8_BRscal_hgg("@0*@0/@1", kgamma, c8_Gscal_tot)')
        if self.doHZg:
            self.modelBuilder.factory_('expr::c8_BRscal_hzg("@0*@0/@1", kZgamma, c8_Gscal_tot)')
        if self.doHInv:
            self.modelBuilder.factory_('expr::c8_BRscal_hinv("@0>=0?@0:-100", BRInv)')



    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c8_XSBRscal_%s_%s" % (production,decay)
        print '[LOFullParametrization::C7]'
        print name, production, decay, energy
        if self.modelBuilder.out.function(name) == None:
            XSscal = "kgluon"
            if production in ["WH","ZH","VH","qqH"]: XSscal = "kV"
            if production == "ttH": XSscal = "ktop"
            BRscal = "hgg"
            if decay in ["hbb", "htt"]: BRscal = decay
            if decay in ["hww", "hzz"]: BRscal = "hvv"
            if self.doHZg and decay == "hzg": BRscal = "hzg"
            if self.doHInv and decay == "hinv": BRscal = "hinv"
            self.modelBuilder.factory_('expr::%s("@0*@0 * @1", %s, c8_BRscal_%s)' % (name, XSscal, BRscal))
        return name


class CWidth(SMLikeHiggsModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
        self.doHZg = False
        self.doHInv = False
        self.doWidth = False
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po == 'doHZg':
                self.doHZg = True
            if po == 'doHInv':
                self.doHInv = True
            if po == 'doWidth':
                self.doWidth = True
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kV[1,0.0,1.0]") # bounded to 1
        self.modelBuilder.doVar("ktau[1,0.0,2.0]")
        self.modelBuilder.doVar("ktop[1,0.0,4.0]")
        self.modelBuilder.doVar("kbottom[1,0.0,3.0]")
        self.modelBuilder.doVar("kgluon[1,0.0,2.0]")
        self.modelBuilder.doVar("kgamma[1,0.0,2.5]")
        if self.doWidth:
		self.modelBuilder.doVar("Gscale[1,0,10]")
            	pois = 'kV,ktau,ktop,kbottom,kgluon,kgamma,Gscale'
	else:
		self.modelBuilder.doVar("BRInvUndet[0,0,1]")
            	pois = 'kV,ktau,ktop,kbottom,kgluon,kgamma,BRInvUndet'
        if self.doHZg:
            self.modelBuilder.doVar("kZgamma[1,0.0,30.0]")
            pois += ",kZgamma"
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
            self.modelBuilder.doSet("POI",pois+',MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
            self.modelBuilder.doSet("POI",pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):

        # SM BR
        for d in [ "htt", "hbb", "hcc", "hww", "hzz", "hgluglu", "htoptop", "hgg", "hzg", "hmm", "hss" ]: self.SMH.makeBR(d)

        ## total witdhs, normalized to the SM one
        self.modelBuilder.factory_('expr::c7_Gscal_Vectors("@0*@0 * (@1+@2)", kV, SM_BR_hzz, SM_BR_hww)')
        self.modelBuilder.factory_('expr::c7_Gscal_tau("@0*@0 * (@1+@2)", ktau, SM_BR_htt, SM_BR_hmm)')
        self.modelBuilder.factory_('expr::c7_Gscal_top("@0*@0 * (@1+@2)", ktop, SM_BR_htoptop, SM_BR_hcc)')
        self.modelBuilder.factory_('expr::c7_Gscal_bottom("@0*@0 * (@1+@2)", kbottom, SM_BR_hbb, SM_BR_hss)')
        self.modelBuilder.factory_('expr::c7_Gscal_gluon("@0*@0 * @1", kgluon, SM_BR_hgluglu)')
        if not self.doHZg:
            self.modelBuilder.factory_('expr::c7_Gscal_gamma("@0*@0 * (@1+@2)", kgamma, SM_BR_hgg, SM_BR_hzg)')
        else:
            self.modelBuilder.factory_('expr::c7_Gscal_gamma("@0*@0 *@1+@2*@2*@3", kgamma, SM_BR_hgg, kZgamma, SM_BR_hzg)')

	# 
#RooFormulaVar::c7_Gscal_tot[ actualVars=(Gscale,c7_Gscal_Vectors,c7_Gscal_tau,c7_Gscal_top,c7_Gscal_bottom,c7_Gscal_gluon,c7_Gscal_gamma) formula="(@1+@2+@3+@4+@5+@6)<=@0?@0:0.001" ] = 0.001
#root [8] w->function("c7_Gscal_tau")->Print()
#RooFormulaVar::c7_Gscal_tau[ actualVars=(ktau,SM_BR_htt,SM_BR_hmm) formula="@0*@0*(@1+@2)" ] = 0.063419
#root [9] w->function("c7_Gscal_top")->Print()
#RooFormulaVar::c7_Gscal_top[ actualVars=(ktop,SM_BR_htoptop,SM_BR_hcc) formula="@0*@0*(@1+@2)" ] = 0.0291
#root [10] w->function("c7_Gscal_bottom")->Print()
#RooFormulaVar::c7_Gscal_bottom[ actualVars=(kbottom,SM_BR_hbb,SM_BR_hss) formula="@0*@0*(@1+@2)" ] = 0.577246    #### The spline interpolation has small difference w.r.t YR3 number 5.769462E-01
#root [11] w->function("c7_Gscal_gluon")->Print()
#RooFormulaVar::c7_Gscal_gluon[ actualVars=(kgluon,SM_BR_hgluglu) formula="@0*@0*@1" ] = 0.0857
#root [12] w->function("c7_Gscal_gamma")->Print()
#RooFormulaVar::c7_Gscal_gamma[ actualVars=(kgamma,SM_BR_hgg,SM_BR_hzg) formula="@0*@0*(@1+@2)" ] = 0.00382
#root [13] w->function("c7_Gscal_Vectors")->Print()
#RooFormulaVar::c7_Gscal_Vectors[ actualVars=(kV,SM_BR_hzz,SM_BR_hww) formula="@0*@0*(@1+@2)" ] = 0.2414
#############################Then the above sum = 1.000685 , so that the asimov data set will be generated with  c7_Gscal_tot = 0.0001, and BRInv = -100 
#root [14] w->var("ktau")->Print()
#RooRealVar::ktau = 1  L(0 - 2) 
#root [15] w->var("ktop")->Print()
#RooRealVar::ktop = 1  L(0 - 4) 
#root [16] w->var("kgluon")->Print()
#RooRealVar::kgluon = 1  L(0 - 2) 
#root [17] w->var("kgamma")->Print()
#RooRealVar::kgamma = 1  L(0 - 2.5) 
############################either generate asimov dataset with kbottom slightly less 1    i.e.  kb=0.999 is fine


	if self.doWidth:
            self.modelBuilder.factory_('expr::c7_Gscal_tot("(@1+@2+@3+@4+@5+@6)<=@0?@0:0.001",Gscale,c7_Gscal_Vectors, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma)')
            self.modelBuilder.factory_('expr::BRInvUndet("1 - (@1+@2+@3+@4+@5+@6)/@0",c7_Gscal_tot,c7_Gscal_Vectors, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma)')
	else:
            self.modelBuilder.factory_('expr::c7_Gscal_tot("(@1+@2+@3+@4+@5+@6)/(1-@0)",BRInvUndet,c7_Gscal_Vectors, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma)')

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM)^2 / (total/total_SM)^2
        self.modelBuilder.factory_('expr::c7_BRscal_hvv("@0*@0/@1", kV, c7_Gscal_tot)')
        self.modelBuilder.factory_('expr::c7_BRscal_htt("@0*@0/@1", ktau, c7_Gscal_tot)')
        self.modelBuilder.factory_('expr::c7_BRscal_hbb("@0*@0/@1", kbottom, c7_Gscal_tot)')
        self.modelBuilder.factory_('expr::c7_BRscal_hgg("@0*@0/@1", kgamma, c7_Gscal_tot)')
        if self.doHZg:
            self.modelBuilder.factory_('expr::c7_BRscal_hzg("@0*@0/@1", kZgamma, c7_Gscal_tot)')
        if self.doHInv:
            self.modelBuilder.factory_('expr::c7_BRscal_hinv("@0>=0?@0:-100", BRInvUndet)')



    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "c7_XSBRscal_%s_%s" % (production,decay)
        print '[LOFullParametrization::C7]'
        print name, production, decay, energy
        if self.modelBuilder.out.function(name) == None:
            XSscal = "kgluon"
            if production in ["WH","ZH","VH","qqH"]: XSscal = "kV"
            if production == "ttH": XSscal = "ktop"
            BRscal = "hgg"
            if decay in ["hbb", "htt"]: BRscal = decay
            if decay in ["hww", "hzz"]: BRscal = "hvv"
            if self.doHZg and decay == "hzg": BRscal = "hzg"
            if self.doHInv and decay == "hinv": BRscal = "hinv"
            self.modelBuilder.factory_('expr::%s("@0*@0 * @1", %s, c7_BRscal_%s)' % (name, XSscal, BRscal))
        return name


class TwoHDM(SMLikeHiggsModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
	self.thdmtype = ['1']  
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("thdmtype="):
                self.thdmtype= po.replace("thdmtype=","")
                if len(self.thdmtype) != 1:
                    raise RuntimeError, "2HDM type requires one value"
                elif int(self.thdmtype[0]) != 1 and int(self.thdmtype[0]) !=2 and int(self.thdmtype[0]) !=3 and int(self.thdmtype[0]) !=4:
                    raise RuntimeError, "2HDM type must be 1 (default) or 2 or 3 or 4 "
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("cosbma[0,-1,1]")
        self.modelBuilder.doVar("tanbeta[0,0.0,10]")
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",'cosbma,tanbeta,MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",'cosbma,tanbeta')
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):

	self.modelBuilder.factory_('expr::kV("sqrt(1-@0*@0)",cosbma)')
	self.modelBuilder.factory_('expr::tana("(@0*@1-@2)/(@1-@0*@2)", tanbeta, cosbma, kV)')
	self.modelBuilder.factory_('expr::cosa("1/sqrt(1+@0*@0)",tana)')
	self.modelBuilder.factory_('expr::sinb("tanbeta/sqrt(1+@0*@0)",tanbeta)')
	self.modelBuilder.factory_('expr::ku("@0/@1", cosa, sinb)')
	if int(self.thdmtype[0]) == 1: 
		self.modelBuilder.factory_('expr::kd("@0", ku)')
		self.modelBuilder.factory_('expr::kl("@0", ku)')
	elif int(self.thdmtype[0]) == 2: 
		self.modelBuilder.factory_('expr::cosb("1/sqrt(1+@0*@0)",tanbeta)')
		self.modelBuilder.factory_('expr::sina("tana/sqrt(1+@0*@0)",tana)')
		self.modelBuilder.factory_('expr::kd("-@0/@1", sina,cosb)')
		self.modelBuilder.factory_('expr::kl("@0", kd)')
	elif int(self.thdmtype[0]) == 3: 
		self.modelBuilder.factory_('expr::cosb("1/sqrt(1+@0*@0)",tanbeta)')
		self.modelBuilder.factory_('expr::sina("tana/sqrt(1+@0*@0)",tana)')
		self.modelBuilder.factory_('expr::kd("@0", ku)')
		self.modelBuilder.factory_('expr::kl("-@0/@1", sina,cosb)')
	elif int(self.thdmtype[0]) == 4: 
		self.modelBuilder.factory_('expr::cosb("1/sqrt(1+@0*@0)",tanbeta)')
		self.modelBuilder.factory_('expr::sina("tana/sqrt(1+@0*@0)",tana)')
		self.modelBuilder.factory_('expr::kd("-@0/@1", sina,cosb)')
		self.modelBuilder.factory_('expr::kl("@0", ku)')
        self.decayScaling = {
            'hgg':'hgg',
            'hzg':'hzg',
            'hww':'hvv',
            'hzz':'hvv',
            'hbb':'hdd',
            'htt':'hll',
            'hss':'hdd',
            'hmm':'hll',
            'hcc':'huu',
            'hgluglu':'hgluglu',
            }
        self.productionScaling = {
            'ttH':'ku',
            'qqH':'kV',
            'WH':'kV',
            'ZH':'kV',
            'VH':'kV',
            }
        
        # scalings of the loops
        self.SMH.makeScaling('ggH', Cb='kd', Ctop='ku')
        self.SMH.makeScaling('hgg', Cb='kd', Ctop='ku', CW='kV', Ctau='kl')
        self.SMH.makeScaling('hzg', Cb='kd', Ctop='ku', CW='kV', Ctau='kl')
        self.SMH.makeScaling('hgluglu', Cb='kd', Ctop='ku')

        # SM BR
        for d in [ "htt", "hbb", "hcc", "hww", "hzz", "hgluglu", "htoptop", "hgg", "hzg", "hmm", "hss" ]:
            self.SMH.makeBR(d)

        ## total witdhs, normalized to the SM one
        self.modelBuilder.factory_('expr::twohdm_Gscal_Vectors("@0*@0 * (@1+@2)", kV, SM_BR_hzz, SM_BR_hww)') 
        self.modelBuilder.factory_('expr::twohdm_Gscal_up("@0*@0 * (@1+@2)", ku, SM_BR_hcc, SM_BR_htoptop)') 
        self.modelBuilder.factory_('expr::twohdm_Gscal_down("@0*@0 * (@1+@2)", kd, SM_BR_hbb, SM_BR_hss)')
        self.modelBuilder.factory_('expr::twohdm_Gscal_leptons("@0*@0 * (@1+@2)", kl, SM_BR_htt, SM_BR_hmm)')
        self.modelBuilder.factory_('expr::twohdm_Gscal_gg("@0 * @1", Scaling_hgg, SM_BR_hgg)') 
        self.modelBuilder.factory_('expr::twohdm_Gscal_Zg("@0 * @1", Scaling_hzg, SM_BR_hzg)')
        self.modelBuilder.factory_('expr::twohdm_Gscal_gluglu("@0 * @1", Scaling_hgluglu, SM_BR_hgluglu)')
        self.modelBuilder.factory_('sum::twohdm_Gscal_tot(twohdm_Gscal_Vectors, twohdm_Gscal_up, twohdm_Gscal_down, twohdm_Gscal_leptons, twohdm_Gscal_gg, twohdm_Gscal_Zg, twohdm_Gscal_gluglu)')

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM)^2 / (total/total_SM)^2 
        self.modelBuilder.factory_('expr::twohdm_BRscal_hvv("@0*@0/@1", kV, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_huu("@0*@0/@1", ku, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_hdd("@0*@0/@1", kd, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_hll("@0*@0/@1", kl, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_hgg("@0/@1", Scaling_hgg, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_hgluglu("@0/@1", Scaling_hgluglu, twohdm_Gscal_tot)')
        self.modelBuilder.factory_('expr::twohdm_BRscal_hzg("@0/@1", Scaling_hzg, twohdm_Gscal_tot)')

        # verbosity
        #self.modelBuilder.out.Print()

    def getHiggsSignalYieldScale(self,production,decay,energy):

        name = 'twohdm_XSBRscal_%(production)s_%(decay)s' % locals()

        #Special case that depends on Energy
        if production == 'ggH':
            self.productionScaling[production] = 'Scaling_ggH_' + energy
            name += '_%(energy)s' % locals()
            
        if self.modelBuilder.out.function(name):
            return name

        XSscal = self.productionScaling[production]
        BRscal = self.decayScaling[decay]
        if 'Scaling_' in XSscal: # it's a Scaling, which means it's already squared
            self.modelBuilder.factory_('expr::%s("@0 * @1", %s, twohdm_BRscal_%s)' % (name, XSscal, BRscal))
        else: # It's a kappa, so it's linear and I must square it
            self.modelBuilder.factory_('expr::%s("@0*@0 * @1", %s, twohdm_BRscal_%s)' % (name, XSscal, BRscal))
        return name


c8 = C8()
cwidth = CWidth()
twohdm = TwoHDM()
