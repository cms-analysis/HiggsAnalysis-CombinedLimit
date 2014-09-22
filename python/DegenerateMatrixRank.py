from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class AllMuiLambdaHiggs(SMLikeHiggsModel):
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.decays = [ "hbb", "htt", "hgg", "hww", "hzz" ]
        self.productions = ["ggH", "qqH", "VH", "ttH"]
        self.fixDecays = []
        self.mHRange = []
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("decays="): self.decays = po.replace("decays=","").split(",")
            if po.startswith("productions="): self.productions = po.replace("productions=","").split(",")
            if po.startswith("fixDecays="): self.fixDecays = po.replace("fixDecays=","").split(",") 
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI ---
        #self.modelBuilder.doVar("Rvf[1,-5,20]")
	self.modelBuilder.doVar("lambda[1,0,10]")
	self.modelBuilder.doVar("lambdav[1,0,10]")
	self.modelBuilder.doVar("lambdat[1,0,10]")
        # --- put a higgs mass, just because the code might need it
        #self.modelBuilder.doVar("MH[125]")
        poi=[]
	poi.append('lambda');
	poi.append('lambdav');
	poi.append('lambdat');
	for decay in self.decays:
            #poi.append('mu_'+decay);
            if decay in self.fixDecays:
                print "It seems that you set signal strength ggH->"+decay+" as missing, I will set that as POI. Please set it to unity when running HybridNew."
                self.modelBuilder.doVar("mu_%s[1,0,5]" % decay)
                poi.append('mu_'+decay)
            else:
                self.modelBuilder.doVar("mu_%s[1,0,5]" % decay)
            self.modelBuilder.factory_("expr::lambdamu_%s(\"@0*@1\",lambda,mu_%s)" % (decay,decay))
            self.modelBuilder.factory_("expr::lambdavmu_%s(\"@0*@1\",lambdav,mu_%s)" % (decay,decay))
            self.modelBuilder.factory_("expr::lambdatmu_%s(\"@0*@1\",lambdat,mu_%s)" % (decay,decay))
        self.modelBuilder.doSet("POI", ",".join(poi) )
        if self.modelBuilder.out.var("MH"):
                if len(self.mHRange):
                        print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                        self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                        self.modelBuilder.out.var("MH").setConstant(False)
                        poi.append('MH');
                        self.modelBuilder.doSet("POI", ",".join(poi) )
                        #self.modelBuilder.doSet('POI','poi,MH')
                else:
                     	print 'MH will be assumed to be', self.options.mass
                        self.modelBuilder.out.var("MH").removeRange()
                        self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
             	if len(self.mHRange):
                        print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                        self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                        poi.append('MH');
                        self.modelBuilder.doSet("POI", ",".join(poi) )
                        #self.modelBuilder.doSet('POI','poi,MH')
                else:
                     	print 'MH (not there before) will be assumed to be', self.options.mass
                        self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        for item in poi:
           print item
        #the next part is for the special case, where the ratio (look for your papers) is constant
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if decay=="hmm":
        	decay="htt"
	if decay=="hzg":
		decay="hgg"
	if decay in [ "hgluglu", "hcc", "hss", "huu", "hdd" ]:
		decay="hbb"
        if production == "ggH":
                return 'mu_'+decay
        elif production == "qqH":
                return 'lambdamu_'+decay
        elif production == "VH" or production=="ZH" or production=="WH":
                return 'lambdavmu_'+decay
        elif production == "ttH":
                return 'lambdatmu_'+decay
        else: raise RuntimeError, "Unknown production mode '%s'" %production



class AllMuiLambdasHiggs(SMLikeHiggsModel):
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.decays = [ "hbb", "htt", "hgg", "hww", "hzz" ]
        self.productions = ["ggH", "qqH", "VH", "ttH"]
        self.fixDecays = []
	self.mHRange = []
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("decays="): self.decays = po.replace("decays=","").split(",")
            if po.startswith("productions="): self.productions = po.replace("productions=","").split(",")
            if po.startswith("fixDecays="): self.fixDecays = po.replace("fixDecays=","").split(",")
	    if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI ---
        #self.modelBuilder.doVar("Rvf[1,-5,20]")
        # --- put a higgs mass, just because the code might need it
        #self.modelBuilder.doVar("MH[125]")
        poi=[]
        for decay in self.decays:
            #poi.append('mu_'+decay);
            for production in self.productions:
                if production=="qqH":
                        poi.append('lambda_'+decay);
                elif production == "VH":
                        poi.append('lambdav_'+decay);
                elif production == "ttH":
                        poi.append('lambdat_'+decay);
            if decay in self.fixDecays:
                print "It seems that you set signal strength ggH->"+decay+" as missing, I will set that equal to unity."
                self.modelBuilder.doVar("mu_%s[1.0]" % decay)
            else:
                self.modelBuilder.doVar("mu_%s[1,0,5]" % decay)
	    self.modelBuilder.doVar("lambda_%s[1,0,10]" % decay)
	    self.modelBuilder.doVar("lambdav_%s[1,0,10]" % decay)
            self.modelBuilder.doVar("lambdat_%s[1,0,10]" % decay)
            self.modelBuilder.factory_("expr::lambda_%smu_%s(\"@0*@1\",lambda_%s,mu_%s)" % (decay,decay,decay,decay))
            self.modelBuilder.factory_("expr::lambdav_%smu_%s(\"@0*@1\",lambdav_%s,mu_%s)" % (decay,decay,decay,decay))
            self.modelBuilder.factory_("expr::lambdat_%smu_%s(\"@0*@1\",lambdat_%s,mu_%s)" % (decay,decay,decay,decay))
        self.modelBuilder.doSet("POI", ",".join(poi) )
	if self.modelBuilder.out.var("MH"):
                if len(self.mHRange):
                        print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                        self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                        self.modelBuilder.out.var("MH").setConstant(False)
			poi.append('MH');
			self.modelBuilder.doSet("POI", ",".join(poi) )
                        #self.modelBuilder.doSet('POI','poi,MH')
                else:
                     	print 'MH will be assumed to be', self.options.mass
                        self.modelBuilder.out.var("MH").removeRange()
                        self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
             	if len(self.mHRange):
                        print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                        self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
			poi.append('MH');
                        self.modelBuilder.doSet("POI", ",".join(poi) )
                        #self.modelBuilder.doSet('POI','poi,MH')
                else:
                     	print 'MH (not there before) will be assumed to be', self.options.mass
                        self.modelBuilder.doVar("MH[%g]" % self.options.mass)
	for item in poi:
           print item
    #the next part is for the special case, where the ratio (look for your papers) is constant
    def getHiggsSignalYieldScale(self,production,decay, energy):
	#print 'Printing out:', production, decay, energy, 'mu_'+decay
	#return 'mu_'+decay
	if decay=="hmm":
                decay="htt"
        if decay=="hzg":
                decay="hgg"
        if decay in [ "hgluglu", "hcc", "hss", "huu", "hdd" ]:
                decay="hbb"
	if production == "ggH":
                return 'mu_'+decay
        elif production == "qqH":
                return 'lambda_'+decay+'mu_'+decay
        elif production == "VH" or production=="ZH" or production=="WH":
                return 'lambdav_'+decay+'mu_'+decay
        elif production == "ttH":
		return 'lambdat_'+decay+'mu_'+decay	
        else: raise RuntimeError, "Unknown production mode '%s'" %production

allmuilambda= AllMuiLambdaHiggs()
allmuilambdas = AllMuiLambdasHiggs()

