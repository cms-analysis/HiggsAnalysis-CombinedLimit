from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class CommonMatrix(SMLikeHiggsModel):
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.decays = [ "hbb", "htt", "hgg", "hww", "hzz" ]
        self.productions = ["ggH", "qqH", "ZH", "WH", "ttH"]
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
        self.modelBuilder.doVar("lambda[1,0,10]")
        self.modelBuilder.doVar("lambdaw[1,0,10]")
        self.modelBuilder.doVar("lambdaz[1,0,10]")
        self.modelBuilder.doVar("lambdat[1,0,10]")
        poi.append('lambda');
        poi.append('lambdaw');
        poi.append('lambdaz');
        poi.append('lambdat');
        for decay in self.decays:
            #poi.append('mu_'+decay);
            for production in self.productions:
                if production=="qqH":
                        poi.append('lambda_'+decay);
                elif production == "ZH":
                        poi.append('lambdaz_'+decay);
                elif production == "WH":
                        poi.append('lambdaw_'+decay);
                elif production == "ttH":
                        poi.append('lambdat_'+decay);
            if decay in self.fixDecays:
                print "It seems that you set signal strength ggH->"+decay+" as missing, I will set that as a POI. Take this into account when running HybridNew."
                self.modelBuilder.doVar("mu_%s[1,0,5]" % decay)
                poi.append('mu_'+decay)
            else:
                self.modelBuilder.doVar("mu_%s[1,0,5]" % decay)
            self.modelBuilder.doVar("lambda_%s[1,0,10]" % decay)
            self.modelBuilder.doVar("lambdaz_%s[1,0,10]" % decay)
            self.modelBuilder.doVar("lambdaw_%s[1,0,10]" % decay)
            self.modelBuilder.doVar("lambdat_%s[1,0,10]" % decay)
            self.modelBuilder.factory_("expr::lambda_%smu_%slambda(\"@0*@1*@2\",lambda_%s,mu_%s,lambda)" % (decay,decay,decay,decay))
            self.modelBuilder.factory_("expr::lambdaz_%smu_%slambdaz(\"@0*@1*@2\",lambdaz_%s,mu_%s,lambdaz)" % (decay,decay,decay,decay))
            self.modelBuilder.factory_("expr::lambdaw_%smu_%slambdaw(\"@0*@1*@2\",lambdaw_%s,mu_%s,lambdaw)" % (decay,decay,decay,decay))
            self.modelBuilder.factory_("expr::lambdat_%smu_%slambdat(\"@0*@1*@2\",lambdat_%s,mu_%s,lambdat)" % (decay,decay,decay,decay))
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
        if production == "ggH" or production == "bbH":
                return 'mu_'+decay
        elif production == "qqH":
                return 'lambda_'+decay+'mu_'+decay+'lambda'
        elif production == "ZH" or production == "ggZH" or production == "qqZH":
                return 'lambdaz_'+decay+'mu_'+decay+'lambdaz'
        elif production == "WH":
                return 'lambdaw_'+decay+'mu_'+decay+'lambdaw'
        elif production == "ttH" or production == "tH":
                return 'lambdat_'+decay+'mu_'+decay+'lambdat'	
        else: raise RuntimeError, "Unknown production mode '%s'" %production

commonMatrix= CommonMatrix()

