from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import fnmatch 

ALL_STXS_PROCS = {
	"Stage0": {
	 "ggH*": 	  "r_ggH"
	,"bbH*": 	  "r_ggH"
	,"qqH*": 	  "r_qqH"
	,"ttH*":    	  "r_ttH"
	,"tH[Wq]*":    	  "r_ttH"
	,"ZH_lep*": 	  "r_QQ2HLL"
	,"WH_lep*": 	  "r_QQ2HLNU"
	,"[VWZ]H_had*":   "r_VH2HQQ"
	}
}

def getSTXSProdDecMode(bin,process,options):
    """Return a triple of (production)"""
    processSource = process
    decaySource   = options.fileName+":"+bin # by default, decay comes from the datacard name or bin label
    if "_" in process: 
        (processSource, decaySource) = "_".join(process.split("_")[0:-1]),process.split("_")[-1]
    foundDecay = None
    for D in ALL_HIGGS_DECAYS:
        if D in decaySource:
            if foundDecay: raise RuntimeError, "Validation Error: decay string %s contains multiple known decay names" % decaySource
            foundDecay = D
    if not foundDecay: raise RuntimeError, "Validation Error: decay string %s does not contain any known decay name" % decaySource
    #
    foundEnergy = None
    for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
        if D in decaySource:
            if foundEnergy: raise RuntimeError, "Validation Error: decay string %s contains multiple known energies" % decaySource
            foundEnergy = D
    if not foundEnergy:
        for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
            if D in options.fileName+":"+bin:
                if foundEnergy: raise RuntimeError, "Validation Error: decay string %s contains multiple known energies" % decaySource
                foundEnergy = D
    if not foundEnergy:
        eoundEnergy = '13TeV' ## if using 81x, chances are its 13 TeV
        print "Warning: decay string %s does not contain any known energy, assuming %s" % (decaySource, foundEnergy)
    #
    return (processSource, foundDecay, foundEnergy)


class STXSBaseModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False
    def preProcessNuisances(self,nuisances):
    	# add here any pre-processed nuisances such as constraint terms for the mass profiling?
    	return 
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
	    print po
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doMH(self):
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
	    self.POIs+=",MH"
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        if not self.DC.isSignal[process]: return 1
        (processSource, foundDecay, foundEnergy) = getSTXSProdDecMode(bin,process,self.options)
        return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)


class StageZero(STXSBaseModel):
    "Allow different signal strength fits for the stage-0 model"
    def __init__(self):
        STXSBaseModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.POIs = ""
    def doVar(self,x,constant=True):
        self.modelBuilder.doVar(x)
        vname = re.sub(r"\[.*","",x)
        #self.modelBuilder.out.var(vname).setConstant(constant)
        print "SignalStrengths:: declaring %s as %s" % (vname,x)
    def doParametersOfInterest(self):
        """Create POI out of signal strengths (and MH)"""
	pois = []
        for X in [ "qqH", "ggH", "ttH", "QQ2HLNU", "QQ2HLL", "VH2HQQ"]:
            self.doVar("r_%s[1,0,10]" % X)
	    pois.append("r_%s"%X)
	self.POIs=",".join(pois)
        self.doMH()
        print "Default parameters of interest: ", self.POIs
        self.modelBuilder.doSet("POI",self.POIs)
        #self.SMH = SMHiggsBuilder(self.modelBuilder)
        #self.setup()
    def setup(self):
	return 

    def getHiggsSignalYieldScale(self,production,decay,energy):
        name = "%s_%s_%s" % (production,decay,energy)
	for regproc in ALL_STXS_PROCS["Stage0"].keys():
	    if	fnmatch.fnmatch(production, regproc): 
	    	retpoi = "%s"%(ALL_STXS_PROCS["Stage0"][regproc])
	    	print "Will scale %s with POI %s"%(name,retpoi)
	    	return retpoi 
        raise RuntimeError, "No production process matching %s for Stage0 found !"%production

stage0 = StageZero()
