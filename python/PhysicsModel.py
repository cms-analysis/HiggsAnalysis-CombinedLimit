from abc import ABCMeta, abstractmethod
import re

### Class that takes care of building a physics model by combining individual channels and processes together
### Things that it can do:
###   - define the parameters of interest (in the default implementation , "r")
###   - define other constant model parameters (e.g., "MH")
###   - yields a scaling factor for each pair of bin and process (by default, constant for background and linear in "r" for signal)
###   - possibly modifies the systematical uncertainties (does nothing by default)

class PhysicsModelBase(object):
    __metaclass__ = ABCMeta
    def __init__(self):
        pass
    def setModelBuilder(self, modelBuilder):
        "Connect to the ModelBuilder to get workspace, datacard and options. Should not be overloaded."
        self.modelBuilder = modelBuilder
        self.DC = modelBuilder.DC
        self.options = modelBuilder.options
    def setPhysicsOptions(self,physOptions):
        "Receive a list of strings with the physics options from command line"
    @abstractmethod
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
    def preProcessNuisances(self,nuisances):
        "receive the usual list of (name,nofloat,pdf,args,errline) to be edited"
        pass # do nothing by default
    def getYieldScale(self,bin,process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        return "r" if self.DC.isSignal[process] else 1;
    def getChannelMask(self, bin):
        "Return the name of a RooAbsReal to mask the given bin (args != 0 => masked)"
        name = 'mask_%s' % bin
        # Check that the mask expression does't exist already, it might do
        # if it was already defined in the datacard
        if not self.modelBuilder.out.arg(name):
            self.modelBuilder.doVar('%s[0]' % name)
        return name
    def done(self):
        "Called after creating the model, except for the ModelConfigs"
        pass

class PhysicsModelBase_NiceSubclasses(PhysicsModelBase):
  """Subclass this so that subclasses work nicer"""
  def doParametersOfInterest(self):
    """
    do not override this if you want subclasses to work nicely.
    put everything that would have gone here into getPOIList instead.
    """
    self.modelBuilder.doSet("POI", ",".join(self.getPOIList()))

  @abstractmethod
  def getPOIList(self):
    """
    Create POI and other parameters, and return a list of POI variable names.
    Make sure to include:
      pois += super([classname], self).getPOIList()
    !!!!
    """
    return []

  def setPhysicsOptions(self,physOptions):
    """
    Better error checking: instead of overriding this one, override processPhysicsOptions.
    It should remove each physicsOption from the list after processing it.
    """
    processed = self.processPhysicsOptions(physOptions)
    processed = set(processed) #remove duplicates
    for _ in processed:
        physOptions.remove(_)
    if physOptions:
        raise ValueError("Unknown physicsOptions:\n{}".format(physOptions))

  @abstractmethod
  def processPhysicsOptions(self, physOptions):
    """
    Process physics options.  Make sure to return a list of physicsOptions processed,
    and to include:
      processed += super([classname], self).processPhysicsOptions(physOptions)
    !!!!
    """
    return []

class PhysicsModel(PhysicsModelBase):
    """Example class with signal strength as only POI"""
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("r[1,0,20]");
        self.modelBuilder.doSet("POI","r")
        # --- Higgs Mass as other parameter ----
        if self.options.mass != 0:
            if self.modelBuilder.out.var("MH"):
              self.modelBuilder.out.var("MH").removeRange()
              self.modelBuilder.out.var("MH").setVal(self.options.mass)
            else:
              self.modelBuilder.doVar("MH[%g]" % self.options.mass);

class MultiSignalModelBase(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        self.poiMap  = []
        self.pois    = {}
        self.verbose = False
        self.factories = []
        super(MultiSignalModelBase, self).__init__()

    def setPhysicsOptions(self,physOptions):
        for po in physOptions[:]:
            if po.startswith("turnoff="):   #shorthand:  turnoff=process1,process2,process3 --> map=.*/(process1|process2|process3):0
                physOptions.remove(po)
                turnoff = po.replace("turnoff=", "").split(",")
                physOptions.append("map=.*/({}):0".format("|".join(turnoff)))
        super(MultiSignalModelBase, self).setPhysicsOptions(physOptions)

    def processPhysicsOptions(self,physOptions):
        processed = []
        physOptions.sort(key=lambda x: x.startswith("verbose"), reverse=True) #put verbose at the beginning
        for po in physOptions:
            if po == "verbose":
                self.verbose = True
                processed.append(po)
            if po.startswith("map="):
                (maplist,poi) = po.replace("map=","").split(":",1)
                maps = maplist.split(",")
                poiname = re.sub("\[.*","", poi)
                if poi == "super":
                    pass
                elif "=" in poi:
                    poiname,expr = poi.split("=")
                    poi = expr.replace(";",":")
                    if self.verbose: print "Will create expression ",poiname," with factory ",poi
                    self.factories.append(poi)
                elif poiname not in self.pois and poi not in [ "1", "0"]: 
                    if self.verbose: print "Will create a POI ",poiname," with factory ",poi
                    self.pois[poiname] = poi
                if self.verbose:
                    if poi == "super":
                        print "Using super method to get scaling for ", maps, " patterns"
                    else:
                        print "Mapping ",poiname," to ",maps," patterns"
                self.poiMap.append((poiname, maps))
                processed.append(po)
        return processed + super(MultiSignalModelBase, self).processPhysicsOptions(physOptions)
    def getPOIList(self):
        """Create POI and other parameters, and define the POI set."""
        # --- Higgs Mass as other parameter ----
        poiNames = []
        poiNames += super(MultiSignalModelBase, self).getPOIList()
        # first do all non-factory statements, so all params are defined
        for pn,pf in self.pois.items():
            poiNames.append(pn)
            self.modelBuilder.doVar(pf)
        # then do all factory statements (so vars are already defined)
        for pf in self.factories:
            self.modelBuilder.factory_(pf)
        return poiNames

    def getYieldScale(self,bin,process):
        string = "%s/%s" % (bin,process)
        poi = None
        for p, list in self.poiMap:
            for l in list:
                if re.match(l, string): poi = p
        if poi == "super":
            poi = super(MultiSignalModelBase, self).getYieldScale(bin,process)
        if poi is None:
            poi = "1"
        print "Will scale ", string, " by ", poi
        if poi in ["1","0"]: return int(poi)
        return poi

class CanTurnOffBkgModel(PhysicsModelBase_NiceSubclasses):
    """
    Generally this should be the FIRST superclass given in any subclass
    If not it might work anyway if the other class's getYieldScale calls super properly
       but no guarantees

    If --PO nobkg is given, bkg yields will be set to 0
    """
    def __init__(self, *args, **kwargs):
        self.usebkg = True
        super(CanTurnOffBkgModel, self).__init__(*args, **kwargs)
    def processPhysicsOptions(self,physOptions):
        processed = super(CanTurnOffBkgModel, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.lower() == "nobkg":
                self.usebkg = False
                print "turning off all background"
                processed.append(po)
        return processed
    def getYieldScale(self,bin,process):
        result = super(CanTurnOffBkgModel, self).getYieldScale(bin, process)
        if not self.usebkg and not self.DC.isSignal[process]:
            print "turning off {}".format(process)
            return 0
        return result

class HiggsMassRangeModel(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        self.mHRange = []
        super(HiggsMassRangeModel, self).__init__()

    def processPhysicsOptions(self,physOptions):
        processed = super(HiggsMassRangeModel, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
                processed.append(po)
        return processed
    def getPOIList(self):
        poiNames = []
        poiNames += super(HiggsMassRangeModel, self).getPOIList()
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        return poiNames

class MultiSignalModel(MultiSignalModelBase, HiggsMassRangeModel):
    pass

### This base class implements signal yields by production and decay mode
### Specific models can be obtained redefining getHiggsSignalYieldScale
SM_HIGG_DECAYS   = [ "hww", "hzz", "hgg", "htt", "hbb", 'hzg', 'hmm', 'hcc', 'hgluglu' ]
SM_HIGG_PROD     = [ "ggH", "qqH", "VH", "WH", "ZH", "ttH", "tHq", "tHW", "ggZH", "bbH" ,"WPlusH","WMinusH"]
BSM_HIGGS_DECAYS = [ "hinv" ]
ALL_HIGGS_DECAYS = SM_HIGG_DECAYS + BSM_HIGGS_DECAYS
ALL_HIGGS_PROD   = SM_HIGG_PROD
def getHiggsProdDecMode(bin,process,options):
    """Return a triple of (production, decay, energy)"""
    processSource = process
    decaySource   = options.fileName+":"+bin # by default, decay comes from the datacard name or bin label
    if "_" in process: 
        (processSource, decaySource) = (process.split("_")[0],process.split("_")[-1]) # ignore anything in the middle for SM-like higgs
        if decaySource not in ALL_HIGGS_DECAYS:
            print "ERROR", "Validation Error in bin %r: signal process %s has a postfix %s which is not one recognized higgs decay modes (%s)" % (bin,process,decaySource,ALL_HIGGS_DECAYS)
    if processSource not in ALL_HIGGS_PROD :
        raise RuntimeError, "Validation Error in bin %r, process %r: signal process %s not among the allowed ones." % (bin,process,decaySource)
    #
    foundDecay = None
    for D in ALL_HIGGS_DECAYS:
        if D in decaySource:
            if foundDecay: raise RuntimeError, "Validation Error in bin %r, process %r: decay string %s contains multiple known decay names" % (bin,process,decaySource)
            foundDecay = D
    if not foundDecay: raise RuntimeError, "Validation Error in bin %r, process %r: decay string %s does not contain any known decay name" % (bin,process,decaySource)
    #
    foundEnergy = None
    for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
        if D in decaySource:
            if foundEnergy: raise RuntimeError, "Validation Error in bin %r, process %r: decay string %s contains multiple known energies" % (bin,process,decaySource)
            foundEnergy = D
    if not foundEnergy:
        for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
            if D in options.fileName+":"+bin:
                if foundEnergy: raise RuntimeError, "Validation Error in bin %r, process %r: decay string %s contains multiple known energies" % (bin,process,decaySource)
                foundEnergy = D
    if not foundEnergy:
        foundEnergy = '13TeV' ## if using 81x, chances are its 13 TeV
        print "Warning: decay string %s does not contain any known energy, assuming %s" % (decaySource, foundEnergy)
    if (processSource=="WPlusH" or processSource=="WMinusH"): processSource = "WH" # treat them the same for now
    return (processSource, foundDecay, foundEnergy)

class SMLikeHiggsModel(PhysicsModel):
    @abstractmethod
    def getHiggsSignalYieldScale(self, production, decay, energy):
        pass
    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        if not self.DC.isSignal[process]: return 1
        (processSource, foundDecay, foundEnergy) = getHiggsProdDecMode(bin,process,self.options)
        return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)

class StrictSMLikeHiggsModel(SMLikeHiggsModel):
    "Doesn't do anything more, but validates that the signal process names are correct"
    def getHiggsSignalYieldScale(self,production,decay, energy):
            if production == "VH": print "WARNING: VH production is deprecated and not supported in coupling fits"
            return "r"

class FloatingHiggsMass(SMLikeHiggsModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.mHRange = ['115','135'] # default
        self.rMode   = 'poi'
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("signalStrengthMode="): 
                self.rMode = po.replace("signalStrengthMode=","")
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI --- 
        POIs="MH"
        if self.rMode.startswith("fixed,"):
            self.modelBuilder.doVar("r[%s]" % self.rMode.replace("fixed,",""))
        else:
            self.modelBuilder.doVar("r[1,0,10]")
            if   self.rMode == "poi": POIs = "r,MH"
            elif self.rMode == "nuisance":  self.modelBuilder.out.var("r").setAttribute("flatParam")
            else: raise RuntimeError, "FloatingHiggsMass: the signal strength must be set to 'poi'(default), 'nuisance' or 'fixed,<value>'"
        if self.modelBuilder.out.var("MH"):
            self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
            self.modelBuilder.out.var("MH").setConstant(False)
        else:
            self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
        self.modelBuilder.doSet("POI",POIs)
    def getHiggsSignalYieldScale(self,production,decay, energy):
            return "r"


class FloatingXSHiggs(SMLikeHiggsModel):
    "Float independently ggH and qqH cross sections"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.modes = SM_HIGG_PROD
        self.mHRange  = []
        self.ggHRange = ['0', '4']
        self.qqHRange = ['0','10']
        self.VHRange  = ['0','20']
        self.WHRange  = ['0','20']
        self.ZHRange  = ['0','20']
        self.ttHRange = ['0','20']
        self.ttHasggH = False
        self.pois     = None
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("modes="): self.modes = po.replace("modes=","").split(",")
            if po.startswith("ttH=ggH"): 
                self.ttHasggH = True
            if po.startswith("poi="):
                self.pois = ",".join(["r_%s" % X for X in po.replace("poi=","").split(",")])
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Higgs mass range: Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("ggHRange="):
                self.ggHRange = po.replace("ggHRange=","").split(":")
                if len(self.ggHRange) != 2:
                    raise RuntimeError, "ggH signal strength range requires minimal and maximal value"
                elif float(self.ggHRange[0]) >= float(self.ggHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"
            if po.startswith("qqHRange="):
                self.qqHRange = po.replace("qqHRange=","").split(":")
                if len(self.qqHRange) != 2:
                    raise RuntimeError, "qqH signal strength range requires minimal and maximal value"
                elif float(self.qqHRange[0]) >= float(self.qqHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"                
            if po.startswith("VHRange="):
                self.VHRange = po.replace("VHRange=","").split(":")
                if len(self.VHRange) != 2:
                    raise RuntimeError, "VH signal strength range requires minimal and maximal value"
                elif float(self.VHRange[0]) >= float(self.VHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"
            if po.startswith("WHRange="):
                self.WHRange = po.replace("WHRange=","").split(":")
                if len(self.WHRange) != 2:
                    raise RuntimeError, "WH signal strength range requires minimal and maximal value"
                elif float(self.WHRange[0]) >= float(self.WHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"
            if po.startswith("ZHRange="):
                self.ZHRange = po.replace("ZHRange=","").split(":")
                if len(self.ZHRange) != 2:
                    raise RuntimeError, "ZH signal strength range requires minimal and maximal value"
                elif float(self.ZHRange[0]) >= float(self.ZHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"                
            if po.startswith("ttHRange="):
                self.ttHRange = po.replace("ttHRange=","").split(":")
                if len(self.ttHRange) != 2:
                    raise RuntimeError, "ttH signal strength range requires minimal and maximal value"
                elif float(self.ttHRange[0]) >= float(self.ttHRange[1]):
                    raise RuntimeError, "minimal and maximal range swapped. Second value must be larger first one"
        if self.ttHasggH:
            if "ggH" not in self.modes: raise RuntimeError, "Cannot set ttH=ggH if ggH is not an allowed mode"
            if "ttH" in self.modes: self.modes.remove("ttH")
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- Signal Strength as only POI ---
        if "ggH" in self.modes: self.modelBuilder.doVar("r_ggH[1,%s,%s]" % (self.ggHRange[0], self.ggHRange[1]))
        if "qqH" in self.modes: self.modelBuilder.doVar("r_qqH[1,%s,%s]" % (self.qqHRange[0], self.qqHRange[1]))
        if "VH"  in self.modes: self.modelBuilder.doVar("r_VH[1,%s,%s]"  % (self.VHRange [0], self.VHRange [1]))
        if "WH"  in self.modes: self.modelBuilder.doVar("r_WH[1,%s,%s]"  % (self.WHRange [0], self.WHRange [1]))
        if "ZH"  in self.modes: self.modelBuilder.doVar("r_ZH[1,%s,%s]"  % (self.ZHRange [0], self.ZHRange [1]))
        if "ttH" in self.modes: self.modelBuilder.doVar("r_ttH[1,%s,%s]" % (self.ttHRange[0], self.ttHRange[1]))
        poi = ",".join(["r_"+m for m in self.modes])
        if self.pois: poi = self.pois
        # --- Higgs Mass as other parameter ----
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poi+=',MH'
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poi+=',MH'
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if production in ["ggH","bbH"]: return ("r_ggH" if "ggH" in self.modes else 1)
        if production == "qqH": return ("r_qqH" if "qqH" in self.modes else 1)
        if production in ["ttH","tHq","tHW"]: return ("r_ttH" if "ttH" in self.modes else ("r_ggH" if self.ttHasggH else 1))
        if production in [ "WPlusH", "WMinusH", "WH", "ZH", "VH", "ggZH" ]: return ("r_VH" if "VH" in self.modes else 1)
        raise RuntimeError, "Unknown production mode '%s'" % production

class RvRfXSHiggs(SMLikeHiggsModel):
    "Float ggH and ttH together and VH and qqH together"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False        

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema."
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first."
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI --- 
        self.modelBuilder.doVar("RV[1,-5,15]")
        self.modelBuilder.doVar("RF[1,-4,8]")
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",'RV,RF,MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",'RV,RF')

    def getHiggsSignalYieldScale(self,production,decay, energy):
        if production in ['ggH', 'ttH', 'bbH', 'tHq', 'tHW']:
            return 'RF'
        if production in ['qqH', 'WH', 'ZH', 'VH', 'ggZH','WPlusH','WMinusH']:
            return 'RV'
        raise RuntimeError, "Unknown production mode '%s'" % production

class FloatingBRHiggs(SMLikeHiggsModel):
    "Float independently branching ratios"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.modes = SM_HIGG_DECAYS   #[ "hbb", "htt", "hgg", "hww", "hzz" ]
        self.modemap = {}
        self.mHRange = []
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("modes="): self.modes = po.replace("modes=","").split(",")
            if po.startswith("map="): 
                (mfrom,mto) = po.replace("map=","").split(":")
                self.modemap[mfrom] = mto
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- Signal Strength as only POI --- 
        for m in self.modes: 
            self.modelBuilder.doVar("r_%s[1,0,10]" % m);
        poi = ",".join(["r_"+m for m in self.modes])
        # --- Higgs Mass as other parameter ----
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poi+=',MH'
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poi+=',MH'
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if decay in self.modes: 
            return "r_"+decay
        if decay in self.modemap:
            if self.modemap[decay] in [ "1", "0" ]:
                return int(self.modemap[decay])
            else:
                return "r_"+self.modemap[decay]
        raise RuntimeError, "Unknown decay mode '%s'" % decay

class RvfBRHiggs(SMLikeHiggsModel):
    "Float ratio of (VH+qqH)/(ggH+ttH) and BR's"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False        
        self.modes = SM_HIGG_DECAYS #[ "hbb", "htt", "hgg", "hww", "hzz" ]
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("modes="): self.modes = po.replace("modes=","").split(",")
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema."
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first."
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI --- 
        self.modelBuilder.doVar("Rvf[1,-5,20]")
        poi = "Rvf"
        for mode in self.modes:
            poi += ',r_'+mode;
            self.modelBuilder.doVar("r_%s[1,0,5]" % mode)
            self.modelBuilder.factory_("expr::rv_%s(\"@0*@1\",Rvf,r_%s)" % (mode,mode))
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",poi+',MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if production in ['ggH', 'ttH', "bbH", "tHq", "tHW"]:
            return 'r_'+decay
        if production in ['qqH', 'WH','WPlusH','WMinusH', 'ZH', 'VH', "ggZH"]:
            return 'rv_'+decay
        raise RuntimeError, "Unknown production mode '%s'" % production

class ThetaVFBRHiggs(SMLikeHiggsModel):
    "Float ratio of (VH+qqH)/(ggH+ttH) and BR's"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False        
        self.modes = SM_HIGG_DECAYS #[ "hbb", "htt", "hgg", "hww", "hzz" ]
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("modes="): self.modes = po.replace("modes=","").split(",")
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema."
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first."
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        # --- Signal Strength as only POI --- 
        self.modelBuilder.doVar("thetaVF[0.78539816339744828,-1.5707963267948966,3.1415926535897931]")
        #self.modelBuilder.doVar("thetaVF[0.78539816339744828,0,1.5707963267948966]")
        poi = "thetaVF"
        for mode in self.modes:
            poi += ',r_'+mode;
            self.modelBuilder.doVar("r_%s[1,0,5]" % mode)
            self.modelBuilder.factory_("expr::rv_%s(\"sin(@0)*@1\",thetaVF,r_%s)" % (mode,mode))
            self.modelBuilder.factory_("expr::rf_%s(\"cos(@0)*@1\",thetaVF,r_%s)" % (mode,mode))
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",poi+',MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if production in ['ggH', 'ttH', 'bbH', 'tHq', 'tHW']:
            return 'rf_'+decay
        if production in ['qqH', 'WH','WPlusH','WMinusH', 'ZH', 'VH', 'ggZH']:
            return 'rv_'+decay
        raise RuntimeError, "Unknown production mode '%s'" % production


class FloatingXSBRHiggs(SMLikeHiggsModel):
    "Float independently cross sections and branching ratios"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.mHRange = []
        self.poiNames = []
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # --- Higgs Mass as other parameter ----
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                self.poiNames += [ 'MH' ]
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                self.poiNames += [ 'MH' ]
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        prod = 'VH' if production in [ 'VH','WH','WPlusH','WMinusH', 'ZH', 'ggZH' ] else production
        name = "r_%s_%s" % (prod,decay)
        if name not in self.poiNames: 
            self.poiNames += [ name ]
            self.modelBuilder.doVar(name+"[1,0,10]")
        return name
    def done(self):
        self.modelBuilder.doSet("POI",",".join(self.poiNames))
        
class DoubleRatioHiggs(SMLikeHiggsModel):
    "Measure the ratio of two BR's profiling mu_V/mu_F"
    def __init__(self):
        SMLikeHiggsModel.__init__(self) # not using 'super(x,self).__init__' since I don't understand it
        self.floatMass = False        
        self.modes = [ ]
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("modes="): self.modes = po.replace("modes=","").split(",")
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema."
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first."
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        if len(self.modes) != 2: raise RuntimeError, "must profide --PO modes=decay1,decay2"
        # --- Signal Strength as only POI --- 
        self.modelBuilder.doVar("rho[1,0,4]")
        self.modelBuilder.doVar("Rvf[1,0,4]")
        self.modelBuilder.doVar("rf_%s[1,0,4]" % self.modes[0])
        self.modelBuilder.factory_("prod::rf_%s(    rho, rf_%s)" % (self.modes[1], self.modes[0]))
        self.modelBuilder.factory_("prod::rv_%s(Rvf,rho, rf_%s)" % (self.modes[1], self.modes[0]))
        self.modelBuilder.factory_("prod::rv_%s(Rvf,     rf_%s)" % (self.modes[0], self.modes[0]))
        poi = "rho,Rvf,rf_%s" % self.modes[0]
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",poi+',MH')
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",poi)
    def getHiggsSignalYieldScale(self,production,decay, energy):
        if decay not in self.modes:
            print "Warning: BR of extra decay %s will be kept to SM value."
            return 1 if production in ['ggH', 'ttH', 'bbH', 'tHq', 'tHW'] else "Rvf"
        if production in ['ggH', 'ttH', 'bbH', 'tHq', 'tHW']:
            return 'rf_'+decay
        if production in ['qqH', 'WH', 'ZH', 'VH', 'ggZH']:
            return 'rv_'+decay
        raise RuntimeError, "Unknown production mode '%s'" % production

class RatioBRSMHiggs(SMLikeHiggsModel): 
    "Measure the ratio of BR's for two decay modes" 
    def __init__(self): 
        SMLikeHiggsModel.__init__(self)  
        self.floatMass = False        
        self.modes = SM_HIGG_DECAYS  #set( ("hbb", "htt", "hgg", "hzz", "hww") ) 
	self.denominator = "hww" 

    def setPhysicsOptions(self,physOptions): 
        for po in physOptions: 
            if po.startswith("denominator="):
		self.denominator = po.replace("denominator=","") 
            if po.startswith("higgsMassRange="): 
                self.floatMass = True 
                self.mHRange = po.replace("higgsMassRange=","").split(",") 
                print 'The Higgs mass range:', self.mHRange 
                if len(self.mHRange) != 2: 
                    raise RuntimeError, "Higgs mass range definition requires two extrema." 
                elif float(self.mHRange[0]) >= float(self.mHRange[1]): 
                    raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first."      
	self.numerators = tuple(self.modes - set((self.denominator,)))
	print 'denominator: ',self.denominator
	print 'numerators: ',self.numerators
	
	
    def doParametersOfInterest(self): 
        """Create POI out of signal strength, MH and BR's""" 
        
	den = self.denominator
        self.modelBuilder.doVar("r_VF[1,-5,5]")
        self.modelBuilder.doVar("r_F_%(den)s[1,0,5]" % locals())
	self.modelBuilder.factory_("prod::r_V_%(den)s(r_VF, r_F_%(den)s)" % locals())
	
	pois = []
	for numerator in self.numerators:
		names = {'num':numerator,'den':self.denominator}
		pois.append("r_%(num)s_%(den)s" % names )
	        self.modelBuilder.doVar("r_%(num)s_%(den)s[1,-5,5]" % names)
	        self.modelBuilder.factory_("prod::r_F_%(num)s(r_F_%(den)s, r_%(num)s_%(den)s)" % names)
        	self.modelBuilder.factory_("prod::r_V_%(num)s(r_VF, r_F_%(num)s)" % names)

	poi = ','.join(pois)
	
        # --- Higgs Mass as other parameter ---- 
        if self.floatMass: 
            if self.modelBuilder.out.var("MH"): 
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1])) 
                self.modelBuilder.out.var("MH").setConstant(False) 
            else: 
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
            self.modelBuilder.doSet("POI",poi+',MH') 
        else: 
            if self.modelBuilder.out.var("MH"): 
                self.modelBuilder.out.var("MH").setVal(self.options.mass) 
                self.modelBuilder.out.var("MH").setConstant(True) 
            else: 
                self.modelBuilder.doVar("MH[%g]" % self.options.mass) 
            self.modelBuilder.doSet("POI",poi)     


    def getHiggsSignalYieldScale(self,production,decay, energy): 
#        if decay not in self.numerators and not in self.denominator:
        if production in ['ggH', 'ttH', 'bbH', 'tHq', 'tHW']:
	    print '%(production)s/%(decay)s scaled by r_F_%(decay)s'%locals()
            return 'r_F_'+decay 
        if production in ['qqH', 'WH',"WPlusH","WMinusH", 'ZH', 'VH', 'ggZH']: 
	    print '%(production)s/%(decay)s scaled by r_V_%(decay)s'%locals()
            return 'r_V_'+decay 
        raise RuntimeError, "Unknown production mode '%s'" % production 

defaultModel = PhysicsModel()
multiSignalModel = MultiSignalModel()
strictSMLikeHiggs = StrictSMLikeHiggsModel()
floatingXSHiggs = FloatingXSHiggs()
rVrFXSHiggs = RvRfXSHiggs()
floatingBRHiggs = FloatingBRHiggs()
rVFBRHiggs = RvfBRHiggs()
thetaVFBRHiggs = ThetaVFBRHiggs()
floatingXSBRHiggs = FloatingXSBRHiggs()
floatingHiggsMass = FloatingHiggsMass()
doubleRatioHiggs = DoubleRatioHiggs()
ratioBRSMHiggs = RatioBRSMHiggs()
