from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from HiggsAnalysis.CombinedLimit.stage1 import stage1_procs # UPDATE THIS TO PUT DICT IN DATA
import ROOT, os, re

# NEED TO CONFIGURE ALSO FOR STAGE 0 and 1.1. Currently just working on stage 1

#List of all stage 1 processes
PROCESSES = [x for v in stage1_procs.itervalues() for x in v]
#List of decays which are defined in model
DECAYS = ['hzz','hbb','hww','hgg','hcc','tot']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to extract STXS production, decay mode and energy from process name (extracted in datacard)
def getSTXSProdDecMode(bin,process,options):
  processSource = process
  decaySource = options.fileName+":"+bin # by default, decay comes from datacard name or bin label
  if "_" in process:
    #decay at end of process name after final _: join all previous parts to give processSource 
    (processSource, decaySource) = "_".join(process.split("_")[0:-1]),process.split("_")[-1]
  foundDecay = None
  for D in ALL_HIGGS_DECAYS:
    if D in decaySource:
      if foundDecay: raise RuntimeError, "Validation Error: decay string %s contains multiple known decay names" % decaySource
      foundDecay = D
  if not foundDecay: raise RuntimeError, "Validation Error: decay string %s does not contain any known decay name" % decaySource

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
    foundEnergy = '13TeV' ## if using 81x, chances are its 13 TeV
    print "Warning: decay string %s does not contain any known energy, assuming %s" % (decaySource, foundEnergy)
  return (processSource, foundDecay, foundEnergy)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# STXS to EFT abstract base class: inherited classes for different stages
class STXStoEFTHiggsModel(SMLikeHiggsModel):

  # initialisation: include options for STXS bin and BR uncertainties (FOR NOW FALSE)
  # TO DO: add function for STXS uncertainties, like partial width unc in SMHiggsBuilder
  def __init__(self,STXSU=False,BRU=False,fixTHandBBH=True):
    SMLikeHiggsModel.__init__(self)
    self.floatMass = False #Initally false, require external option to float mass
    self.doSTXSU = STXSU
    self.doBRU = BRU
    self.fixTHandBBH = fixTHandBBH
  def setPhysicsOption(self,physOptions):
    for po in physOptions:
      if po.startswith("higgsMassRange="):
        self.floatMass = True
        self.mHRange = po.replace("higgsMassRange=","").split(",")
        if len(self.mHRange) != 2:
          raise RuntimeError, "Higgs mass range definition requires two extrema"
        elif float(self.mHRange[0]) >= float(self.mHRange[1]):
          raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
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

  #overwrite getYieldScale to use getSTXSProdDecMode
  def getYieldScale(self,bin,process):
    "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
    if not self.DC.isSignal[process]: return 1
    (processSource, foundDecay, foundEnergy) = getSTXSProdDecMode(bin,process,self.options)
    return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Functions for extracting the pois. scaling functions etc from .txt file defined in datadir+"/eft/XX.txt
  def textToPOIList( self, filename, skipRows=1 ):
    self.pois = {} #initiate empty dict
    #print "[DEBUG] filename = %s"%filename
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    for line in lines[skipRows:]:
      if( len(line.strip())== 0 ): continue
      self.pois['%s'%line.split()[0]] = "[%s,%s,%s]"%(line.split()[1],line.split()[2],line.split()[3])
    #print "[DEBUG] self.pois =", self.pois
    file_in.close()

  def textToSTXSScalingFunction( self, filename, skipRows=0 ):
    stxs_id = {} #dict for integer id to STXS process name
    self.STXSScalingFunction = {} #dict for STXS process name to function (in terms of EFT coefficients)
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    #Initially loop over first block of text until empty line to get integer id of processes
    emptyLinePosition = 0
    for line in lines[skipRows:]:
      emptyLinePosition += 1
      if( len(line.strip())==0 ): break
      proc, procid = line.split()
      #parse process to desired convention
      proc = re.sub('r_','',proc)
      proc = proc.upper()
      proc = re.sub('GG2H','ggH',proc)
      proc = re.sub('VBF_QQ2HQQ','qqH',proc)
      proc = re.sub('QQ2HLNU','WH_lep',proc)
      proc = re.sub('QQ2HLL','ZH_lep',proc)
      proc = re.sub('WH_QQ2HQQ','WH_had',proc)
      proc = re.sub('ZH_QQ2HQQ','ZH_had',proc)
      proc = re.sub('TTH','ttH',proc) 
      #Add to dictionary
      stxs_id[ procid ] = proc
    #print "[DEBUG] stxs_id =", stxs_id
    #print "[DEBUG] emptyLinePosition = %g"%emptyLinePosition
    #Loop over second block of text (starting at empty line) to get scaling functions
    _id = None
    for line in lines[emptyLinePosition:]:
      if( len(line.strip())==0 ): continue
      # procid and scaling function on subsequent lines 
      # save proc_id of line starting with bin number
      if line.startswith("Bin number"): _id = line.split()[2]
      # extract scaling function
      elif line.startswith("1"): #perturbation theory: safer than using else statement
        if _id == None:
          raise ValueError("[ERROR] Process ID not set. Cannot link function to STXS bin")
        function = line.replace('\n','')
        function = "".join(function.split(" "))
        self.STXSScalingFunction[ stxs_id[ _id ] ] = function 
        _id = None #Reset
    #print "[DEBUG] self.STXSScalingFunction =", self.STXSScalingFunction
    file_in.close()
  
  def textToDecayScalingFunction( self, filename, skipRows=0 ):
    self.DecayScalingFunction = {}
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    #Decay definition and scaling function on subsequent line
    decay_ = None
    for line in lines[skipRows:]:
      if( len(line.strip())==0 ): continue
      #save decay definion of lines starting w/ Bin number
      if line.startswith("Bin number"):
        decay_ = line.split()[-1]
        if( decay_ == 'hzzto4l' ): decay_ = re.sub("hzzto4l","hzz",decay_)
        elif( decay_ == 'dec_full' ): decay_ = re.sub("dec_full","tot",decay_)
      elif line.startswith("1"): #perturbation theory: safer than using else statement
        if decay_ == None:
          raise ValueError("[ERROR] Decay not set. Cannot link function to decay channel")
        function = line.replace('\n','')
        function = "".join(function.split(" ")) #remove spaces
        self.DecayScalingFunction[ decay_ ] = function
        _decay = None #reset
    #print "[DEBUG] self.DecayScalingFunction =", self.DecayScalingFunction
    file_in.close() 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Functions for making scaling
  def makeScalingFunction( self, what ): 

    #if in processes/decays extract formula from corresponding dict
    if what in PROCESSES: formula = self.STXSScalingFunction[ what ]
    elif what in DECAYS: formula = self.DecayScalingFunction[ what ]
    else:
      raise ValueError("[ERROR] Scaling function for %s does not exist"%what)

    #turn formula into list: splitting by +/*/- (keeping delimiters)
    formula = re.split('([+,*,-])', formula)
    #dict to store integer for poi: iterate integer when using new poi
    poi_id = {}
    poi_ = 0
    #strings to hold formula + list of parameters of interest
    formulaStr = 'expr::scaling_%s(\"'%what
    poiStr = '\",'
    #loop over elements in formula
    for element in formula:
      if(element.startswith("c"))|(element.startswith("t")): #if element is a poi
        if( element not in self.pois )&( element not in ['cWW','cB'] ):
          raise ValueError("[ERROR] %s not defined in POI list"%element)
        #check if poi has already been used in formula
        if element in poi_id: formulaStr += "@%g"%poi_id[ element ] #if yes: add @(id) to formula string
        else: 
          #add poi to dictionary and poiStr and add one to iterator
          poi_id[ element ] = poi_
          poi_ += 1
          #add @(id) to formulaStr and poi to correct position in poiStr
          formulaStr += "@%g"%poi_id[ element ]
          poiStr += "%s,"%element
      else:
        #element is not a poi: simply add to formula Str
        formulaStr += element
    #combine removing final comma of poiStr
    functionStr = formulaStr+poiStr[:-1]+")"
    #print "[DEBUG] Expression for %s: %s\n"%(what,functionStr)
    self.modelBuilder.factory_( functionStr )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Function extracting the STXS bin uncertainties
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def doParametersOfInterest(self):
    if self.floatMass: print "[WARNING] Floating Higgs mass selected. STXStoEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameter list from file using dedicated: textToPOIList function, setting self.pois    
    # self.pois directionary of pois name and string of allowed range: e.g {'cG':'[0,-3.2e-4,1.1e-4]',...}
    self.textToPOIList( os.path.join(self.SMH.datadir,'eft/pois_plus_undefined.txt'))
    POIs = ','.join(self.pois.keys())
    for poi, poi_range in self.pois.iteritems(): 
      self.modelBuilder.doVar("%s%s"%(poi,poi_range))
    self.modelBuilder.doSet("POI",POIs)      
    #POIs for cWW and cB defined in terms of constraints on cWW+cB and cWW-cB: define expression for individual coefficient
    self.modelBuilder.factory_("expr::cWW(\"0.5*(@0+@1)\",cWWPluscB,cWWMinuscB)")
    self.modelBuilder.factory_("expr::cB(\"0.5*(@0-@1)\",cWWPluscB,cWWMinuscB)")

    #set up model
    self.setup()



  def setup(self):
    #For additional options e.g. STXS/BR uncertainties
    # ...

    #Read scaling functions for STXS bins and decays from text files
    self.textToSTXSScalingFunction( os.path.join(self.SMH.datadir, 'eft/crosssections.txt') )
    self.textToDecayScalingFunction( os.path.join(self.SMH.datadir, 'eft/decay.txt' ))
    #Make scaling for each STXS process and decay:

    for proc in PROCESSES: 
      #if fixTHandBBH...
      if( self.fixTHandBBH ): 
        if proc in ['tHq','tHW','bbH']: self.modelBuilder.factory_("expr::scaling_%s(\"@0\",1.)"%proc)
        else: self.makeScalingFunction( proc )
      #else: scaling function must be defined in text file
      else: self.makeScalingFunction( proc )
    for dec in DECAYS: self.makeScalingFunction( dec ) 

  def getHiggsSignalYieldScale(self,production,decay,energy):
    name = "stxs2eft_scaling_%s_%s_%s"%(production,decay,energy)
    if self.modelBuilder.out.function(name) == None:
      #Check production scaling exists:
      XSscal = None
      if production in PROCESSES: XSscal = "scaling_%s"%(production)
      else:
        raise ValueError("[ERROR] Process %s is not supported in Stage 1 Model"%production)
      
      #Check decay scaling exists
      BRscal = None
      if decay in DECAYS:
        BRscal = "scaling_%s"%decay
      else:
        raise ValueError("[ERROR] Decay %d is not defined"%decay)
      
      #Combine XS and BR scaling
      self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))

      print '[STXStoEFT Stage 1]', name, production, decay, energy, ":", self.modelBuilder.out.function(name).Print("")
    return name

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

STXStoEFT = STXStoEFTHiggsModel()

#For debugging
#if __name__ == "__main__":
#    print " ------------------- "
#    STXStoEFT.textToPOIList( '/afs/cern.ch/work/j/jlangfor/combine/STXStoEFT/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/eft/pois_plus_undefined.txt')
#    STXStoEFT.textToSTXSScalingFunction( '/afs/cern.ch/work/j/jlangfor/combine/STXStoEFT/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/eft/crosssections.txt')
#    STXStoEFT.textToDecayScalingFunction( '/afs/cern.ch/work/j/jlangfor/combine/STXStoEFT/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/eft/decay.txt')
#    for proc in STAGE1_PROCESSES: STXStoEFT.makeScalingFunction( proc )
#    for dec in STAGE1_DECAYS: STXStoEFT.makeScalingFunction( dec )
#    print "--------------------"

# Author: Jonathon Langford (ICL)
# Date: 06/02/2019
# Description: Model to describe how bins in STXS (0,1,1.1) scale using full set of dimension-6 EFT parameters
#              Equations calculated using Higgs Effective Lagrangian (Madgraph)
#              Model encompasses S0,S1 and S1.1
#              TO DO: Ask about sub-bins (these will be a separate gen-level process?). If so likely to only have 
#                     scaling for combined bin. Therefore bin name should have commonality to allow merging
#                     Also in response matrices: this will completely blow up (is it just for theory uncertainties?)

# TO DO LIST
#    * Configure for bbh, thq, thW: What to do here (fix to SM i.e. scaling = 1), copy off previous models



