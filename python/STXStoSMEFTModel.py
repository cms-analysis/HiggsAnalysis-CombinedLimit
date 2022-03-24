# Author: Jonathon Langford (ICL)
# Date: 3/2022
# Description: Model to describe how bins in STXS stage 1.2 scale using full set of dimension-6 EFT parameters
#              Equations calculated using nanoAOD reweighting using Madgraph reweighting modules
#              SMEFTsim

from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from math import exp
import ROOT, os, re, sys
import json
import yaml
from collections import OrderedDict as od

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordered dicts crucial: e.g. choose WH_had before WH

MAP_HIGGS_DECAY_SMEFT = od()
MAP_HIGGS_DECAY_SMEFT["hgg"] = "gamgam"
MAP_HIGGS_DECAY_SMEFT["hzz"] = "ZZ"
MAP_HIGGS_DECAY_SMEFT["hbb"] = "bb"
MAP_HIGGS_DECAY_SMEFT["hww"] = "WW"
MAP_HIGGS_DECAY_SMEFT["htt"] = "tautau"
MAP_HIGGS_DECAY_SMEFT["hmm"] = "mumu"
MAP_HIGGS_DECAY_SMEFT["hzg"] = "Zgam"

MAP_HIGGS_PROD_SMEFT = od()
MAP_HIGGS_PROD_SMEFT["ggH"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["qqH"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["WH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["ZH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["ggZH_had"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["ggZH_qq"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["WH_lep"] = "QQ2HLNU"
MAP_HIGGS_PROD_SMEFT["ZH_lep"] = "QQ2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_lep"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_ll"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_nunu"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ttH"] = "TTH"
MAP_HIGGS_PROD_SMEFT["tHq"] = "THQ"
MAP_HIGGS_PROD_SMEFT["tHW"] = "THW"
MAP_HIGGS_PROD_SMEFT["bbH"] = "BBH"
# If only specify VH: use leptonic equations as most likely to enter VH leptonic tag?
#MAP_HIGGS_PROD_SMEFT["WH"] = "QQ2HLNU"
#MAP_HIGGS_PROD_SMEFT["ZH"] = "QQ2HLL"
#MAP_HIGGS_PROD_SMEFT["ggZH"] = "GG2HLL"

MAP_HIGGS_PROD_SMEFT_ATLAS = od()
MAP_HIGGS_PROD_SMEFT_ATLAS["ggH"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["qqH"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["WH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["ZH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_had"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_qq"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["WH_lep"] = "QQ2HLNU"
MAP_HIGGS_PROD_SMEFT_ATLAS["ZH_lep"] = "QQ2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_lep"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_ll"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_nunu"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ttH"] = "TTH"
MAP_HIGGS_PROD_SMEFT_ATLAS["tHq"] = "TH"
MAP_HIGGS_PROD_SMEFT_ATLAS["tHW"] = "TH"
MAP_HIGGS_PROD_SMEFT_ATLAS["bbH"] = "BBH"
# If only specify VH: use leptonic equations as most likely to enter VH leptonic tag?
#MAP_HIGGS_PROD_SMEFT_ATLAS["WH"] = "QQ2HLNU"
#MAP_HIGGS_PROD_SMEFT_ATLAS["ZH"] = "QQ2HLL"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Global function to extract reco category, STXS bin, decay mode and energy from process name
def getProcessInfo(bin,process):
  foundRecoCategory = bin
  foundSTXSBin = process
  foundDecay = None
  foundEnergy = "13TeV"
  #Iterate over Higgs decays
  matchedDecayString = False 
  for D in ALL_HIGGS_DECAYS:
    if matchedDecayString: continue
    if "_%s"%D in foundSTXSBin:
      foundSTXSBin = re.sub('_%s'%D,'',foundSTXSBin)
      foundDecay = D
      matchedDecayString = True
  # Also drop year tag in STXS bin name if present
  for Y in ['2016','2017','2018']:
    if "_%s"%Y in foundSTXSBin:
      foundSTXSBin = re.sub('_%s'%Y,'',foundSTXSBin)

  # Catch for H->Zgam
  if( foundDecay == "hzg" )|( "bkg" in foundSTXSBin ): foundSTXSBin = foundSTXSBin.split("_")[0]

  if not matchedDecayString: raise RuntimeError, "Validation Error: no supported decay found in process"

  return (foundRecoCategory, foundSTXSBin, foundDecay, foundEnergy)

#################################################################################################################
# STXS to EFT abstract base class: inherited classes for different stages
class STXStoSMEFTBaseModel(SMLikeHiggsModel):

  def __init__(self,fixProcesses=[]):
    SMLikeHiggsModel.__init__(self)
    self.PROCESSES = None
    self.DECAYS = None
    # Dicts to store pois + scaling functions
    self.pois = None
    self.poiNameMap = {} # To account for exponents in poi names, which aren't in input json files
    self.STXSScalingTerms = None
    self.DecayScalingTerms = None
    self.map_prod = MAP_HIGGS_PROD_SMEFT
    self.map_decay = MAP_HIGGS_DECAY_SMEFT
    # Options
    self.floatMass = False
    self.fixProcesses = fixProcesses #Option to fix certain STXS bins: comma separated list of STXS bins
    self.linearOnly=False
    self.stage0=False
    self.parametrisation="CMS-prelim-SMEFT-topU3l_22_03_21"
 
  def setPhysicsOptionsBase(self,physOptions):
    for po in physOptions:
      if po.startswith("higgsMassRange="):
        self.floatMass = True
        self.mHRange = po.replace("higgsMassRange=","").split(",")
        if len(self.mHRange) != 2:
          raise RuntimeError, "Higgs mass range definition requires two extrema"
        elif float(self.mHRange[0]) >= float(self.mHRange[1]):
          raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
      if po.startswith("fixProcesses="): 
        self.fixProcesses = (po.replace("fixProcesses=","")).split(",")
      if po.startswith("linearOnly="): 
        self.linearOnly = (po.replace("linearOnly=","") in ["yes","1","Yes","True","true"])
      if po.startswith("stage0="): 
        self.stage0 = (po.replace("stage0=","") in ["yes","1","Yes","True","true"])
      if po.startswith("parametrisation="): 
        self.parametrisation = po.replace("parametrisation=","")
        if "ATLAS" in self.parametrisation: self.map_prod = MAP_HIGGS_PROD_SMEFT_ATLAS

    #Output options to screen
    print " --> [STXStoSMEFT] Using (%s) parametrisation"%self.parametrisation
    if( len( self.fixProcesses ) > 0 ): print " --> [STXStoSMEFT] Fixing following processes to SM: %s"%self.fixProcesses
    if self.linearOnly: print " --> [STXStoSMEFT] Only linear terms (Aj)"

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

  # Overwrite getYieldScale to extract (RECO-category,STXS bin,decay,energy)
  def getYieldScale(self,bin,process):
    if not self.DC.isSignal[process]: 
      return 1.

    # Extract process line info
    (recocat, stxsbin, decay, energy) = getProcessInfo(bin,process)

    # Return 1 (no scaling) for fixed processes and scaling for non-fixed
    if stxsbin in self.fixProcesses: 
      return 1. 
    else: 
      procStr = stxsbin
      return self.getHiggsSignalYieldScale(procStr, decay, energy)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract pois from yaml file
  def extractPOIs( self, filename ):
    with open( filename, 'r' ) as fpois:
      try:
        self.pois = yaml.safe_load(fpois)
      except yaml.YAMLERROR as exc:
        print exc

  #Function to extract STXS scaling terms from json file
  def extractSTXSScalingTerms( self, filename="" ):
    if filename != "":
      with open(filename,"r") as jf: self.STXSScalingTerms = json.load(jf) 
    else: 
      self.STXSScalingTerms = {}
  
  #Function to extract decay scaling functions from file
  def extractDecayScalingTerms( self, filename="" ):
    if filename != "":
      with open(filename,"r") as jf: self.DecayScalingTerms = json.load(jf)
    else:
      self.DecayScalingTerms = {}
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Function to make scaling function in workspace from terms
  def makeScalingFunction( self, what, isDecay = False ):

    # Apply mapping of production mode/decay to match inputs in json file
    k = what
    if isDecay:
      for D in self.map_decay.keys():
        if what == D:
          k = self.map_decay[D]
    else:
      for P in self.map_prod.keys():
        if P in what: 
          k = re.sub(P,self.map_prod[P],what)  

    # Fix for ttH multilepton: missing lep label in VH bins
    if "WH_PTV" in k: k = re.sub("WH","QQ2HLNU",k)
    if "ZH_PTV" in k: k = re.sub("ZH","QQ2HLL",k)
    if "ggZH_PTV" in k: k = re.sub("ggZH","GG2HLL",k)

    # Fix for ttH multilepton: duplicate of proc names
    if not isDecay:
      for P in self.map_prod.values():
        if "%s_%s"%(P,P) in k: k = re.sub("%s_%s"%(P,P),P,k)

    # Fix for VH procs without had/lep label: use leptonic scaling function. Is this accurate?
    if k == "WH": k = "QQ2HLNU"
    if k == "ZH": k = "QQ2HLL"
    if k == "ggZH": k = "GG2HLL"

    # Extract terms for dict
    if k in self.STXSScalingTerms: terms = self.STXSScalingTerms[k] 
    elif k in self.DecayScalingTerms: terms = self.DecayScalingTerms[k]
    else:
      print " --> [WARNING] Scaling terms for %s do not exist in input json. Setting to 1"%k
      terms = {}
      #raise ValueError("[ERROR] Scaling terms for %s do not exist"%what)

    # Loop over pois and extract the terms from scaling function, stored in C++ map
    coeffs = ROOT.std.map("string","double")()
    A, B = od(), od()
    for jpoi in self.pois:
      # Interference terms: Aj
      e_jpoi = 10**(-1*self.pois[jpoi]['exponent'])
      jpoi_name = self.poiNameMap[jpoi]
      if "A_%s"%jpoi in terms: coeffs[jpoi_name] = e_jpoi*terms["A_%s"%jpoi]
      # BSM-only terms: Bjk
      if not self.linearOnly:
	if "B_%s_2"%jpoi in terms: coeffs['%s_2'%jpoi_name] = e_jpoi*e_jpoi*terms["B_%s_2"%jpoi]
	# Cross terms
	for kpoi in self.pois:
          e_kpoi = 10**(-1*self.pois[kpoi]['exponent'])
          kpoi_name = self.poiNameMap[kpoi]
	  if "B_%s_%s"%(jpoi,kpoi) in terms: coeffs['%s_%s'%(jpoi_name,kpoi_name)] = e_jpoi*e_kpoi*terms["B_%s_%s"%(jpoi,kpoi)]
        
    # Make RooEFTScalingFunction
    if( isDecay )&( what != "tot" ): name = "scaling_partial_%s"%what 
    else: name = "scaling_%s"%what 
    eft_scaling = ROOT.RooEFTScalingFunction(name,name,coeffs,self.POIs)
          
    #Add scaling function as RooAddition into model
    self.modelBuilder.out._import(eft_scaling)

  #Function to make BR scaling functions: partial width/total width
  def makeBRScalingFunction( self, what ): self.modelBuilder.factory_( 'expr::scaling_BR_%s("@0/@1", scaling_partial_%s, scaling_tot)'%(what,what) )
    
#################################################################################################################
# Combination of different stages
class STXSToSMEFTModel(STXStoSMEFTBaseModel):
  def __init__(self):
    STXStoSMEFTBaseModel.__init__(self)

  def setPhysicsOptions(self,physOptions):
    self.setPhysicsOptionsBase(physOptions)
  
  def doParametersOfInterest(self):
    if self.floatMass: print " --> [WARNING] Floating Higgs mass selected. STXStoSMEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameters of interest from yaml file
    self.extractPOIs("%s/src/HiggsAnalysis/CombinedLimit/data/eft/EFTScalingEquations/equations/%s/pois.yaml"%(os.environ['CMSSW_BASE'],self.parametrisation))

    # Create list of pois and build RooRealVars
    POIs = []
    for poi in self.pois: 
      poi_name = poi if self.pois[poi]['exponent'] == 0 else "%sXE%g"%(poi,self.pois[poi]['exponent'])
      self.poiNameMap[poi] = poi_name
      POIs.append(poi_name)
      self.modelBuilder.doVar("%s[%g,%g,%g]"%(poi_name,self.pois[poi]['val'],self.pois[poi]['min'],self.pois[poi]['max']))
      self.modelBuilder.out.var(poi_name).setConstant(True)

    self.modelBuilder.doSet("POI",",".join(POIs))
    self.POIs = ROOT.RooArgList( self.modelBuilder.out.set("POI") )

    #set up model
    self.setup()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def setup(self):
 
    # Extract scaling terms from json files: inclusive vs reco-level
    self.extractSTXSScalingTerms(filename="%s/src/HiggsAnalysis/CombinedLimit/data/eft/EFTScalingEquations/equations/%s/prod.json"%(os.environ['CMSSW_BASE'],self.parametrisation))    
    self.extractDecayScalingTerms(filename="%s/src/HiggsAnalysis/CombinedLimit/data/eft/EFTScalingEquations/equations/%s/decay.json"%(os.environ['CMSSW_BASE'],self.parametrisation))    

    # Make total scaling function for decay side
    self.makeScalingFunction("tot", isDecay=True)
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getHiggsSignalYieldScale(self,production,decay,energy):

    # Function to convert troublesome procs into viable one for HC combination
    production = convert_to_STXS(production,decay)

    # Stage 0 option: use stage 0 bin scaling
    if self.stage0:
      for P in self.map_prod.keys():
        if P in production: production = P

    name = "stxstosmeft_scaling_%s_%s_%s"%(production,decay,energy)

    if self.modelBuilder.out.function(name) == None:

      # Build scaling functions if they do not exist
      if self.modelBuilder.out.function("scaling_%s"%production) == None:
        print " --> [STXStoSMEFT] Making scaling function for STXS bin: %s"%production
        self.makeScalingFunction(production)
      XSscal = "scaling_%s"%production

      if self.modelBuilder.out.function("scaling_BR_%s"%decay) == None:
        print " --> [STXStoSMEFT] Making scaling function for decay: %s"%decay
        self.makeScalingFunction(decay, isDecay=True)
        self.makeBRScalingFunction(decay)
      BRscal = "scaling_BR_%s"%decay
  
      #Combine XS and BR scaling: incuding theory unc if option selected
      self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))
      
    return name
 

#################################################################################################################
# Function to convert troublesome procs to the names in the json files
def convert_to_STXS( _production, _decay ):
  # Add string replace functions
  return _production


#################################################################################################################
# Instantiation of STXStoSMEFT model
STXStoSMEFT = STXSToSMEFTModel()

