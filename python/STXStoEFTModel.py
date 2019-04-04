# Author: Jonathon Langford (ICL)
# Date: 06/02/2019
# Description: Model to describe how bins in STXS (0,1,1.1) scale using full set of dimension-6 EFT parameters
#              Equations calculated using Higgs Effective Lagrangian (Madgraph)

from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from math import exp
import ROOT, os, re

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Global function to extract STXS production, decay mode and energy from process name
#   * this has changed with STXS under naming convention
def getSTXSProdDecMode(bin,process,options):
  matchedDecayString = False#Boolean 
  processSource = process
  decaySource = options.fileName+":"+bin # by default, decay comes from datacard name or bin label
  foundDecay = None
  foundEnergy = "13TeV"
  #Iterate over Higgs decays
  for D in ALL_HIGGS_DECAYS:
    if matchedDecayString: continue
    # if decay in process name: set process and decay
    if "_%s"%D in process:
      processSource = re.sub('_%s'%D,'',process)
      foundDecay = D
      matchedDecayString = True
  if not matchedDecayString: raise RuntimeError, "Validation Error: no supported decay found in process"
  return (processSource, foundDecay, foundEnergy)

#################################################################################################################
# STXS to EFT abstract base class: inherited classes for different stages
class STXStoEFTBaseModel(SMLikeHiggsModel):

  # initialisation: include options for STXS bin and BR uncertainties
  #    * note: STXS bin uncertainties are defined in data/lhc-hxswg/eft/stageX/BinUncertainties.txt. Needs updating!
  def __init__(self,STXSU=False,BRU=False,fixTHandBBH=True,freezeOtherParameters=True,fixProcesses=[]):
    SMLikeHiggsModel.__init__(self)
    self.PROCESSES = []
    self.DECAYS = []
    self.floatMass = False #Initally false, require external option to float mass
    self.doSTXSU = STXSU
    self.doBRU = BRU
    self.fixTHandBBH = fixTHandBBH #if false scaling function for tH, bbH MUST be defined in input text file: data/lhc-hxswg/eft/stageX/crosssections.txt
    self.freezeOtherParameters = freezeOtherParameters #Option to freeze majority of parameters in model. Leaving those used in LHCHXSWG-INT-2017-001 fit to float
    self.fixProcesses = fixProcesses #Option to fix certain STXS bins: comma separated list of STXS bins

  def setPhysicsOptionsBase(self,physOptions):
    for po in physOptions:
      if po.startswith("higgsMassRange="):
        self.floatMass = True
        self.mHRange = po.replace("higgsMassRange=","").split(",")
        if len(self.mHRange) != 2:
          raise RuntimeError, "Higgs mass range definition requires two extrema"
        elif float(self.mHRange[0]) >= float(self.mHRange[1]):
          raise RuntimeError, "Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"
      if po.startswith("BRU="):
        self.doBRU = (po.replace("BRU=","") in ["yes","1","Yes","True","true"])
      if po.startswith("STXSU="):
        self.doSTXSU = (po.replace("STXSU=","") in ["yes","1","Yes","True","true"])
      if po.startswith("freezeOtherParameters="):
        self.freezeOtherParameters = (po.replace("freezeOtherParameters=","") in ["yes","1","Yes","True","true"])
      if po.startswith("fixProcesses="): 
        self.fixProcesses = (po.replace("fixProcesses=","")).split(",")
    #Output option used to terminal
    print "[STXStoEFT] Theory uncertainties in partial widths: %s"%self.doBRU
    print "[STXStoEFT] Theory uncertainties in STXS bins: %s"%self.doSTXSU
    if( self.doSTXSU ): print "   [WARNING]: theory uncertainties in STXS bins are currently incorrect. Need to update: data/lhc-hxswg/eft/stageX/BinUncertainties.txt"
    if( self.freezeOtherParameters ): print "[STXStoEFT] Freezing all but [cG,cA,cu,cHW,cWWMinuscB] to 0"
    else: print "[STXStoEFT] Allowing all HEL parameters to float"
    if( len( self.fixProcesses ) > 0 ): print "[STXStoEFT] Fixing following processes to SM: %s"%self.fixProcesses

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
    #Return 1 for fixed processes and scaling for non-fixed
    if processSource in self.fixProcesses: return 1 
    else: return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # READING FROM TEXT FILES
  #Function for extracting the pois: save as dictionary of pois
  def textToPOIList( self, filename, skipRows=1 ):
    self.pois = {} #initiate empty dict
    self.poi_scaling = {} #separate dictionary for how POIs scale
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    for line in lines[skipRows:]:
      if( len(line.strip())== 0 ): continue
      self.pois['%s'%line.split()[0]] = "[%s,%s,%s]"%(line.split()[1],line.split()[2],line.split()[3])
      #extract unscaled poi and scaling factor in list
      poi_scale = line.split()[0].split("_")
      #save poi scaling to dictionary: used when making scaling function
      if len(poi_scale) == 1: self.poi_scaling[ poi_scale[0] ] = "1.*%s"%line.split()[0]
      else:
        if poi_scale[1] == "x05": self.poi_scaling[ poi_scale[0] ] = "0.00001*%s"%line.split()[0]
        elif poi_scale[1] == "x04": self.poi_scaling[ poi_scale[0] ] = "0.0001*%s"%line.split()[0]
        elif poi_scale[1] == "x03": self.poi_scaling[ poi_scale[0] ] = "0.001*%s"%line.split()[0]
        elif poi_scale[1] == "x02": self.poi_scaling[ poi_scale[0] ] = "0.01*%s"%line.split()[0]
        elif poi_scale[1] == "x01": self.poi_scaling[ poi_scale[0] ] = "0.1*%s"%line.split()[0]
        else: raise ValueError("[ERROR] Parameter of interest scaling not in allowed range [0.1-0.0001]")
      #if poi: cG, cA, tcG, tcA add factor of 16pi2
      if poi_scale[0] in ['cG','tcG','cA','tcA']: self.poi_scaling[ poi_scale[0] ] = "157.914*%s"%(self.poi_scaling[ poi_scale[0] ])
    file_in.close()

  #Function to extract STXS scaling functions
  def textToSTXSScalingFunctions( self, filename, skipRows=0 ):
    stxs_id = {} #dict for integer id to STXS process name
    self.STXSScalingFunctions = {} #dict for STXS process name to function (in terms of EFT coefficients)
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

    #Loop over second block of text (starting at empty line) to get scaling functions
    _id = None
    for line in lines[emptyLinePosition:]:
      if( len(line.strip())==0 ): continue
      # procid and scaling function on subsequent lines 
      # save proc_id of line beginning with "Bin number"
      if line.startswith("Bin number"): _id = line.split()[2]
      # extract scaling function
      elif line.startswith("1"): #perturbation theory: safer than using else statement
        if _id == None:
          raise ValueError("[ERROR] Process ID not set. Cannot link function to STXS bin")
        function = line.replace('\n','')
        function = "".join(function.split(" ")) #remove space
        self.STXSScalingFunctions[ stxs_id[ _id ] ] = function 
        _id = None #Reset
    file_in.close()
  
  #Function to extract decay scaling functions from file
  def textToDecayScalingFunctions( self, filename, skipRows=0 ):
    self.DecayScalingFunctions = {}
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
        self.DecayScalingFunctions[ decay_ ] = function
        _decay = None #reset
    file_in.close() 

  # ** NEW **
  #Function to make scaling function from string
  #   > defines each terms as RooProduct
  #   > sum of terms as RooAddition
  def makeScalingFunction( self, what, STXSstage="1" ):
  
    #if in processes/decays extract formula from corresponding dict
    if what in self.STXSScalingFunctions: formula = self.STXSScalingFunctions[ what ]
    elif what in self.DecayScalingFunctions: formula = self.DecayScalingFunctions[ what ]
    else:
      raise ValueError("[ERROR] Scaling function for %s does not exist for STXS Stage %s"%(what,STXSstage))

    #replace "-" in formula string by "+-" and then turn into list, splitting by delimeter "+"
    formula = re.sub("-","+-",formula).split("+")
    
    #define list to hold name of terms 
    formulaTerms = []
    _termIdx = 0
 
    #loop over terms in formula (ignoring first term which is just 1 [perturbation theory])
    for term in formula[1:]:
      constituents = term.split("*")
      #replace constituents by POIs (i.e scaled)
      scaled_constituents = []
      for c in constituents:
        if c in self.poi_scaling: scaled_constituents.extend( self.poi_scaling[c].split("*") )
        else: scaled_constituents.append( c )
      #define string for term: RooProduct (factory)
      termStr = "prod::term_%s_%g("%(what,_termIdx)
      for sc in scaled_constituents: termStr += "%s,"%sc
      termStr = termStr[:-1]+")"
      # add term to model and dictionary
      self.modelBuilder.factory_( termStr )
      formulaTerms.append("term_%s_%g"%(what,_termIdx))
      # add one to the iterator
      _termIdx += 1

    #split up sums into sizeable chunks: 15 terms max
    sumTerms = {}
    _sumIdx = -1 #start at -1
    _termIdx = 0    
    for _termIdx in range(len(formulaTerms)): 
      if _termIdx % 25 == 0: _sumIdx += 1
      #Add terms to sum
      if( what in self.DecayScalingFunctions )&( what != "tot" ): sumString = "scaling_partial_%s_%s"%(what,_sumIdx)
      else: sumString = "scaling_%s_%s"%(what,_sumIdx)
      if sumString in sumTerms: sumTerms[ sumString ] += "%s,"%formulaTerms[_termIdx]
      else: sumTerms[ sumString ] = "%s,"%formulaTerms[_termIdx ]
    #Add sizeable sums as RooAdditions
    for key, value in sumTerms.iteritems(): self.modelBuilder.factory_( "sum::%s(%s)"%(key,value[:-1]) )

    #Define string for total: 1 + sizeable sums
    if( what in self.DecayScalingFunctions )&( what != "tot" ): totalStr = "sum::scaling_partial_%s(1,"%what
    else: totalStr = "sum::scaling_%s(1,"%what
    for key in sumTerms: totalStr += "%s,"%key
    totalStr = totalStr[:-1]+")"
          
    #Add scaling function as RooAddition into model
    self.modelBuilder.factory_( totalStr )

  #Function to make BR scaling functions: partial width/total width
  def makeBRScalingFunction( self, what ): self.modelBuilder.factory_( 'expr::scaling_BR_%s("@0/@1", scaling_partial_%s, scaling_tot)'%(what,what) )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Function extracting the STXS bin uncertainties
  def makeSTXSBinUncertainties( self, STXSstage="1" ):
    stxsUncertainties = {}; stxsUncertaintiesKeys = []
    for line in open( os.path.join(self.SMH.datadir, 'eft/stage%s/BinUncertainties.txt'%STXSstage) ):
      if stxsUncertaintiesKeys == []:
        stxsUncertaintiesKeys = line.split()[1:]
      else:
        fields = line.split()
        stxsUncertainties[fields[0]] = dict([(k,0.01*float(v)) for (k,v) in zip(stxsUncertaintiesKeys, fields[1:])])
    for _u in stxsUncertaintiesKeys: self.modelBuilder.doVar("param_%s[-7,7]"%_u)
    for proc in self.PROCESSES:
      #if STXS process not defined in text file: build uncertainty scaling but set equal to 1
      if proc not in stxsUncertainties: 
        self.modelBuilder.doVar("STXS%s_UncertaintyScaling_%s[1]"%(STXSstage,proc))
        continue
      else:
        #if defined...
        pnorm = ROOT.ProcessNormalization("STXS%s_UncertaintyScaling_%s"%(STXSstage,proc), "")
        for _u in stxsUncertaintiesKeys:
          var = self.modelBuilder.out.var("param_%s" %_u)
          pnorm.addLogNormal(exp(stxsUncertainties[proc][_u]),var)
        self.modelBuilder.out._import(pnorm)


#################################################################################################################
# Define inherited classes for different STXS Stages

# _________________________________________________________________________________________________________________
# STAGE 0...




# _________________________________________________________________________________________________________________
# STAGE 1...
class Stage1ToEFTModel(STXStoEFTBaseModel):
  def __init__(self):
     STXStoEFTBaseModel.__init__(self)
     from HiggsAnalysis.CombinedLimit.STXS import stage1_procs 
     self.PROCESSES = [x for v in stage1_procs.itervalues() for x in v]
     self.DECAYS = ['hzz','hbb','hww','hgg','hgluglu','hcc','tot']

  def setPhysicsOptions(self,physOptions):
    self.setPhysicsOptionsBase(physOptions)
  
  def doParametersOfInterest(self):
    if self.floatMass: print "[WARNING] Floating Higgs mass selected. STXStoEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameter list from file using textToPOIList function
    self.textToPOIList( os.path.join(self.SMH.datadir,'eft/stage1/pois.txt') )

    POIs = ','.join(self.pois.keys())
    for poi, poi_range in self.pois.iteritems(): 
      self.modelBuilder.doVar("%s%s"%(poi,poi_range))
    self.modelBuilder.doSet("POI",POIs)      
    #POIs for cWW and cB defined in terms of constraints on cWW+cB and cWW-cB: define expression for individual coefficient
    self.modelBuilder.factory_("expr::cWW_x02(\"0.5*(@0+@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.modelBuilder.factory_("expr::cB_x02(\"0.5*(@0-@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.poi_scaling['cWW'] = "0.01*cWW_x02"
    self.poi_scaling['cB'] = "0.01*cB_x02"

    #If specified in options: fix all parameters not used in LHCHXSWG-INT-2017-001 fit to 0 (freeze)
    if( self.freezeOtherParameters ):
      for poi in self.pois:
        if poi not in ['cG_x05','cA_x04','cu_x01','cHW_x02','cWWMinuscB_x02']: self.modelBuilder.out.var( poi ).setConstant(True)

    #set up model
    self.setup()

  def setup(self):
    #For additional options e.g. STXS/BR uncertainties: defined in base class
    if self.doBRU: self.SMH.makePartialWidthUncertainties()
    if self.doSTXSU: self.makeSTXSBinUncertainties( STXSstage="1" )
 
    #Read scaling functions for STXS bins and decays from text files
    self.textToSTXSScalingFunctions( os.path.join(self.SMH.datadir, 'eft/stage1/crosssections.txt') )
    self.textToDecayScalingFunctions( os.path.join(self.SMH.datadir, 'eft/stage1/decay.txt' ))
    #Make scaling for each STXS process and decay:
    for proc in self.PROCESSES: 
      if( self.fixTHandBBH ): 
        if proc in ['tHq','tHW','bbH']: self.modelBuilder.factory_("expr::scaling_%s(\"@0\",1.)"%proc)
        else: self.makeScalingFunction( proc, STXSstage="1" )
      else: self.makeScalingFunction( proc, STXSstage="1" )
    #First make partial width + total width scaling functions
    for dec in self.DECAYS: self.makeScalingFunction( dec, STXSstage="1" ) 
    #loop over decays again and make BR scaling functions: partial width/total width
    for dec in self.DECAYS:
      if dec != "tot": self.makeBRScalingFunction( dec )

  def getHiggsSignalYieldScale(self,production,decay,energy):
    name = "stxs1toeft_scaling_%s_%s_%s"%(production,decay,energy)
    if self.modelBuilder.out.function(name) == None:

      #Check production scaling exists:
      XSscal = None
      if production in self.PROCESSES: 
        XSscal = "scaling_%s"%(production)
      else:
        raise ValueError("[ERROR] Process %s is not supported in Stage 1 Model"%production)
      #Check decay scaling exists
      BRscal = None
      if decay in self.DECAYS:
        BRscal = "scaling_BR_%s"%decay
      else:
        raise ValueError("[ERROR] Decay %d is not supported"%decay)

      #For including BR and STXS bin uncertainties
      if( self.doSTXSU )&( self.doBRU ): 
        THUscaler = "uncertainty_scaling_%s_%s"%(production,decay)
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s_%s(\"@0*@1\",STXS1_UncertaintyScaling_%s,HiggsDecayWidth_UncertaintyScaling_%s)'%(production,decay,production,decay))
      elif( self.doSTXSU ): 
        THUscaler = "uncertainty_scaling_%s"%production
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",STXS1_UncertaintyScaling_%s)'%(production,production))
      elif( self.doBRU ): 
        THUscaler = "uncertainty_scaling_%s"%decay
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",HiggsDecayWidth_UncertaintyScaling_%s)'%(decay,decay))

      #Combine XS and BR scaling: including theory uncertainties if option selected
      if( self.doSTXSU )|( self.doBRU ): self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal,THUscaler])))
      else: self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))

    return name

# _________________________________________________________________________________________________________________
# STAGE 1.1...


#################################################################################################################
Stage1toEFT = Stage1ToEFTModel()

