# Author: Jonathon Langford (ICL)
# Date: 06/02/2019
# Description: Model to describe how bins in STXS (0,1,1.1) scale using full set of dimension-6 EFT parameters
#              Equations calculated using Higgs Effective Lagrangian (Madgraph)

from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from math import exp
import ROOT, os, re, sys

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
    self.PROCESSES = {}
    self.DECAYS = []
    # Dicts to store pois + scaling functions
    self.pois = {}
    self.poi_scaling = {}
    self.STXSScalingFunctions = {}
    self.DecayScalingFunctions = {}
    # Options
    self.floatMass = False #Initally false, require external option to float mass
    self.doSTXSU = STXSU
    self.doBRU = BRU
    self.fixTHandBBH = fixTHandBBH #if false scaling function for tH, bbH MUST be defined in input txt file
    self.freezeOtherParameters = freezeOtherParameters #Option to freeze majority of parameters in model. Leaving those used in LHCHXSWG-INT-2017-001 fit to float
    self.fixProcesses = fixProcesses #Option to fix certain STXS bins: comma separated list of STXS bins
    self.useLHCHXSWGNumbers=False
    self.useExtendedVBFScheme=False
    self.linearOnly=False
    if self.freezeOtherParameters: 
      self.parametersOfInterest = ['cG','cA','cWWMinuscB','cWWPluscB','cHW','cHB','cu','cd','cl'] # note cWW+cB is frozen, but required to define cWW and cB
      self.distinctParametersOfInterest = set([])
      for p in self.parametersOfInterest:
        if "Plus" in p: self.distinctParametersOfInterest = self.distinctParametersOfInterest | set(p.split("Plus"))
        elif "Minus" in p: self.distinctParametersOfInterest = self.distinctParametersOfInterest | set(p.split("Minus"))
        else: self.distinctParametersOfInterest.add(p)

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
      if po.startswith("useExtendedVBFScheme="): 
        self.useExtendedVBFScheme = (po.replace("useExtendedVBFScheme=","") in ["yes","1","Yes","True","true"])
      if po.startswith("useLHCHXSWGStage1="): 
        self.useLHCHXSWGNumbers = (po.replace("useLHCHXSWGStage1=","") in ["yes","1","Yes","True","true"])
      if po.startswith("linearOnly="): 
        self.linearOnly = (po.replace("linearOnly=","") in ["yes","1","Yes","True","true"])

    #Output options to screen
    print " --> [STXStoEFT] Theory uncertainties in partial widths: %s"%self.doBRU
    print " --> [STXStoEFT] Theory uncertainties in STXS bins: %s"%self.doSTXSU
    if( self.doSTXSU ): print " --> [WARNING]: theory uncertainties in STXS bins are currently incorrect. Need to update: data/lhc-hxswg/eft/HEL/*_binuncertainties.txt"
    if( self.freezeOtherParameters ): print " --> [STXStoEFT] Freezing all but [cG,cu,cd,cHW,cHB,cA,cWWMinuscB,cl] (distinct: %s) to 0"%(self.distinctParametersOfInterest)
    else: print " --> [STXStoEFT] Allowing all HEL parameters to float"
    if( len( self.fixProcesses ) > 0 ): print " --> [STXStoEFT] Fixing following processes to SM: %s"%self.fixProcesses
    if self.useExtendedVBFScheme: print " --> [STXStoEFT] Use extended VBF scheme (different scaling for pure + V(qq)H)"
    if self.useLHCHXSWGNumbers: print " --> [STXStoEFT] Using LHCHXSWG numbers for stage 1 scaling functions"
    if self.linearOnly: print " --> [STXStoEFT] Only linear terms (Ai)"

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
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    for line in lines[skipRows:]:
      if( len(line.strip())== 0 ): continue

      # If param is frozen (except cWW+cB) then skip poi
      skip_poi = False
      if self.freezeOtherParameters:
        if line.split()[0].split("_")[0] not in self.parametersOfInterest: skip_poi = True
      if not skip_poi:

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
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    # Loop over lines in txt file
    for line in lines[skipRows:]:      
      # Extract STXS bin and scaling function
      stxs_bin_, function = line.split(":")
      # Format function: remove spaces
      function = function.replace('\n','')
      function = "".join(function.split(" "))
      # Add to dict
      self.STXSScalingFunctions[stxs_bin_] = function
    # Close file
    file_in.close()
  
  #Function to extract decay scaling functions from file
  def textToDecayScalingFunctions( self, filename, skipRows=0 ):
    file_in = open(filename,"r")
    lines = [l for l in file_in]
    # Loop over lines in txt file
    for line in lines[skipRows:]:
      # Extract decay mode and scaling function
      decay_, function = line.split(":")
      # Format function: remove spaces
      function = function.replace('\n','')
      function = "".join(function.split(" "))
      # Add to dict
      self.DecayScalingFunctions[decay_] = function
    # Close file
    file_in.close()
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Function to make scaling function in workspace from string
  #   > defines each terms as RooProduct
  #   > sum of terms as RooAddition
  def makeScalingFunction( self, what, STXSstage="" ):
  
    #if in processes/decays extract formula from corresponding dict
    if what in self.STXSScalingFunctions: formula = self.STXSScalingFunctions[ what ]
    elif what in self.DecayScalingFunctions: formula = self.DecayScalingFunctions[ what ]
    else:
      if STXSstage != "": raise ValueError("[ERROR] Scaling function for %s does not exist for STXS Stage %s"%(what,STXSstage))
      else: raise ValueError("[ERROR] Scaling function for %s does not exist"%what)

    #replace "-" in formula string by "+-" and then turn into list, splitting by delimeter "+"
    formula = re.sub("-","+-",formula).split("+")
    
    #define list to hold name of terms 
    formulaTerms = []
    _termIdx = 0

    #loop over terms in formula (ignoring first term which is just 1 [perturbation theory])
    for term in formula[1:]:
      constituents = term.split("*")
      # If param is frozen then skip term
      skipTerm = False
      if self.freezeOtherParameters:
        for c in constituents[1:]: #Ignore pre-factor
          if c not in self.distinctParametersOfInterest: skipTerm = True
      # If set: only include linear terms
      if self.linearOnly:
        if len(constituents) > 2: skipTerm = True
      if not skipTerm:
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
      if _termIdx % 15 == 0: _sumIdx += 1
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
    if STXSstage == "all": stage_list = ['0','1','1_1','fixedproc']
    else: stage_list = [STXSstage, 'fixedproc']
    storedUncertainties = []
    for s in stage_list:
      if s == "fixedproc": key=s
      else: key = "stage%s"%s
      stxsUncertainties = {}; stxsUncertaintiesKeys = []
      for line in open( os.path.join(self.SMH.datadir, 'eft/HEL/%s_binuncertainties.txt'%key) ):
        if stxsUncertaintiesKeys == []:
          stxsUncertaintiesKeys = line.split()[1:]
        else:
          fields = line.split()
          stxsUncertainties[fields[0]] = dict([(k,0.01*float(v)) for (k,v) in zip(stxsUncertaintiesKeys, fields[1:])])
      for _u in stxsUncertaintiesKeys: 
        if _u not in storedUncertainties:
          self.modelBuilder.doVar("param_%s[-7,7]"%_u)
          storedUncertainties.append(_u)
      for proc in self.PROCESSES[key]:
        # if uncertainties not defined for given process: build scaling but set = 1
        if proc not in stxsUncertainties:
          self.modelBuilder.doVar("%s_UncertaintyScaling_%s[1]"%(key,proc))
          continue
        else:
          pnorm = ROOT.ProcessNormalization("%s_UncertaintyScaling_%s"%(key,proc), "")
          for _u in stxsUncertaintiesKeys:
            var = self.modelBuilder.out.var("param_%s" %_u)
            pnorm.addLogNormal(exp(stxsUncertainties[proc][_u]),var)
          self.modelBuilder.out._import(pnorm)
    
#################################################################################################################
# Define inherited classes: AllStageToEFT and StageXToEFT

# Combination of different stages
class AllStagesToEFTModel(STXStoEFTBaseModel):
  def __init__(self):
    STXStoEFTBaseModel.__init__(self)
    from HiggsAnalysis.CombinedLimit.STXS import fixed_procs, stage0_procs, stage1_procs, stage1_1_procs
    #self.PROCESSES = {}
    for s in ['0','1','1_1']:
      if s=='0': self.PROCESSES["stage%s"%s] = [x for v in stage0_procs.itervalues() for x in v] 
      elif s=='1': self.PROCESSES["stage%s"%s] = [x for v in stage1_procs.itervalues() for x in v] 
      else: self.PROCESSES["stage%s"%s] = [x for v in stage1_1_procs.itervalues() for x in v] 
    self.PROCESSES["fixedproc"] = fixed_procs
    self.DECAYS = ['hzz','hbb','htt','hww','hgg','hgluglu','hcc','hzg','hmm','tot']

  def setPhysicsOptions(self,physOptions):
    self.setPhysicsOptionsBase(physOptions)
  
  def doParametersOfInterest(self):
    if self.floatMass: print " --> [WARNING] Floating Higgs mass selected. STXStoEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameter list from file using textToPOIList function
    self.textToPOIList( os.path.join(self.SMH.datadir,'eft/HEL/pois.txt') )
    POIs = ','.join(self.pois.keys())
    for poi, poi_range in self.pois.iteritems(): self.modelBuilder.doVar("%s%s"%(poi,poi_range))
    # Remove cWW+cB from POI list if freezing other parameters
    if self.freezeOtherParameters: POIs = re.sub("cWWPluscB_x02,","",POIs)
    self.modelBuilder.doSet("POI",POIs)
    #POIs for cWW and cB defined in terms of constraints on cWW+cB and cWW-cB: define expression for individual coefficient
    self.modelBuilder.factory_("expr::cWW_x02(\"0.5*(@0+@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.modelBuilder.factory_("expr::cB_x02(\"0.5*(@0-@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.poi_scaling['cWW'] = "0.01*cWW_x02"
    self.poi_scaling['cB'] = "0.01*cB_x02"

    # Freeze cWW+cB if freezeOtherParameters
    if self.freezeOtherParameters: 
      for poi in self.pois: 
        if 'cWWPluscB' in poi: self.modelBuilder.out.var( poi ).setConstant(True)

    #set up model
    self.setup()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def setup(self):
    #For additional options e.g. STXS/BR uncertainties: defined in base class
    if self.doBRU: self.SMH.makePartialWidthUncertainties()
    if self.doSTXSU: self.makeSTXSBinUncertainties( STXSstage="all" )
 
    # Make scaling functions for each STXS process
    stored = []
    for s in ['1','0','1_1']:
      stage = "stage%s"%s
      # Read scaling functions for STXS bins from txt file
      if s=="1" and self.useLHCHXSWGNumbers: self.textToSTXSScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/%s_xs_LHCHXSWG-INT-2017-001.txt'%stage) )
      else: 
        if self.useExtendedVBFScheme: self.textToSTXSScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/%s_xs_extended_vbf.txt'%stage) )
        else: self.textToSTXSScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/%s_xs.txt'%stage) )
      for proc in self.PROCESSES[stage]: 
        if proc not in stored:
          self.makeScalingFunction( proc, STXSstage=s )
          stored.append( proc )

    # Make dummy scaling (=1) for fixed procs
    for proc in self.PROCESSES["fixedproc"]: self.modelBuilder.factory_("expr::scaling_%s(\"@0\",1.)"%proc)

    # Read scaling functions for decays from txt files
    self.textToDecayScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/decay.txt' ) )

    # Make partial width + total width scaling functions
    for dec in self.DECAYS: self.makeScalingFunction( dec )
    # Make BR scaling functions: partial width/total width
    for dec in self.DECAYS:
      if dec != "tot": self.makeBRScalingFunction( dec )
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getHiggsSignalYieldScale(self,production,decay,energy):

    # Function to convert troublesome procs into viable one for HC combination
    production = convert_to_STXS(production,decay)

    name = "stxstoeft_scaling_%s_%s_%s"%(production,decay,energy)
    if self.modelBuilder.out.function(name) == None:

      XSscal = None
      BRscal = None

      # Extract STXS stage process belongs to: in descreasing order as want most recent th. unc
      if production in self.PROCESSES['stage1_1']: key = "stage1_1"
      elif production in self.PROCESSES['stage1']: key = "stage1"
      elif production in self.PROCESSES['stage0']: key = "stage0"
      elif production in self.PROCESSES['fixedproc']: key = "fixedproc" 
      # For fwd production: set to SM
      elif "fwd" in production:
        self.modelBuilder.factory_("expr::scaling_%s(\"@0\",1.)"%production)
        key = None
      else:
        print "[WARNING] Process %s is not supported in STXStoEFT Model, Setting to 1"%production
        return 1
        #raise ValueError("[ERROR] Process %s is not supported in STXStoEFT Model"%production)
      # Give production correct scaling
      XSscal = "scaling_%s"%production
     
      # Check decay scaling exists:
      if decay in self.DECAYS: BRscal = "scaling_BR_%s"%decay
      else:
        print "[WARNING] Decay %s is not supported in STXStoEFT Model, setting to 1"%decay
        return 1
        #raise ValueError("[ERROR] Decay %s is not supported in STXStoEFT Model"%decay)

      # Uncertainty scaling: BR and STXS bin uncertainties
      if( self.doSTXSU )&( self.doBRU ):
        if key==None:
          THUscaler = "uncertainty_scaling_%s"%(decay)
          self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",HiggsDecayWidth_UncertaintyScaling_%s)'%(decay,decay))
        else:
          THUscaler = "uncertainty_scaling_%s_%s"%(production,decay)
          self.modelBuilder.factory_('expr::uncertainty_scaling_%s_%s(\"@0*@1\",%s_UncertaintyScaling_%s,HiggsDecayWidth_UncertaintyScaling_%s)'%(production,decay,key,production,decay))
      elif( self.doSTXSU ):
        if key==None:
          THUscaler = "uncertainty_scaling_dummy"
          self.modelBuilder.factory_('expr::uncertainty_scaling_dummy(\"@0\",1.)')
        else:
          THUscaler = "uncertainty_scaling_%s"%production
          self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",%s_UncertaintyScaling_%s)'%(production,key,production))
      elif( self.doBRU ):
        THUscaler = "uncertainty_scaling_%s"%decay
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",HiggsDecayWidth_UncertaintyScaling_%s)'%(decay,decay))
      
      #Combine XS and BR scaling: incuding theory unc if option selected
      if( self.doSTXSU )|( self.doBRU ): self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal,THUscaler])))
      else: self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))
      
    return name
 
# _________________________________________________________________________________________________________________
# Single stage to EFT model
class StageXToEFTModel(STXStoEFTBaseModel):
  def __init__(self,stage):
    STXStoEFTBaseModel.__init__(self)
    self.stage = stage
    import HiggsAnalysis.CombinedLimit.STXS as STXS
    self.PROCESSES['stage%s'%self.stage] = [x for v in getattr(STXS,"stage%s_procs"%self.stage).itervalues() for x in v] 
    self.PROCESSES["fixedproc"] = STXS.fixed_procs
    self.DECAYS = ['hzz','hbb','htt','hww','hgg','hgluglu','hcc','hzg','hmm','tot']

  def setPhysicsOptions(self,physOptions):
    self.setPhysicsOptionsBase(physOptions)
  
  def doParametersOfInterest(self):
    if self.floatMass: print " --> [WARNING] Floating Higgs mass selected. STXStoEFT model assumes MH=125.0 GeV"
    self.doMH()
    self.SMH = SMHiggsBuilder(self.modelBuilder)
    
    #Read in parameter list from file using textToPOIList function
    self.textToPOIList( os.path.join(self.SMH.datadir,'eft/HEL/pois.txt') )
    POIs = ','.join(self.pois.keys())
    for poi, poi_range in self.pois.iteritems(): 
      self.modelBuilder.doVar("%s%s"%(poi,poi_range))
    self.modelBuilder.doSet("POI",POIs)
    #POIs for cWW and cB defined in terms of constraints on cWW+cB and cWW-cB: define expression for individual coefficient
    self.modelBuilder.factory_("expr::cWW_x02(\"0.5*(@0+@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.modelBuilder.factory_("expr::cB_x02(\"0.5*(@0-@1)\",cWWPluscB_x02,cWWMinuscB_x02)")
    self.poi_scaling['cWW'] = "0.01*cWW_x02"
    self.poi_scaling['cB'] = "0.01*cB_x02"
    
    # Freeze cWW+cB if freezeOtherParameters
    if self.freezeOtherParameters: 
      for poi in self.pois: 
        if 'cWWPluscB' in poi: self.modelBuilder.out.var( poi ).setConstant(True)

    #set up model
    self.setup()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def setup(self):
    #For additional options e.g. STXS/BR uncertainties: defined in base class
    if self.doBRU: self.SMH.makePartialWidthUncertainties()
    if self.doSTXSU: self.makeSTXSBinUncertainties( STXSstage=self.stage )

    # Read scaling functions for STXS bins and decays from txt files
    self.textToSTXSScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/stage%s_xs.txt'%self.stage) )
    self.textToDecayScalingFunctions( os.path.join(self.SMH.datadir, 'eft/HEL/decay.txt' ) )

    # Make scaling functions for STXS processes
    for proc in self.PROCESSES["stage%s"%self.stage]: self.makeScalingFunction( proc, STXSstage=self.stage )

    # Make dummy scaling (=1) for fixed procs
    for proc in self.PROCESSES["fixedproc"]: self.modelBuilder.factory_("expr::scaling_%s(\"@0\",1.)"%proc)

    # Make partial width + total width scaling functions
    for dec in self.DECAYS: self.makeScalingFunction( dec )
    # Make BR scaling functions: partial width/total width
    for dec in self.DECAYS:
      if dec != "tot": self.makeBRScalingFunction( dec )
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getHiggsSignalYieldScale(self,production,decay,energy):
    name = "stxstoeft_scaling_%s_%s_%s"%(production,decay,energy)
    if self.modelBuilder.out.function(name) == None:

      XSscal = None
      BRscal = None

      # Extract STXS stage process belongs to: in descreasing order as want most recent th. unc
      if production in self.PROCESSES['stage%s'%self.stage]: key = "stage%s"%self.stage
      elif production in self.PROCESSES['fixedproc']: key = "fixedproc" 
      else:
        raise ValueError("[ERROR] Process %s is not supported in STXStoEFT Model (stage %s)"%(production,self.stage))
      # Give production correct scaling
      XSscal = "scaling_%s"%production
     
      # Check decay scaling exists:
      if decay in self.DECAYS: BRscal = "scaling_BR_%s"%decay
      else:
        raise ValueError("[ERROR] Decay %s is not supported in STXStoEFT Model"%decay)

      # Uncertainty scaling: BR and STXS bin uncertainties
      if( self.doSTXSU )&( self.doBRU ):
        THUscaler = "uncertainty_scaling_%s_%s"%(production,decay)
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s_%s(\"@0*@1\",%s_UncertaintyScaling_%s,HiggsDecayWidth_UncertaintyScaling_%s)'%(production,decay,key,production,decay))
      elif( self.doSTXSU ):
        THUscaler = "uncertainty_scaling_%s"%production
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",%s_UncertaintyScaling_%s)'%(production,key,production))
      elif( self.doBRU ):
        THUscaler = "uncertainty_scaling_%s"%decay
        self.modelBuilder.factory_('expr::uncertainty_scaling_%s(\"@0\",HiggsDecayWidth_UncertaintyScaling_%s)'%(decay,decay))
      
      #Combine XS and BR scaling: incuding theory unc if option selected
      if( self.doSTXSU )|( self.doBRU ): self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal,THUscaler])))
      else: self.modelBuilder.factory_("prod::%s(%s)"%(name,",".join([XSscal,BRscal])))
      
    return name

#################################################################################################################
# Function to convert troublesome procs to those defined in STXS
def convert_to_STXS( _production, _decay ):
  # For hmm: split into W+/W- -> remove from tag
  _production = re.sub("Plus","",_production)
  _production = re.sub("Minus","",_production)

  # Some VH signal processes missing had/lep: for now assume all lep
  if _production in ['WH','ZH','ggZH']: _production = "%s_lep"%_production

  # hzz: NIGHTMARE
  if _decay=="hzz": 
    # VBFori: set to qqH
    if "VBFori" in _production: _production = "qqH"

    # Not split VH: assume all WH
    _production = re.sub("VH_","WH_", _production)
    # VH: PTV_GT150 is not an STXS bin, set to PTV_15_250_0J (majority of events)
    _production = re.sub("GT150","150_250_0J", _production)

    # ggH: 0J high PT wrong name
    _production = re.sub("ggH_0J_PTH_10_200","ggH_0J_PTH_GT10", _production)
    # ggH: 2J should be GE2J
    _production = re.sub("ggH_2J","ggH_GE2J", _production)
    # ggH: havent split VBFTOPO, choose one with larger cross section, 3J
    _production = re.sub("ggH_VBF","ggH_VBFTOPO_JET3", _production)

    # qqH: names all wrong, combination of 1 and 1.1
    _production = re.sub("qqH_2j_mjj_350_700_2j","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25", _production)
    _production = re.sub("qqH_2j_mjj_GT350_3j", "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25", _production)
    _production = re.sub("qqH_2j_mjj_GT700_2j", "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25", _production)
    _production = re.sub("qqH_GT200_VHori", "WH_had_PTJET1_GT200", _production)
    _production = re.sub("qqH_GT200", "qqH_PTJET1_GT200", _production)
    _production = re.sub("qqH_Rest_VHori", "WH_had_REST", _production)
    _production = re.sub("qqH_Rest", "qqH_REST", _production)

    # assume tH is all tHQ
    if _production == "tH": _production = "tHq"

  # No binnning in ttH up to stage 1.1 so can remove any X in ttH_decay 
  if "ttH" in _production: _production = "ttH"
    #_production = re.sub("_faketau","",_production)
    #_production = re.sub("_gentau","",_production)
    
  return _production


#################################################################################################################
# Instantiation of STXStoEFT models
Stage0toEFT = StageXToEFTModel("0")
Stage1toEFT = StageXToEFTModel("1")
Stage1_1toEFT = StageXToEFTModel("1_1")
AllStagestoEFT = AllStagesToEFTModel()

