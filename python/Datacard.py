class Datacard():
    """
    Description:

    This is a container class that is filled by the function parseCards in HiggsAnalysis/CombinedLimit/python/DatacardParser.py.
    This function parses a given datacards file and fills Datacard as a data structure that is easily accessible within python
    scripts. To simplify access a set of getter c++ like function have been added to the class. Some more access function for
    systematics will be added later.
    """
    def __init__(self):
        ## list of [bins in datacard]
        self.bins = []
        ## dict of {bin : number of observed events}
        self.obs  = {}
        ## list of [processes]
        self.processes = []
        ## list of [signal processes]
        self.signals = []
        ## dict of {processes : boolean to indicate whether process is signal or not}
        self.isSignal = {}
        ## list of [(bin, process, boolean to indicate whether process is signal or not)]
        self.keyline = []
        ## dict of {bin : {process : yield}}
        self.exp     = {}
        ## list of [(name of uncert, boolean to indicate whether to float this nuisance or not, type, list of what additional arguments (e.g. for gmN), keyline element)]
        self.systs   = []
        ## list of [{bin : {process : [input file, path to shape, path to shape for uncertainty]}}]
        self.shapeMap = {}
        ## boolean that indicates whether the datacard contains shapes or not
        self.hasShapes = False
        ## dirct of {name of uncert, boolean to indicate whether it is a flat parametric uncertainty or not}
        self.flatParamNuisances = {}
	      ## dict of rateParam
        self.rateParams = {}
	      ## dict of extArgs
        self.extArgs = {}
	      ## maintain the order for rate modifiers
        self.rateParamsOrder = set()
        ## dirct of {name of uncert, boolean to indicate whether this nuisance is floating or not}
        self.frozenNuisances = set()

	# Allows for nuisance renaming of "shape" systematics
	self.systematicsShapeMap = {}

	# Allows for nuisance renaming of "param" systematics
	self.systematicsParamMap = {}

	# Allow to pick out entry in self.systs.
	self.systIDMap = {}

        # Keep edits
	self.nuisanceEditLines = []

        # map of which bins should have automated Barlow-Beeston parameters
	self.binParFlags = {}

	self.groups = {}
	self.discretes = []

    def print_structure(self):
	"""
	Print the contents of the -> should allow for direct text2workspace on python config
	"""
	print """
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ModelTools import *
from HiggsAnalysis.CombinedLimit.ShapeTools import *
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

from sys import exit
from optparse import OptionParser
parser = OptionParser()
addDatacardParserOptions(parser)
options,args = parser.parse_args()
options.bin = True # make a binary workspace

DC = Datacard()
MB = None

############## Setup the datacard (must be filled in) ###########################
	"""

	print "DC.bins = 	"		, self.bins			,"#",type(self.bins)
	print "DC.obs = 	"		, self.obs                      ,"#",type(self.obs)
	print "DC.processes = 	"		, self.processes                ,"#",type(self.processes)
	print "DC.signals = 	"		, self.signals                  ,"#",type(self.signals)
	print "DC.isSignal = 	"		, self.isSignal                 ,"#",type(self.isSignal)
	print "DC.keyline = 	"		, self.keyline                  ,"#",type(self.keyline)
	print "DC.exp = 	"		, self.exp                      ,"#",type(self.exp)
	print "DC.systs = 	"		, self.systs                    ,"#",type(self.systs)
	print "DC.shapeMap = 	"		, self.shapeMap                 ,"#",type(self.shapeMap)
	print "DC.hasShapes = 	"		, self.hasShapes                ,"#",type(self.hasShapes)
	print "DC.flatParamNuisances = "	, self.flatParamNuisances       ,"#",type(self.flatParamNuisances)
	print "DC.rateParams = "		, self.rateParams               ,"#",type(self.rateParams)
	print "DC.extArgs = 	"		, self.extArgs                  ,"#",type(self.extArgs)
	print "DC.rateParamsOrder 	= "	, self.rateParamsOrder          ,"#",type(self.rateParamsOrder)
	print "DC.frozenNuisances 	= "	, self.frozenNuisances          ,"#",type(self.frozenNuisances)
	print "DC.systematicsShapeMap = "	, self.systematicsShapeMap      ,"#",type(self.systematicsShapeMap)
	print "DC.systematicsParamMap = "	, self.systematicsParamMap      ,'#',type(self.systematicsParamMap)
	print "DC.systIDMap = "			, self.systIDMap		,'#',type(self.systIDMap)
	print "DC.nuisanceEditLines 	= "	, self.nuisanceEditLines        ,"#",type(self.nuisanceEditLines)
  	print "DC.binParFlags 	= "	  	, self.binParFlags        	,"#",type(self.binParFlags)
	print "DC.groups 	= "		, self.groups        		,"#",type(self.groups)
	print "DC.discretes 	= "		, self.discretes        	,"#",type(self.discretes)

	print """

###### User defined options #############################################

options.out 	 = "combine_workspace.root"  	# Output workspace name
options.fileName = "./" 			# Path to input ROOT files
options.verbose  = "1" 				# Verbosity

##########################################################################

if DC.hasShapes:
    MB = ShapeBuilder(DC, options)
else:
    MB = CountingModelBuilder(DC, options)

# Set physics models
MB.setPhysics(defaultModel)
MB.doModel()
	"""

	# map of which bins should have automated Barlow-Beeston parameters
	self.binParFlags = {}

    def list_of_bins(self) :
        """
        Return the list of all bins in the datacard.
        """
        return self.bins

    def list_of_procs(self, type='') :
        """
        Return the list of all processes in the datacard. For type='' all processes in the datacards are returned. For
        type='s' only signal processes will be returned. For type='b' only background processes will be returned.
        """
        if type == 's' :
            return self.signals
        elif type == 'b' :
            bgs = []
            for (proc, sig) in self.isSignal.iteritems() :
                if not sig : bgs.append(proc)
            return bgs
        else :
            return self.processes

    def list_of_signals(self) :
        """
        Return the list of signal processes in the datacards.
        """
        return self.list_of_procs('s')

    def list_of_backgrounds(self) :
        """
        Return the list of background processes in the datacards.
        """
        return self.list_of_procs('b')

    def barcode(self, bin, proc, idx) :
        """
        Return the barcode where to find corresponding shapes for a given bin and process. The barcode is a list of
        three elements: [0] path to the file that contains the the shapes, [1] path to the shape histogram in the
        input file, [2] path to the uncertainty histograms in the input file. The elements of the list are returned
        as obtained from the datacard. The function actually does not return the full list but the first second or
        thrid element as indicated by idx. If there is no entry for a given bin, process and idx an empty list is
        returned.
        """
        path = ''
        if not bin in self.shapeMap.keys() :
            if '*' in self.shapeMap.keys() :
                if not proc in self.shapeMap['*'] :
                    if '*' in self.shapeMap['*'].keys() :
                        path = self.shapeMap['*']['*'][idx]
                else :
                    path = self.shapeMap['*'][proc][idx]
        else :
            path = self.shapeMap[bin][proc][idx]
        return path

    def path_to_file(self, bin, proc) :
        """
        Return the path to the root input files for shape analyses for a given bin and process. If there is no such
        path an empty string will be returned.
        """
        return self.barcode(bin, proc, 0)

    def path_to_shape(self, bin, proc, resolve=True) :
        """
        Return the path to the histogram within the root input file for shape analyses for a given bin and process. If
        there is no such path an empty string will be returned. The path will be returned as obtained from the datacard.
        If the option resolve is True the keywords $CHANNEL correp. to bin and $PROCESS corresp. to proc are replaced
        by their corresponding real values. Other potential keywords remain unchanged.
        """
        return self.barcode(bin, proc, 1).replace('$CHANNEL',bin).replace('$PROCESS',proc) if resolve else self.barcode(bin, proc, 1)

    def shape(self, bin, proc, resolve) :
        """
        This is a safe way to determine the histogram name of the shape template for a given bin and proc.
        """
        path =  self.path_to_shape(bin, proc, resolve)
        if '/' in path :
            return path[path.rfind('/')+1:]
        else :
            return path

    def path_to_uncert(self, bin, proc, resolve=True) :
        """
        Return the path to the uncertainty histograms within the root input file for shape analyses for a given bin and
        process. If there is no such path an empty string will be returned. The path will be returned as obtained from
        the datacard. If option resolve is True the keywords $CHANNEL correp. to bin and $PROCESS corresp. to proc are
        replaced by their corresponding real values. Other potential keywords remain unchanged.
        """
        return self.barcode(bin, proc, 2).replace('$CHANNEL',bin).replace('$PROCESS',proc) if resolve else self.barcode(bin, proc, 2)

    def uncert(self, bin, proc, resolve) :
        """
        This is a safe way to determine the histogram name of the uncert template for a given bin and proc.
        """
        path =  self.path_to_uncert(bin, proc, contact)
        if '/' in path :
            return path[path.rfind('/')+1:]
        else :
            return path

    def obs(self, bin) :
        """
        Return the number of observed events for a given bin.
        """
        return self.obs[bin]

    def rate(self, bin, proc) :
        """
        Return thenumber of expected events for a given bin and process.
        """
        return self.exp[bin][proc]

    def getAllVariables(self):
    	"""
	Return all variables defined in the datacard
	"""
	allVars = tuple([syst[0] for syst in self.systs]+self.flatParamNuisances.keys()+self.extArgs.keys()+self.discretes)
	for rp in self.rateParams:
	  modifiers = self.rateParams[rp]
	  for p in modifiers : allVars+=tuple([p[0][0]])

	return list(set(allVars))

    def add_syst_id(self,lsyst):
        if lsyst in self.systIDMap.keys(): self.systIDMap[lsyst].append(len(self.systs)-1)
        else: self.systIDMap[lsyst] = [len(self.systs)-1]
    
    def renameNuisanceParameter(self,oldname,newname,process_list=[],channel_list=[]): 
	"""
	Rename nuisance parameter from oldname to newname

	This will by default rename the parameter globally. However, 
	if you only want to modify the name of the nuisance parameter 
	for specific channels/processes, then you should specify a 
	process (list or leave empty for all) and channel (list or leave empty for all)
	"""
	newnameExists = False
	existingclashes = {}
        for lsyst,nofloat,pdf0,args0,errline0 in (self.systs[:]):
	  if lsyst == newname : # found the nuisance exists
	    newnameExists = True 
	    existingclashes[lsyst]=(nofloat,pdf0,args0,errline0)

	found=False 
        nuisanceID = i = -1
	
	if process_list == [] : process_list = self.processes 
	if channel_list == [] : channel_list = self.bins

        for lsyst,nofloat,pdf0,args0,errline0 in (self.systs[:]):
	  i+=1
	  if lsyst == oldname : # found the nuisance
	    nuisanceID = i
	    found = True
	    # check if the new name exists 
	    if lsyst in existingclashes.keys(): 
	      nofloat1,pdf1,args1,errline1 = existingclashes[lsyst]
		
	    else: pass
	    datacard.systs[nuisanceID][0]=newname 
	    if "shape" in pdf0 :       
	      for b in channel_list: 
	        for p in process_list :  
		  datacard.systematicsShapeMap[newname,b,p]=oldname
	    if "param" in pdf0 : datacard.systematicsParamMap[oldname]=newname

        if not found: raise RuntimeError, "Error: no parameter found with = %s\n" % oldname
	return 0 
