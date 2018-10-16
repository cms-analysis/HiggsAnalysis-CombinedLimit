#Script to plot XS*BR as a function of POIs (kappas/lambdas) 
# Run with "python test/scanScalefactorsFromWS.py -M c7 -o outfolder -m 125.09 
#   Note that most LHCHCG models use the c7 base. This is just the name prefix of the scaling function 
#   you may need to use eg CvCf, c6 etc for older modelts in the CMS HCG repo - if in doublt look for the XSBR functions in the WS
# for A1/B1 (mu based) models, the output of t2w will show how the model scales
import os,sys,numpy,array
import itertools

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] file \nrun with --help to get list of options")
parser.add_option("-M","--Model",dest="model",default="DEFAULT",type='str',help="Name of model used in model builder (note its the name of the scaling functions, eg the L1/K models use 'c7'")
parser.add_option("-m","--mh",dest="mh",default=125.09,type='float',help="Lightest Higgs Mass")
#parser.add_option("-e","--energydependance",dest="energydependant",action='store_true')
parser.add_option("-s","--step",dest="stepsize",type='float',default=0.5)
parser.add_option("","--scaling_prefix",dest="scaling_prefix",type='str',default="XSBRscal", help="Prefix for scaler, eg usually this is XSBRscal")
parser.add_option("","--slice",dest="sliceval",type='str', default = "")
parser.add_option("","--energies",dest="energies",type='str', default = "13TeV", help= "A comma separated list of the energies to show expressions for")
parser.add_option("","--skipTree", default = False, action='store_true', help= "Skip filling the tree of points in the ND space (slow)")
parser.add_option("-o","--out",dest="out",type='str', default = "", help= "Output folder for the plots/trees etc")
parser.add_option("-P","--prod",dest="prod",type='str', default = "", help= "Overwrite default production channels (comma separated list in a string)")
parser.add_option("-D","--decay",dest="decay",type='str', default = "", help= "Overwrite default decay channels (comma separated list in a string)")
(options,args)=parser.parse_args()

import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gROOT.SetBatch(1)

if not options.out: options.out = "."
else: os.system("mkdir -p %s"%options.out)

#ROOT.gStyle.SetNumberContours(255)
# Can do full matrix of prod*decay
# In the future we should really think anout just importing these from the SMHiggs Builder!
energies = list(options.energies.split(","))
production_channels 	= ["ggH","qqH","WH","ZH","ttH","ggZH","tHq","tHW","bbH"]
decay_channels 		= ['hww','hzz','hgg','hbb','hcc','htt','hmm','hzg','hgluglu','hinv']
decay_modes 		= ["hgg","hvv","hff"]
if options.prod: 
	production_channels = options.prod.split(",")
	print "Will look for productions ", production_channels
if options.decay: 
	decay_channels      = options.decay.split(",")
	print "Will look for decays ", decay_channels

# Stepping of parameters in the model 
step = options.stepsize

# A way to keep x,y,z parameters for slice plots
config ={}

# Coloured plots
def set_palette(ncontours=999):

    # default palette, looks cool
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [0.00, 0.00, 0.87, 1.00, 0.51]
    green = [0.00, 0.81, 1.00, 0.20, 0.00]
    blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array.array('d', stops)
    r = array.array('d', red)
    g = array.array('d', green)
    b = array.array('d', blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)


# Function which add a Double_t branch for each entry in a list
def createBranches(tree, params):
    for p in params.keys():
	tree.Branch(p,params[p],"%s/Double_t"%(p))

def resetVals():

  # set the parameters 
  it = params.createIterator()
  for j in range(nparams):
	p = it.Next()
	if p==None : break		
	p.setVal(default_parameter_vals[p.GetName()])

def fillGrid(func,graph,txtfile,tree,c_vals):

	point = 0
	for val in c_vals:

	  # set the parameters 
	  it = params.createIterator()
	  for j in range(nparams):
		p = it.Next()
		if p==None : break		
		p.setVal(val[j])
		parameter_vals[p.GetName()][0]=val[j]
		txtfile.write("%1.2f   "%(val[j]))

	  mu = func.getVal()

	  if nparams == 2: 
	 	graph.SetPoint(point,val[0],val[1],mu)
		point+=1
	  elif options.sliceval:
	    if abs(config["sliceval"] - val[config["fixparameter"]]) < 0.001 and nparams < 4:
	 	graph.SetPoint(point,val[config["xparameter"]],val[config["yparameter"]],mu)
		point+=1

	  txtfile.write("%2.4f\n"%(mu))
	  parameter_vals["mu"][0]=mu
	  tree.Fill()

def fillOutputs(func,graph,txtfile,tree):

	full_grid = []

	# create a range for each parameter	
	it = params.createIterator()
	for j in range(nparams):
		p = it.Next()
		vals = numpy.arange(p.getMin(),p.getMax()+step,step)
		full_grid.append(vals)
	
	c_vals = itertools.product(*full_grid)

        fillGrid(func,graph,txtfile,tree,c_vals)

	return 0

def makePlot(name,tgraph):


	if not options.sliceval: return 
	it = params.createIterator()
	if "xparameter_name" not in config.keys(): 
	  for pcounter in range(nparams):
		p = it.Next()
		#if options.sliceval:
			 #if pcounter == config["fixparameter"]: config["fixparameter_name"] = p.GetName()
  		if pcounter == config["xparameter"] : config["xparameter_name"] = p.GetName()
  		if pcounter == config["yparameter"] : config["yparameter_name"] = p.GetName()

	xparam = params.find(config["xparameter_name"])
	yparam = params.find(config["yparameter_name"])

	if options.sliceval : zparam = params.find(config["fixparameter_name"])
	
	c = ROOT.TCanvas("c","c",600,600)

	minxval = xparam.getMin()
	maxxval = xparam.getMax()
	minyval = yparam.getMin()
	maxyval = yparam.getMax()

	thist = tgraph.GetHistogram()

	if options.sliceval: name+="_%s_%.3f"%(zparam.GetName(),config["sliceval"])
	thist.SetTitle(name)
	thist.GetXaxis().SetTitle(xparam.GetName())
	thist.GetYaxis().SetTitle(yparam.GetName())
	#thist.GetZaxis().SetTitle("#sigma(%s)*BR(%s)/#sigma_{SM}"%(prod,decay))
	
	thist.GetXaxis().SetRangeUser(minxval,maxxval)
	thist.GetYaxis().SetRangeUser(minyval,maxyval)

	# Contour for mu=1 (no scaling applied)
	thist.Draw("COLZ")
	thistContour = thist.Clone()
	thistContour.SetContour(2)
	thistContour.SetContourLevel(1,1)
	thistContour.SetLineColor(10)
	thistContour.Draw("CONT3same")
	
	# lines for the SM values
	
	lvert = ROOT.TLine(1,minyval,1,maxyval)
	lhorz = ROOT.TLine(minxval,1,maxxval,1)
	lvert.SetLineStyle(2)
	lhorz.SetLineStyle(2)
	lvert.Draw()
	lhorz.Draw()
	
	c.SaveAs("%s/%s.pdf"%(options.out,name))
	c.SaveAs("%s/%s.png"%(options.out,name))
	#c.SetLogz()
	#c.SaveAs("%s_logscale.pdf"%(name))

def produceScan(modelname,extname,proddecaystring,work,energy=""):

	# Get the appropriate mapping 
	proddecay = proddecaystring+"_%s"%energy
	name = "*%s_%s*"%(extname,proddecay)
	func = work.allFunctions().selectByName(name).first()
	print "Looking for ", name 
	if not func: 
		# mostly modesl will NOT be energy dependant
	        #proddecay+="_%s"%energy	
		proddecay = proddecaystring
		name = "*%s_%s*"%(extname,proddecay)
	        print "Nope!, Looking for ", name 
	 	work.allFunctions().selectByName(name).Print()
		func = work.allFunctions().selectByName(name).first()
	if not func: #Give up!
		return 

	# Produce a plot of the scaling parameters
	tgraph = ROOT.TGraph2D()

	# And a .txt file of the numbers:
	txtfile = open("%s_%s.txt"%(modelname,proddecay),"w")
	txtfile.write("%s - %s scaling factors\n"%(modelname, proddecay))
	
	it = params.createIterator()
	for j in range(nparams):
		p = it.Next()
		if p==None : break		
		txtfile.write("%s    "%(p.GetName()))
	txtfile.write("mu\n")
	
	# And a TTree 
	tr = ROOT.TTree("%s"%(proddecay),"%s"%(proddecay))
	createBranches(tr,parameter_vals)

	# This is the loop over points in the model 		
	if not options.skipTree: fillOutputs(func,tgraph,txtfile,tr)

	# make a 2D plot 
	if nparams == 2 or (options.sliceval and nparams <4): makePlot("%s_%s"%(modelname,proddecay),tgraph)
	else: print "Skipping 2D plots (nparams != 2, and nparams!=3 with a slice value given)"

	# make 1D scans 
	it = params.createIterator()
	allplots = []
	for j in range(nparams):
		p = it.Next()
		if p==None : break		
	        plot = p.frame()
		resetVals()
		func.plotOn(plot)
		allplots.append(plot)

	cc = ROOT.TCanvas("c_PARAM_%s_%s"%(modelname,proddecay),"",1080,600)
	if not nparams == 1: cc.Divide((nparams+1)//2,2)
	for j,p in enumerate(allplots):
		cc.cd(j+1)
		p.Draw()

	cc.SaveAs("%s/%s_%s.pdf"%(options.out,modelname,proddecay))
	cc.SaveAs("%s/%s_%s.png"%(options.out,modelname,proddecay))

	# Write The Tree:
	tr.Write()

	# close the txtfile
	txtfile.close()

# MAIN FUNCTION 

wsfile = args[0]

# Select a mass for the (light) Higgs boson 
mHval = options.mh

tfile = ROOT.TFile.Open(wsfile)
work  = tfile.Get("w")

mH = work.var("MH")
mH.setVal(mHval)

# model config defines the parameters of interest
mc_s   = work.genobj("ModelConfig")
params  = mc_s.GetParametersOfInterest()
nparams = params.getSize()
print "Number of parameters in Model: ", nparams
params.Print()

parameter_vals	= {"mu":array.array('d',[0])}

# create the map of default values 
default_parameter_vals = {}
iter = params.createIterator()
while 1: 
  p = iter.Next()
  if p == None: break
  name = p.GetName()
  # should ideally sort these out, rather than hard code/set to defaults in the workspace
  #work.var(name).setMin(-2)
  #work.var(name).setMax(2)   
  default_parameter_vals[name] = p.getVal() 

# create a dictionary object which containts name, empty array
iter = params.createIterator()
pcounter = 0
index_x = 0
index_y = 1
index_z = 2

doslice = False
if options.sliceval: 
	myFixed,sliceval = options.sliceval.split(",")
	sliceval = float(sliceval)
	doslice = True
	config["sliceval"]=sliceval
	config["fixparameter_name"]=myFixed

while 1: 
  p = iter.Next()
  if p == None: break
  name = p.GetName()
  if doslice and name == myFixed:
	print "putting z index as ", pcounter 
	index_z = pcounter
	index_x = (pcounter+1)%3
	index_y = (pcounter+2)%3
  parameter_vals[name] = array.array('d',[0])
  pcounter+=1

# config will say which index is x and y (and z)

config["xparameter"]=index_x
config["yparameter"]=index_y
config["fixparameter"]=index_z


print config

#for p in params: print p.GetName()
# Output file for ROOT TTrees
tree_output = ROOT.TFile("%s/%s_trees.root"%(options.out,options.model),"RECREATE")
set_palette(ncontours=255) # For colored plots

"""
# Simple production x-section and branching ratios
for decay in decay_modes:

	ext = "BRscal"
	produceScan(options.model,ext,decay,work)
"""

# Full Matrix sigma*BR available in Model:
for prod in production_channels:
   for decay in decay_channels:
	
	ext = options.scaling_prefix
	proddecay = "%s_%s"%(prod,decay)
	for e in energies:
	  produceScan(options.model,ext,proddecay,work,energy=e)
	
#	if options.energydependant: 
#	else: 
#			produceScan(options.model,ext,proddecay,work)

	

tree_output.Close()

