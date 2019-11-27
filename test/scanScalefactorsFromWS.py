# Massive simplification of the scaling function plotter which is more generic/much better
import os,sys,numpy,array
import itertools

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] file \nrun with --help to get list of options")
parser.add_option("-m","--mh",dest="mh",default=125.09,type='float',help="Lightest Higgs Mass")
parser.add_option("","--scaling_prefix",dest="scaling_prefix",type='str',default="XSBRscal", help="Prefix for scaler, eg usually this is XSBRscal")
parser.add_option("","--snapshotName",dest="snapshotName",type='str',default="", help="load a snapshot in the workspace at which the other parameters will be set to in each scan")
parser.add_option("","--setDefaultParameterValues",dest="setDefaultParameterValues",type='str',default="", help="set parameter values when not scanning (default is to use the values in the workspace/loaded from the snapshot). use p1=x,p2=y,p3=...")
parser.add_option("","--out",dest="out",type='str',default="", help="output folder")
parser.add_option("","--loadScan",dest="loadScan",type='str',default="", help="This is a bit special, but if you already have a likelihood scan (with a bunch of other parameters being profiled), then put --loadScan higgsCombineBlah.root:X. Then when it comes to scan X, the scan will use this grid of points with other parameters set to the saved values!")

(options,args)=parser.parse_args()

import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gROOT.SetBatch(1)

inFile = ROOT.TFile.Open(args[0])
work = inFile.Get("w")

mH = work.var("MH")
mH.setVal(options.mh)

mc_s   = work.genobj("ModelConfig")
params  = mc_s.GetParametersOfInterest()

if options.snapshotName: work.loadSnapshot(snapshotName)

if len(options.setDefaultParameterValues) : paramCList=options.setDefaultParameterValues.split(",")
else: paramCList=[]
print paramCList
# create the map of default values 
default_parameter_vals = {}
iter = params.createIterator()
while 1: 
  p = iter.Next()
  if p == None: break
  name = p.GetName()
  #p.setVal(1)
  if len(paramCList): 
   for parameterc in paramCList: 
    pn,v = parameterc.split("=")
    if pn==name: 
          p.setVal(float(v))
          print "setting default of parameter %s to %g"%(name,float(v))
  default_parameter_vals[name] = p.getVal() 

nparams = params.getSize()
print "Number of parameters in Model: ", nparams
params.Print("V")

def resetVals():
  # set the parameters 
  it = params.createIterator()
  for j in range(nparams):
	p = it.Next()
	if p==None : break		
	work.var(p.GetName()).setVal(default_parameter_vals[p.GetName()])

def runScan(function,tfilename,param): 
  print "Plotting function values %s for profiled parameters saved in %s for parameter %s as POI"%(function,tfilename,param)
  tfi = ROOT.TFile.Open(tfilename)
  tree = tfi.Get("limit")
  g = ROOT.TGraph()
  g.SetLineColor(4)
  g.SetMarkerColor(4)
  g.SetMarkerStyle(21)
  g.SetMarkerSize(1)
  g.SetLineWidth(2)
  for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    xv = getattr(tree,param)
    iter = params.createIterator()
    while 1: 
      p = iter.Next()
      if p == None: break
      name = p.GetName()
      p.setVal(getattr(tree,name))
      #print "Set parameter %s to %g"%(p.GetName(),p.getVal())
    g.SetPoint(i,xv,function.getVal())
  
  g.GetXaxis().SetTitle(param)
  g.GetYaxis().SetTitle(function.GetName())
  g.SetTitle("Profiled from %s"%tfilename)
  return g

  
allScalingFunctions = work.allFunctions().selectByName("*%s*"%(options.scaling_prefix))

iter = params.createIterator()
while 1: 
  p = iter.Next()
  if p == None: break
  name = p.GetName()
  # now go through the scaling functions and plot on parameter 
  resetVals()
  iter_f = allScalingFunctions.createIterator()
  for i in range(allScalingFunctions.getSize()):
    f = iter_f.Next()
    print  work.function(f.GetName()).getParameters(ROOT.RooArgSet()).find(p.GetName()) 
    if ( work.function(f.GetName()).getParameters(ROOT.RooArgSet()).find(p.GetName()) ): 
      resetVals()
      pl=p.frame()
      c = ROOT.TCanvas()
      param=""
      if options.loadScan: 
        fil,param=options.loadScan.split(":")
      if param==name: 
        g = runScan(f,fil,param)
	g.Draw("AP")
      else: 
        f.plotOn(pl)
        pl.GetYaxis().SetTitle(f.GetName())
	pl.Draw()
	resetVals()
	ymax=pl.GetYaxis().GetXmax()
	ymin=pl.GetYaxis().GetXmin()
	xmax=pl.GetXaxis().GetXmax()
	xmin=pl.GetXaxis().GetXmin()
	l1 = ROOT.TLine(p.getVal(),ymin,p.getVal(),ymax)
	l2 = ROOT.TLine(xmin,1,xmax,1)
	l1.SetLineColor(2)
	l2.SetLineColor(2)
	l1.Draw()
	l2.Draw()
      c.SetGridy()
      c.SetGridx()
      c.SaveAs("%s/%s_vs_%s.pdf"%(options.out,f.GetName(),p.GetName()))
      c.SaveAs("%s/%s_vs_%s.png"%(options.out,f.GetName(),p.GetName()))
print "Done!"
