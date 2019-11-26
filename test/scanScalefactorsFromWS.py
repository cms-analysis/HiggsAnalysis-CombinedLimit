# Massive simplification of the scaling function plotter which is more generic/much better
import os,sys,numpy,array
import itertools

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] file \nrun with --help to get list of options")
parser.add_option("-m","--mh",dest="mh",default=125.09,type='float',help="Lightest Higgs Mass")
parser.add_option("","--scaling_prefix",dest="scaling_prefix",type='str',default="XSBRscal", help="Prefix for scaler, eg usually this is XSBRscal")
parser.add_option("","--snapshotName",dest="snapshotName",type='str',default="", help="load a snapshot in the workspace at which the other parameters will be set to in each scan")

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


allScalingFunctions = work.allFunctions().selectByName("*%s*"%(options.scaling_prefix))

if options.snapshotName: work.loadSnapshot(snapshotName)

# create the map of default values 
default_parameter_vals = {}
iter = params.createIterator()
while 1: 
  p = iter.Next()
  if p == None: break
  name = p.GetName()
  #p.setVal(1)
  default_parameter_vals[name] = p.getVal() 

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
      c.SaveAs("%s_vs_%s.pdf"%(f.GetName(),p.GetName()))
      c.SaveAs("%s_vs_%s.png"%(f.GetName(),p.GetName()))
print "Done!"
