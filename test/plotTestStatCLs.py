#!/usr/bin/env python 
# plotTestStatCLs.py -- Nicholas Wardle 

# Plot test-statistics distributions from HybridNew Grids
# Currently only supported for 1D (eg limits) although running 
# --val all will work, labels will not be correct.

import sys
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-i","--input",type='str',help="Input Grid file")
parser.add_option("-I","--invert",default=False,action='store_true',help="Invert Alt for Null in Plotting (+ numbers)")
parser.add_option("-o","--output",default="cls_qmu_distributions.root",type='str',help="outputfile for plots")
parser.add_option("-p","--poi",default="r",type='str',help="poi name")
parser.add_option("-v","--val",default="",type='str',help="poi values, comma separated (type all to make a plot for every value found in the file?!)")
parser.add_option("-r","--rebin",default=0,type='int',help="Rebin the histos by this rebin factor (nbins)")
parser.add_option("-m","--mass",default=[],action='append',help="mass value(s) (same as -m for combine)")
parser.add_option("-P","--Print",default=False,action='store_true',help="Just print the toys directory to see whats there")
parser.add_option("","--doublesided",action='store_true',default=False,help="If 2 sided (i.e LEP style or non-nested hypos e.g for spin)")
(options,args)=parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(1)
ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx")
ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/qmuPlot.cxx")
from ROOT import hypoTestResultTree
from ROOT import qmuPlot

def findPOIvals(fi,m):
 retvals = []
 ti = fi.Get("toys")
 for k in ti.GetListOfKeys():
  obj = k
  nam = obj.GetName()
  if not "HypoTestResult" in nam: continue 
  if not m in nam: continue 
  if not options.poi in nam: continue
  lhs = nam[nam.find(options.poi)+len(options.poi):-1]
  val = float(lhs[0:lhs.find('_')])
  retvals.append(val)
 return retvals

ifile = ROOT.TFile.Open(options.input)
ifile.cd()

canvases = []
if not len(options.mass): sys.exit("Need at least one mass (or put any number if you know it will work)")
if options.Print:  
  tdir = ifile.Get("toys")
  print "Contents of toys directory in file ", ifile.GetName()
  tdir.ls()
  sys.exit()

for m in options.mass: 
 if options.val == "all": 
  poivals = findPOIvals(ifile,m)
  poivals = list(set(poivals))
 else: poivals = [float(v) for v in options.val.split(",")]
 for i,pv in enumerate(poivals):
  print "Point %d/%d"%(i+1,len(poivals))
  ifile = ROOT.TFile.Open(options.input)
  print "Converting HypoTestResults from %s to Tree (%s=%g,mH=%s) ..."%(ifile.GetName(),options.poi,pv,m)
  #hypoTestResultTree("tmp_out%s%d.root"%(m,i),float(m),pv,options.poi)
  hypoTestResultTree("tmp_out.root",float(m),pv,options.poi)
  ifile.Close()
  #ft = ROOT.TFile.Open("tmp_out%s%d.root"%(m,i))
  ft = ROOT.TFile.Open("tmp_out.root")
  print "Plotting ... "
  can = qmuPlot(float(m),options.poi,float(pv),int(options.doublesided),int(options.invert),options.rebin)
  canvases.append(can.Clone())
  ft.Close()
ofile = ROOT.TFile(options.output,"RECREATE")
for can in canvases:
  ofile.WriteTObject(can)
print "All Plots saved to ", ofile.GetName()
ofile.Close()
