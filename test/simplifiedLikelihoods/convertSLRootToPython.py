import sys,optparse
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-s', '--signal'	,help='.root file:histogram for signal expectation')
parser.add_option('-d', '--data'	,help='.root file:histogram for data')
parser.add_option('-b', '--background'  ,help='.root file:histogram for nominal background (eg from fit_b)')
parser.add_option('-c', '--covariance'  ,help='.root file:2D histogram for covariance (should be same fit as background)')
parser.add_option('-O', '--outname'	,default="test.py"  ,help='outname for python config')
options,args = parser.parse_args()

outname = options.outname
object_cache = []

import ROOT

def getHisto1D(gr):
  nb =  gr.GetN()
  hist = ROOT.TH1F("hist_%s"%gr.GetName(),"",nb,0,nb)
  for i in range(nb): 
    hist.SetBinContent(i+1,gr.GetY()[i])
  object_cache.append(hist)
  return hist
    

def getObject(string):
  fi,loc = string.split(":")
  print("getting object", loc, " from", fi)
  tfi=ROOT.TFile.Open(fi)
  obj = tfi.Get(loc)
  if obj.InheritsFrom(ROOT.TGraph.Class()): obj = getHisto1D(obj)
  object_cache.append(obj)
  object_cache.append(tfi)
  return obj

def printObject(obj,out,name,is1D=True):
  out.write("%s=array.array('d',["%name)
  nb  = obj.GetNbinsX()
  numbers = []
  if is1D: 
    numbers = ["%g"%(obj.GetBinContent(i)) for i in range(1,nb+1)]
  else: 
    for i in range(1,nb+1):
     for j in range(1,nb+1):
        numbers.append("%g"%(obj.GetBinContent(i,j)))
  out.write(",".join(numbers))
  out.write("])"+'\n')
  

outf = open(outname,"w")
outf.write("import numpy"+"\n")
outf.write("import array"+"\n")
outf.write("name='dummy name'"+'\n')

signal = getObject(options.signal)
nb = signal.GetNbinsX()
outf.write("nbins=%d"%nb+'\n')

printObject(getObject(options.data),outf,"data",1)
printObject(getObject(options.signal),outf,"signal",1)
printObject(getObject(options.background),outf,"background",1)
printObject(getObject(options.covariance),outf,"covariance",0)

outf.close()
