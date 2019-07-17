"""
Python wrapper around RooAddPdfFixer in utils.
Use this e.g. if you're debugging a RooWorkspace by using python interactively.
Be careful about using it in scripts, because it loads things globally in root.

usage:
from HiggsAnalysis.CombinedLimit.RooAddPdfFixer import FixAll
f = ROOT.TFile("yourworkspace.root")
w = f.Get("w")
FixAll(w)
"""

import ROOT

__inited = False

def __init():
  global __inited
  if __inited: return
  __inited = True
  ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
  ROOT.gROOT.ProcessLine("#include <HiggsAnalysis/CombinedLimit/interface/utils.h>")

def FixAll(workspace):
  __init()
  fixer = ROOT.utils.RooAddPdfFixer()
  fixer.FixAll(workspace)
