from __future__ import absolute_import, print_function

import sys

import ROOT

fin = ROOT.TFile(sys.argv[1])
for name in ROOT.gDirectory.GetListOfKeys():
    if name.GetClassName() == "TCanvas":
        print(name)
        fin.Get(name.GetName()).Print(".pdf")
        fin.Get(name.GetName()).Print(".png")
