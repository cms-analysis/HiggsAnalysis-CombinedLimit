#!/usr/bin/env python3
"""Print the content and error of every bin of the shapes saved by
combine -M FitDiagnostics --saveShapes, for comparisons with reference files.
Usage: python3 printFitDiagnosticsShapes.py fitDiagnostics.root"""

import sys

import ROOT

ROOT.gROOT.SetBatch(True)

f = ROOT.TFile.Open(sys.argv[1])
if not f or f.IsZombie():
    sys.exit("Could not open %s" % sys.argv[1])

for fit in ["shapes_prefit", "shapes_fit_b", "shapes_fit_s"]:
    fitdir = f.GetDirectory(fit)
    if not fitdir:
        print("Missing directory %s" % fit)
        continue
    for chkey in fitdir.GetListOfKeys():
        chdir = fitdir.GetDirectory(chkey.GetName())
        if not chdir:
            continue
        for key in chdir.GetListOfKeys():
            h = key.ReadObj()
            # Skip the data graph and the 2D covariance histogram
            if not h.InheritsFrom("TH1") or h.InheritsFrom("TH2"):
                continue
            print("%s/%s/%s" % (fit, chdir.GetName(), h.GetName()))
            for i in range(1, h.GetNbinsX() + 1):
                print("  bin %3d: %.6e +/- %.6e" % (i, h.GetBinContent(i), h.GetBinError(i)))

f.Close()
