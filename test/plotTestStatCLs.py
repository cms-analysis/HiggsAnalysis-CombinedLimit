#!/usr/bin/env python3
# plotTestStatCLs.py -- Nicholas Wardle

# Plot test-statistics distributions from HybridNew Grids
# Currently only supported for 1D (eg limits) although running
# --val all will work, labels will not be correct.

from __future__ import absolute_import, print_function

import sys
from optparse import OptionParser

import ROOT

ROOT.gROOT.SetBatch(1)
ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx")
ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/qmuPlot.cxx")
ROOT.gStyle.SetOptStat(0)

from ROOT import hypoTestResultTree, q0Plot, qmuPlot

parser = OptionParser()
parser.add_option("-i", "--input", type="str", help="Input Grid file")
parser.add_option(
    "-I",
    "--invert",
    default=False,
    action="store_true",
    help="Invert Alt for Null in Plotting (+ numbers). Necessary for some test-statistics (but not typical for limits/significance)",
)
parser.add_option(
    "-o",
    "--output",
    default="test_stat_distributions.root",
    type="str",
    help="outputfile for plots",
)
parser.add_option("-p", "--poi", default="r", type="str", help="poi name")
parser.add_option(
    "-v",
    "--val",
    default="all",
    type="str",
    help="POI values, comma separated (type all to make a plot for every value found in the file)",
)
parser.add_option(
    "",
    "--sub-label",
    default="",
    type="str",
    dest="sub_label",
    help="Change sub-label from q to q_sub_label (doesn't apply if using --signif option",
)
parser.add_option(
    "-r",
    "--rebin",
    default=0,
    type="int",
    help="Rebin the histos by this rebin factor (nbins)",
)
parser.add_option(
    "-m",
    "--mass",
    default=[],
    action="append",
    help="Mass value(s) (same as -m for combine)",
)
parser.add_option(
    "-P",
    "--Print",
    default=False,
    action="store_true",
    help="Just print the toys directory to see whats there",
)
parser.add_option(
    "-E",
    "--expected",
    default=False,
    action="store_true",
    help="Replace observation with expected (under alt hyp - useful for limits)",
)
parser.add_option(
    "-q",
    "--quantileExpected",
    default=-1,
    type="float",
    help="Replace observation with expected quantile (under alt hyp, i.e 1-Pb=quantileExpected)",
)
parser.add_option(
    "",
    "--doublesided",
    action="store_true",
    default=False,
    help="If 2 sided (i.e LEP style or non-nested hypos e.g for spin)",
)
parser.add_option(
    "",
    "--signif",
    action="store_true",
    default=False,
    help="If significance, plot the distribution of q0 by default for the background only hypothesis",
)

parser.add_option(
    "",
    "--plot_both",
    action="store_true",
    default=False,
    help="If --signif, Force plotting of both distributions (including for signal injected).",
)
parser.add_option(
    "",
    "--save-as-pdf",
    dest="save_as_pdf",
    action="store_true",
    default=False,
    help="Save plots as pdfs as well as Canvases in root files.",
)
(options, args) = parser.parse_args()

if options.quantileExpected >= 0:
    options.expected = True
if options.expected and options.quantileExpected < 0:
    sys.exit("You should specify a quantile for the expected result, eg 0.5 for median expected")
if options.plot_both and not options.signif:
    sys.exit("plot_both only works with significance test-stat distribution")


def findPOIvals(fi, m):
    retvals = []
    ti = fi.Get("toys")
    for k in ti.GetListOfKeys():
        obj = k
        nam = obj.GetName()
        if "HypoTestResult" not in nam:
            continue
        if m not in nam:
            continue
        if options.poi not in nam:
            continue
        lhs = nam[nam.find(options.poi) + len(options.poi) : -1]
        val = float(lhs[0 : lhs.find("_")])
        retvals.append(val)
    return retvals


ifile = ROOT.TFile.Open(options.input)
ifile.cd()

canvases = []
if not len(options.mass):
    sys.exit("Need at least one mass (or put any number if you know it will work)")
if options.Print:
    tdir = ifile.Get("toys")
    print("Contents of toys directory in file ", ifile.GetName())
    tdir.ls()
    sys.exit()

for m in options.mass:
    if options.val == "all":
        poivals = findPOIvals(ifile, m)
        poivals = list(set(poivals))
    else:
        poivals = [float(v) for v in options.val.split(",")]
    for i, pv in enumerate(poivals):
        print("Point %d/%d" % (i + 1, len(poivals)))
        ifile = ROOT.TFile.Open(options.input)
        print("Converting HypoTestResults from %s to Tree (%s=%g,mH=%s) ..." % (ifile.GetName(), options.poi, pv, m))
        # hypoTestResultTree("tmp_out%s%d.root"%(m,i),float(m),pv,options.poi)
        hypoTestResultTree("tmp_out.root", float(m), pv, options.poi)
        ifile.Close()
        # ft = ROOT.TFile.Open("tmp_out%s%d.root"%(m,i))
        ft = ROOT.TFile.Open("tmp_out.root")
        print("Plotting ... ")
        if options.signif:
            can = q0Plot(float(m), options.poi, float(pv), options.rebin, options.plot_both)
        else:
            can = qmuPlot(
                float(m),
                options.poi,
                float(pv),
                int(options.doublesided),
                int(options.invert),
                options.rebin,
                options.expected,
                options.quantileExpected,
                options.sub_label,
            )
        canvases.append(can.Clone())
        if options.save_as_pdf:
            can.SaveAs(options.output + "_" + can.GetName() + ".pdf")
        ft.Close()
ofile = ROOT.TFile(options.output, "RECREATE")
for can in canvases:
    ofile.WriteTObject(can)
print("All Plots saved to ", ofile.GetName())
ofile.Close()
