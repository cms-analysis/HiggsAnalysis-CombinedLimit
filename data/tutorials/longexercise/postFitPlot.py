from __future__ import absolute_import

import HiggsAnalysis.CombinedLimit.util.plotting as plot
import ROOT
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle()

parser = argparse.ArgumentParser(description="Arguments parser")
parser.add_argument("--input_file", dest="input_file", help="Input root file from fitDiagnostics", required=True)
parser.add_argument(
    "--shape_type", dest="shape_type", help="Shape directory in input root file from fitDiagnostics (e.g. shapes_fit_b, shapes_fit_s)", required=True
)
parser.add_argument("--region", dest="cards_dir", help="Region of interest (e.g. signal_region, ch1, ch2, etc.)", required=True)
args = parser.parse_args()


canvas = ROOT.TCanvas()

fin = ROOT.TFile(args.input_file)

first_dir = args.shape_type
second_dir = args.cards_dir

h_bkg = fin.Get(first_dir + "/" + second_dir + "/total_background")
h_sig = fin.Get(first_dir + "/" + second_dir + "/total_signal")
h_dat = fin.Get(first_dir + "/" + second_dir + "/data")  # This is a TGraphAsymmErrors, not a TH1F


h_bkg.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
h_bkg.Draw("HIST")

h_err = h_bkg.Clone()
h_err.SetFillColorAlpha(12, 0.3)  # Set grey colour (12) and alpha (0.3)
h_err.SetMarkerSize(0)
h_err.Draw("E2SAME")

h_sig.SetLineColor(ROOT.kRed)
h_sig.Draw("HISTSAME")

h_dat.Draw("PSAME")

h_bkg.SetMaximum(h_bkg.GetMaximum() * 1.4)

legend = ROOT.TLegend(0.60, 0.70, 0.90, 0.91, "", "NBNDC")
legend.AddEntry(h_bkg, "Background", "F")
legend.AddEntry(h_sig, "Signal", "L")
legend.AddEntry(h_err, "Background uncertainty", "F")
legend.Draw()


canvas.SaveAs("plot_%s_%s.png" % (first_dir, second_dir))
canvas.SaveAs("plot_%s_%s.pdf" % (first_dir, second_dir))
