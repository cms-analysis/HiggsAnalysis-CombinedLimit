from __future__ import absolute_import

import plotting as plot
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle()

canvas = ROOT.TCanvas()

fin = ROOT.TFile("fitDiagnosticsTest.mH120.root")

first_dir = "shapes_fit_s"
second_dir = "ch4"

h_bkg = fin.Get(first_dir + "/" + second_dir + "/total_background")
h_sig = fin.Get(first_dir + "/" + second_dir + "/total_signal")
h_dat = fin.Get(first_dir + "/" + second_dir + "/data")  # This is a TGraphAsymmErrors, not a TH1F


#scale by 100 the signal
h_sig.Scale(100)

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
legend.AddEntry(h_sig, "Signal*100", "L")
legend.AddEntry(h_err, "Background uncertainty", "F")
legend.Draw()

#set log scale y axis
canvas.SetLogy()

canvas.SaveAs("plot.pdf")
canvas.SaveAs("plot.png")