import ROOT
import glob
import numpy as np
import array
import json
import argparse
import math
import sys

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleYOffset(1.1)
ROOT.gStyle.SetTitleXOffset(1.2)

import CMS_tdrStyle_lumi

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
CMS_tdrStyle_lumi.extraText = "Combine Tutorial"
ROOT.gStyle.SetOptStat(0)
CMS_tdrStyle_lumi.setTDRStyle()
ROOT.gStyle.SetErrorX(0.5)

bin_labels_ZH = [
    "ZH 75-150 GeV",
    "ZH 150-250 GeV0J",
    "ZH 150-250 GeV GE1J",
    "ZH 250-400 GeV",
    "ZH >400 GeV",
]
sm_xs_input = ROOT.TFile("sm_xs.root", "r")
sm_xs = sm_xs_input.Get("hist_ZH")
unfolded = sm_xs.Clone()
pois = [
    "r_zh_75_150",
    "r_zh_150_250noj",
    "r_zh_150_250wj",
    "r_zh_250_400",
    "r_zh_gt400",
]

with open(sys.argv[1]) as json_file:
    meas_dict = json.load(json_file)["STXS"]

for i, p in enumerate(pois):
    symm_err = (float(meas_dict[p]["ErrorHi"]) + float(-1 * meas_dict[p]["ErrorLo"])) / 2.0
    unfolded.SetBinContent(i + 2, meas_dict[p]["Val"] * sm_xs.GetBinContent(i + 2))
    unfolded.SetBinError(i + 2, symm_err * sm_xs.GetBinContent(i + 2))
    sm_xs.GetXaxis().SetBinLabel(i + 2, bin_labels_ZH[i])
sm_xs.GetXaxis().SetBinLabel(1, "")
sm_xs.GetYaxis().SetTitle("#sigma #times BR [fb]")
sm_xs.SetBinContent(1, 0)
unfolded.SetBinContent(1, 0)

leg = ROOT.TLegend(0.6, 0.77, 0.8, 0.85)
leg.SetTextSize(0.03)
leg.AddEntry(0, "Material for Combine tutorial", "")
leg.AddEntry(sm_xs, "Expected")
leg.AddEntry(unfolded, "Observed")

cB = ROOT.TCanvas()
sm_xs.SetLineColor(1)
sm_xs.Draw()
unfolded.SetLineColor(9)
unfolded.Draw("same e1")
sm_xs.Draw("same e1")
leg.Draw("same")
cB.SetLogy()
cB.SaveAs("stxs_zh.png")
cB.SaveAs("stxs_zh.pdf")
