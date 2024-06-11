#!/usr/bin/env python

import ROOT
import argparse
import HiggsAnalysis.CombinedLimit.util.plotting as plot

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.2f")


translate = {
    "r_zh_250_400": ["ZH, 250<p_{T}(V)<400"],
    "r_zh_gt400": ["ZH, p_{T}(V)>400"],
    "r_zh_150_250noj": ["ZH, 150<p_{T}(V)<=250,=0J"],
    "r_zh_150_250wj": ["ZH, 150<p_{T}(V)<=250,#geq1J"],
    "r_zh_75_150": ["ZH, 75<p_{T}(V)<=150"],
}


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input mlfit file")
parser.add_argument("-p", "--parameters", help="parameters for listing correlations")
parser.add_argument("-o", "--outname", default="correlationMatrix", help="output name")

args = parser.parse_args()

fin = ROOT.TFile(args.input.split(":")[0])
rfr = fin.Get(args.input.split(":")[1])

parameters = args.parameters.split(",")

parameters_stxs = parameters
outname = args.outname

plot.SetCorrMatrixPalette()


correlation_matrix_pruned = ROOT.TH2F(
    "correlation_matrix_pruned",
    "",
    len(parameters_stxs),
    0,
    len(parameters_stxs),
    len(parameters_stxs),
    0,
    len(parameters_stxs),
)

for p1 in range(0, len(parameters_stxs)):
    for p2 in range(0, len(parameters_stxs)):
        param1 = parameters_stxs[p1]
        param2 = parameters_stxs[p2]
        correlation_matrix_pruned.SetBinContent(p1 + 1, p2 + 1, rfr.correlation(param1, param2))
    correlation_matrix_pruned.GetXaxis().SetBinLabel(p1 + 1, parameters_stxs[p1])
    correlation_matrix_pruned.GetYaxis().SetBinLabel(p1 + 1, parameters_stxs[p1])

correlation_matrix_pruned.SetStats(0)
correlation_matrix_pruned.GetZaxis().SetRangeUser(-1, 1)
correlation_matrix_pruned.GetYaxis().SetLabelSize(0.04)
correlation_matrix_pruned.GetXaxis().SetLabelSize(0.04)
correlation_matrix_pruned.GetZaxis().SetLabelSize(0.03)

outfile = ROOT.TFile("%(outname)s.root" % vars(), "RECREATE")
correlation_matrix_pruned.Write()
outfile.Close()

c2 = ROOT.TCanvas()
multiplier = 1
c2.cd()
pads = plot.OnePad()
pads[0].SetLeftMargin(0.15)
pads[0].SetRightMargin(0.10)
pads[0].SetBottomMargin(0.10)
pads[0].SetTopMargin(0.1)
ROOT.gStyle.SetTitleXOffset(1.2)
c2.SetFixedAspectRatio()
pads[0].cd()
correlation_matrix_pruned.Draw("COLZTEXT")
plot.DrawCMSLogo(pads[0], "CMS", "Material for Combine tutorial", 0.0, 0.15, 0.04, 1.0, "", 0.6)
plot.DrawTitle(pads[0], "", 3, textSize=0.4)

c2.SaveAs("%(outname)s_pois.pdf" % vars())
c2.SaveAs("%(outname)s_pois.png" % vars())
