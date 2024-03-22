import CombineHarvester.CombineTools.ch as ch
import HiggsAnalysis.CombinedLimit.util.plotting as plot
import ROOT
import sys

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.2f")


plot.SetCorrMatrixPalette()
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

cb = ch.CombineHarvester()
cb.SetVerbosity(0)

cb.ParseDatacard(sys.argv[1], analysis="vhbb")
signals = [
    "ZH_lep_PTV_75_150_hbb",
    "ZH_lep_PTV_150_250_0J_hbb",
    "ZH_lep_PTV_150_250_GE1J_hbb",
    "ZH_lep_PTV_250_400_hbb",
    "ZH_lep_PTV_GT400_hbb",
]
bins = [
    "vhbb_Zee_75_150_13TeV",
    "vhbb_Zee_150_250_noj_13TeV",
    "vhbb_Zee_150_250_wj_13TeV",
    "vhbb_Zee_250_400_13TeV",
    "vhbb_Zee_gt400_13TeV",
]
hist = ROOT.TH2F("migration", "", len(signals), 0, len(signals), len(bins), 0, len(bins))
for i, b in enumerate(bins):
    total_signal = 0
    for j, p in enumerate(signals):
        procs = cb.cp().bin([b]).process([p]).GetRate()
        procs += cb.cp().bin([b]).process([p.replace("ZH", "ggZH")]).GetRate()
        procs += cb.cp().bin([b.replace("Zee", "Zmm")]).process([p]).GetRate()
        procs += cb.cp().bin([b.replace("Zee", "Zmm")]).process([p.replace("ZH", "ggZH")]).GetRate()
        total_signal += procs
        hist.SetBinContent(j + 1, i + 1, procs)
        hist.GetXaxis().SetBinLabel(j + 1, p)
    hist.GetYaxis().SetBinLabel(i + 1, b.replace("Zee", "Zll"))
    for j, p in enumerate(signals):
        hist.SetBinContent(j + 1, i + 1, hist.GetBinContent(j + 1, i + 1) / total_signal)

c1 = ROOT.TCanvas()
c1.cd()
pads = plot.OnePad()
pads[0].SetLeftMargin(0.25)
pads[0].SetBottomMargin(0.15)
hist.Draw("colz text")
c1.SaveAs("migration_matrix_zh.pdf")
c1.SaveAs("migration_matrix_zh.png")
