import ROOT
from config import plot_dir

ROOT.gROOT.SetBatch(True)

# Open file with fits
f = ROOT.TFile("higgsCombine.scan.multidimfit.MultiDimFit.mH125.root")
t = f.Get("limit")

r, pdfindex = [], []

for ev in t:
    r.append(getattr(ev, "r"))
    pdfindex.append(getattr(ev, "pdfindex_Tag0"))

gr = ROOT.TGraph()
for i in range(len(r)):
    gr.SetPoint(gr.GetN(), r[i], pdfindex[i])

gr.GetXaxis().SetTitle("r")
gr.GetYaxis().SetTitle("pdfindex_Tag0")

gr.SetMarkerStyle(20)
gr.SetMarkerSize(1.5)
gr.SetLineWidth(0)

canv = ROOT.TCanvas()
gr.Draw()

canv.Update()

canv.SaveAs("%s/part5_r_vs_pdfindex.png" % plot_dir)
