import ROOT
import numpy as np
from config import plot_dir

ROOT.gROOT.SetBatch(True)

N_toys = 1000
r_truth = 1

truth_function = "exp"
fit_function = "poly"

name = "truth_%s_fit_%s" % (truth_function, fit_function)

# Open file with fits
f = ROOT.TFile("higgsCombine.bias_%s.MultiDimFit.mH125.123456.root" % name)
t = f.Get("limit")

hist_pull = ROOT.TH1F("pull_%s" % name,
                      "Pull distribution: truth=%s, fit=%s" %
                      (truth_function, fit_function),
                      80, -4, 4
                      )
hist_pull.GetXaxis().SetTitle("Pull = (r_{truth}-r_{fit})/#sigma_{fit}")
hist_pull.GetYaxis().SetTitle("Entries")

sigma_values = np.array([])

for i_toy in range(N_toys):
    # Best-fit value
    t.GetEntry(i_toy*3)
    r_fit = getattr(t, "r")

    # -1 sigma value
    t.GetEntry(i_toy*3+1)
    r_lo = getattr(t, "r")

    # +1 sigma value
    t.GetEntry(i_toy*3+2)
    r_hi = getattr(t, "r")

    diff = r_truth-r_fit
    # Use uncertainty depending on where mu_truth is relative to mu_fit
    if diff > 0:
        sigma = abs(r_hi-r_fit)
    else:
        sigma = abs(r_lo-r_fit)

    if sigma != 0:
        sigma_values = np.append(sigma_values, sigma)
    else:
        sigma = sigma_values.mean()

    if sigma != 0:
        hist_pull.Fill(diff/sigma)

canv = ROOT.TCanvas()
hist_pull.Draw()

# Fit Gaussian to pull distribution
ROOT.gStyle.SetOptFit(111)
hist_pull.Fit("gaus")

canv.SaveAs("%s/part4_pull_%s.png" % (plot_dir, name))
