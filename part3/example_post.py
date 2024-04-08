from __future__ import absolute_import, print_function

from six.moves import range

import ROOT

rmin = 0
rmax = 30
nbins = 100
CL = 0.95
chains = "higgsCombineTest.MarkovChainMC.blahblahblah.root"


def findSmallestInterval(hist, CL):
    bins = hist.GetNbinsX()
    best_i = 1
    best_j = 1
    bd = bins + 1
    val = 0
    for i in range(1, bins + 1):
        integral = hist.GetBinContent(i)
        for j in range(i + 1, bins + 2):
            integral += hist.GetBinContent(j)
            if integral > CL:
                val = integral
                break
        if integral > CL and j - i < bd:
            bd = j - i
            best_j = j + 1
            best_i = i
            val = integral
    return hist.GetBinLowEdge(best_i), hist.GetBinLowEdge(best_j), val


fi_MCMC = ROOT.TFile.Open(chains)
# Sum up all of the chains (or we could take the average limit)
mychain = 0
for k in fi_MCMC.Get("toys").GetListOfKeys():
    obj = k.ReadObj
    if mychain == 0:
        mychain = k.ReadObj().GetAsDataSet()
    else:
        mychain.append(k.ReadObj().GetAsDataSet())
hist = ROOT.TH1F("h_post", ";r;posterior probability", nbins, rmin, rmax)
for i in range(mychain.numEntries()):
    mychain.get(i)
    hist.Fill(mychain.get(i).getRealValue("r"), mychain.weight())
hist.Scale(1.0 / hist.Integral())
hist.SetLineColor(1)
vl, vu, trueCL = findSmallestInterval(hist, CL)
histCL = hist.Clone()
for b in range(nbins):
    if histCL.GetBinLowEdge(b + 1) < vl or histCL.GetBinLowEdge(b + 2) > vu:
        histCL.SetBinContent(b + 1, 0)
c6a = ROOT.TCanvas()
histCL.SetFillColor(ROOT.kAzure - 3)
histCL.SetFillStyle(1001)
hist.Draw()
histCL.Draw("histFsame")
hist.Draw("histsame")
ll = ROOT.TLine(vl, 0, vl, 2 * hist.GetBinContent(hist.FindBin(vl)))
ll.SetLineColor(2)
ll.SetLineWidth(2)
lu = ROOT.TLine(vu, 0, vu, 2 * hist.GetBinContent(hist.FindBin(vu)))
lu.SetLineColor(2)
lu.SetLineWidth(2)
ll.Draw()
lu.Draw()

print(" %g %% (%g %%) interval (target)  = %g < r < %g " % (trueCL, CL, vl, vu))
