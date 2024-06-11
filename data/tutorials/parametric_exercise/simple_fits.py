import ROOT
from config import plot_dir

ROOT.gROOT.SetBatch(True)
print("Plotting directory: %s" % plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Open the best-fit file and load the postfit workspace
# f = ROOT.TFile("higgsCombine.bestfit.MultiDimFit.mH125.root")
# w = f.Get("w")
# w.Print("v")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Lets plot the post-fit and prefit model to the data
# n_bins = 80
# binning = ROOT.RooFit.Binning(n_bins,100,180)
#
# can = ROOT.TCanvas()
# plot = w.var("CMS_hgg_mass").frame()
# w.data("data_obs").plotOn( plot, binning )
#
# # Load the S+B model
# sb_model = w.pdf("model_s").getPdf("Tag0")
#
# # Prefit
# sb_model.plotOn(plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit") )
#
# # Postfit
# w.loadSnapshot("MultiDimFit")
# sb_model.plotOn(plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit") )
# r_bestfit = w.var("r").getVal()
#
# plot.Draw()
#
# leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
# leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
# leg.AddEntry("postfit", "Postfit S+B model (r=%.2f)"%r_bestfit, "L")
# leg.Draw("Same")
#
# can.Update()
# can.SaveAs("%s/part2_sb_model.png"%plot_dir)
