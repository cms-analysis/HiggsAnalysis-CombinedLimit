import ROOT
from config import plot_dir

ROOT.gROOT.SetBatch(True)
print("Plotting directory: %s" % plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Signal modelling
# f = ROOT.TFile("mc_part1.root","r")
# # Load TTree
# t = f.Get("ggH_Tag0")
#
# # Define mass and weight variables
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight","weight",0,0,1)
#
# # Convert to RooDataSet
# mc = ROOT.RooDataSet("ggH_Tag0","ggH_Tag0",
#     t, ROOT.RooArgSet(mass,weight), "", "weight"
#     )
#
# # Lets plot the signal mass distribution
# can = ROOT.TCanvas()
# plot = mass.frame()
# mc.plotOn(plot)
# plot.Draw()
# can.Update()
# can.SaveAs("%s/part1_signal_mass.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Introduce a RooRealVar into the workspace for the Higgs mass
# MH = ROOT.RooRealVar("MH", "MH", 125, 120, 130 )
# MH.setConstant(True)
#
# # Signal peak width
# sigma = ROOT.RooRealVar("sigma_ggH_Tag0", "sigma_ggH_Tag0", 2, 1, 5)
#
# # Define the Gaussian with mean=MH and width=sigma
# model = ROOT.RooGaussian("model_ggH_Tag0", "model_ggH_Tag0",
#     mass, MH, sigma
#     )
#
# # Fit Gaussian to MC events and plot
# model.fitTo(mc,ROOT.RooFit.SumW2Error(True))
#
# can = ROOT.TCanvas()
# plot = mass.frame()
# mc.plotOn(plot)
# model.plotOn( plot, ROOT.RooFit.LineColor(2) )
# plot.Draw()
# can.Update()
# can.Draw()
# can.SaveAs("%s/part1_signal_model_v0.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define a model with a variable mean
# dMH = ROOT.RooRealVar("dMH_ggH_Tag0", "dMH_ggH_Tag0", 0, -1, 1 )
# mean = ROOT.RooFormulaVar("mean_ggH_Tag0", "mean_ggH_Tag0",
#     "(@0+@1)", ROOT.RooArgList(MH,dMH)
#     )
# model = ROOT.RooGaussian( "model_ggH_Tag0", "model_ggH_Tag0",
#     mass, mean, sigma
#     )
#
# # Fit the new model with a variable mean
# model.fitTo(mc,ROOT.RooFit.SumW2Error(True),ROOT.RooFit.PrintLevel(-1))
#
# # Show model for different values of MH
# can = ROOT.TCanvas()
# plot = mass.frame()
# MH.setVal(120)
# model.plotOn( plot, ROOT.RooFit.LineColor(2) )
# MH.setVal(125)
# model.plotOn( plot, ROOT.RooFit.LineColor(3) )
# MH.setVal(130)
# model.plotOn( plot, ROOT.RooFit.LineColor(4) )
# plot.Draw()
# can.Update()
# can.Draw()
# can.SaveAs("%s/part1_signal_model_v1.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Save the model to a workspace
# MH.setVal(125)
# dMH.setConstant(True)
# sigma.setConstant(True)
# f_out = ROOT.TFile("workspace_sig.root", "RECREATE")
# w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
# getattr(w_sig, "import")(model)
# w_sig.Print()
# w_sig.Write()
# f_out.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define SM cross section and branching fraction values
# xs_ggH = 48.58 #in [pb]
# br_gamgam = 2.7e-3
#
# # Calculate the efficiency and print output
# sumw = mc.sumEntries()
# eff = sumw/(xs_ggH*br_gamgam)
# print("Efficiency of ggH events landing in Tag0 is: %.2f%%"%(eff*100))
#
# # Calculate the total yield and print output
# lumi = 138000
# N = xs_ggH*br_gamgam*eff*lumi
# print("For 138fb^-1, total normalisation of signal is: \
#     N = xs * br * eff * lumi = %.2f events"%N
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Background modelling
# f = ROOT.TFile("data_part1.root","r")
# t = f.Get("data_Tag0")
#
# # Convert TTree to a RooDataSet
# data = ROOT.RooDataSet("data_Tag0", "data_Tag0",
#     t, ROOT.RooArgSet(mass, weight), "", "weight"
#     )
#
# # Define mass sideband ranges on the mass variable
# n_bins = 80
# binning = ROOT.RooFit.Binning(n_bins,100,180)
# mass.setRange("loSB", 100, 115 )
# mass.setRange("hiSB", 135, 180 )
# mass.setRange("full", 100, 180 )
# fit_range = "loSB,hiSB"
#
# # Plot the data in the mass sidebands
# can = ROOT.TCanvas()
# plot = mass.frame()
# data.plotOn( plot, ROOT.RooFit.CutRange(fit_range), binning )
# plot.Draw()
# can.Update()
# can.Draw()
# can.SaveAs("%s/part1_data_sidebands.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define the background model slope parameter
# alpha = ROOT.RooRealVar("alpha", "alpha", -0.05, -0.2, 0 )
# model_bkg = ROOT.RooExponential("model_bkg_Tag0",
#     "model_bkg_Tag0", mass, alpha
#     )
#
# # Fit model to data sidebands
# model_bkg.fitTo( data, ROOT.RooFit.Range(fit_range),
#     ROOT.RooFit.PrintLevel(-1)
#     )
#
# # Let's plot the model fit to the data
# can = ROOT.TCanvas()
# plot = mass.frame()
# # We have to be careful with the normalisation as we only fit over sidebands
# # First do an invisible plot of the full data set
# data.plotOn( plot, binning,
#     ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0)
#     )
# model_bkg.plotOn( plot, ROOT.RooFit.NormRange(fit_range),
#     ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(2)
#     )
# data.plotOn( plot, ROOT.RooFit.CutRange(fit_range), binning )
# plot.Draw()
# can.Update()
# can.Draw()
# can.SaveAs("%s/part1_bkg_model.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define background model normalisation term, which is freely floating
# norm = ROOT.RooRealVar("model_bkg_Tag0_norm",
#     "Number of background events in Tag0",
#     data.numEntries(), 0, 3*data.numEntries()
#     )
# alpha.setConstant(False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Save the background model and data set to a RooWorkspace
# f_out = ROOT.TFile("workspace_bkg.root", "RECREATE")
# w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
# getattr(w_bkg, "import")(data)
# getattr(w_bkg, "import")(norm)
# getattr(w_bkg, "import")(model_bkg)
# w_bkg.Print()
# w_bkg.Write()
# f_out.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Extension: signal model normalisation object
# # Build signal model normalisation components
# xs_ggH = ROOT.RooRealVar("xs_ggH", "Cross section of ggH in [pb]", 48.58 )
# br_gamgam = ROOT.RooRealVar("BR_gamgam",
#     "Branching ratio of Higgs to gamma gamma", 0.0027
#     )
# eff_ggH_Tag0 = ROOT.RooRealVar("eff_ggH_Tag0",
#     "Efficiency for ggH events to land in Tag0", eff
#     )
# # Set values to be constant
# xs_ggH.setConstant(True)
# br_gamgam.setConstant(True)
# eff_ggH_Tag0.setConstant(True)
# # Define normalisation component as product of these three variables
# norm_sig = ROOT.RooProduct("model_ggH_Tag0_norm",
#     "Normalisation term for ggH in Tag 0",
#     ROOT.RooArgList(xs_ggH,br_gamgam,eff_ggH_Tag0)
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Build new signal model workspace with signal normalisation term.
# f_out = ROOT.TFile("workspace_sig_with_norm.root", "RECREATE")
# w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
# getattr(w_sig, "import")(model)
# getattr(w_sig, "import")(norm_sig)
# w_sig.Print()
# w_sig.Write()
# f_out.Close()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Extension: binned likelihood model
# f = ROOT.TFile("data_part1.root","r")
# t = f.Get("data_Tag0")
#
# # Convert TTree to a RooDataSet
# data = ROOT.RooDataSet("data_Tag0", "data_Tag0",
#     t, ROOT.RooArgSet(mass,weight), "", "weight"
#     )
#
# # Set bin number for mass variables
# mass.setBins(320)
# data_hist = ROOT.RooDataHist("data_hist_Tag0", "data_hist_Tag0", mass, data)
#
# # Save the background model with the RooDataHist instead
# f_out = ROOT.TFile("workspace_bkg_binned.root", "RECREATE")
# w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
# getattr(w_bkg, "import")(data_hist)
# getattr(w_bkg, "import")(norm)
# getattr(w_bkg, "import")(model_bkg)
# w_bkg.Print()
# w_bkg.Write()
# f_out.Close()
