import ROOT
from config import plot_dir

ROOT.gROOT.SetBatch(True)
print("Plotting directory: %s" % plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define mass and weight variables
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight","weight",0,0,1)
#
# # Load the data
# f = ROOT.TFile("data_part1.root","r")
# t = f.Get("data_Tag0")
#
# # Convert to RooDataSet
# data = ROOT.RooDataSet("data_Tag0", "data_Tag0",
#     t, ROOT.RooArgSet(mass, weight), "", "weight"
#     )
#
# # Define ranges to fit
# n_bins = 80
# binning = ROOT.RooFit.Binning(n_bins,100,180)
# mass.setRange("loSB", 100, 115 )
# mass.setRange("hiSB", 135, 180 )
# mass.setRange("full", 100, 180 )
# fit_range = "loSB,hiSB"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define the different background model pdf choices
# # And fit to the data mass sidebands
# # RooExponential
# alpha = ROOT.RooRealVar("alpha", "alpha", -0.05, -0.2, 0 )
# model_exp_bkg = ROOT.RooExponential("model_exp_bkg_Tag0",
#     "model_exp_bkg_Tag0", mass, alpha
#     )
# # Fit model to data sidebands
# model_exp_bkg.fitTo( data, ROOT.RooFit.Range(fit_range),
#     ROOT.RooFit.Minimizer("Minuit2","minimize"),
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
#
# # RooChebychev polynomial: 4th order
# poly_1 = ROOT.RooRealVar("poly_1","T1 of chebychev polynomial", 0.01, -4, 4)
# poly_2 = ROOT.RooRealVar("poly_2","T2 of chebychev polynomial", 0.01, -4, 4)
# poly_3 = ROOT.RooRealVar("poly_3","T3 of chebychev polynomial", 0.01, -4, 4)
# poly_4 = ROOT.RooRealVar("poly_4","T4 of chebychev polynomial", 0.01, -4, 4)
# model_poly_bkg = ROOT.RooChebychev("model_poly_bkg_Tag0",
#     "model_poly_bkg_Tag0", mass,
#     ROOT.RooArgList(poly_1,poly_2,poly_3,poly_4)
#     )
# # Fit model to data sidebands
# model_poly_bkg.fitTo( data, ROOT.RooFit.Range(fit_range),
#     ROOT.RooFit.Minimizer("Minuit2","minimize"),
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
#
# # Power law function: using RooGenericPdf functionality
# pow_1 = ROOT.RooRealVar("pow_1","Exponent of power law", -3, -10, -0.0001)
# model_pow_bkg = ROOT.RooGenericPdf("model_pow_bkg_Tag0",
#     "TMath::Power(@0,@1)", ROOT.RooArgList(mass,pow_1)
#     )
# # Fit model to data sidebands
# model_pow_bkg.fitTo( data, ROOT.RooFit.Range(fit_range),
#     ROOT.RooFit.Minimizer("Minuit2","minimize"),
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Lets plot the model fit to the data
# can = ROOT.TCanvas()
# plot = mass.frame()
# # We have to be careful with the normalisation as we only fit over sidebands
# data.plotOn( plot, binning,
#     ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0)
#     )
# model_exp_bkg.plotOn(plot, ROOT.RooFit.NormRange(fit_range),
#     ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(2),
#     ROOT.RooFit.Name("Exponential")
#     )
# model_poly_bkg.plotOn(plot, ROOT.RooFit.NormRange(fit_range),
#     ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(3),
#     ROOT.RooFit.Name("Polynomial")
#     )
# model_pow_bkg.plotOn(plot, ROOT.RooFit.NormRange(fit_range),
#     ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(4),
#     ROOT.RooFit.Name("PowerLaw")
#     )
# data.plotOn( plot, ROOT.RooFit.CutRange(fit_range), binning )
# plot.Draw()
#
# leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
# leg.AddEntry("Exponential", "Exponential", "L")
# leg.AddEntry("Polynomial", "Chebychev polynomial (4th order)", "L")
# leg.AddEntry("PowerLaw", "Power law", "L")
# leg.Draw("Same")
#
# can.Update()
# can.SaveAs("%s/part4_data_sidebands.png"%plot_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define normalisation objects
# # Data-driven fit: bkg model to have freely floating yield in final fit
# norm_exp = ROOT.RooRealVar("model_exp_bkg_Tag0_norm",
#     "Number of background events in Tag0 (exponential)",
#     data.numEntries(), 0, 3*data.numEntries()
#     )
# norm_poly = ROOT.RooRealVar("model_poly_bkg_Tag0_norm",
#     "Number of background events in Tag0 (polynomial)",
#     data.numEntries(), 0, 3*data.numEntries()
#     )
# norm_pow = ROOT.RooRealVar("model_pow_bkg_Tag0_norm",
#     "Number of background events in Tag0 (power law)",
#     data.numEntries(), 0, 3*data.numEntries()
#     )
#
# # Also we want parameters of models to be floating in final fit to data
# # Therefore no need to set the shape parameters of the model to be constant

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Lets save the data as a RooDataHist with 320 bins between 100 and 180
# mass.setBins(320)
# data_hist = ROOT.RooDataHist("data_hist_Tag0", "data_hist_Tag0", mass, data )
#
# # Save the background model and data set to a RooWorkspace
# w_bkg = {}
# for pdf in ['exp','poly','pow']:
#     f_out = ROOT.TFile("workspace_bkg_%s.root"%pdf, "RECREATE")
#     w_bkg[pdf] = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
#     getattr(w_bkg[pdf], "import")(data_hist)
#     if pdf == 'exp':
#         getattr(w_bkg[pdf], "import")(norm_exp)
#         getattr(w_bkg[pdf], "import")(model_exp_bkg)
#     elif pdf == 'poly':
#         getattr(w_bkg[pdf], "import")(norm_poly)
#         getattr(w_bkg[pdf], "import")(model_poly_bkg)
#     elif pdf == 'pow':
#         getattr(w_bkg[pdf], "import")(norm_pow)
#         getattr(w_bkg[pdf], "import")(model_pow_bkg)
#     w_bkg[pdf].Print()
#     w_bkg[pdf].Write()
#     f_out.Close()
