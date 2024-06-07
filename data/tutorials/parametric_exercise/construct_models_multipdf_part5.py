import ROOT

ROOT.gROOT.SetBatch(True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define mass and weight variables
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight","weight",0,0,1)
#
# # Load the data
# f = ROOT.TFile("data_part1.root","r")
# t = f.Get("data_Tag0")
#
# data = ROOT.RooDataSet("data_Tag0", "data_Tag0",
#     t, ROOT.RooArgSet(mass,weight), "", "weight"
#     )
#
# # Define ranges to fit for initial parameter values
# mass.setRange("loSB", 100, 115 )
# mass.setRange("hiSB", 135, 180 )
# mass.setRange("full", 100, 180 )
# fit_range = "loSB,hiSB"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define the different background model pdf choices
# # And fit to the data mass sidebands
# # RooExponential
# alpha = ROOT.RooRealVar("alpha", "alpha", -0.05, -0.2, 0 )
# model_exp_bkg = ROOT.RooExponential("model_exp_bkg_Tag0",
#     "model_exp_bkg_Tag0", mass, alpha
#     )
# # Fit model to data sidebands
# model_exp_bkg.fitTo(data, ROOT.RooFit.Range(fit_range),
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
# # Make a RooCategory object: this will control which PDF is "active"
# cat = ROOT.RooCategory("pdfindex_Tag0",
#     "Index of Pdf which is active for Tag0"
#     )
#
# # Make a RooArgList of the models
# models = ROOT.RooArgList()
# models.add(model_exp_bkg)
# models.add(model_poly_bkg)
# models.add(model_pow_bkg)
#
# # Build the RooMultiPdf object
# multipdf = ROOT.RooMultiPdf("multipdf_Tag0",
#     "MultiPdf for Tag0", cat, models
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Define normalisation object
# # Data-driven fit, want the background model to have a freely floating yield
# norm = ROOT.RooRealVar("multipdf_Tag0_norm",
#     "Number of background events in Tag0",
#     data.numEntries(), 0, 3*data.numEntries()
#     )
#
# # Lets save the data as a RooDataHist with 320 bins between 100 and 180
# mass.setBins(320)
# data_hist = ROOT.RooDataHist("data_hist_Tag0", "data_hist_Tag0", mass, data )
#
# # Save the background model and data set to a RooWorkspace
# f_out = ROOT.TFile("workspace_bkg_multipdf.root", "RECREATE")
# w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
# getattr(w_bkg, "import")(data_hist)
# getattr(w_bkg, "import")(cat)
# getattr(w_bkg, "import")(norm)
# getattr(w_bkg, "import")(multipdf)
# w_bkg.Print()
# w_bkg.Write()
# f_out.Close()
