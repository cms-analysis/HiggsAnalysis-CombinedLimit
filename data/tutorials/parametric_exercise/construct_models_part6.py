import ROOT

print("Using ROOT version: %s" % ROOT.__version__)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Signal modelling
# # Open ROOT file containing trees for both ggH and VBF
# f = ROOT.TFile("mc_part6.root")
#
# # Define mass and weight variables
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight","weight",0,0,1)
#
# # Load the MC as RooDataSets
# mc = {}
#
# procs = ['ggH','VBF']
# cats = ['Tag0','Tag1']
#
# systs = ['scale','smear','JEC','photonID']
#
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         t = f.Get(key)
#         mc[key] = ROOT.RooDataSet(key, key, t,
#                                   ROOT.RooArgSet(mass,weight), "", "weight"
#                                   )
#         for syst in systs:
#             for direction in ['Up','Down']:
#                 key_syst = "%s_%s%s01Sigma"%(key,syst,direction)
#                 t = f.Get(key_syst)
#                 mc[key_syst] = ROOT.RooDataSet(key_syst, key_syst,
#                     t, ROOT.RooArgSet(mass,weight), "", "weight"
#                     )
#
# # Calculate the variations in mean/sigma for the systematic-varied MC
# # Let's first build the model to fit
# mean = ROOT.RooRealVar("mean", "mean", 125, 124, 126)
# sigma = ROOT.RooRealVar("sigma", "sigma", 2, 1, 3)
# gaus = ROOT.RooGaussian("model", "model", mass, mean, sigma)
#
# mean_values, sigma_values = {}, {}
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         gaus.fitTo(mc[key], ROOT.RooFit.SumW2Error(True) )
#         gaus.fitTo(mc[key], ROOT.RooFit.SumW2Error(True) )
#         mean_values[key] = [mean.getVal(), mean.getError()]
#         sigma_values[key] = [sigma.getVal(), sigma.getError()]
#         for syst in ['scale','smear']:
#             key_up = "%s_%sUp01Sigma"%(key,syst)
#             key_down = "%s_%sDown01Sigma"%(key,syst)
#             gaus.fitTo(mc[key_up], ROOT.RooFit.SumW2Error(True) )
#             gaus.fitTo(mc[key_up], ROOT.RooFit.SumW2Error(True) )
#             mean_values[key_up] = [mean.getVal(), mean.getError()]
#             sigma_values[key_up] = [sigma.getVal(), sigma.getError()]
#             gaus.fitTo(mc[key_down], ROOT.RooFit.SumW2Error(True) )
#             gaus.fitTo(mc[key_down], ROOT.RooFit.SumW2Error(True) )
#             mean_values[key_down] = [mean.getVal(), mean.getError()]
#             sigma_values[key_down] = [sigma.getVal(), sigma.getError()]
#
# # Store variations to bake into model
# scale_variations, smear_variations = {}, {}
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         # Scale
#         syst = "scale"
#         key_up = "%s_%sUp01Sigma"%(key,syst)
#         key_down = "%s_%sDown01Sigma"%(key,syst)
#         scale_variations[key] = 0.5*(
#             abs(mean_values[key_up][0]/mean_values[key][0]-1) +
#             abs(mean_values[key_down][0]/mean_values[key][0]-1)
#             )
#         # Smear
#         syst = "smear"
#         key_up = "%s_%sUp01Sigma"%(key,syst)
#         key_down = "%s_%sDown01Sigma"%(key,syst)
#         smear_variations[key] = 0.5*(
#             abs(sigma_values[key_up][0]/sigma_values[key][0]-1) +
#             abs(sigma_values[key_down][0]/sigma_values[key][0]-1)
#             )
#
# # Now lets construct the signal models to save in the workspace
# MH = ROOT.RooRealVar("MH", "MH", 125, 120, 130 )
# MH.setConstant(True)
#
# # Dicts to store the shape parameters
# dMH, mean_formula, sigma, sigma_formula = {}, {}, {}, {}
# models = {}
#
# # Also build nuisance parameters for scale and smearing effects
# eta = ROOT.RooRealVar("nuisance_scale", "nuisance_scale", 0, -5, 5)
# eta.setConstant(True)
# chi = ROOT.RooRealVar("nuisance_smear", "nuisance_smear", 0, -5, 5)
# chi.setConstant(True)
#
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         dMH[key] = ROOT.RooRealVar("dMH_%s"%key, "dMH_%s"%key, 0, -2, 2 )
#         mean_formula[key] = ROOT.RooFormulaVar("mean_%s"%key, "mean_%s"%key,
#             "(@0+@1)*(1+%.4f*@2)"%scale_variations[key],
#             ROOT.RooArgList(MH,dMH[key],eta)
#             )
#
#         sigma[key] = ROOT.RooRealVar("sigma_%s_nominal"%key,
#             "sigma_%s_nominal"%key, 2, 1, 5
#             )
#         sigma_formula[key] = ROOT.RooFormulaVar("sigma_%s"%key,
#             "sigma_%s"%key, "@0*(1+%.4f*@1)"%smear_variations[key],
#             ROOT.RooArgList(sigma[key],chi)
#             )
#
#         models[key] = ROOT.RooGaussian("model_%s"%key, "model_%s"%key, mass,
#             mean_formula[key], sigma_formula[key]
#             )
#
#         # Fit model to MC
#         models[key].fitTo( mc[key], ROOT.RooFit.SumW2Error(True) )
#
#         # Set shape parameters of model to be constant
#         # i.e. fixed in fit to data
#         dMH[key].setConstant(True)
#         sigma[key].setConstant(True)
#
# # Now let's build the normalisation objects
#
# # From LHCHWG twiki:
# # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
# xs = {}
# xs['ggH'] = ROOT.RooRealVar("xs_ggH", "Cross section of ggH in [pb]", 48.58 )
# xs['ggH'].setConstant(True)
# xs['VBF'] = ROOT.RooRealVar("xs_VBF", "Cross section of VBF in [pb]", 3.782 )
# xs['VBF'].setConstant(True)
#
# br_gamgam = ROOT.RooRealVar("BR_gamgam",
#     "Branching ratio of Higgs to gamma gamma", 2.7e-3
#     )
# br_gamgam.setConstant(True)
#
# # Calculate selection efficiency for the different proc/cat combinations
# # And define norm objects
# eff, norms = {}, {}
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         sumw = mc[key].sumEntries()
#         e = sumw/(xs[proc].getVal()*br_gamgam.getVal())
#         eff[key] = ROOT.RooRealVar("eff_%s"%key,
#             "Efficiency for %s events to land in %s"%(proc,cat), e
#             )
#         # Set constant
#         eff[key].setConstant(True)
#         # Define normalisation objects
#         norms[key] = ROOT.RooProduct("model_%s_norm"%key,
#             "Normalisation term for %s in %s"%(proc,cat),
#             ROOT.RooArgList(xs[proc],br_gamgam,eff[key])
#             )
#
# # Build new signal model workspace with all models
# f_out = ROOT.TFile("workspace_sig_part6.root", "RECREATE")
# w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
# for model in models.values(): getattr(w_sig, "import")(model)
# for norm in norms.values(): getattr(w_sig, "import")(norm)
# w_sig.Print()
# w_sig.Write()
# f_out.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Background model construction
# # Build a RooMultiPdf to fit the background distribution in each category
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight","weight",0,0,1)
#
# # Define mass ranges to find initial background model parameter values
# mass.setRange("loSB", 100, 115 )
# mass.setRange("hiSB", 135, 180 )
# mass.setRange("full", 100, 180 )
# fit_range = "loSB,hiSB"
#
# # Define dicts to store the background model objects
# alpha, poly_1, poly_2, poly_3, poly_4, pow_1 = {}, {}, {}, {}, {}, {}
# pdfs = {}
# index = {}
# models_bkg = {}
# multipdfs = {}
# norms_bkg = {}
#
# # Open the file
# f = ROOT.TFile("data_part6.root","r")
#
# data = {}
# for cat in cats:
#
#     # Build a RooDataSet for the category
#     t = f.Get("data_%s"%cat)
#     data[cat] = ROOT.RooDataSet("data_%s"%cat, "data_%s"%cat,
#         t, ROOT.RooArgSet(mass, weight), "", "weight"
#         )
#
#     alpha[cat] = ROOT.RooRealVar("alpha_%s"%cat, "alpha_%s"%cat,
#         -0.05, -0.2, 0
#         )
#     pdfs["%s_exp"%cat] = ROOT.RooExponential("model_bkg_exp_%s"%cat,
#         "model_bkg_exp_%s"%cat, mass, alpha[cat]
#         )
#     pdfs["%s_exp"%cat].fitTo( data[cat], ROOT.RooFit.Range(fit_range),
#         ROOT.RooFit.Minimizer("Minuit2","minimize"),
#         ROOT.RooFit.SumW2Error(True)
#         )
#
#     poly_1[cat] = ROOT.RooRealVar("poly_1_%s"%cat,"poly_1_%s"%cat,
#          0.01, -4, 4
#          )
#     poly_2[cat] = ROOT.RooRealVar("poly_2_%s"%cat,"poly_2_%s"%cat,
#         0.01, -4, 4
#         )
#     poly_3[cat] = ROOT.RooRealVar("poly_3_%s"%cat,"poly_3_%s"%cat,
#         0.01, -4, 4
#         )
#     poly_4[cat] = ROOT.RooRealVar("poly_4_%s"%cat,"poly_4_%s"%cat,
#         0.01, -4, 4
#         )
#     pdfs["%s_poly"%cat] = ROOT.RooChebychev("model_bkg_poly_%s"%cat,
#         "model_bkg_poly_%s"%cat, mass,
#         ROOT.RooArgList(poly_1[cat],poly_2[cat],poly_3[cat],poly_4[cat])
#         )
#     pdfs["%s_poly"%cat].fitTo(data[cat], ROOT.RooFit.Range(fit_range),
#         ROOT.RooFit.Minimizer("Minuit2","minimize"),
#         ROOT.RooFit.SumW2Error(True)
#         )
#
#     pow_1[cat] = ROOT.RooRealVar("pow_1_%s"%cat,"pow_1_%s"%cat,
#         -3, -10, -0.0001
#     )
#     pdfs["%s_pow"%cat] = ROOT.RooGenericPdf("model_bkg_pow_%s"%cat,
#         "TMath::Power(@0,@1)", ROOT.RooArgList(mass,pow_1[cat])
#         )
#     pdfs["%s_pow"%cat].fitTo( data[cat], ROOT.RooFit.Range(fit_range),
#         ROOT.RooFit.Minimizer("Minuit2","minimize"),
#         ROOT.RooFit.SumW2Error(True)
#         )
#
#     # Make RooCategory index object to control which pdf is acitve
#     index[cat] = ROOT.RooCategory("pdfindex_%s"%cat,
#         "Index of Pdf which is active for %s"%cat
#         )
#
#     # Make ArgList of models
#     models_bkg[cat] = ROOT.RooArgList()
#     models_bkg[cat].add( pdfs["%s_exp"%cat] )
#     models_bkg[cat].add( pdfs["%s_poly"%cat] )
#     models_bkg[cat].add( pdfs["%s_pow"%cat] )
#
#     # Build the RooMultiPdf object
#     multipdfs[cat] = ROOT.RooMultiPdf("multipdf_%s"%cat,
#         "Multipdf for %s"%cat, index[cat], models_bkg[cat]
#         )
#
#     # As usual for the data driven fit
#     # We want the background to have a freely floating yield
#     norms_bkg[cat] = ROOT.RooRealVar("multipdf_%s_norm"%cat,
#         "Number of background events in %s"%cat,
#         data[cat].numEntries(), 0, 3*data[cat].numEntries()
#         )
#
# # Lets again save the data as RooDataHist with 320 bins
# mass.setBins(320)
# data_hist = {}
# for cat in cats:
#     data_hist[cat] = ROOT.RooDataHist("data_hist_%s"%cat,
#         "data_hist_%s"%cat, mass, data[cat]
#         )
#
# # Now save the background models and data sets to a workspace
# f_out = ROOT.TFile("workspace_bkg_part6.root", "RECREATE")
# w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
# for cat in cats:
#     getattr(w_bkg, "import")(data_hist[cat])
#     getattr(w_bkg, "import")(index[cat])
#     getattr(w_bkg, "import")(norms_bkg[cat])
#     getattr(w_bkg, "import")(multipdfs[cat])
# w_bkg.Print()
# w_bkg.Write()
# f_out.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Finally lets calculate the systematic yield variations
# # And add to the datacard
# print(" --> Yield variations to add to the datacard...")
# for proc in procs:
#     for cat in cats:
#         key = "%s_%s"%(proc,cat)
#         for syst in ['JEC','photonID']:
#             key_up = "%s_%sUp01Sigma"%(key,syst)
#             key_down = "%s_%sDown01Sigma"%(key,syst)
#             yield_variation_up = mc[key_up].sumEntries() / \
#                                  mc[key].sumEntries()
#             yield_variation_down = mc[key_down].sumEntries() / \
#                                    mc[key].sumEntries()
#             print("Yield variation for %s syst for (%s,%s): %.3f/%.3f" %
#                 (syst,proc,cat,yield_variation_down,yield_variation_up)
#                 )
