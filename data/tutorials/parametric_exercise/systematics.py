import ROOT
print("Using ROOT version: %s" % ROOT.__version__)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Open the ROOT file containing the systematic-varied trees
# f = ROOT.TFile("mc_part3.root")
# f.ls()
#
# Define mass and weight variables
# mass = ROOT.RooRealVar("CMS_hgg_mass", "CMS_hgg_mass", 125, 100, 180)
# weight = ROOT.RooRealVar("weight", "weight", 0, 0, 1)
#
# mc = {}
#
# # Load the nominal dataset
# t = f.Get("ggH_Tag0")
# mc['nominal'] = ROOT.RooDataSet(
#     "ggH_Tag0", "ggH_Tag0", t, ROOT.RooArgSet(mass, weight), "", "weight"
#     )

# # Load the systematic-varied datasets
# for syst in ['JEC', 'photonID', 'scale', 'smear']:
#     for direction in ['Up', 'Down']:
#         key = "%s%s01Sigma" % (syst, direction)
#         name = "ggH_Tag0_%s" % (key)
#         t = f.Get(name)
#         mc[key] = ROOT.RooDataSet(
#             name, name, t, ROOT.RooArgSet(mass, weight), "", "weight" )
#             )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for syst in ['JEC','photonID']:
#     for direction in ['Up','Down']:
#         key = "%s%s01Sigma" % (syst, direction)
#         yield_variation = mc[key].sumEntries()/mc['nominal'].sumEntries()
#         print("Systematic varied yield (%s,%s): %.3f" %
#             (syst,direction,yield_variation)
#         )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Parametric shape uncertainties
# # Build the model to fit the systematic-varied datasets
# mean = ROOT.RooRealVar("mean", "mean", 125, 124, 126)
# sigma = ROOT.RooRealVar("sigma", "sigma", 2, 1.5, 2.5)
# gaus = ROOT.RooGaussian("model", "model", mass, mean, sigma)
#
# # Run the fits twice to obtain more reliable results
# gaus.fitTo(mc['scaleUp01Sigma'],
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
# gaus.fitTo(mc['scaleUp01Sigma'],
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
# print("Mean = %.3f +- %.3f GeV, Sigma = %.3f +- %.3f GeV" %
#     (mean.getVal(),mean.getError(),sigma.getVal(),sigma.getError())
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # First fit the nominal dataset
# gaus.fitTo(mc['nominal'],
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
# gaus.fitTo(mc['nominal'],
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
# # Save the mean and sigma values and errors to python dicts
# mean_values, sigma_values = {}, {}
# mean_values['nominal'] = [mean.getVal(),mean.getError()]
# sigma_values['nominal'] = [sigma.getVal(),sigma.getError()]
#
# # Next for the systematic varied datasets
# for syst in ['scale','smear']:
#     for direction in ['Up','Down']:
#         key = "%s%s01Sigma"%(syst,direction)
#         gaus.fitTo(mc[key],
#             ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#         )
#         gaus.fitTo(mc[key],
#             ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#         )
#         mean_values[key] = [mean.getVal(), mean.getError()]
#         sigma_values[key] = [sigma.getVal(), sigma.getError()]
#
# # Print the variations in mean and sigma
# for key in mean_values.keys():
#     print("%s: mean = %.3f +- %.3f GeV, sigma = %.3f +- %.3f GeV" %
#         (key,mean_values[key][0],mean_values[key][1],
#         sigma_values[key][0],sigma_values[key][1])
#         )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # Building the workspace with systematic variations
# MH = ROOT.RooRealVar("MH", "MH", 125, 120, 130 )
# MH.setConstant(True)
#
# # Define formula for mean of Gaussian
# dMH = ROOT.RooRealVar("dMH_ggH_Tag0", "dMH_ggH_Tag0", 0, -5, 5 )
# eta = ROOT.RooRealVar("nuisance_scale", "nuisance_scale", 0, -5, 5)
# eta.setConstant(True)
# mean_formula = ROOT.RooFormulaVar("mean_ggH_Tag0", "mean_ggH_Tag0",
#     "(@0+@1)*(1+0.003*@2)", ROOT.RooArgList(MH,dMH,eta)
#     )
#
# # Define formula width of Gaussian
# sigma = ROOT.RooRealVar("sigma_ggH_Tag0_nominal", "sigma_ggH_Tag0_nominal",
#     2, 1, 5
#     )
# chi = ROOT.RooRealVar("nuisance_smear", "nuisance_smear", 0, -5, 5)
# chi.setConstant(True)
# sigma_formula = ROOT.RooFormulaVar("sigma_ggH_Tag0", "sigma_ggH_Tag0",
#     "@0*(1+0.045*@1)", ROOT.RooArgList(sigma,chi)
#     )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # Define Gaussian
# model = ROOT.RooGaussian( "model_ggH_Tag0", "model_ggH_Tag0",
#     mass, mean_formula, sigma_formula
#     )
#
# # Fit model to MC
# model.fitTo( mc['nominal'],
#     ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)
#     )
#
# # Build signal model normalisation object
# xs_ggH = ROOT.RooRealVar("xs_ggH", "Cross section of ggH in [pb]", 48.58 )
# br_gamgam = ROOT.RooRealVar("BR_gamgam",
#     "Branching ratio of Higgs to gamma gamma", 0.0027
#     )
# eff = mc['nominal'].sumEntries()/(xs_ggH.getVal()*br_gamgam.getVal())
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
#
# # Set shape parameters of model to be constant (i.e. fixed in fit to data)
# dMH.setConstant(True)
# sigma.setConstant(True)
#
# # Build new signal model workspace with signal normalisation term.
# f_out = ROOT.TFile("workspace_sig_with_syst.root", "RECREATE")
# w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
# getattr(w_sig, "import")(model)
# getattr(w_sig, "import")(norm_sig)
# w_sig.Print()
# w_sig.Write()
# f_out.Close()
