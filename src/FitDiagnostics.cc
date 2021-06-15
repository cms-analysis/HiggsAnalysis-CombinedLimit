#include "HiggsAnalysis/CombinedLimit/interface/FitDiagnostics.h"
#include "RooMinimizer.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFit.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooRealSumPdf.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooTrace.h"
#include <RooMinimizer.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TFile.h"
#include <RooStats/ModelConfig.h>
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"
#include "HiggsAnalysis/CombinedLimit/interface/Significance.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logger.h"
#include "HiggsAnalysis/CombinedLimit/interface/RobustHesse.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"

#include <Math/MinimizerOptions.h>

#include <iomanip>
using namespace RooStats;

std::string FitDiagnostics::name_ = "";
std::string FitDiagnostics::massName_ = "";
std::string FitDiagnostics::toyName_ = "";
std::string FitDiagnostics::minos_ = "poi";
std::string FitDiagnostics::out_ = ".";
bool        FitDiagnostics::makePlots_ = false;
float       FitDiagnostics::rebinFactor_ = 1.0;
int         FitDiagnostics::numToysForShapes_ = 200;
std::string FitDiagnostics::signalPdfNames_     = "shapeSig*";
std::string FitDiagnostics::filterString_     = "";
std::string FitDiagnostics::backgroundPdfNames_ = "shapeBkg*";
bool        FitDiagnostics::saveNormalizations_ = false;
bool        FitDiagnostics::savePredictionsPerToy_ = false;
bool        FitDiagnostics::oldNormNames_ = false;
bool        FitDiagnostics::saveShapes_ = false;
bool        FitDiagnostics::saveOverallShapes_ = false;
bool        FitDiagnostics::saveWithUncertainties_ = false;
bool        FitDiagnostics::justFit_ = false;
bool        FitDiagnostics::skipBOnlyFit_ = false;
bool        FitDiagnostics::skipSBFit_ = false;
bool        FitDiagnostics::noErrors_ = false;
bool        FitDiagnostics::reuseParams_ = false;
bool        FitDiagnostics::customStartingPoint_ = false;
bool        FitDiagnostics::robustHesse_ = false;
bool        FitDiagnostics::saveWithUncertsRequested_=false;
bool        FitDiagnostics::ignoreCovWarning_=false;


FitDiagnostics::FitDiagnostics() :
    FitterAlgoBase("FitDiagnostics specific options"),
    globalObservables_(0),
    nuisanceParameters_(0),
    processNormalizations_(0),
    processNormalizationsShapes_(0),
    t_fit_b_(nullptr),
    t_fit_sb_(nullptr),
    t_prefit_(nullptr)
{
    options_.add_options()
        ("minos",              	boost::program_options::value<std::string>(&minos_)->default_value(minos_), "Compute MINOS errors for: 'none', 'poi', 'all'")
        ("noErrors",  	       	"Don't compute uncertainties on the best fit value. Best if using toys (-t N) to evaluate distributions of results")
        ("out",                	boost::program_options::value<std::string>(&out_)->default_value(out_), "Directory to put the diagnostics output file in")
        ("plots",              	"Make pre/post-fit RooPlots of 1D distributions of observables and fitted models")
        ("rebinFactor",        	boost::program_options::value<float>(&rebinFactor_)->default_value(rebinFactor_), "Rebin by this factor before plotting (does not affect fitting!)")
        ("signalPdfNames",     	boost::program_options::value<std::string>(&signalPdfNames_)->default_value(signalPdfNames_), "Names of signal pdfs in plots (separated by ,)")
        ("backgroundPdfNames", 	boost::program_options::value<std::string>(&backgroundPdfNames_)->default_value(backgroundPdfNames_), "Names of background pdfs in plots (separated by ',')")
        ("saveNormalizations",  "Save post-fit normalizations RooArgSet (single toy only)")
        ("savePredictionsPerToy",  "Save post-fit normalizations and shapes per toy")
        ("oldNormNames",  	"Name the normalizations as in the workspace, and not as channel/process")
        ("saveShapes",  	"Save pre and post-fit distributions as TH1 in fitDiagnostics.root")
        ("saveWithUncertainties",  "Save also pre/post-fit uncertainties on the shapes and normalizations (from resampling the covariance matrix)")
        ("saveOverallShapes",  "Save total shapes (and covariance if used with --saveWithUncertainties), ie will produce TH1 (TH2) merging bins across all channels")
        ("numToysForShapes", 	boost::program_options::value<int>(&numToysForShapes_)->default_value(numToysForShapes_),  "Choose number of toys for re-sampling of the covariance (for shapes with uncertainties)")
        ("filterString",	boost::program_options::value<std::string>(&filterString_)->default_value(filterString_), "Filter to search for when making covariance and shapes")
        ("justFit",  		"Just do the S+B fit, don't do the B-only one, don't save output file")
        ("robustHesse",  boost::program_options::value<bool>(&robustHesse_)->default_value(robustHesse_),  "Use a more robust calculation of the hessian/covariance matrix")
        ("skipBOnlyFit",  	"Skip the B-only fit (do only the S+B fit)")
        ("skipSBFit",  	"Skip the S+B fit (do only the B-only fit)")
        ("initFromBonly",  	"Use the values of the nuisance parameters from the background only fit as the starting point for the s+b fit. Can help fit convergence")
        ("customStartingPoint", "Don't set the first POI to 0 for the background-only fit. Instead if using this option, the parameter will be fixed to its default value, which can be set with the --setParameters option.")
        ("ignoreCovWarning",    "Override the default behaviour of saveWithUncertainties being ignored if the covariance matrix is not accurate.")
   ;

    // setup a few defaults
    currentToy_=0; nToys=0; fitStatus_=0; mu_=0; muLoErr_=0; muHiErr_=0; numbadnll_=-1; nll_nll0_=-1; nll_bonly_=-1; nll_sb_=-1;
    overallBins_=0; overallNorms_=0,overallNuis_=0,overallCons_=0;
}

FitDiagnostics::~FitDiagnostics(){
   // delete the Arrays used to fill the trees;
   delete globalObservables_;
   delete nuisanceParameters_;
   delete processNormalizations_;
   delete processNormalizationsShapes_;
}

void FitDiagnostics::setToyNumber(const int iToy){
	currentToy_ = iToy;
}
void FitDiagnostics::setNToys(const int iToy){
	nToys = iToy;
}
void FitDiagnostics::applyOptions(const boost::program_options::variables_map &vm) 
{
    applyOptionsBase(vm);
    makePlots_ = vm.count("plots");
    name_ = vm["name"].as<std::string>();
    massName_ = vm["massName"].as<std::string>();
    toyName_ = vm["toyName"].as<std::string>();
    saveOverallShapes_  = vm.count("saveOverallShapes");
    saveShapes_  = saveOverallShapes_ || vm.count("saveShapes");
    saveNormalizations_  = saveShapes_ || vm.count("saveNormalizations");
    savePredictionsPerToy_ = vm.count("savePredictionsPerToy");
    oldNormNames_  = vm.count("oldNormNames");
    saveWithUncertainties_  = vm.count("saveWithUncertainties");
    saveWithUncertsRequested_ = saveWithUncertainties_;
    justFit_  = vm.count("justFit");
    skipBOnlyFit_ = vm.count("skipBOnlyFit");
    skipSBFit_ = vm.count("skipSBFit");
    noErrors_ = vm.count("noErrors");
    reuseParams_ = vm.count("initFromBonly");
    customStartingPoint_ = vm.count("customStartingPoint");
    ignoreCovWarning_ = vm.count("ignoreCovWarning");
     
    if (justFit_) { out_ = "none"; makePlots_ = false; savePredictionsPerToy_ = false; saveNormalizations_ = false; reuseParams_ = false, skipBOnlyFit_ = true; skipSBFit_ = false; }
}

bool FitDiagnostics::runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {

  if (reuseParams_ && minos_!="none"){
	std::cout << "Cannot reuse b-only fit params when running minos. Parameters will be reset when running S+B fit"<<std::endl;
	reuseParams_=false;
  }

  if (!justFit_ && out_ != "none"){
	if (currentToy_ < 1){
		const bool longName = runtimedef::get(std::string("longName"));
		std::string fdname(out_+"/fitDiagnostics"+name_);
		if (longName)
			fdname += "."+massName_+toyName_+"root";
		else
			fdname += ".root";
		fitOut.reset(TFile::Open(fdname.c_str(), "RECREATE")); 
		createFitResultTrees(*mc_s,withSystematics,savePredictionsPerToy_);
	}
  }

  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
  TCanvas *c1 = 0;
  if (makePlots_) {
      utils::tdrStyle();
      c1 = new TCanvas("c1","c1");
  }


  if ( currentToy_ > 1 && (saveShapes_) ) {
      std::cerr << " ERROR, cannot use saveShapes  with > 1 toy dataset, \n you should run multiple times with -t 1 using random seeds (-s -1) or remove those options." << std::endl;
      if ( verbose > 0 ) Logger::instance().log(std::string(Form("FitDiagnostics.cc: %d -- cannot use saveShapes with > 1 toy dataset, \n you should run multiple times with -t 1 using random seeds (-s -1) or remove those options",__LINE__)),Logger::kLogLevelError,__func__);
      assert(0);
  }

  // Determine pre-fit values of nuisance parameters
  // if (currentToy_ < 1){
    const RooArgSet *nuis      = mc_s->GetNuisanceParameters();
    const RooArgSet *globalObs = mc_s->GetGlobalObservables();
    if (!justFit_ && nuis && globalObs ) {
      std::auto_ptr<RooAbsPdf> nuisancePdf(utils::makeNuisancePdf(*mc_s));
      w->loadSnapshot("toyGenSnapshot");
      r->setVal(preFitValue_);
      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_prefit");
          ToySampler sampler(&*nuisancePdf, nuis);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_prefit",data);
          delete norms;
      }
      if (withSystematics)	{
	  setFitResultTrees(mc_s->GetNuisanceParameters(),nuisanceParameters_);
	  setFitResultTrees(mc_s->GetGlobalObservables(),globalObservables_);
      }
      if (savePredictionsPerToy_){
	   RooArgSet *norms = new RooArgSet();
	   norms->setName("norm_prefit");
	   getNormalizationsSimple(mc_s->GetPdf(), *mc_s->GetObservables(), *norms);
	   setNormsFitResultTrees(norms,processNormalizations_);
	   std::map<std::string,ShapeAndNorm> snm;
	   getShapesAndNorms(mc_s->GetPdf(),*mc_s->GetObservables(), snm, "");
	   setShapesFitResultTrees(snm,processNormalizationsShapes_);
	   delete norms;
      }
      w->loadSnapshot("clean");
      r->setVal(preFitValue_);
      RooCategory dummyCat("dummyCat", "");
      RooSimultaneousOpt simNuisancePdf("simNuisancePdf", "", dummyCat);
      simNuisancePdf.addExtraConstraints(((RooProdPdf*)(nuisancePdf.get()))->pdfList());
      std::auto_ptr<RooDataSet> globalData(new RooDataSet("globalData","globalData", RooArgSet(dummyCat)));
      std::auto_ptr<RooAbsReal> nuisanceNLL(simNuisancePdf.RooAbsPdf::createNLL(*globalData, RooFit::Constrain(*nuis)));
      RooFitResult *res_prefit = 0;
       // Fit to nuisance pdf to get fitRes for sampling
      {
            CloseCoutSentry sentry(verbose < 2);
            CascadeMinimizer minim(*nuisanceNLL, CascadeMinimizer::Constrained);
            minim.minimize();
            minim.hesse();
            if (minos_ == "all") minim.minos(*nuis);
            res_prefit = minim.save();
      }
      if (fitOut.get() && currentToy_ < 1) fitOut->WriteTObject(res_prefit, "nuisances_prefit_res");
      if (fitOut.get() && currentToy_ < 1) fitOut->WriteTObject(nuis->snapshot(), "nuisances_prefit");

      nuisancePdf.reset();
      globalData.reset();

      if (makePlots_ && currentToy_ < 1) {
	std::vector<RooPlot *> plots = utils::makePlots(*mc_s->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_,res_prefit);
	for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
	    (*it)->Draw(); 
	    c1->Print((out_+"/"+(*it)->GetName()+"_prefit.png").c_str());
	    c1->SetLogy();c1->Print((out_+"/"+(*it)->GetName()+"_prefit_logy.png").c_str()); c1->SetLogy(false);
	    if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_prefit").c_str());
	}
      }
      delete res_prefit;
    } else if (nuis) {
      if (fitOut.get() ) fitOut->WriteTObject(nuis->snapshot(), "nuisances_prefit");
    }
  // }
  if (t_prefit_) {
      t_prefit_->Fill();
      resetFitResultTrees(withSystematics);
  }
 
  RooFitResult *res_b = 0, *res_s = 0;
  const RooCmdArg &constCmdArg_s = withSystematics  ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooFit::NumCPU(1); // use something dummy 
  //const RooCmdArg &minosCmdArg = minos_ == "poi" ?  RooFit::Minos(*mc_s->GetParametersOfInterest())   : RooFit::Minos(minos_ != "none");  //--> dont use fitTo!
  w->loadSnapshot("clean");


  /* Background only fit (first POI set to customStartingPoint or 0) ****************************************************************/

  if (!customStartingPoint_) r->setVal(0.0);
  else std::cout << "customStartingPoint set to true, background-only fit will correspond to " << r->GetName() << " = " << r->getVal() << std::endl;
  r->setConstant(true);

  // Setup Nll before calling fits;
  if (currentToy_<1) nll.reset(mc_s->GetPdf()->createNLL(data,constCmdArg_s,RooFit::Extended(mc_s->GetPdf()->canBeExtended())));
  // Get the nll value on the prefit
  double nll0 = nll->getVal();

  if (justFit_ || skipBOnlyFit_ ) { 
    // skip b-only fit
  } else if (minos_ != "all") {
    RooArgList minos; 
    res_b = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/true,/*ndim*/1,/*reuseNLL*/ true); 
    nll_bonly_=nll->getVal()-nll0;   
  } else {
    CloseCoutSentry sentry(verbose < 2);
    RooArgList minos = (*mc_s->GetNuisanceParameters()); 
    res_b = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/true,/*ndim*/1,/*reuseNLL*/ true); 

    if (res_b) nll_bonly_ = nll->getVal() - nll0;

  }

  if (res_b && robustHesse_) {
    RobustHesse robustHesse(*nll, verbose - 1);
    robustHesse.ProtectArgSet(*mc_s->GetParametersOfInterest());
    robustHesse.hesse();
    auto res_b_new = robustHesse.GetRooFitResult(res_b);
    delete res_b;
    res_b = res_b_new;
  }

  if (res_b) { 
      if (verbose > 1) res_b->Print("V");
      if (fitOut.get()) {
        if (currentToy_< 1)	fitOut->WriteTObject(res_b,"fit_b");
        if (withSystematics)	{
          setFitResultTrees(mc_s->GetNuisanceParameters(),nuisanceParameters_);
          setFitResultTrees(mc_s->GetGlobalObservables(),globalObservables_);
        }
        fitStatus_ = res_b->status();
      }
      numbadnll_=res_b->numInvalidNLL();

      if (!robustHesse_ && res_b->covQual() < 3){
          if(!saveWithUncertainties_){
              std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
              Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
          } else if (ignoreCovWarning_) {
              std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. Caution: by passing --ignoreCovWarning the shapes and uncertainties will be stored as configured via the command line despite this issue. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
              Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. Caution: by passing --ignoreCovWarning the shapes and uncertainties will be stored as configured via the command line despite this issue. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
          } else {
              saveWithUncertainties_=false;
              std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. The option --saveWithUncertainties will be ignored as it would lead to incorrect results. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
              Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in b-only fit. The option --saveWithUncertainties will be ignored as it would lead to incorrect results. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
          }
      }
      if ( verbose > 0 ) Logger::instance().log(std::string(Form("FitDiagnostics.cc: %d -- Fit B-only, status = %d, numBadNLL = %d, covariance quality = %d",__LINE__,fitStatus_,numbadnll_,res_b->covQual())),Logger::kLogLevelDebug,__func__);

      if (makePlots_ && currentToy_<1) {
          std::vector<RooPlot *> plots = utils::makePlots(*mc_b->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_,res_b);
          for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
              c1->cd(); (*it)->Draw(); 
              c1->Print((out_+"/"+(*it)->GetName()+"_fit_b.png").c_str());
              c1->SetLogy(); c1->Print((out_+"/"+(*it)->GetName()+"_fit_b_logy.png").c_str()); c1->SetLogy(false);
              if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_fit_b").c_str());
          }
      }
      if (savePredictionsPerToy_){
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_b");
          getNormalizationsSimple(mc_s->GetPdf(), *mc_s->GetObservables(), *norms);
          setNormsFitResultTrees(norms,processNormalizations_);
          std::map<std::string,ShapeAndNorm> snm;
          getShapesAndNorms(mc_s->GetPdf(),*mc_s->GetObservables(), snm, "");
          setShapesFitResultTrees(snm,processNormalizationsShapes_);
          delete norms;
      }
      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_b");
          CovarianceReSampler sampler(res_b);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_fit_b",data);
          delete norms;
      }

      if (makePlots_ && currentToy_<1)  {
          TH2 *corr = res_b->correlationHist();
          c1->SetLeftMargin(0.25);  c1->SetBottomMargin(0.25);
          corr->SetTitle("Correlation matrix of fit parameters");
          gStyle->SetPaintTextFormat(res_b->floatParsFinal().getSize() > 10 ? ".1f" : ".2f");
          gStyle->SetOptStat(0);
          corr->SetMarkerSize(res_b->floatParsFinal().getSize() > 10 ? 2 : 1);
          corr->Draw("COLZ TEXT");
          c1->Print((out_+"/covariance_fit_b.png").c_str());
          c1->SetLeftMargin(0.16);  c1->SetBottomMargin(0.13);
          if (fitOut.get()) fitOut->WriteTObject(corr, "covariance_fit_b");
      }

      //take the limit value from "b-only" when skipping s+b
      if (skipSBFit_) {
        Combine::toggleGlobalFillTree(true);
        limit    = r->getVal();
        limitErr = r->getError();
        Combine::commitPoint(/*expected=*/false, /*quantile=*/-1.);
        Combine::toggleGlobalFillTree(false);
      }
  }
  else {
    fitStatus_=-1;
    numbadnll_=-1;
    //if just doing b-only fit, we care if it failed
    if(!justFit_ and !skipBOnlyFit_ and skipSBFit_){
      std::cout << "\n --- FitDiagnostics ---" << std::endl;
      std::cout << "Fit failed."  << std::endl;
    }
  }
  mu_=r->getVal();
  if (t_fit_b_) {
      t_fit_b_->Fill();
      resetFitResultTrees(withSystematics);}
  // no longer need res_b
  delete res_b;

  /**********************************************************************************************************************************/
  /* S+B fit (Signal parameters free to float ) *************************************************************************************/

  if (!reuseParams_) w->loadSnapshot("clean"); // Reset, also ensures nll_prefit is same in call to doFit for b and s+b
  r->setVal(preFitValue_); r->setConstant(false); 
  if (skipSBFit_) {
    // skip s+b fit
  }
  else if (minos_ != "all") {
    RooArgList minos; if (minos_ == "poi") minos.add(*r);
    res_s = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/!noErrors_,/*ndim*/1,/*reuseNLL*/ true); 
    nll_sb_ = nll->getVal()-nll0;
  } else {
    CloseCoutSentry sentry(verbose < 2);
    RooArgList minos = (*mc_s->GetNuisanceParameters()); 
    minos.add((*mc_s->GetParametersOfInterest()));  // Add POI this time 
    res_s = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/true,/*ndim*/1,/*reuseNLL*/ true); 
    if (res_s) nll_sb_= nll->getVal()-nll0;

  }

  if (res_s && robustHesse_) {
    RobustHesse robustHesse(*nll, verbose - 1);
    robustHesse.ProtectArgSet(*mc_s->GetParametersOfInterest());
    robustHesse.hesse();
    auto res_s_new = robustHesse.GetRooFitResult(res_s);
    delete res_s;
    res_s = res_s_new;
  }
  /**********************************************************************************************************************************/
  if (res_s) { 
      limit    = r->getVal();
      limitErr = r->getError();
      if (verbose > 1) res_s->Print("V");
      if (fitOut.get()){
	 if (currentToy_<1) fitOut->WriteTObject(res_s, "fit_s");

	 if (withSystematics)	{
	   setFitResultTrees(mc_s->GetNuisanceParameters(),nuisanceParameters_);
	   setFitResultTrees(mc_s->GetGlobalObservables(),globalObservables_);
	 }
	 fitStatus_ = res_s->status();
         numbadnll_ = res_s->numInvalidNLL();

         saveWithUncertainties_=saveWithUncertsRequested_; //Reset saveWithUncertainties flag to original value in case it has been set to false due to covariance matrix issues in the b-only fit.
         if (!robustHesse_ && res_s->covQual() < 3){
             if(!saveWithUncertainties_){
                 std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
                 Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
             } else if (ignoreCovWarning_) {
                  std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. Caution: by passing --ignoreCovWarning the shapes and uncertainties will be stored as configured via the command line despite this issue. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
                  Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. Caution: by passing --ignoreCovWarning the shapes and uncertainties will be stored as configured via the command line despite this issue. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
             } else {
                  saveWithUncertainties_=false;
                  std::cerr<<"[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. The option --saveWithUncertainties will be ignored as it would lead to incorrect results. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."<<std::endl;
                  Logger::instance().log(std::string("[WARNING]: Unable to determine uncertainties on all fit parameters in s+b fit. The option --saveWithUncertainties will be ignored as it would lead to incorrect results. Have a look at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fit-parameter-uncertainties for more information."),Logger::kLogLevelError,__func__);
             }
         }
         if ( verbose > 0 ) Logger::instance().log(std::string(Form("FitDiagnostics.cc: %d -- Fit S+B, status = %d, numBadNLL = %d, covariance quality = %d",__LINE__,fitStatus_,numbadnll_,res_s->covQual())),Logger::kLogLevelDebug,__func__);
	 // Additionally store the nll_sb - nll_bonly (=0.5*q0)
	 nll_nll0_ =  nll_sb_ -  nll_bonly_;
      }

      if (makePlots_ && currentToy_<1) {
          std::vector<RooPlot *> plots = utils::makePlots(*mc_s->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_,res_s);
          for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
              c1->cd(); (*it)->Draw(); 
              c1->Print((out_+"/"+(*it)->GetName()+"_fit_s.png").c_str());
              c1->SetLogy(); c1->Print((out_+"/"+(*it)->GetName()+"_fit_s_logy.png").c_str()); c1->SetLogy(false);
              if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_fit_s").c_str());
          }
      }

      if (savePredictionsPerToy_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_s");
          getNormalizationsSimple(mc_s->GetPdf(), *mc_s->GetObservables(), *norms);
          setNormsFitResultTrees(norms,processNormalizations_);
	  std::map<std::string,ShapeAndNorm> snm;
	  getShapesAndNorms(mc_s->GetPdf(),*mc_s->GetObservables(), snm, "");
          setShapesFitResultTrees(snm,processNormalizationsShapes_);
	  delete norms;
      }
      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_s");
          CovarianceReSampler sampler(res_s);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_fit_s",data);
	  delete norms;
      }

      if (makePlots_&& currentToy_< 1)  {
          TH2 *corr = res_s->correlationHist();
          c1->SetLeftMargin(0.25);  c1->SetBottomMargin(0.25);
          corr->SetTitle("Correlation matrix of fit parameters");
          gStyle->SetPaintTextFormat(res_s->floatParsFinal().getSize() > 10 ? ".1f" : ".2f");
          gStyle->SetOptStat(0);
          corr->SetMarkerSize(res_s->floatParsFinal().getSize() > 10 ? 2 : 1);
          corr->Draw("COLZ TEXT");
          c1->Print((out_+"/covariance_fit_s.png").c_str());
          c1->SetLeftMargin(0.16);  c1->SetBottomMargin(0.13);
          if (fitOut.get() ) fitOut->WriteTObject(corr, "covariance_fit_s");
      }
  }  else {
	fitStatus_=-1;
	numbadnll_=-1;
  	nll_nll0_ = -1;
  }
  mu_=r->getVal();

  if (res_s) {
      RooRealVar *rf = dynamic_cast<RooRealVar*>(res_s->floatParsFinal().find(r->GetName()));
      if (rf){
	double bestFitVal = rf->getVal();

	double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal : rf->getAsymErrorHi());
	double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal : rf->getAsymErrorLo());
	double maxError = std::max<double>(std::max<double>(hiErr, loErr), rf->getError());

	if (fabs(hiErr) < 0.001*maxError) hiErr = -bestFitVal + rf->getMax();
	if (fabs(loErr) < 0.001*maxError) loErr = +bestFitVal - rf->getMin();

	muLoErr_=loErr;
	muHiErr_=hiErr;
	muErr_  =rf->getError();

	double hiErr95 = +(do95_ && rf->hasRange("err95") ? rf->getMax("err95") - bestFitVal : 0);
	double loErr95 = -(do95_ && rf->hasRange("err95") ? rf->getMin("err95") - bestFitVal : 0);

	limit = bestFitVal;  limitErr = 0;
	if (!noErrors_) Combine::commitPoint(/*expected=*/true, /*quantile=*/0.5);
	limit = bestFitVal - loErr; limitErr = 0;
	if (!noErrors_) Combine::commitPoint(/*expected=*/true, /*quantile=*/0.16);
	limit = bestFitVal + hiErr; limitErr = 0;
	if (!noErrors_) Combine::commitPoint(/*expected=*/true, /*quantile=*/0.84);
	if (do95_ && rf->hasRange("err95") && !noErrors_) {
	  limit = rf->getMax("err95"); Combine::commitPoint(/*expected=*/true, /*quantile=*/0.975);
	  limit = rf->getMin("err95"); Combine::commitPoint(/*expected=*/true, /*quantile=*/0.025);
	}

	limit = bestFitVal;
	limitErr = maxError;
	std::cout << "\n --- FitDiagnostics ---" << std::endl;
	std::cout << "Best fit " << r->GetName() << ": " << rf->getVal() << "  "<<  -loErr << "/+" << +hiErr << "  (68% CL)" << std::endl;
	if (do95_) {
	  std::cout << "         " << r->GetName() << ": " << rf->getVal() << "  "<<  -loErr95 << "/+" << +hiErr95 << "  (95% CL)" << std::endl;
	}
      } else {
      	std::cout << "\n --- FitDiagnostics ---" << std::endl;
      	std::cout << "Done! fit s+b and fit b should be identical" << std::endl;
      }
  } else if (!skipSBFit_) {
      std::cout << "\n --- FitDiagnostics ---" << std::endl;
      std::cout << "Fit failed."  << std::endl;
  }
  if (t_fit_sb_) {
      t_fit_sb_->Fill();
      resetFitResultTrees(withSystematics);}

  if (currentToy_==nToys-1 || nToys==0 ) {
        
        if (fitOut.get()) {	
		fitOut->cd();
		t_fit_sb_->Write(); t_fit_b_->Write(); t_prefit_->Write();
		fitOut.release()->Close();
	}

  } 
  bool fitreturn = (res_s!=0);
  delete res_s;

  /*
  if(saveWorkspace_){
	  RooWorkspace *ws = new RooWorkspace("FitDiagnosticsResult");
	  ws->import(*mc_s->GetPdf());
	  ws->import(data);
	  std::cout << "Saving pdfs and data to FitDiagnosticsResult.root" << std::endl;
	  ws->writeToFile("FitDiagnosticsResult.root");
  }
  */

  return fitreturn;
}

void FitDiagnostics::getNormalizationsSimple(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out) {
    std::map<std::string,ShapeAndNorm> snm;
    getShapesAndNorms(pdf,obs, snm, "");
    std::map<std::string,ShapeAndNorm>::const_iterator bg = snm.begin(), ed = snm.end(), pair; int i;
    for (pair = bg, i = 0; pair != ed; ++pair, ++i) {
        RooRealVar *val = new RooRealVar((oldNormNames_ ? pair->first : pair->second.channel+"/"+pair->second.process).c_str(), "",pair->second.norm->getVal());
        // val->setError(sumx2[i]);
        out.addOwned(*val); 
    }
    return;
}
void FitDiagnostics::getShapesAndNorms(RooAbsPdf *pdf, const RooArgSet &obs, std::map<std::string,ShapeAndNorm> &out, const std::string &channel) {
    RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(pdf);
    if (sim != 0) {
        RooAbsCategoryLValue &cat = const_cast<RooAbsCategoryLValue &>(sim->indexCat());
        for (int i = 0, n = cat.numBins((const char *)0); i < n; ++i) {
            cat.setBin(i);
            RooAbsPdf *pdfi = sim->getPdf(cat.getLabel());
            if (pdfi) getShapesAndNorms(pdfi, obs, out, cat.getLabel());
        }        
        return;
    }
    RooProdPdf *prod = dynamic_cast<RooProdPdf *>(pdf);
    if (prod != 0) {
        RooArgList list(prod->pdfList());
        for (int i = 0, n = list.getSize(); i < n; ++i) {
            RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
            if (pdfi->dependsOn(obs)) getShapesAndNorms(pdfi, obs, out, channel);
        }
        return;
    }
    RooAddPdf *add = dynamic_cast<RooAddPdf *>(pdf);
    if (add != 0) {
        RooArgList clist(add->coefList());
        RooArgList plist(add->pdfList());
        for (int i = 0, n = clist.getSize(); i < n; ++i) {
            RooAbsReal *coeff = (RooAbsReal *) clist.at(i);
	    std::string coeffName = coeff->GetName();
	    if (coeffName.find(filterString_) == std::string::npos) continue; 
            ShapeAndNorm &ns = out[coeffName];
            ns.norm = coeff;
            ns.pdf = (RooAbsPdf*) plist.at(i);
            ns.channel = (coeff->getStringAttribute("combine.channel") ? coeff->getStringAttribute("combine.channel") : channel.c_str());
            ns.process = (coeff->getStringAttribute("combine.process") ? coeff->getStringAttribute("combine.process") : ns.norm->GetName());
            ns.signal  = (coeff->getStringAttribute("combine.process") ? coeff->getAttribute("combine.signal") : (strstr(ns.norm->GetName(),"shapeSig") != 0));
            std::auto_ptr<RooArgSet> myobs(ns.pdf->getObservables(obs));
            ns.obs.add(*myobs);
            ns.isfunc = false;
        }
        return;
    }
    RooRealSumPdf *sum = dynamic_cast<RooRealSumPdf *>(pdf);
    if (sum != 0) {
      RooArgList clist(sum->coefList());
      RooArgList plist(sum->funcList());
      if (plist.getSize() == 1) {
         CMSHistErrorPropagator *err = dynamic_cast<CMSHistErrorPropagator*>(plist.at(0));
         if (err) {
           clist.removeAll();
           plist.removeAll();
           clist.add(err->coefList());
           plist.add(err->wrapperList());
         }
      }
      for (int i = 0, n = clist.getSize(); i < n; ++i) {
        RooAbsReal *coeff = (RooAbsReal *) clist.at(i);
        std::string coeffName = coeff->GetName();
        if (coeffName.find(filterString_) == std::string::npos) continue;
        ShapeAndNorm &ns = out[coeffName];
        RooAbsReal* shape = (RooAbsReal*)plist.at(i);
        std::auto_ptr<RooArgSet> myobs(shape->getObservables(obs));
        ns.pdf = shape;
        TString normProdName = TString::Format("%s_mlf_fullnorm", coeff->GetName());
        RooAbsReal * normProd = nullptr;
        if (coeff->ownedComponents()) {
          normProd = dynamic_cast<RooAbsReal*>(coeff->ownedComponents()->find(normProdName));
        }
        if (!normProd) {
          RooAbsReal* integral = shape->createIntegral(*myobs);
          normProd = new RooProduct(normProdName, "", RooArgList(*integral, *coeff));
          normProd->addOwnedComponents(RooArgSet(*integral));
          coeff->addOwnedComponents(RooArgSet(*normProd));
        }

        ns.norm = normProd;
        ns.channel = (coeff->getStringAttribute("combine.channel") ? coeff->getStringAttribute("combine.channel") : channel.c_str());
        ns.process = (coeff->getStringAttribute("combine.process") ? coeff->getStringAttribute("combine.process") : ns.norm->GetName());
        ns.signal  = (coeff->getStringAttribute("combine.process") ? coeff->getAttribute("combine.signal") : (strstr(ns.norm->GetName(),"shapeSig") != 0));
        ns.obs.add(*myobs);
        ns.isfunc = true;
      }
    }
}

void FitDiagnostics::getNormalizations(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out, NuisanceSampler & sampler, TDirectory *fOut, const std::string &postfix,RooAbsData &data) {
    // fill in a map
    std::map<std::string,ShapeAndNorm> snm;
    getShapesAndNorms(pdf,obs, snm, "");
    typedef std::map<std::string,ShapeAndNorm>::const_iterator IT;
    typedef std::map<std::string,TH1*>::const_iterator IH;
    typedef std::map<std::string,TGraphAsymmErrors*>::const_iterator IG;
    typedef std::map<std::string,TH2*>::const_iterator IH2;
    // create directory structure for shapes
    TDirectory *shapeDir = fOut && saveShapes_ ? fOut->mkdir((std::string("shapes")+postfix).c_str()) : 0;
    std::map<std::string,TDirectory*> shapesByChannel;
    std::map<std::string,TGraphAsymmErrors*> datByCh;
    if (shapeDir) {
        for (IT it = snm.begin(), ed = snm.end(); it != ed; ++it) {
            TDirectory *& sub = shapesByChannel[it->second.channel];
            if (sub == 0) {
		sub = shapeDir->mkdir(it->second.channel.c_str());
		RooRealVar *x = (RooRealVar*)it->second.obs.at(0);
		RooDataHist *dataCut = (RooDataHist*)data.reduce(RooFit::Cut(TString("CMS_channel==CMS_channel::"+it->second.channel)));
		TH1* dH = dataCut->createHistogram("data",*x);
		TGraphAsymmErrors * dGraph = utils::makeDataGraph(dH,/*asDensity*/true);
		datByCh[it->second.channel] = (TGraphAsymmErrors *) dGraph->Clone();
		datByCh[it->second.channel]->SetNameTitle("data", (it->second.process+" in "+it->second.channel).c_str());
		delete dH;
	    }
        }
    }
    // now let's start with the central values
    std::vector<double> vals(snm.size(), 0.), sumx2(snm.size(), 0.);
    std::vector<TH1*>   shapes(snm.size(), 0), shapes2(snm.size(), 0);
    std::vector<int>    bins(snm.size(), 0), sig(snm.size(), 0);
    std::map<std::string,TH1*> totByCh, totByCh2, sigByCh, sigByCh2, bkgByCh, bkgByCh2, widthByCh;
    std::map<std::string,TH2*> totByCh2Covar;
    IT bg = snm.begin(), ed = snm.end(), pair; int i;
    for (pair = bg, i = 0; pair != ed; ++pair, ++i) {  
        vals[i] = pair->second.norm->getVal();
        //out.addOwned(*(new RooConstVar(pair->first.c_str(), "", pair->second.norm->getVal())));
        if (fOut != 0 && saveShapes_ && pair->second.obs.getSize() == 1) {
            RooRealVar *x = (RooRealVar*)pair->second.obs.at(0);
            TH1* hist = pair->second.pdf->createHistogram("", *x, pair->second.isfunc ? RooFit::Extended(false) : RooCmdArg::none());
            for (int binN = 1; binN <= hist->GetNbinsX(); ++binN){
                hist->SetBinError(binN,0);
            }
            hist->SetNameTitle(pair->second.process.c_str(), (pair->second.process+" in "+pair->second.channel).c_str());
            hist->Scale(vals[i] / hist->Integral("width"));
            hist->SetDirectory(shapesByChannel[pair->second.channel]);
            shapes[i] = hist;
            //if (saveWithUncertainties_) {
            shapes2[i] = (TH1*) hist->Clone();
            shapes2[i]->SetDirectory(0);
            shapes2[i]->Reset();
            bins[i] = hist->GetNbinsX();
            TH1 *&htot = totByCh[pair->second.channel];
            if (htot == 0) {
                    htot = (TH1*) hist->Clone();
                    htot->SetName("total");
		    htot->SetTitle(Form("Total signal+background in %s", pair->second.channel.c_str()));
                    htot->SetDirectory(shapesByChannel[pair->second.channel]);
                    TH1 *htot2 = (TH1*) hist->Clone(); htot2->Reset();
                    htot2->SetDirectory(0);
                    totByCh2[pair->second.channel] = htot2;
		    TH2F *htot2covar = new TH2F("total_covar","Covariance signal+background",bins[i],0,bins[i],bins[i],0,bins[i]);
		    htot2covar->GetXaxis()->SetTitle("Bin number");
		    htot2covar->GetYaxis()->SetTitle("Bin number");
		    htot2covar->GetZaxis()->SetTitle(Form("covar (%s)",hist->GetYaxis()->GetTitle()));
		    htot2covar->SetDirectory(0);
		    totByCh2Covar[pair->second.channel] = htot2covar;


                // The bin width is stored so that
                // one can later reverse bin width normalization
                // if desired
                TH1 * hwidth = (TH1*) htot->Clone("width");
                for(int i=1; i < hwidth->GetNbinsX()+1; i++) {
                    float width = hwidth->GetBinWidth(i);
                    hwidth->SetBinContent(i, width);
                }
                widthByCh[pair->second.channel] = hwidth;

            } else {
                    htot->Add(hist);
            }
            sig[i] = pair->second.signal;
            TH1 *&hpart = (sig[i] ? sigByCh : bkgByCh)[pair->second.channel];
            if (hpart == 0) {
                    hpart = (TH1*) hist->Clone();
                    hpart->SetName((sig[i] ? "total_signal" : "total_background"));
		    hpart->SetTitle(Form((sig[i] ? "Total signal in %s" : "Total background in %s"),pair->second.channel.c_str()));
                    hpart->SetDirectory(shapesByChannel[pair->second.channel]);
                    TH1 *hpart2 = (TH1*) hist->Clone(); hpart2->Reset();
                    hpart2->SetDirectory(0);
                    (sig[i] ? sigByCh2 : bkgByCh2)[pair->second.channel] = hpart2;
            } else {
                    hpart->Add(hist);
            }
            //}
        }
    }

    int totalBins = 0;

    for (IH h = totByCh.begin(), eh = totByCh.end(); h != eh; ++h) {
	totalBins +=  h->second->GetNbinsX();
    }
    if (totalBins==0) totalBins=1; // cover the case, where there is only a counting experiment without a "fake" shape, which is essentially just 1 bin. 

    //Total covariance
    TH2D* totOverall2Covar = new TH2D("overall_total_covar","Covariance signal+background",totalBins,0,totalBins,totalBins,0,totalBins);
    totOverall2Covar->GetXaxis()->SetTitle("Bin number");
    totOverall2Covar->GetYaxis()->SetTitle("Bin number");
    totOverall2Covar->SetDirectory(0);
    //Total background
    TH1D* totOverall = new TH1D("total_overall","signal+background",totalBins,0,totalBins);
    totOverall->SetDirectory(0);
    TH1D* wdtOverall = new TH1D("total_bin_width","Bin widths",totalBins,0,totalBins);
    wdtOverall->SetDirectory(0);
    TH1D* bkgOverall = new TH1D("total_background","Total background",totalBins,0,totalBins);
    bkgOverall->SetDirectory(0);
    TH1D* sigOverall = new TH1D("total_signal","Total signal",totalBins,0,totalBins);
    sigOverall->SetDirectory(0);
    TH1D* datOverallHist = new TH1D("total_data","Total data",totalBins,0,totalBins);
    datOverallHist->SetDirectory(0);

    int iBinOverall = 1;
    //Map to hold info on bins across channels
    std::map<TString,int> binMap;
    for (IH h = totByCh.begin(), eh = totByCh.end(); h != eh; ++h){
        for (int iBin = 0; iBin < h->second->GetNbinsX(); iBin++,iBinOverall++){
            TString label = Form("%s_%d",h->first.c_str(),iBin);
            binMap[label] = iBinOverall;
            totOverall->GetXaxis()->SetBinLabel(iBinOverall,label);
            totOverall->SetBinContent(iBinOverall,h->second->GetBinContent(iBin+1));
            wdtOverall->GetXaxis()->SetBinLabel(iBinOverall,label);
            wdtOverall->SetBinContent(iBinOverall,widthByCh[h->first]->GetBinContent(iBin+1));

            datOverallHist->GetXaxis()->SetBinLabel(iBinOverall,label);
            double x,y;
            datByCh[h->first]->GetPoint(iBin,x,y);
            datOverallHist->SetBinContent(iBinOverall,y);

            sigOverall->GetXaxis()->SetBinLabel(iBinOverall,label);
            //For signal have to deal with empty channels
            std::map<std::string,TH1*>::iterator iH = sigByCh.find(h->first);
            if (iH != sigByCh.end()){
            sigOverall->SetBinContent(iBinOverall,iH->second->GetBinContent(iBin+1));
            }

            bkgOverall->GetXaxis()->SetBinLabel(iBinOverall,label);
            bkgOverall->SetBinContent(iBinOverall,bkgByCh[h->first]->GetBinContent(iBin+1));

            totOverall2Covar->GetXaxis()->SetBinLabel(iBinOverall,label);
            totOverall2Covar->GetYaxis()->SetBinLabel(iBinOverall,label);
        }
    }

    TGraphAsymmErrors * datOverall = utils::makeDataGraph(datOverallHist);
    datOverall->SetNameTitle("total_data","Total data");
    delete datOverallHist;

    if (saveWithUncertainties_) {
        int ntoys = numToysForShapes_;

        if ( verbose > 0 ) Logger::instance().log(std::string(Form("FitDiagnostics.cc: %d -- Generating toy data for evaluating per-bin uncertainties and covariances with post-fit nuisance parameters with %d toys",__LINE__,ntoys)),Logger::kLogLevelInfo,__func__);

        sampler.generate(ntoys);
        std::auto_ptr<RooArgSet> params(pdf->getParameters(obs));
        // prepare histograms for running sums
        std::map<std::string,TH1*> totByCh1, sigByCh1, bkgByCh1;
        for (IH h = totByCh.begin(), eh = totByCh.end(); h != eh; ++h) totByCh1[h->first] = (TH1*) h->second->Clone();
        for (IH h = sigByCh.begin(), eh = sigByCh.end(); h != eh; ++h) sigByCh1[h->first] = (TH1*) h->second->Clone();
        for (IH h = bkgByCh.begin(), eh = bkgByCh.end(); h != eh; ++h) bkgByCh1[h->first] = (TH1*) h->second->Clone();
        for (int t = 0; t < ntoys; ++t) {
            // zero out partial sums
            for (IH h = totByCh1.begin(), eh = totByCh1.end(); h != eh; ++h) h->second->Reset();
            for (IH h = sigByCh1.begin(), eh = sigByCh1.end(); h != eh; ++h) h->second->Reset();
            for (IH h = bkgByCh1.begin(), eh = bkgByCh1.end(); h != eh; ++h) h->second->Reset();
            // randomize numbers
            params->assignValueOnly( sampler.get(t) );
            for (pair = bg, i = 0; pair != ed; ++pair, ++i) { 
                // add up deviations in numbers for each channel
                sumx2[i] += std::pow(pair->second.norm->getVal() - vals[i], 2);  
                if (saveShapes_ && pair->second.obs.getSize() == 1) {
                    // and also deviations in the shapes
                    RooRealVar *x = (RooRealVar*)pair->second.obs.at(0);
                    std::auto_ptr<TH1> hist(pair->second.pdf->createHistogram(pair->second.pdf->GetName(), *x,
                      pair->second.isfunc ? RooFit::Extended(false) : RooCmdArg::none()));
                    hist->Scale(pair->second.norm->getVal() / hist->Integral("width"));
                    for (int b = 1; b <= bins[i]; ++b) {
                        shapes2[i]->AddBinContent(b, std::pow(hist->GetBinContent(b) - shapes[i]->GetBinContent(b), 2));
                    }
                    // and cumulate in the total for this toy as well
                    totByCh1[pair->second.channel]->Add(&*hist);
                    (sig[i] ? sigByCh1 : bkgByCh1)[pair->second.channel]->Add(&*hist);
                }
            }
            // now add up the deviations within channels in this toy
            for (IH h = totByCh1.begin(), eh = totByCh1.end(); h != eh; ++h) {
                TH1 *target = totByCh2[h->first], *reference = totByCh[h->first];
                TH2 *targetCovar = totByCh2Covar[h->first];
                for (int b = 1, nb = target->GetNbinsX(); b <= nb; ++b) {
		    double deltaBi = h->second->GetBinContent(b) - reference->GetBinContent(b);
		    target->AddBinContent(b, std::pow(deltaBi, 2));
		    int binX = b - 1;
		    TString xLabel = Form("%s_%d",h->first.c_str(),binX);
		    for (int bj = 1;bj <= b; bj++) {
			int binY = bj - 1;
			TString yLabel = Form("%s_%d",h->first.c_str(),binY);
			double deltaBj = h->second->GetBinContent(bj) - reference->GetBinContent(bj);
			targetCovar->AddBinContent(targetCovar->GetBin(b,bj),deltaBj*deltaBi);  // covariance
			totOverall2Covar->Fill(xLabel,yLabel,deltaBj*deltaBi);
			if (b != bj) {
			    targetCovar->AddBinContent(targetCovar->GetBin(bj,b),deltaBj*deltaBi);  // covariance
			    totOverall2Covar->Fill(yLabel,xLabel,deltaBj*deltaBi);
			}
		    }
		}
	    }
            // deviations across channels in this toy
	    if (saveOverallShapes_){
		for (IH h = totByCh1.begin(), eh = totByCh1.end(); h != eh; ++h) {
		    TH1 *reference = totByCh[h->first];
		    for (IH h2 = totByCh1.begin();h2 != h; ++h2) {
			TH1 *reference2 = totByCh[h2->first];
			for (int b = 1, nb = reference->GetNbinsX(); b <= nb; ++b) {
			    double deltaBi = h->second->GetBinContent(b) - reference->GetBinContent(b);
			    int binX = b - 1;
			    TString xLabel = Form("%s_%d",h->first.c_str(),binX);
			    for (int bj = 1, nb2 = reference2->GetNbinsX(); bj <= nb2; ++bj) {
				int binY = bj - 1;
				double deltaBj = h2->second->GetBinContent(bj) - reference2->GetBinContent(bj);
				TString yLabel = Form("%s_%d",h2->first.c_str(),binY);
				totOverall2Covar->Fill(xLabel,yLabel,deltaBj*deltaBi);
				totOverall2Covar->Fill(yLabel,xLabel,deltaBj*deltaBi);
			    }
			}
		    }
		}
	    }
            for (IH h = sigByCh1.begin(), eh = sigByCh1.end(); h != eh; ++h) {
                TH1 *target = sigByCh2[h->first], *reference = sigByCh[h->first];
                for (int b = 1, nb = target->GetNbinsX(); b <= nb; ++b) {
                    target->AddBinContent(b, std::pow(h->second->GetBinContent(b) - reference->GetBinContent(b), 2));
                }
            }           
            for (IH h = bkgByCh1.begin(), eh = bkgByCh1.end(); h != eh; ++h) {
                TH1 *target = bkgByCh2[h->first], *reference = bkgByCh[h->first];
                for (int b = 1, nb = target->GetNbinsX(); b <= nb; ++b) {
                    target->AddBinContent(b, std::pow(h->second->GetBinContent(b) - reference->GetBinContent(b), 2));
                }
            }           
        } // end of the toy loop
        // now take square roots and such
        for (pair = bg, i = 0; pair != ed; ++pair, ++i) {
            sumx2[i] = sqrt(sumx2[i]/ntoys);
            if (shapes2[i]) {
                for (int b = 1; b <= bins[i]; ++b) {
                    shapes[i]->SetBinError(b, std::sqrt(shapes2[i]->GetBinContent(b)/ntoys));
                }
                delete shapes2[i]; shapes2[i] = 0;
            }

        }
        // and the same for the total histograms
        for (IH h = totByCh.begin(), eh = totByCh.end(); h != eh; ++h) {
            TH1 *sum2   = totByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
		TString xLabel = Form("%s_%d",h->first.c_str(),b-1);
		totOverall->SetBinError(binMap[xLabel],std::sqrt(sum2->GetBinContent(b)/ntoys));
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete totByCh1[h->first];
	}
	// same for covariance matrix 
	for (int b = 1, nb = totOverall2Covar->GetNbinsX(); b <= nb; ++b) {
	    for (int bj = 1, nbj = totOverall2Covar->GetNbinsY(); bj <= nbj; ++bj) {    
		totOverall2Covar->SetBinContent(b,bj, (totOverall2Covar->GetBinContent(b,bj)/ntoys));
	    }
	}
        for (IH2 h = totByCh2Covar.begin(), eh = totByCh2Covar.end(); h != eh; ++h) {
            TH2 *covar2 = h->second;
            for (int b = 1, nb = covar2->GetNbinsX(); b <= nb; ++b) {
              for (int bj = 1, nbj = covar2->GetNbinsY(); bj <= nbj; ++bj) {    
		  h->second->SetBinContent(b,bj, (covar2->GetBinContent(b,bj)/ntoys));
	      }
	    }
	}

        for (IH h = sigByCh.begin(), eh = sigByCh.end(); h != eh; ++h) {
            TH1 *sum2 = sigByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
		TString xLabel = Form("%s_%d",h->first.c_str(),b-1);
		sigOverall->SetBinError(binMap[xLabel],std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete sigByCh1[h->first];
        }
        for (IH h = bkgByCh.begin(), eh = bkgByCh.end(); h != eh; ++h) {
            TH1 *sum2 = bkgByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
		TString xLabel = Form("%s_%d",h->first.c_str(),b-1);
		bkgOverall->SetBinError(binMap[xLabel],std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete bkgByCh1[h->first];
        }
        totByCh1.clear(); totByCh2.clear(); sigByCh1.clear(); sigByCh2.clear(); bkgByCh1.clear(); bkgByCh2.clear();
	totByCh2.clear();
        // finally reset parameters
        params->assignValueOnly( sampler.centralValues() );
    }

    for (IG h = datByCh.begin(), eh = datByCh.end(); h != eh; ++h) {
	shapesByChannel[h->first]->WriteTObject(h->second);
    }
    for (pair = bg, i = 0; pair != ed; ++pair, ++i) {
        RooRealVar *val = new RooRealVar((oldNormNames_ ? pair->first : pair->second.channel+"/"+pair->second.process).c_str(), "", vals[i]);
        val->setError(sumx2[i]);
        out.addOwned(*val); 
        if (shapes[i]) shapesByChannel[pair->second.channel]->WriteTObject(shapes[i]);
    }
    if (fOut) {
        fOut->WriteTObject(&out, (std::string("norm")+postfix).c_str());
        for (IH h = totByCh.begin(), eh = totByCh.end(); h != eh; ++h) { shapesByChannel[h->first]->WriteTObject(h->second); }
        for (IH h = sigByCh.begin(), eh = sigByCh.end(); h != eh; ++h) { shapesByChannel[h->first]->WriteTObject(h->second); }
        for (IH h = bkgByCh.begin(), eh = bkgByCh.end(); h != eh; ++h) { shapesByChannel[h->first]->WriteTObject(h->second); }
        for (IH2 h = totByCh2Covar.begin(), eh = totByCh2Covar.end(); h != eh; ++h) { shapesByChannel[h->first]->WriteTObject(h->second); }
	//Save total shapes or clean up if not keeping
	if (saveShapes_) shapeDir->cd();

	if (saveShapes_ && saveOverallShapes_){
	    wdtOverall->Write();
	    totOverall->Write();
	    sigOverall->Write();
	    datOverall->Write();
	    bkgOverall->Write();
	}
	else{
	    delete totOverall;
	    delete sigOverall;
	    delete datOverall;
	    delete bkgOverall;
	}
	if (saveWithUncertainties_ && saveOverallShapes_){
	    totOverall2Covar->Write();
	}
	else{
	    delete totOverall2Covar;
	}
    }
}


//void FitDiagnostics::setFitResultTrees(const RooArgSet *args, std::vector<double> *vals){
void FitDiagnostics::setFitResultTrees(const RooArgSet *args, double * vals){
	
         TIterator* iter(args->createIterator());
	 int count=0;
	 
         for (TObject *a = iter->Next(); a != 0; a = iter->Next()) { 
                 RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);        
		 //std::string name = rrv->GetName();
		 vals[count]=rrv->getVal();
		 count++;
         }
	 delete iter;
	 return;
}

void FitDiagnostics::setShapesFitResultTrees(std::map<std::string,ShapeAndNorm> &snm, double * vals){
    int iBinOverall = 1;	
    std::map<std::string,ShapeAndNorm>::const_iterator bg = snm.begin(), ed = snm.end(), pair; int i;
    for (pair = bg, i = 0; pair != ed; ++pair, ++i) {  
        if (pair->second.obs.getSize() == 1) {
            RooRealVar *x = (RooRealVar*)pair->second.obs.at(0);
            TH1* hist = pair->second.pdf->createHistogram(Form("%d",iBinOverall), *x, pair->second.isfunc ? RooFit::Extended(false) : RooCmdArg::none());
            hist->Scale(pair->second.norm->getVal() / hist->Integral("width"));
	    for (int iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
		vals[iBinOverall-1] = hist->GetBinContent(iBin);
		iBinOverall++;
	    }
	    delete hist;
	}
    }
	 return;
}
void FitDiagnostics::resetFitResultTrees(bool withSys){
	
	 for (int count = 0; count < overallNorms_; count ++){
	     processNormalizations_[count] = -999;
         }
	 for (int count = 0; count < overallCons_; count ++){
	     globalObservables_[count] = -999;
         }
	 for (int count = 0; count < overallNuis_; count ++){
	     nuisanceParameters_[count] = -999;
         }
	 for (int count = 0; count < overallBins_; count ++){
	     processNormalizationsShapes_[count] = -999;
         }
	 return;
}
void FitDiagnostics::setNormsFitResultTrees(const RooArgSet *args, double * vals){
	
         TIterator* iter(args->createIterator());
	 int count=0;
	 
         for (TObject *a = iter->Next(); a != 0; a = iter->Next()) { 
                 RooRealVar *rcv = dynamic_cast<RooRealVar *>(a);   
		 // std::cout << "index " << count << ", Name " << rcv->GetName() << ", val " <<  rcv->getVal() << std::endl;
		 //std::string name = rcv->GetName();
		 vals[count]=rcv->getVal();
		 count++;
         }
	 delete iter;
	 return;
}

void FitDiagnostics::createFitResultTrees(const RooStats::ModelConfig &mc, bool withSys, bool savePerToy){

	 // Initiate the arrays to store parameters

	 std::string poiName = (dynamic_cast<RooRealVar *>(mc.GetParametersOfInterest()->first()))->GetName();

	 // create TTrees to store fit results:
	 t_prefit_  = new TTree("tree_prefit","tree_prefit");
	 t_fit_b_  = new TTree("tree_fit_b","tree_fit_b");
	 t_fit_sb_ = new TTree("tree_fit_sb","tree_fit_sb");

    	 t_fit_b_->Branch("fit_status",&fitStatus_,"fit_status/I");
   	 t_fit_sb_->Branch("fit_status",&fitStatus_,"fit_status/I");

	 t_fit_b_->Branch(poiName.c_str(),&mu_,Form("%s/D",poiName.c_str()));
	 t_fit_sb_->Branch(poiName.c_str(),&mu_,Form("%s/D",poiName.c_str()));
	 
	 t_fit_b_->Branch(Form("%sErr",poiName.c_str()),&muErr_,Form("%sErr/D",poiName.c_str()));
	 t_fit_sb_->Branch(Form("%sErr",poiName.c_str()),&muErr_,Form("%sErr/D",poiName.c_str()));

	 t_fit_b_->Branch(Form("%sLoErr",poiName.c_str()),&muLoErr_,Form("%sLoErr/D",poiName.c_str()));
	 t_fit_sb_->Branch(Form("%sLoErr",poiName.c_str()),&muLoErr_,Form("%sLoErr/D",poiName.c_str()));

	 t_fit_b_->Branch(Form("%sHiErr",poiName.c_str()),&muHiErr_,Form("%sHiErr/D",poiName.c_str()));
	 t_fit_sb_->Branch(Form("%sHiErr",poiName.c_str()),&muHiErr_,Form("%sHiErr/D",poiName.c_str()));

	 t_fit_b_->Branch("numbadnll",&numbadnll_,"numbadnll/I");
	 t_fit_sb_->Branch("numbadnll",&numbadnll_,"numbadnll/I");

	 t_fit_b_->Branch("nll_min",&nll_bonly_,"nll_min/D");
	 t_fit_sb_->Branch("nll_min",&nll_sb_,"nll_min/D");

	 t_fit_sb_->Branch("nll_nll0",&nll_nll0_,"nll_nll0/D");

	 int count=0; 
         // fill the maps for the nuisances, and global observables
         RooArgSet *norms= new RooArgSet();
         //getNormalizationsSimple(mc.GetPdf(), *mc.GetObservables(), *norms);  <-- This is useless as the order is messed up !

         std::map<std::string,ShapeAndNorm> snm;
         int totalBins = 0;
         typedef std::map<std::string,ShapeAndNorm>::const_iterator IT;
         if(savePerToy){
           getShapesAndNorms(mc.GetPdf(),*mc.GetObservables(), snm, "");
           IT bg = snm.begin(), ed = snm.end(), pair; int i;
           for (pair = bg, i = 0; pair != ed; ++pair, ++i) {
             RooRealVar *val = new RooRealVar(pair->first.c_str(), "", 0.);
             //val->setError(sumx2[i]);
             norms->addOwned(*val);
             RooRealVar *x = (RooRealVar*)pair->second.obs.at(0);
             if (pair->second.obs.getSize() == 1) {
                TH1* hist = pair->second.pdf->createHistogram(Form("%d",totalBins), *x, pair->second.isfunc ? RooFit::Extended(false) : RooCmdArg::none());
	        totalBins += hist->GetNbinsX();
                delete hist;
             }
           }
         }

         overallBins_ = totalBins;
         overallNorms_ = norms->getSize();

         processNormalizations_ = new double[norms->getSize()];
         processNormalizationsShapes_ = new double[totalBins];

	 // If no systematic (-S 0), then don't make nuisance trees
	 if (withSys){
          const RooArgSet *cons = mc.GetGlobalObservables();
          const RooArgSet *nuis = mc.GetNuisanceParameters();
 	  globalObservables_ = new double[cons->getSize()];
	  overallCons_ = cons->getSize();
	  nuisanceParameters_= new double[nuis->getSize()];
	  overallNuis_ = nuis->getSize();

          TIterator* iter_c(cons->createIterator());
          for (TObject *a = iter_c->Next(); a != 0; a = iter_c->Next()) { 
                 RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rrv->GetName();
		 globalObservables_[count]=0;
		 t_fit_sb_->Branch(name.c_str(),&(globalObservables_[count]),Form("%s/D",name.c_str()));
		 t_fit_b_->Branch(name.c_str(),&(globalObservables_[count]),Form("%s/D",name.c_str()));
		 t_prefit_->Branch(name.c_str(),&(globalObservables_[count]),Form("%s/D",name.c_str()));
		 count++;
	  }         
	  count = 0;
          TIterator* iter_n(nuis->createIterator());
          for (TObject *a = iter_n->Next(); a != 0; a = iter_n->Next()) { 
                 RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rrv->GetName();
		 nuisanceParameters_[count] = 0;
		 t_fit_sb_->Branch(name.c_str(),&(nuisanceParameters_[count]),Form("%s/D",name.c_str()));
		 t_fit_b_->Branch(name.c_str(),&(nuisanceParameters_[count]),Form("%s/D",name.c_str()));
		 t_prefit_->Branch(name.c_str(),&(nuisanceParameters_[count]),Form("%s/D",name.c_str()));
		 count++;
          }

         }

	 count = 0;
         TIterator* iter_no(norms->createIterator());
         for (TObject *a = iter_no->Next(); a != 0; a = iter_no->Next()) { 
                 RooRealVar *rcv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rcv->GetName();
		 processNormalizations_[count] = -999;
		 //std::cout << " Creating the TREE -- " << count << ", Branch Name =  " << name << ", Param name " << rcv->GetName() << std::endl; 
		 t_fit_sb_->Branch(name.c_str(),&(processNormalizations_[count]),Form("%s/D",name.c_str()));
		 t_fit_b_->Branch(name.c_str(),&(processNormalizations_[count]),Form("%s/D",name.c_str()));
		 t_prefit_->Branch(name.c_str(),&(processNormalizations_[count]),Form("%s/D",name.c_str()));
		 count++;
         }

         if(savePerToy){
           count = 0;
           IT itb = snm.begin(), ite = snm.end(), npair; int j;
           itb = snm.begin(), ite = snm.end();
           for (npair = itb, j = 0; npair != ite; ++npair, ++j) {
             if (npair->second.obs.getSize() == 1) {
               RooRealVar *x = (RooRealVar*)npair->second.obs.at(0);
               TH1* hist = npair->second.pdf->createHistogram(Form("%d",count), *x, npair->second.isfunc ? RooFit::Extended(false) : RooCmdArg::none());
               int bins = hist->GetNbinsX();
                for (int iBin = 1; iBin <= bins; iBin++){
                     processNormalizationsShapes_[count] = -999;
                     TString label = Form("%s_%d",npair->first.c_str(),iBin);
                     t_fit_sb_->Branch(label,&(processNormalizationsShapes_[count]),label+"/D");
                     t_fit_b_->Branch(label,&(processNormalizationsShapes_[count]),label+"/D");
                     t_prefit_->Branch(label,&(processNormalizationsShapes_[count]),label+"/D");
                     count++;
                }
               delete hist;
             }
           }
         }

         delete norms;

	 //std::cout << "Created Branches for toy diagnostics" <<std::endl;
         return;	
}

FitDiagnostics::ToySampler::ToySampler(RooAbsPdf *pdf, const RooArgSet *nuisances) :
    pdf_(pdf),
    data_(0)
{
    nuisances->snapshot(snapshot_);
}

FitDiagnostics::ToySampler::~ToySampler() 
{
    delete data_;
}

void FitDiagnostics::ToySampler::generate(int ntoys) {
    delete data_;
    data_ = pdf_->generate(snapshot_, ntoys);
}

const RooAbsCollection & FitDiagnostics::ToySampler::get(int itoy) { 
    return *data_->get(itoy); 
}

const RooAbsCollection & FitDiagnostics::ToySampler::centralValues() {
    return snapshot_;
}
