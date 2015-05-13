#include "HiggsAnalysis/CombinedLimit/interface/MaxLikelihoodFit.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
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
#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfileLikelihood.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"


#include <Math/MinimizerOptions.h>

using namespace RooStats;

std::string MaxLikelihoodFit::name_ = "";
std::string MaxLikelihoodFit::minos_ = "poi";
std::string MaxLikelihoodFit::out_ = ".";
bool        MaxLikelihoodFit::makePlots_ = false;
bool        MaxLikelihoodFit::saveWorkspace_ = false;
float       MaxLikelihoodFit::rebinFactor_ = 1.0;
std::string MaxLikelihoodFit::signalPdfNames_     = "shapeSig*";
std::string MaxLikelihoodFit::backgroundPdfNames_ = "shapeBkg*";
bool        MaxLikelihoodFit::saveNormalizations_ = false;
bool        MaxLikelihoodFit::oldNormNames_ = false;
bool        MaxLikelihoodFit::saveShapes_ = false;
bool        MaxLikelihoodFit::saveWithUncertainties_ = false;
bool        MaxLikelihoodFit::justFit_ = false;
bool        MaxLikelihoodFit::skipBOnlyFit_ = false;
bool        MaxLikelihoodFit::noErrors_ = false;
bool        MaxLikelihoodFit::reuseParams_ = false;
bool        MaxLikelihoodFit::customStartingPoint_ = false;


MaxLikelihoodFit::MaxLikelihoodFit() :
    FitterAlgoBase("MaxLikelihoodFit specific options"),
    globalObservables_(0),
    nuisanceParameters_(0),
    processNormalizations_(0)
{
    options_.add_options()
        ("minos",              boost::program_options::value<std::string>(&minos_)->default_value(minos_), "Compute MINOS errors for: 'none', 'poi', 'all'")
        ("out",                boost::program_options::value<std::string>(&out_)->default_value(out_), "Directory to put output in")
        ("plots",              "Make plots")
        ("rebinFactor",        boost::program_options::value<float>(&rebinFactor_)->default_value(rebinFactor_), "Rebin by this factor before plotting (does not affect fitting!)")
        ("signalPdfNames",     boost::program_options::value<std::string>(&signalPdfNames_)->default_value(signalPdfNames_), "Names of signal pdfs in plots (separated by ,)")
        ("backgroundPdfNames", boost::program_options::value<std::string>(&backgroundPdfNames_)->default_value(backgroundPdfNames_), "Names of background pdfs in plots (separated by ',')")
        ("saveNormalizations",  "Save post-fit normalizations of all components of the pdfs")
        ("oldNormNames",  "Name the normalizations as in the workspace, and not as channel/process")
//        ("saveWorkspace",       "Save post-fit pdfs and data to MaxLikelihoodFitResults.root")
        ("saveShapes",  "Save post-fit binned shapes")
        ("saveWithUncertainties",  "Save also post-fit uncertainties on the shapes and normalizations (from resampling the covariance matrix)")
        ("justFit",  "Just do the S+B fit, don't do the B-only one, don't save output file")
        ("skipBOnlyFit",  "Skip the B-only fit (do only the S+B fit)")
        ("noErrors",  "Don't compute uncertainties on the best fit value")
        ("initFromBonly",  "Use the values of the nuisance parameters from the background only fit as the starting point for the s+b fit")
        ("customStartingPoint",  "Don't set the signal model parameters to zero before the fit")
   ;

    // setup a few defaults
    nToys=0; fitStatus_=0; mu_=0; muLoErr_=0; muHiErr_=0; numbadnll_=-1; nll_nll0_=-1; nll_bonly_=-1; nll_sb_=-1;
}

MaxLikelihoodFit::~MaxLikelihoodFit(){
   // delete the Arrays used to fill the trees;
   delete globalObservables_;
   delete nuisanceParameters_;
   delete processNormalizations_;
}

void MaxLikelihoodFit::setToyNumber(const int iToy){
	currentToy_ = iToy;
}
void MaxLikelihoodFit::setNToys(const int iToy){
	nToys = iToy;
}
void MaxLikelihoodFit::applyOptions(const boost::program_options::variables_map &vm) 
{
    applyOptionsBase(vm);
    makePlots_ = vm.count("plots");
    name_ = vm["name"].defaulted() ?  std::string() : vm["name"].as<std::string>();
    saveShapes_  = vm.count("saveShapes");
    saveNormalizations_  = saveShapes_ || vm.count("saveNormalizations");
    oldNormNames_  = vm.count("oldNormNames");
    saveWithUncertainties_  = vm.count("saveWithUncertainties");
    saveWorkspace_ = vm.count("saveWorkspace");
    justFit_  = vm.count("justFit");
    skipBOnlyFit_ = vm.count("skipBOnlyFit");
    noErrors_ = vm.count("noErrors");
    reuseParams_ = vm.count("initFromBonly");
    customStartingPoint_ = vm.count("customStartingPoint");
     
    if (justFit_) { out_ = "none"; makePlots_ = false; saveNormalizations_ = false; reuseParams_ = false;}
    // For now default this to true;
}

bool MaxLikelihoodFit::runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {

  if (reuseParams_ && minos_!="none"){
	std::cout << "Cannot reuse b-only fit params when running minos. Parameters will be reset when running S+B fit"<<std::endl;
	reuseParams_=false;
  }


  if (!justFit_ && out_ != "none"){
	if (currentToy_ < 1){
		fitOut.reset(TFile::Open((out_+"/mlfit"+name_+".root").c_str(), "RECREATE")); 
		createFitResultTrees(*mc_s,withSystematics);
	}
  }

  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());

  TCanvas *c1 = 0;
  if (makePlots_) {
      utils::tdrStyle();
      c1 = new TCanvas("c1","c1");
  }

  // Make pre-plots before the fit
  r->setVal(preFitValue_);
  if (makePlots_) {
      std::vector<RooPlot *> plots = utils::makePlots(*mc_s->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_);
      for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
          (*it)->Draw(); 
          c1->Print((out_+"/"+(*it)->GetName()+"_prefit.png").c_str());
          if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_prefit").c_str());
      }
  }


  // Determine pre-fit values of nuisance parameters
  if (currentToy_ < 1){
    const RooArgSet *nuis      = mc_s->GetNuisanceParameters();
    const RooArgSet *globalObs = mc_s->GetGlobalObservables();
    if (!justFit_ && nuis && globalObs ) {
      std::auto_ptr<RooAbsPdf> nuisancePdf(utils::makeNuisancePdf(*mc_s));
      std::auto_ptr<RooDataSet> globalData(new RooDataSet("globalData","globalData", *globalObs));
      globalData->add(*globalObs);
      RooFitResult *res_prefit = 0;
      {     
            CloseCoutSentry sentry(verbose < 2);
            res_prefit = nuisancePdf->fitTo(*globalData,
            RooFit::Save(1),
            RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()),
            RooFit::Strategy(minimizerStrategy_),
            RooFit::Minos(minos_ == "all")
            );
      }
      if (fitOut.get() ) fitOut->WriteTObject(res_prefit, "nuisances_prefit_res");
      if (fitOut.get() ) fitOut->WriteTObject(nuis->snapshot(), "nuisances_prefit");

      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_prefit");
          ToySampler sampler(&*nuisancePdf, nuis);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_prefit");
          delete norms;
      }

      nuisancePdf.reset();
      globalData.reset();
      delete res_prefit;
    } else if (nuis) {
      if (fitOut.get() ) fitOut->WriteTObject(nuis->snapshot(), "nuisances_prefit");
    }
  }
  
  RooFitResult *res_b = 0, *res_s = 0;
  const RooCmdArg &constCmdArg_s = withSystematics  ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooFit::NumCPU(1); // use something dummy 
  const RooCmdArg &minosCmdArg = minos_ == "poi" ?  RooFit::Minos(*mc_s->GetParametersOfInterest())   : RooFit::Minos(minos_ != "none"); 
  w->loadSnapshot("clean");
  if (!customStartingPoint_) r->setVal(0.0); 
  r->setConstant(true);

  // Setup Nll before calling fits;
  if (currentToy_<1) nll.reset(mc_s->GetPdf()->createNLL(data,constCmdArg_s,RooFit::Extended(mc_s->GetPdf()->canBeExtended())));
  // Get the nll value on the prefit
  double nll0 = nll->getVal();

  if (justFit_ || skipBOnlyFit_ ) { 
    // skip b-only fit
  } else if (minos_ != "all") {
    RooArgList minos; 
    res_b = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/true,/*reuseNLL*/ true); 
    nll_bonly_=nll->getVal()-nll0;   
  } else {
    CloseCoutSentry sentry(verbose < 2);
    res_b = mc_s->GetPdf()->fitTo(data, 
            RooFit::Save(1), 
            RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()), 
            RooFit::Strategy(minimizerStrategy_),
            RooFit::Extended(mc_s->GetPdf()->canBeExtended()), 
            constCmdArg_s, minosCmdArg
            );
    if (res_b) nll_bonly_ = nll->getVal() - nll0;

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

      if (makePlots_) {
          std::vector<RooPlot *> plots = utils::makePlots(*mc_b->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_);
          for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
              c1->cd(); (*it)->Draw(); 
              c1->Print((out_+"/"+(*it)->GetName()+"_fit_b.png").c_str());
              if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_fit_b").c_str());
          }
      }

      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_b");
          CovarianceReSampler sampler(res_b);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_fit_b");
          setNormsFitResultTrees(norms,processNormalizations_);
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
  }
  else {
	fitStatus_=-1;
	numbadnll_=-1;	
  }
  mu_=r->getVal();
  if (t_fit_b_) t_fit_b_->Fill();
  // no longer need res_b
  delete res_b;

  if (!reuseParams_) w->loadSnapshot("clean"); // Reset, also ensures nll_prefit is same in call to doFit for b and s+b
  r->setVal(preFitValue_); r->setConstant(false); 
  if (minos_ != "all") {
    RooArgList minos; if (minos_ == "poi") minos.add(*r);
    res_s = doFit(*mc_s->GetPdf(), data, minos, constCmdArg_s, /*hesse=*/!noErrors_,/*reuseNLL*/ true); 
    nll_sb_ = nll->getVal()-nll0;
  } else {
    CloseCoutSentry sentry(verbose < 2);
    res_s = mc_s->GetPdf()->fitTo(data, 
            RooFit::Save(1), 
            RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()), 
            RooFit::Strategy(minimizerStrategy_),
            RooFit::Extended(mc_s->GetPdf()->canBeExtended()), 
            constCmdArg_s, minosCmdArg
            );
    if (res_s) nll_sb_= nll->getVal()-nll0;

  }
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

	 // Additionally store the nll_sb - nll_bonly (=0.5*q0)
	 nll_nll0_ =  nll_sb_ -  nll_bonly_;
      }

      if (makePlots_) {
          std::vector<RooPlot *> plots = utils::makePlots(*mc_s->GetPdf(), data, signalPdfNames_.c_str(), backgroundPdfNames_.c_str(), rebinFactor_);
          for (std::vector<RooPlot *>::iterator it = plots.begin(), ed = plots.end(); it != ed; ++it) {
              c1->cd(); (*it)->Draw(); 
              c1->Print((out_+"/"+(*it)->GetName()+"_fit_s.png").c_str());
              if (fitOut.get() && currentToy_< 1) fitOut->WriteTObject(*it, (std::string((*it)->GetName())+"_fit_s").c_str());
          }
      }

      if (saveNormalizations_) {
          RooArgSet *norms = new RooArgSet();
          norms->setName("norm_fit_s");
          CovarianceReSampler sampler(res_s);
          getNormalizations(mc_s->GetPdf(), *mc_s->GetObservables(), *norms, sampler, currentToy_<1 ? fitOut.get() : 0, "_fit_s");
          setNormsFitResultTrees(norms,processNormalizations_);
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
      double bestFitVal = rf->getVal();

      double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal : rf->getAsymErrorHi());
      double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal : rf->getAsymErrorLo());
      double maxError = std::max<double>(std::max<double>(hiErr, loErr), rf->getError());

      if (fabs(hiErr) < 0.001*maxError) hiErr = -bestFitVal + rf->getMax();
      if (fabs(loErr) < 0.001*maxError) loErr = +bestFitVal - rf->getMin();

      muLoErr_=loErr;
      muHiErr_=hiErr;

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
      std::cout << "\n --- MaxLikelihoodFit ---" << std::endl;
      std::cout << "Best fit " << r->GetName() << ": " << rf->getVal() << "  "<<  -loErr << "/+" << +hiErr << "  (68% CL)" << std::endl;
      if (do95_) {
        std::cout << "         " << r->GetName() << ": " << rf->getVal() << "  "<<  -loErr95 << "/+" << +hiErr95 << "  (95% CL)" << std::endl;
      }
  } else {
      std::cout << "\n --- MaxLikelihoodFit ---" << std::endl;
      std::cout << "Fit failed."  << std::endl;
  }
  if (t_fit_sb_) t_fit_sb_->Fill();

  if (currentToy_==nToys-1 || nToys==0 ) {
        
        if (fitOut.get()) {	
		fitOut->cd();
		t_fit_sb_->Write(); t_fit_b_->Write();
		fitOut.release()->Close();
	}

  } 
  bool fitreturn = (res_s!=0);
  delete res_s;

  if(saveWorkspace_){
	  RooWorkspace *ws = new RooWorkspace("MaxLikelihoodFitResult");
	  ws->import(*mc_s->GetPdf());
	  ws->import(data);
	  std::cout << "Saving pdfs and data to MaxLikelihoodFitResult.root" << std::endl;
	  ws->writeToFile("MaxLikelihoodFitResult.root");
  }
  std::cout << "nll S+B -> "<<nll_sb_ << "  nll B -> " << nll_bonly_ <<std::endl;
  return fitreturn;
}

void MaxLikelihoodFit::getNormalizationsSimple(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out) {
    RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(pdf);
    if (sim != 0) {
        RooAbsCategoryLValue &cat = const_cast<RooAbsCategoryLValue &>(sim->indexCat());
        for (int i = 0, n = cat.numBins((const char *)0); i < n; ++i) {
            cat.setBin(i);
            RooAbsPdf *pdfi = sim->getPdf(cat.getLabel());
            if (pdfi) getNormalizationsSimple(pdfi, obs, out);
        }        
        return;
    }
    RooProdPdf *prod = dynamic_cast<RooProdPdf *>(pdf);
    if (prod != 0) {
        RooArgList list(prod->pdfList());
        for (int i = 0, n = list.getSize(); i < n; ++i) {
            RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
            if (pdfi->dependsOn(obs)) getNormalizationsSimple(pdfi, obs, out);
        }
        return;
    }
    RooAddPdf *add = dynamic_cast<RooAddPdf *>(pdf);
    if (add != 0) {
        RooArgList list(add->coefList());
        for (int i = 0, n = list.getSize(); i < n; ++i) {
            RooAbsReal *coeff = (RooAbsReal *) list.at(i);
            out.addOwned(*(new RooRealVar(coeff->GetName(), "", coeff->getVal())));
        }
        return;
    }
}
void MaxLikelihoodFit::getShapesAndNorms(RooAbsPdf *pdf, const RooArgSet &obs, std::map<std::string,ShapeAndNorm> &out, const std::string &channel) {
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
            ShapeAndNorm &ns = out[coeff->GetName()];
            ns.norm = coeff;
            ns.pdf = (RooAbsPdf*) plist.at(i);
            ns.channel = (coeff->getStringAttribute("combine.channel") ? coeff->getStringAttribute("combine.channel") : channel.c_str());
            ns.process = (coeff->getStringAttribute("combine.process") ? coeff->getStringAttribute("combine.process") : ns.norm->GetName());
            ns.signal  = (coeff->getStringAttribute("combine.process") ? coeff->getAttribute("combine.signal") : (strstr(ns.norm->GetName(),"shapeSig") != 0));
            std::auto_ptr<RooArgSet> myobs(ns.pdf->getObservables(obs));
            ns.obs.add(*myobs);
        }
        return;
    }
}

void MaxLikelihoodFit::getNormalizations(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out, NuisanceSampler & sampler, TDirectory *fOut, const std::string &postfix) {
    // fill in a map
    std::map<std::string,ShapeAndNorm> snm;
    getShapesAndNorms(pdf,obs, snm, "");
    typedef std::map<std::string,ShapeAndNorm>::const_iterator IT;
    typedef std::map<std::string,TH1*>::const_iterator IH;
    // create directory structure for shapes
    TDirectory *shapeDir = fOut && saveShapes_ ? fOut->mkdir((std::string("shapes")+postfix).c_str()) : 0;
    std::map<std::string,TDirectory*> shapesByChannel;
    if (shapeDir) {
        for (IT it = snm.begin(), ed = snm.end(); it != ed; ++it) {
            TDirectory *& sub = shapesByChannel[it->second.channel];
            if (sub == 0) sub = shapeDir->mkdir(it->second.channel.c_str());
        }
    }
    // now let's start with the central values
    std::vector<double> vals(snm.size(), 0.), sumx2(snm.size(), 0.);
    std::vector<TH1*>   shapes(snm.size(), 0), shapes2(snm.size(), 0);
    std::vector<int>    bins(snm.size(), 0), sig(snm.size(), 0);
    std::map<std::string,TH1*> totByCh, totByCh2, sigByCh, sigByCh2, bkgByCh, bkgByCh2;
    IT bg = snm.begin(), ed = snm.end(), pair; int i;
    for (pair = bg, i = 0; pair != ed; ++pair, ++i) {  
        vals[i] = pair->second.norm->getVal();
        //out.addOwned(*(new RooConstVar(pair->first.c_str(), "", pair->second.norm->getVal())));
        if (fOut != 0 && saveShapes_ && pair->second.obs.getSize() == 1) {
            RooRealVar *x = (RooRealVar*)pair->second.obs.at(0);
            TH1* hist = pair->second.pdf->createHistogram("", *x);
            hist->SetNameTitle(pair->second.process.c_str(), (pair->second.process+" in "+pair->second.channel).c_str());
            hist->Scale(vals[i] / hist->Integral("width"));
            hist->SetDirectory(shapesByChannel[pair->second.channel]);
            shapes[i] = hist;
            if (saveWithUncertainties_) {
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
            }
        }
    }
    if (saveWithUncertainties_) {
        int ntoys = 200;
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
                    std::auto_ptr<TH1> hist(pair->second.pdf->createHistogram(pair->second.pdf->GetName(), *x));
                    hist->Scale(pair->second.norm->getVal() / hist->Integral("width"));
                    for (int b = 1; b <= bins[i]; ++b) {
                        shapes2[i]->AddBinContent(b, std::pow(hist->GetBinContent(b) - shapes[i]->GetBinContent(b), 2));
                    }
                    // and cumulate in the total for this toy as well
                    totByCh1[pair->second.channel]->Add(&*hist);
                    (sig[i] ? sigByCh1 : bkgByCh1)[pair->second.channel]->Add(&*hist);
                }
            }
            // now add up the deviations in this toy
            for (IH h = totByCh1.begin(), eh = totByCh1.end(); h != eh; ++h) {
                TH1 *target = totByCh2[h->first], *reference = totByCh[h->first];
                for (int b = 1, nb = target->GetNbinsX(); b <= nb; ++b) {
                    target->AddBinContent(b, std::pow(h->second->GetBinContent(b) - reference->GetBinContent(b), 2));
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
            TH1 *sum2 = totByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete totByCh1[h->first];
        }
        for (IH h = sigByCh.begin(), eh = sigByCh.end(); h != eh; ++h) {
            TH1 *sum2 = sigByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete sigByCh1[h->first];
        }
        for (IH h = bkgByCh.begin(), eh = bkgByCh.end(); h != eh; ++h) {
            TH1 *sum2 = bkgByCh2[h->first];
            for (int b = 1, nb = sum2->GetNbinsX(); b <= nb; ++b) {
                h->second->SetBinError(b, std::sqrt(sum2->GetBinContent(b)/ntoys));
            }
            delete sum2; delete bkgByCh1[h->first];
        }
        totByCh1.clear(); totByCh2.clear(); sigByCh1.clear(); sigByCh2.clear(); bkgByCh1.clear(); bkgByCh2.clear();
        // finally reset parameters
        params->assignValueOnly( sampler.centralValues() );
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
    }
}


//void MaxLikelihoodFit::setFitResultTrees(const RooArgSet *args, std::vector<double> *vals){
void MaxLikelihoodFit::setFitResultTrees(const RooArgSet *args, double * vals){
	
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

void MaxLikelihoodFit::setNormsFitResultTrees(const RooArgSet *args, double * vals){
	
         TIterator* iter(args->createIterator());
	 int count=0;
	 
         for (TObject *a = iter->Next(); a != 0; a = iter->Next()) { 
                 RooRealVar *rcv = dynamic_cast<RooRealVar *>(a);        
		 //std::string name = rcv->GetName();
		 vals[count]=rcv->getVal();
		 count++;
         }
	 delete iter;
	 return;
}

void MaxLikelihoodFit::createFitResultTrees(const RooStats::ModelConfig &mc, bool withSys){

	 // Initiate the arrays to store parameters

	 // create TTrees to store fit results:
	 t_fit_b_  = new TTree("tree_fit_b","tree_fit_b");
	 t_fit_sb_ = new TTree("tree_fit_sb","tree_fit_sb");

    	 t_fit_b_->Branch("fit_status",&fitStatus_,"fit_status/Int_t");
   	 t_fit_sb_->Branch("fit_status",&fitStatus_,"fit_status/Int_t");

	 t_fit_b_->Branch("mu",&mu_,"mu/Double_t");
	 t_fit_sb_->Branch("mu",&mu_,"mu/Double_t");

	 t_fit_b_->Branch("muLoErr",&muLoErr_,"muLoErr/Double_t");
	 t_fit_sb_->Branch("muLoErr",&muLoErr_,"muLoErr/Double_t");

	 t_fit_b_->Branch("muHiErr",&muHiErr_,"muHiErr/Double_t");
	 t_fit_sb_->Branch("muHiErr",&muHiErr_,"muHiErr/Double_t");

	 t_fit_b_->Branch("numbadnll",&numbadnll_,"numbadnll/Int_t");
	 t_fit_sb_->Branch("numbadnll",&numbadnll_,"numbadnll/Int_t");

	 t_fit_b_->Branch("nll_min",&nll_bonly_,"nll_min/Double_t");
	 t_fit_sb_->Branch("nll_min",&nll_sb_,"nll_min/Double_t");

	 t_fit_sb_->Branch("nll_nll0",&nll_nll0_,"nll_nll0/Double_t");

	 int count=0; 
         // fill the maps for the nuisances, and global observables
         RooArgSet *norms = new RooArgSet();
         getNormalizationsSimple(mc.GetPdf(), *mc.GetObservables(), *norms);
 
         processNormalizations_ = new double[norms->getSize()];

	 // If no systematic (-S 0), then don't make nuisance trees
	 if (withSys){
          const RooArgSet *cons = mc.GetGlobalObservables();
          const RooArgSet *nuis = mc.GetNuisanceParameters();
 	  globalObservables_ = new double[cons->getSize()];
	  nuisanceParameters_= new double[nuis->getSize()];

          TIterator* iter_c(cons->createIterator());
          for (TObject *a = iter_c->Next(); a != 0; a = iter_c->Next()) { 
                 RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rrv->GetName();
		 globalObservables_[count]=0;
		 t_fit_sb_->Branch(name.c_str(),&(globalObservables_[count]),Form("%s/Double_t",name.c_str()));
		 t_fit_b_->Branch(name.c_str(),&(globalObservables_[count]),Form("%s/Double_t",name.c_str()));
		 count++;
	  }         
	  count = 0;
          TIterator* iter_n(nuis->createIterator());
          for (TObject *a = iter_n->Next(); a != 0; a = iter_n->Next()) { 
                 RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rrv->GetName();
		 nuisanceParameters_[count] = 0;
		 t_fit_sb_->Branch(name.c_str(),&(nuisanceParameters_[count])),Form("%s/Double_t",name.c_str());
		 t_fit_b_->Branch(name.c_str(),&(nuisanceParameters_[count]),Form("%s/Double_t",name.c_str()));
		 count++;
          }

         }

	 count = 0;
         TIterator* iter_no(norms->createIterator());
         for (TObject *a = iter_no->Next(); a != 0; a = iter_no->Next()) { 
                 RooRealVar *rcv = dynamic_cast<RooRealVar *>(a);        
		 std::string name = rcv->GetName();
		 processNormalizations_[count] = 0;
		 t_fit_sb_->Branch(name.c_str(),&(processNormalizations_[count])),Form("%s/Double_t",name.c_str());
		 t_fit_b_->Branch(name.c_str(),&(processNormalizations_[count]),Form("%s/Double_t",name.c_str()));
		 count++;
         }
         delete norms;

	std::cout << "Created Branches" <<std::endl;
         return;	
}

MaxLikelihoodFit::ToySampler::ToySampler(RooAbsPdf *pdf, const RooArgSet *nuisances) :
    pdf_(pdf),
    data_(0)
{
    nuisances->snapshot(snapshot_);
}

MaxLikelihoodFit::ToySampler::~ToySampler() 
{
    delete data_;
}

void MaxLikelihoodFit::ToySampler::generate(int ntoys) {
    delete data_;
    data_ = pdf_->generate(snapshot_, ntoys);
}

const RooAbsCollection & MaxLikelihoodFit::ToySampler::get(int itoy) { 
    return *data_->get(itoy); 
}

const RooAbsCollection & MaxLikelihoodFit::ToySampler::centralValues() {
    return snapshot_;
}
