#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "HiggsAnalysis/CombinedLimit/interface/Significance.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logger.h"

#include <Math/MinimizerOptions.h>
#include <Math/IOptions.h>
#include <RooCategory.h>
#include <RooNumIntConfig.h>
#include <TStopwatch.h>
#include <RooStats/RooStatsUtils.h>

#include <iomanip>

boost::program_options::options_description CascadeMinimizer::options_("Cascade Minimizer options");
std::vector<CascadeMinimizer::Algo> CascadeMinimizer::fallbacks_;
bool CascadeMinimizer::preScan_;
double CascadeMinimizer::approxPreFitTolerance_ = 0;
int CascadeMinimizer::approxPreFitStrategy_ = 0;
int  CascadeMinimizer::preFit_ = 0;
bool CascadeMinimizer::poiOnlyFit_;
bool CascadeMinimizer::singleNuisFit_;
bool CascadeMinimizer::setZeroPoint_ = true;
bool CascadeMinimizer::oldFallback_ = false;
bool CascadeMinimizer::firstHesse_ = false;
bool CascadeMinimizer::lastHesse_ = false;
int  CascadeMinimizer::minuit2StorageLevel_ = 0;
bool CascadeMinimizer::runShortCombinations = true;
float CascadeMinimizer::nuisancePruningThreshold_ = 0;
double CascadeMinimizer::discreteMinTol_ = 0.001;
std::string CascadeMinimizer::defaultMinimizerType_="Minuit2"; // default to minuit2 (not always the default !?)
std::string CascadeMinimizer::defaultMinimizerAlgo_=ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
double CascadeMinimizer::defaultMinimizerTolerance_=1e-1;  
double CascadeMinimizer::defaultMinimizerPrecision_=-1.0;
int  CascadeMinimizer::strategy_=1; 

std::map<std::string,std::vector<std::string> > const CascadeMinimizer::minimizerAlgoMap_{
 {"Minuit"	 ,{"Migrad","Simplex","Combined","Scan"}}
,{"Minuit2" 	 ,{"Migrad","Simplex","Combined","Scan"}}
,{"GSLMultiMin"  ,{"ConjugateFR", "ConjugatePR", "BFGS", "BFGS2", "SteepestDescent"}}
};

CascadeMinimizer::CascadeMinimizer(RooAbsReal &nll, Mode mode, RooRealVar *poi) :
    nll_(nll),
    mode_(mode),
    //strategy_(0),
    poi_(poi),
    nuisances_(0),
    autoBounds_(false),
    poisForAutoBounds_(0),
    poisForAutoMax_(0)
{
    remakeMinimizer();
}

void CascadeMinimizer::remakeMinimizer() {
    cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
    if (simnll) simnll->setHideRooCategories(true);
    minimizer_.reset(); // avoid two copies in memory
    minimizer_.reset(new RooMinimizer(nll_));
    if (simnll) simnll->setHideRooCategories(false);
}

bool CascadeMinimizer::freezeDiscParams(const bool freeze)
{
    static bool freezeDisassParams = runtimedef::get(std::string("MINIMIZER_freezeDisassociatedParams"));
    static bool freezeDisassParams_verb = runtimedef::get(std::string("MINIMIZER_freezeDisassociatedParams_verbose"));
    if (freezeDisassParams) {
      if (freezeDisassParams_verb) {
          CascadeMinimizerGlobalConfigs::O().allRooMultiPdfs.Print();
          CascadeMinimizerGlobalConfigs::O().allRooMultiPdfParams.Print();
      }
      bool ret =  utils::freezeAllDisassociatedRooMultiPdfParameters((CascadeMinimizerGlobalConfigs::O().allRooMultiPdfs),(CascadeMinimizerGlobalConfigs::O().allRooMultiPdfParams),freeze);
      return ret;
    } else {
      return false;
    }
}

void CascadeMinimizer::setAutoBounds(const RooArgSet *pois) 
{
    poisForAutoBounds_ = pois;
    autoBounds_ = (poisForAutoBounds_ != 0 || poisForAutoMax_ != 0);
}

void CascadeMinimizer::setAutoMax(const RooArgSet *pois) 
{
    poisForAutoMax_ = pois;
    autoBounds_ = (poisForAutoBounds_ != 0 || poisForAutoMax_ != 0);
}


bool CascadeMinimizer::improve(int verbose, bool cascade, bool forceResetMinimizer) 
{
    cacheutils::CachingSimNLL *simnllbb = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
    if (simnllbb && runtimedef::get(std::string("MINIMIZER_analytic"))) {
      simnllbb->setAnalyticBarlowBeeston(true);
      forceResetMinimizer = true;
    }
    if (forceResetMinimizer || !minimizer_.get()) remakeMinimizer();
    minimizer_->setPrintLevel(verbose-1);
   
    strategy_ = ROOT::Math::MinimizerOptions::DefaultStrategy(); // re-configure 

    minimizer_->setStrategy(strategy_);
    std::string nominalType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
    std::string nominalAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
    float       nominalTol(ROOT::Math::MinimizerOptions::DefaultTolerance());
    minimizer_->setEps(nominalTol);
    if (approxPreFitTolerance_ > 0) {
      double tol = std::max(approxPreFitTolerance_, 10. * nominalTol);
      do {
        if (verbose > 1) std::cout << "Running pre-fit with " << nominalType << "," << nominalAlgo << " and tolerance " << tol << std::endl;
        Significance::MinimizerSentry minimizerConfig(nominalType+","+nominalAlgo, tol);
        minimizer_->setEps(tol);
        minimizer_->setStrategy(approxPreFitStrategy_);
        improveOnce(verbose-1, true);
        if (runtimedef::get("DBG_QUICKEXIT")) {
          exit(0);
        }
        minimizer_->setEps(nominalTol);
        minimizer_->setStrategy(strategy_);
      } while (autoBounds_ && !autoBoundsOk(verbose-1));
    }
    bool outcome;
    do {
      outcome = improveOnce(verbose-1);
      if (cascade && !outcome && !fallbacks_.empty()) {
        int         nominalStrat(strategy_);
        if (verbose > 0) {
		std::cerr << "Failed minimization with " << nominalType << "," << nominalAlgo << " and tolerance " << nominalTol << std::endl;
		Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Failed minimization with %s, %s and tolerance %g",__LINE__,nominalType.c_str(),nominalAlgo.c_str(),nominalTol)),Logger::kLogLevelDebug,__func__);
	}
        for (std::vector<Algo>::const_iterator it = fallbacks_.begin(), ed = fallbacks_.end(); it != ed; ++it) {
            Significance::MinimizerSentry minimizerConfig(it->type + "," + it->algo, it->tolerance != Algo::default_tolerance() ? it->tolerance : nominalTol); // set the global defaults
            int myStrategy = it->strategy; if (myStrategy == Algo::default_strategy()) myStrategy = nominalStrat;
            if (nominalType != ROOT::Math::MinimizerOptions::DefaultMinimizerType() ||
                nominalAlgo != ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo() ||
                nominalTol  != ROOT::Math::MinimizerOptions::DefaultTolerance()     ||
                myStrategy  != nominalStrat) {
                if (verbose > 0) { 
			std::cerr << "Will fallback to minimization using " << it->algo << ", strategy " << myStrategy << " and tolerance " << it->tolerance << std::endl;
			Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Will fallback to minimization using %s, strategy %d and tolerance %g",__LINE__,(it->algo).c_str(),myStrategy,it->tolerance)),Logger::kLogLevelDebug,__func__);
		}
                minimizer_->setEps(ROOT::Math::MinimizerOptions::DefaultTolerance());
                minimizer_->setStrategy(myStrategy); 
                outcome = improveOnce(verbose-2);
                if (outcome) break;
            }
        }
	
      }
    } while (autoBounds_ && !autoBoundsOk(verbose-1));

    if (simnllbb && runtimedef::get(std::string("MINIMIZER_analytic"))) {
      simnllbb->setAnalyticBarlowBeeston(false);
    }
    return outcome;
}

bool CascadeMinimizer::improveOnce(int verbose, bool noHesse) 
{
    static int optConst = runtimedef::get("MINIMIZER_optimizeConst");
    static int rooFitOffset = runtimedef::get("MINIMIZER_rooFitOffset");
    std::string myType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
    std::string myAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
    int myStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    bool outcome = false;
    double tol = ROOT::Math::MinimizerOptions::DefaultTolerance();
    static int maxcalls = runtimedef::get("MINIMIZER_MaxCalls");
    if (!minimizer_.get()) remakeMinimizer();
    if (maxcalls) {
        minimizer_->setMaxFunctionCalls(maxcalls);
        minimizer_->setMaxIterations(maxcalls);
    }
    if (oldFallback_){
        if (optConst) minimizer_->optimizeConst(std::max(0,optConst));
        if (rooFitOffset) minimizer_->setOffsetting(std::max(0,rooFitOffset));
        outcome = nllutils::robustMinimize(nll_, *minimizer_, verbose, setZeroPoint_);
    } else {
        if (verbose+2>0) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimisation configured with Type=%s, Algo=%s, strategy=%d, tolerance=%g",__LINE__,myType.c_str(),myAlgo.c_str(),myStrategy,tol)),Logger::kLogLevelInfo,__func__);
        cacheutils::CachingSimNLL *simnll = setZeroPoint_ ? dynamic_cast<cacheutils::CachingSimNLL *>(&nll_) : 0;
        if (simnll) simnll->setZeroPoint();
        if ((!simnll) && optConst) minimizer_->optimizeConst(std::max(0,optConst));
        if ((!simnll) && rooFitOffset) minimizer_->setOffsetting(std::max(0,rooFitOffset));
        if (firstHesse_ && !noHesse) {
            minimizer_->setPrintLevel(std::max(0,verbose-3)); 
            minimizer_->hesse();
            if (simnll) simnll->updateZeroPoint(); 
            minimizer_->setPrintLevel(verbose-1); 
        }
        int status = minimizer_->minimize(myType.c_str(), myAlgo.c_str());
        if (lastHesse_ && !noHesse) {
            if (simnll) simnll->updateZeroPoint(); 
            minimizer_->setPrintLevel(std::max(0,verbose-3)); 
            status = minimizer_->hesse();
            minimizer_->setPrintLevel(verbose-1); 
    	    if (verbose+2>0 ) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Hesse finished with status=%d",__LINE__,status)),Logger::kLogLevelDebug,__func__);
        }
        if (simnll) simnll->clearZeroPoint();
        outcome = (status == 0 || status == 1);
	if (status==1) std::cerr << "[WARNING] Minimisation finished with status 1 (covariance forced positive definite), this could indicate a problem with the minimim!" << std::endl;
    	if (verbose+2>0 ) {
		Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimisation finished with status=%d",__LINE__,status)),Logger::kLogLevelInfo,__func__);
		if (status==1) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- finished with status 1 (covariance forced positive definite), this could indicate a problem with the minimim.",__LINE__)),Logger::kLogLevelDebug,__func__);
	}
    }
    if (verbose+2>0 ){
     if  (outcome) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimization success! status=0",__LINE__)),Logger::kLogLevelInfo,__func__);
     else Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimization ended with latest status != 0 or 1",__LINE__)),Logger::kLogLevelDebug,__func__);
    }
    return outcome;
}


bool CascadeMinimizer::minos(const RooArgSet & params , int verbose ) {
   
   cacheutils::CachingSimNLL *simnllbb = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
   if (simnllbb && runtimedef::get(std::string("MINIMIZER_analytic"))) {
      // if one of the barlow-beeston params is in "params", we don't actually
      // want to freeze it here. Trick is to set all floating ones constant now,
      // then call setAnalyticBarlowBeeston, which will initiate bb only for the
      // floating ones, before unfreezing params again.
      RooArgSet toFreeze(params);
      RooStats::RemoveConstantParameters(&toFreeze);
      utils::setAllConstant(toFreeze, true);
      simnllbb->setAnalyticBarlowBeeston(true);
      utils::setAllConstant(toFreeze, false);
      remakeMinimizer();
   }
   if (!minimizer_.get()) remakeMinimizer();
   minimizer_->setPrintLevel(verbose-1); // for debugging
   std::string myType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
   std::string myAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());

   if (setZeroPoint_) {
      cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
      if (simnll) { 
         simnll->setZeroPoint();
      }
   }

   //TStopwatch tw;
   // need to re-run Migrad before running minos
   minimizer_->minimize(myType.c_str(), "Migrad");
   int iret = minimizer_->minos(params); 
   if (verbose>0 ) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minos finished with status=%d",__LINE__,iret)),Logger::kLogLevelDebug,__func__);

   //std::cout << "Run Minos in  "; tw.Print(); std::cout << std::endl;

   if (setZeroPoint_) {
      cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
      if (simnll) simnll->clearZeroPoint();
   }

   if (simnllbb && runtimedef::get(std::string("MINIMIZER_analytic"))) {
     simnllbb->setAnalyticBarlowBeeston(false);
   }

   return (iret != 1) ? true : false; 
}

bool CascadeMinimizer::hesse(int verbose ) {
   
   cacheutils::CachingSimNLL *simnllbb = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
   if (simnllbb && runtimedef::get(std::string("MINIMIZER_analytic"))) {
      // Have to reset and minimize again first to get all parameters in
      remakeMinimizer();
      float       nominalTol(ROOT::Math::MinimizerOptions::DefaultTolerance());
      minimizer_->setEps(nominalTol);
      minimizer_->setStrategy(strategy_);
      improveOnce(verbose - 1);
   }
   if (!minimizer_.get()) remakeMinimizer();
   minimizer_->setPrintLevel(verbose-1); // for debugging
   std::string myType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
   std::string myAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());

   if (setZeroPoint_) {
      cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
      if (simnll) { 
         simnll->setZeroPoint();
      }
   }

   int iret = minimizer_->hesse(); 

   if (setZeroPoint_) {
      cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
      if (simnll) simnll->clearZeroPoint();
   }

   return (iret != 1) ? true : false; 
}

bool CascadeMinimizer::iterativeMinimize(double &minimumNLL,int verbose, bool cascade){

   /* 
   If there are discrete parameters, first we cycle through them, 
   fixing all parameters which do not depend on them 

   step 1, the set which is passed contains all of the parameters which 
   are freely floating. We should cut them down to find which ones are
   */

   // Do A reasonable fit if something changed before 
   
   // First freeze all parameters that have nothing to do with the current active pdfs*
   freezeDiscParams(true);

   //std::cout << " Staring in iterativeMinimize and the minimum NLL so far is  " << minimumNLL << std::endl; 
   if ( fabs(minimumNLL - nll_.getVal()) > discreteMinTol_ ) { 
     improve(verbose,cascade);
     //std::cout << " Had to improve further since tolerance is not yet reached   " << nll_.getVal() << std::endl; 
   }

   // Next remove the POIs and constrained nuisances - this is to set up for the fast loop over the Index combinations
   RooArgSet nuisances = CascadeMinimizerGlobalConfigs::O().allFloatingParameters;
   nuisances.remove(CascadeMinimizerGlobalConfigs::O().allRooMultiPdfParams);

   RooArgSet poi = CascadeMinimizerGlobalConfigs::O().parametersOfInterest;
   RooArgSet frozen;

   if (nuisances.getSize() >0) frozen.add(nuisances);
   if (poi.getSize() >0) frozen.add(poi);
   
   RooStats::RemoveConstantParameters(&frozen);
   utils::setAllConstant(frozen,true);

   RooArgSet reallyCleanParameters;
   std::unique_ptr<RooArgSet> nllParams(nll_.getParameters((const RooArgSet*)0));
   nllParams->remove(CascadeMinimizerGlobalConfigs::O().pdfCategories);
   RooStats::RemoveConstantParameters(&*nllParams);
   (nllParams)->snapshot(reallyCleanParameters); 

   // Now cycle and fit
   bool ret=true;
   std::vector<std::vector<bool>> contIndex;
   
   // start from simplest scan, this is the full scan if runShortCombinations is off
   multipleMinimize(reallyCleanParameters,ret,minimumNLL,verbose,cascade,0,contIndex); 
 
   //if (simnll) simnll->clearZeroPoint();

   TStopwatch tw; tw.Start();
   utils::setAllConstant(frozen,false);
   
   // Run one last fully floating fit to maintain RooFitResult
   ret = improve(verbose, cascade); 
   minimumNLL = nll_.getVal();

   // unfreeze from *
   freezeDiscParams(false);

   tw.Stop(); if (verbose > 2) std::cout << "Done the full fit in " << tw.RealTime() << std::endl;

   return ret;
}

bool CascadeMinimizer::minimize(int verbose, bool cascade) 
{
    static int optConst = runtimedef::get("MINIMIZER_optimizeConst");
    static int rooFitOffset = runtimedef::get("MINIMIZER_rooFitOffset");
    if (runtimedef::get("CMIN_CENSURE")) {
        RooMsgService::instance().setStreamStatus(0,kFALSE);
        RooMsgService::instance().setStreamStatus(1,kFALSE);
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    }

    freezeDiscParams(true); // We should do anyway this since there can also be some indeces which are frozen 

    bool doMultipleMini = (CascadeMinimizerGlobalConfigs::O().pdfCategories.getSize()>0);
    if (runtimedef::get(std::string("MINIMIZER_skipDiscreteIterations"))) doMultipleMini=false;
    // if ( doMultipleMini ) preFit_ = 1;
    if (!minimizer_.get()) remakeMinimizer();
    minimizer_->setPrintLevel(verbose-2);  
    minimizer_->setStrategy(strategy_);
    
    RooArgSet nuisances = CascadeMinimizerGlobalConfigs::O().nuisanceParameters;

    if (preFit_ ) {
        RooArgSet frozen(nuisances);
        RooStats::RemoveConstantParameters(&frozen);
        utils::setAllConstant(frozen,true);
        freezeDiscParams(true);

        remakeMinimizer();
        minimizer_->setPrintLevel(verbose-2);
        minimizer_->setStrategy(preFit_-1);
        cacheutils::CachingSimNLL *simnll = setZeroPoint_ ? dynamic_cast<cacheutils::CachingSimNLL *>(&nll_) : 0;
        if (simnll) simnll->setZeroPoint();
        if (optConst) minimizer_->optimizeConst(std::max(0,optConst));
        if (rooFitOffset) minimizer_->setOffsetting(std::max(0,rooFitOffset));
        minimizer_->minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        if (simnll) simnll->clearZeroPoint();
        utils::setAllConstant(frozen,false);
        freezeDiscParams(false);
        remakeMinimizer();
    }
    
    bool ret = true;
    if (!doMultipleMini){
    	if (mode_ == Unconstrained && poiOnlyFit_) {
       	 trivialMinimize(nll_, *poi_, 200);
    	} 

      ret = improve(verbose, cascade);

    }else{
      // Do the discrete nuisance magic

      // clean parameters before minimization but dont include the pdf indeces of course!
      RooArgSet reallyCleanParameters;
      std::unique_ptr<RooArgSet> nllParams(nll_.getParameters((const RooArgSet*)0));
      nllParams->remove(CascadeMinimizerGlobalConfigs::O().pdfCategories);
      (nllParams)->snapshot(reallyCleanParameters); // should remove also the nuisance parameters from here!
      // Before each step, reset the parameters back to their prefit state!
      
      if (runShortCombinations) {
        // Initial fit under current index values
        improve(verbose, cascade);
        double backupApproxPreFitTolerance = approxPreFitTolerance_;
        approxPreFitTolerance_ = 0.;

        double minimumNLL  = nll_.getVal();
        double previousNLL = nll_.getVal();
        int maxIterations = 15; int iterationCounter=0;
        for (;iterationCounter<maxIterations;iterationCounter++){
          ret = iterativeMinimize(minimumNLL,verbose,cascade);
          if ( fabs(previousNLL-minimumNLL) < discreteMinTol_ ) break; // should be minimizer tolerance
          previousNLL = minimumNLL ;
        }
        approxPreFitTolerance_ = backupApproxPreFitTolerance;
      } else {

        double minimumNLL = 10+nll_.getVal();
        std::vector<std::vector<bool>> contIndex;
        multipleMinimize(reallyCleanParameters,ret,minimumNLL,verbose,cascade,0,contIndex);
   
        if (CascadeMinimizerGlobalConfigs::O().pdfCategories.getSize() > 1) {
           multipleMinimize(reallyCleanParameters,ret,minimumNLL,verbose,cascade,1,contIndex);
           multipleMinimize(reallyCleanParameters,ret,minimumNLL,verbose,cascade,2,contIndex);
        }

      }
    }

    // Check boundaries
    std::unique_ptr<RooArgSet> nllParams(nll_.getParameters((const RooArgSet*)0));
    RooStats::RemoveConstantParameters(&*nllParams);
    nllParams->remove(CascadeMinimizerGlobalConfigs::O().pdfCategories);
    nllParams->remove(CascadeMinimizerGlobalConfigs::O().parametersOfInterest);

    bool boundariesNotOk = utils::anyParameterAtBoundaries(*nllParams, verbose);
    if(boundariesNotOk && verbose > 0){
      fprintf(CloseCoutSentry::trueStdOutGlobal(),
        " [WARNING] After the fit some parameters are at their boundary.\n"
        " [WARNING] Are you sure your model is correct?\n");
      Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- After fit, some parameters are found at the boundary (within ~1sigma)",__LINE__)),Logger::kLogLevelInfo,__func__);
    }
    freezeDiscParams(false);
    return ret;
}

bool CascadeMinimizer::multipleMinimize(const RooArgSet &reallyCleanParameters, bool& ret, double& minimumNLL, int verbose, bool cascade,int mode, std::vector<std::vector<bool> >&contributingIndeces){
    static bool freezeDisassParams = runtimedef::get(std::string("MINIMIZER_freezeDisassociatedParams"));
    static bool hideConstants = freezeDisassParams && runtimedef::get(std::string("MINIMIZER_multiMin_hideConstants"));
    static bool maskConstraints = freezeDisassParams && runtimedef::get(std::string("MINIMIZER_multiMin_maskConstraints"));
    static int maskChannels = freezeDisassParams ? runtimedef::get(std::string("MINIMIZER_multiMin_maskChannels")) : 0;
    cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);

    //RooTrace::active(true);
    /* Different modes for minimization 
     Mode 0 -- Generate all combinations but only scan per-index 
     Mode 1 -- Generate only combinations which are orthogonal from the best fit after mode 0
	       Remove functions which cause increase in NLL > 10 (except best fit ones from previous mode)
     Mode 2 -- Full scan over the remaining combinations after mode 1
    */

    //std::cout << " At the start of the looping over the Indeces, minimum NLL is " << minimumNLL << std::endl; 
    // If the barlow-beeston minimisation is being used we can disable it temporarily,
    // saves time if we don't have to call enable/disable on the CMSHistErrorPropagators
    // repeatedly for no purpose
    int currentBarlowBeeston = runtimedef::get(std::string("MINIMIZER_analytic"));
    runtimedef::set("MINIMIZER_analytic", 0);
    
    double backupStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);

    bool newDiscreteMinimum = false;

    RooArgList pdfCategoryIndeces = CascadeMinimizerGlobalConfigs::O().pdfCategories; 
    int numIndeces = pdfCategoryIndeces.getSize();
    
    // create all combinations of indeces 
    std::vector<int> pdfSizes;

    RooCategory *fPdf;

    std::vector<int> bestIndeces(numIndeces,0);

    // Set to the current best indeces
    for (int id=0;id<numIndeces;id++) {
	int c =((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex();
	bestIndeces[id]=c;
    } 

    if (mode==0) { // mode 0 makes the indeces
      contributingIndeces.clear();
      for (int id=0;id<numIndeces;id++){
    	int npdf = ((RooCategory*)(pdfCategoryIndeces.at(id)))->numTypes();
	std::vector<bool> indexFlags(npdf,true);
	contributingIndeces.push_back(indexFlags);
      }
    }

    // now find the number of available pdfs
    for (int id=0;id<numIndeces;id++){
    	int npdf = ((RooCategory*)(pdfCategoryIndeces.at(id)))->numTypes();
	pdfSizes.push_back(npdf);
    }

    // keep hold of best fitted parameters! 
    std::auto_ptr<RooArgSet> params;
    params.reset(nll_.getParameters((const RooArgSet *)0) );
    params->remove(CascadeMinimizerGlobalConfigs::O().pdfCategories);

    //take a snapshot of those parameters
    RooArgSet snap;
    params->snapshot(snap);

    if (maskChannels && simnll) {
        simnll->setMaskNonDiscreteChannels(true);
    }
    if (hideConstants && simnll) {
        simnll->setHideConstants(true);
        if (maskConstraints) simnll->setMaskConstraints(true);
        minimizer_.reset(); // will be recreated when needed by whoever needs it
    }

    std::vector<std::vector<int> > myCombos;

    // Get All Permutations of pdfs
    if ( ( mode==0 ) /*&& runShortCombinations )*/ || mode ==1 ) myCombos = utils::generateOrthogonalCombinations(pdfSizes);
    else myCombos = utils::generateCombinations(pdfSizes);

    // Reorder to start from the "best indeces"
    //if (mode!=0) utils::reorderCombinations(myCombos,pdfSizes,bestIndeces);
    utils::reorderCombinations(myCombos,pdfSizes,bestIndeces);

    int numberOfCombinations = 1;
    if (mode==1 || mode==0) numberOfCombinations=myCombos.size();

    else {
    	for (int i=0;i<numIndeces;i++){
	 int nokpdfs=0;
      	 for (int j=0;j<pdfSizes[i];j++){
	   nokpdfs+=contributingIndeces[i][j];
         }
	 numberOfCombinations*=nokpdfs;
	}
    }

    std::vector<std::vector<int> >::iterator my_it = myCombos.begin();
    if (mode!=0) my_it++; // already did the best fit case
  
    TStopwatch tw; tw.Start();

    int fitCounter = 0;
    for (;my_it!=myCombos.end(); my_it++){

	     bool isValidCombo = true;
	
	     int pdfIndex=0, changedIndex = -1;
	     // Set the current indeces;
	     std::vector<int> cit = *my_it;
	     for (std::vector<int>::iterator it = cit.begin();
	         it!=cit.end(); it++){

		 isValidCombo *= (contributingIndeces)[pdfIndex][*it];
		 if (!isValidCombo ) /*&& runShortCombinations)*/ continue;

	     	 fPdf = (RooCategory*) pdfCategoryIndeces.at(pdfIndex);
                 if (fPdf->getIndex() != *it) changedIndex = pdfIndex;
		 fPdf->setIndex(*it);
		 pdfIndex++;
	     }
	
      if (!isValidCombo )/*&& runShortCombinations)*/ continue;
      
      if (verbose>2) {
	std::cout << "Setting indices := ";
	for (int id=0;id<numIndeces;id++) {
		std::cout << ((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex() << " ";
	}
        std::cout << std::endl;
      }

      if (fitCounter>0) params->assignValueOnly(reallyCleanParameters); // no need to reset from 0'th fit

      if (maskChannels == 2 && simnll) {
        for (int id=0;id<numIndeces;id++)  ((RooCategory*)(pdfCategoryIndeces.at(id)))->setConstant(id != changedIndex && changedIndex != -1);
        simnll->setMaskNonDiscreteChannels(true);
      }
      // Remove parameters which are not associated to the current PDF (only works if using --X-rtd MINIMIZER_freezeDisassociatedParams)
      freezeDiscParams(true);

      // FIXME can be made smarter than this
      if (mode_ == Unconstrained && poiOnlyFit_) {
        trivialMinimize(nll_, *poi_, 200);
      }

      ret =  improve(verbose, cascade, freezeDisassParams);

      if (maskChannels == 2 && simnll) {
        for (int id=0;id<numIndeces;id++)  ((RooCategory*)(pdfCategoryIndeces.at(id)))->setConstant(false);
        simnll->setMaskNonDiscreteChannels(false);
      }
      freezeDiscParams(false);


      fitCounter++;
      double thisNllValue = nll_.getVal();
      
      if ( thisNllValue < minimumNLL ){
		// Now we insert the correction ! 
                if (verbose>2) {
                    std::cout << " .... Found a better fit: new NLL = " << thisNllValue << " (improvement: " << (thisNllValue-minimumNLL) << std::endl;
                }
	        minimumNLL = thisNllValue;	
                //std::cout << " .... Found a better fit! hoorah! " << minimumNLL << std::endl; 
    		snap.assignValueOnly(*params);
		// set the best indeces again
		for (int id=0;id<numIndeces;id++) {
			if (bestIndeces[id] != ((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex() ) newDiscreteMinimum = true;
			bestIndeces[id]=((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex();	
		}
                if (verbose>2 && newDiscreteMinimum) {
                    std::cout << " .... Better fit corresponds to a new set of indices :=" ; 
                    for (int id=0;id<numIndeces;id++) { std::cout << " " << bestIndeces[id]; }
                    std::cout << std::endl;
                }
      }

      // FIXME this should be made configurable!
      double maxDeviation = 5;

      if (mode==1 )/*&& runShortCombinations)*/{

        if (thisNllValue > minimumNLL+maxDeviation){
		// Step 1, find out which index just changed 
		int modid   =0;
		int modcount=0;

      		for (int id=0;id<numIndeces;id++) {
			RooCategory* thisCat = (RooCategory*)(pdfCategoryIndeces.at(id));
			if (thisCat->getIndex()!=bestIndeces[id]){
				modid=id;
				modcount++;
			}
		}
		
		if (modcount==1){
		  // Step 2, remove its current index from the allowed indexes
		  RooCategory* thisCat = (RooCategory*)(pdfCategoryIndeces.at(modid));
		  int cIndex = thisCat->getIndex();
		  if (cIndex!=bestIndeces[modid]){ // don't remove the best pdf for this index!
			(contributingIndeces)[modid][cIndex]=false;
		  }
		}
        }
     }

    }

    // Assign best values ;
    for (int id=0;id<numIndeces;id++) {
	((RooCategory*)(pdfCategoryIndeces.at(id)))->setIndex(bestIndeces[id]);	
    } 
    params->assignValueOnly(snap);

    runtimedef::set("MINIMIZER_analytic", currentBarlowBeeston);
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(backupStrategy);

    tw.Stop(); if (verbose > 2) std::cout << "Done " << myCombos.size() << " combinations in " << tw.RealTime() << " s. New discrete minimum? " << newDiscreteMinimum << std::endl;

    if (maskChannels && simnll) {
        simnll->setMaskNonDiscreteChannels(false);
    }
    if (hideConstants && simnll) {
        simnll->setHideConstants(false);
        if (maskConstraints) simnll->setMaskConstraints(false);
        minimizer_.reset(); // will be recreated when needed by whoever needs it
    }


    return newDiscreteMinimum;
}

void CascadeMinimizer::initOptions() 
{
    options_.add_options()
        ("cminPoiOnlyFit",  "Do first a fit floating only the parameter of interest")
        ("cminPreScan",  "Do a scan before first minimization")
        ("cminPreFit", boost::program_options::value<int>(&preFit_)->default_value(preFit_), "if set to a value N > 0, it will perform a pre-fit with strategy (N-1) with frozen constrained nuisance parameters.")
        ("cminApproxPreFitTolerance", boost::program_options::value<double>(&approxPreFitTolerance_)->default_value(approxPreFitTolerance_), "If non-zero, do first a pre-fit with this tolerance (or 10 times the final tolerance, whichever is largest)")
        ("cminApproxPreFitStrategy", boost::program_options::value<int>(&approxPreFitStrategy_)->default_value(approxPreFitStrategy_), "Strategy to use in the pre-fit")
        ("cminSingleNuisFit", "Do first a minimization of each nuisance parameter individually")
        ("cminFallbackAlgo", boost::program_options::value<std::vector<std::string> >(), "Fallback algorithms if the default minimizer fails (can use multiple ones). Syntax is algo[,subalgo][,strategy][:tolerance]")
        ("cminSetZeroPoint", boost::program_options::value<bool>(&setZeroPoint_)->default_value(setZeroPoint_), "Change the reference point of the NLL to be zero during minimization")
        ("cminOldRobustMinimize", boost::program_options::value<bool>(&oldFallback_)->default_value(oldFallback_), "Use the old 'robustMinimize' logic in addition to the cascade (for debug only)")
        ("cminInitialHesse", boost::program_options::value<bool>(&firstHesse_)->default_value(firstHesse_), "Call Hesse before the minimization")
        ("cminFinalHesse", boost::program_options::value<bool>(&lastHesse_)->default_value(lastHesse_), "Call Hesse after the minimization")
	("cminDefaultMinimizerType",boost::program_options::value<std::string>(&defaultMinimizerType_)->default_value(defaultMinimizerType_), "Set the default minimizer Type")
	("cminDefaultMinimizerAlgo",boost::program_options::value<std::string>(&defaultMinimizerAlgo_)->default_value(defaultMinimizerAlgo_), "Set the default minimizer Algo")
	("cminDefaultMinimizerTolerance",boost::program_options::value<double>(&defaultMinimizerTolerance_)->default_value(defaultMinimizerTolerance_), "Set the default minimizer Tolerance")
	("cminDefaultMinimizerPrecision",boost::program_options::value<double>(&defaultMinimizerPrecision_)->default_value(defaultMinimizerPrecision_), "Set the default minimizer precision")
	("cminDefaultMinimizerStrategy",boost::program_options::value<int>(&strategy_)->default_value(strategy_), "Set the default minimizer (initial) strategy")
        ("cminRunAllDiscreteCombinations",  "Run all combinations for discrete nuisances")
        ("cminDiscreteMinTol", boost::program_options::value<double>(&discreteMinTol_)->default_value(discreteMinTol_), "tolerance on min NLL for discrete combination iterations")
        ("cminM2StorageLevel", boost::program_options::value<int>(&minuit2StorageLevel_)->default_value(minuit2StorageLevel_), "storage level for minuit2 (0 = don't store intermediate covariances, 1 = store them)")
        //("cminNuisancePruning", boost::program_options::value<float>(&nuisancePruningThreshold_)->default_value(nuisancePruningThreshold_), "if non-zero, discard constrained nuisances whose effect on the NLL when changing by 0.2*range is less than the absolute value of the threshold; if threshold is negative, repeat afterwards the fit with these floating")

        //("cminDefaultIntegratorEpsAbs", boost::program_options::value<double>(), "RooAbsReal::defaultIntegratorConfig()->setEpsAbs(x)")
        //("cminDefaultIntegratorEpsRel", boost::program_options::value<double>(), "RooAbsReal::defaultIntegratorConfig()->setEpsRel(x)")
        //("cminDefaultIntegrator1D", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->method1D().setLabel(x)")
        //("cminDefaultIntegrator1DOpen", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->method1DOpen().setLabel(x)")
        //("cminDefaultIntegrator2D", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->method2D().setLabel(x)")
        //("cminDefaultIntegrator2DOpen", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->method2DOpen().setLabel(x)")
        //("cminDefaultIntegratorND", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->methodND().setLabel(x)")
        //("cminDefaultIntegratorNDOpen", boost::program_options::value<std::string>(), "RooAbsReal::defaultIntegratorConfig()->methodNDOpen().setLabel(x)")
        ;
}

bool CascadeMinimizer::checkAlgoInType(std::string type, std::string algo){

    std::map<std::string,std::vector<std::string> >::const_iterator v = minimizerAlgoMap_.find(type);
    if (v != minimizerAlgoMap_.end()) {
      std::vector<std::string>::const_iterator a = (*v).second.end();
      if (std::find((*v).second.begin(), (*v).second.end(), algo) != a){
      	return true;
      }
      return false;
    }
    return false;

}

void CascadeMinimizer::applyOptions(const boost::program_options::variables_map &vm) 
{
    using namespace std;
    preScan_ = vm.count("cminPreScan");
    poiOnlyFit_ = vm.count("cminPoiOnlyFit");
    singleNuisFit_ = vm.count("cminSingleNuisFit");
    setZeroPoint_  = vm.count("cminSetZeroPoint");
    runShortCombinations = !(vm.count("cminRunAllDiscreteCombinations"));

    // check default minimizer type/algo if they are set and make sense
    if (vm.count("cminDefaultMinimizerAlgo")){
      if (! checkAlgoInType(defaultMinimizerType_,defaultMinimizerAlgo_)) {
	std::cerr << Form("The combination of minimizer type/algo %s/%s, is not recognized. Please set these with --cminDefaultMinimizerType and --cminDefaultMinimizerAlgo",defaultMinimizerType_.c_str(),defaultMinimizerAlgo_.c_str());
	//Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- The combination of minimizer type/algo %s/%s, is not recognized. Please set these with --cminDefaultMinimizerType and --cminDefaultMinimizerAlgo",__LINE__,defaultMinimizerType_.c_str(),defaultMinimizerAlgo_.c_str())),Logger::kLogLevelError,__func__);
	exit(0);
      }
    }

    if (vm.count("cminFallbackAlgo")) {
        vector<string> falls(vm["cminFallbackAlgo"].as<vector<string> >());
        for (vector<string>::const_iterator it = falls.begin(), ed = falls.end(); it != ed; ++it) {
            std::string algo = *it;
	    std::string type; 
            float tolerance = Algo::default_tolerance(); 
            int   strategy = Algo::default_strategy(); 
            string::size_type idx = std::min(algo.find(";"), algo.find(":"));
            if (idx != string::npos && idx < algo.length()) {
                 tolerance = atof(algo.substr(idx+1).c_str());
                 algo      = algo.substr(0,idx); // DON'T SWAP THESE TWO LINES
		 type	   = std::string(defaultMinimizerType_);
            }
            idx = algo.find(",");
            if (idx != string::npos && idx < algo.length()) {
                // if after the comma there's a number, then it's a strategy
                if ( '0' <= algo[idx+1] && algo[idx+1] <= '9' ) {
                    strategy = atoi(algo.substr(idx+1).c_str());
                    type     = algo.substr(0,idx); // DON'T SWAP THESE TWO LINES
		    std::map<std::string,std::vector<std::string> >::const_iterator ft = minimizerAlgoMap_.find(type);
		    if (ft!=minimizerAlgoMap_.end()){
		      algo     = (ft->second)[0];
		    } else algo = std::string(defaultMinimizerAlgo_);
		    
                } else {
                    // otherwise, it could be Name,subname,strategy
		    std::vector<std::string> configs;
		    boost::algorithm::split(configs,algo,boost::is_any_of(","));
		    if (configs.size()!=3) {
		    	std::cerr << "The fallback command from --cminFallbackAlgo " << *it << " is malformed. It should be formatted as Type[,Algo],strategy[:tolerance] " << std::endl;
			exit(0);
		    }
		    type = configs[0];
		    algo = configs[1];
		    strategy = atoi(configs[2].c_str());
                }
            }
      	    if (! checkAlgoInType(type,algo)) {
		std::cerr << Form("The fallback combination of minimizer type/algo %s/%s, is not recognized. Please check --cminFallbackAlgo again",type.c_str(),algo.c_str());
		//Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- The fallback combination of minimizer type/algo %s/%s, is not recognized. Please check --cminFallbackAlgo again",__LINE__,defaultMinimizerType_.c_str(),algo.c_str())),Logger::kLogLevelError,__func__);
		exit(0);
     	    }
            fallbacks_.push_back(Algo(type, algo, tolerance, strategy));
            std::cout << "Configured fallback algorithm " << 
	    		    ", type " << fallbacks_.back().type << 
	    		    ", algo " << fallbacks_.back().algo << 
                            ", strategy " << fallbacks_.back().strategy   << 
                            ", tolerance " << fallbacks_.back().tolerance << std::endl;
        }
    }
    
    ROOT::Math::IOptions & options = ROOT::Math::MinimizerOptions::Default("Minuit2");
    options.SetValue("StorageLevel", minuit2StorageLevel_);
    
    // Note that the options are not applied again when recreating a CascadeMinimizer so need to set the global attributes (should we make the modifiable options persistant too?)
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizerType_.c_str(),defaultMinimizerAlgo_.c_str());
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(defaultMinimizerTolerance_);
    if (defaultMinimizerPrecision_ > 0.) {
      ROOT::Math::MinimizerOptions::SetDefaultPrecision(defaultMinimizerPrecision_);
    }
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(strategy_);

    //if (vm.count("cminDefaultIntegratorEpsAbs")) RooAbsReal::defaultIntegratorConfig()->setEpsAbs(vm["cminDefaultIntegratorEpsAbs"].as<double>());
    //if (vm.count("cminDefaultIntegratorEpsRel")) RooAbsReal::defaultIntegratorConfig()->setEpsRel(vm["cminDefaultIntegratorEpsRel"].as<double>());
    //if (vm.count("cminDefaultIntegrator1D")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->method1D(), vm["cminDefaultIntegrator1D"].as<std::string>());
    //if (vm.count("cminDefaultIntegrator1DOpen")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->method1DOpen(), vm["cminDefaultIntegrator1DOpen"].as<std::string>());
    //if (vm.count("cminDefaultIntegrator2D")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->method2D(), vm["cminDefaultIntegrator2D"].as<std::string>());
    //if (vm.count("cminDefaultIntegrator2DOpen")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->method2DOpen(), vm["cminDefaultIntegrator2DOpen"].as<std::string>());
    //if (vm.count("cminDefaultIntegratorND")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->methodND(), vm["cminDefaultIntegratorND"].as<std::string>());
    //if (vm.count("cminDefaultIntegratorNDOpen")) setDefaultIntegrator(RooAbsReal::defaultIntegratorConfig()->methodNDOpen(), vm["cminDefaultIntegratorNDOpen"].as<std::string>());
}

//void CascadeMinimizer::setDefaultIntegrator(RooCategory &cat, const std::string & val) {
//    if (val == "list") {
//        std::cout << "States for " << cat.GetName() << std::endl;
//        int i0 = cat.getBin();
//        for (int i = 0, n = cat.numBins((const char *)0); i < n; ++i) {
//            cat.setBin(i); std::cout << " - " << cat.getLabel() <<  ( i == i0 ? " (current default)" : "") << std::endl;
//        }
//        std::cout << std::endl;
//        cat.setBin(i0);
//    } else {
//        cat.setLabel(val.c_str()); 
//    }
//}


void CascadeMinimizer::trivialMinimize(const RooAbsReal &nll, RooRealVar &r, int points) const {
    double rMin = r.getMin(), rMax = r.getMax(), rStep = (rMax-rMin)/(points-1);
    int iMin = -1; double minnll = 0;
    for (int i = 0; i < points; ++i) {
        double x = rMin + (i+0.5)*rStep;
        r.setVal(x);
        double y = nll.getVal();
        if (iMin == -1 || y < minnll) { minnll = y; iMin = i; }
    }
    r.setVal( rMin + (iMin+0.5)*rStep );
}

//void CascadeMinimizer::collectIrrelevantNuisances(RooAbsCollection &irrelevant) const {
//    if (nuisances_ == 0) return;
//    cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);
//    if (simnll == 0) return;
//    bool do_debug = runtimedef::get("CMIN_REVIEW_NUIS");
//    RooLinkedListIter iter = nuisances_->iterator();
//    const double here = simnll->myEvaluate(false);
//    const double thrsh = std::abs(nuisancePruningThreshold_);
//    for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next()) {
//        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
//        if (rrv->isConstant()) continue;
//        double v0 = rrv->getVal();
//        double rMin = rrv->getMin();
//        double rMax = rrv->getMax();
//        rrv->setVal( std::max(rMin, v0 - 0.2*(rMax-rMin)) );
//        double down = simnll->myEvaluate(false);
//        rrv->setVal( std::min(rMax, v0 + 0.2*(rMax-rMin)) );
//        double up   = simnll->myEvaluate(false);
//        rrv->setVal( v0 );
//        if (do_debug) printf("nuisance %s: %g [%g, %g]; deltaNLLs = %.5f, %.5f\n", rrv->GetName(), v0, rMin, rMax, here-up, here-down);
//        if (std::abs(here-up) < thrsh && std::abs(here-down)  < thrsh) irrelevant.add(*rrv);
//    }
//}

bool CascadeMinimizer::autoBoundsOk(int verbose) {
    bool ok = true;
    for (int bothBounds = 0; bothBounds <= 1; ++bothBounds) {
      const RooArgSet * pois = (bothBounds ? poisForAutoBounds_ : poisForAutoMax_);
      if (!pois) continue;
      RooFIter f = pois->fwdIterator();
      for (RooAbsArg *a = f.next(); a != 0; a = f.next()) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv && !rrv->isConstant() && rrv->hasMax() && rrv->hasMin()) {
            double val = rrv->getVal(), lo = rrv->getMin(), hi = rrv->getMax();
            if (bothBounds && val < (0.9*lo+0.1*hi)) {
                ok = false;
                rrv->setMin(val - (hi-val));
                if (verbose) std::cout << " POI " << rrv->GetName() << " is at " << val << ", within 10% from the low boundary " << lo << ". Will enlarge range to [ " << rrv->getMin() << " , " << hi << " ]" << std::endl;
            } else if (val > (0.9*hi+0.1*lo)) {
                ok = false;
                rrv->setMax(val + (val-lo));
                if (verbose) std::cout << " POI " << rrv->GetName() << " is at " << val << ", within 10% from the high boundary " << hi << ". Will enlarge range to [ " << lo << " , " << rrv->getMax() << " ]" << std::endl;
            }
        }
      }
    }
    if (!ok && verbose) { 
    	std::cout << "At least one of the POIs was close to the boundary, repeating the fit." << std::endl;
	Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- On checking with autoBounds on, At least one of the POIs was close to the boundary, repeating the fit.",__LINE__)),Logger::kLogLevelDebug,__func__);
    }
    return ok;
}
