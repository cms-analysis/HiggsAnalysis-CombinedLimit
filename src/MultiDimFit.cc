#include "HiggsAnalysis/CombinedLimit/interface/MultiDimFit.h"
#include <stdexcept>
#include <cmath>

#include "TMath.h"
#include "TFile.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRandom.h"
#include "RooAbsData.h"
#include "RooCategory.h"
#include "RooFitResult.h"
//#include "HiggsAnalysis/CombinedLimit/interface/RooMinimizerOpt.h"
#include "RooMinimizer.h"
#include <RooStats/ModelConfig.h>
#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"
#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/RobustHesse.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"

#include <Math/Minimizer.h>
#include <Math/MinimizerOptions.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFunc.h>

using namespace RooStats;

std::string MultiDimFit::name_ = "";
std::string MultiDimFit::massName_ = "";
std::string MultiDimFit::toyName_ = "";
std::string MultiDimFit::out_ = ".";
MultiDimFit::Algo MultiDimFit::algo_ = None;
MultiDimFit::GridType MultiDimFit::gridType_ = G1x1;
std::vector<std::string>  MultiDimFit::poi_;
std::vector<RooRealVar *> MultiDimFit::poiVars_;
std::vector<float>        MultiDimFit::poiVals_;
RooArgList                MultiDimFit::poiList_;
float                     MultiDimFit::deltaNLL_ = 0;
unsigned int MultiDimFit::points_ = 50;
unsigned int MultiDimFit::firstPoint_ = 0;
unsigned int MultiDimFit::lastPoint_  = std::numeric_limits<unsigned int>::max();
std::string MultiDimFit::gridPoints_ = "";
bool MultiDimFit::floatOtherPOIs_ = false;
unsigned int MultiDimFit::nOtherFloatingPoi_ = 0;
bool MultiDimFit::fastScan_ = false;
bool MultiDimFit::loadedSnapshot_ = false;
bool MultiDimFit::savingSnapshot_ = false;
bool MultiDimFit::startFromPreFit_ = false;
bool MultiDimFit::alignEdges_ = false;
bool MultiDimFit::hasMaxDeltaNLLForProf_ = false;
bool MultiDimFit::squareDistPoiStep_ = false;
bool MultiDimFit::skipInitialFit_ = false;
bool MultiDimFit::saveFitResult_ = false;
float MultiDimFit::maxDeltaNLLForProf_ = 200;
float MultiDimFit::autoRange_ = -1.0;
std::string MultiDimFit::fixedPointPOIs_ = "";
float MultiDimFit::centeredRange_ = -1.0;
bool        MultiDimFit::robustHesse_ = false;
std::string MultiDimFit::robustHesseLoad_ = "";
std::string MultiDimFit::robustHesseSave_ = "";


std::string MultiDimFit::saveSpecifiedFuncs_;
std::string MultiDimFit::saveSpecifiedIndex_;
std::string MultiDimFit::saveSpecifiedNuis_;
std::string MultiDimFit::setParametersForGrid_;
std::vector<std::string>  MultiDimFit::specifiedFuncNames_;
std::vector<RooAbsReal*> MultiDimFit::specifiedFunc_;
std::vector<float>        MultiDimFit::specifiedFuncVals_;
RooArgList                MultiDimFit::specifiedFuncList_;
std::vector<std::string>  MultiDimFit::specifiedCatNames_;
std::vector<RooCategory*> MultiDimFit::specifiedCat_;
std::vector<int>        MultiDimFit::specifiedCatVals_;
RooArgList                MultiDimFit::specifiedCatList_;
std::vector<std::string>  MultiDimFit::specifiedNuis_;
std::vector<RooRealVar *> MultiDimFit::specifiedVars_;
std::vector<float>        MultiDimFit::specifiedVals_;
RooArgList                MultiDimFit::specifiedList_;
bool MultiDimFit::saveInactivePOI_= false;

MultiDimFit::MultiDimFit() :
    FitterAlgoBase("MultiDimFit specific options")
{
    options_.add_options()
        ("algo",  boost::program_options::value<std::string>()->default_value("none"), "Algorithm to compute uncertainties")
        ("parameters,P",   boost::program_options::value<std::vector<std::string> >(&poi_), "Parameters to fit/scan (default = all parameters of interest)")
        ("floatOtherPOIs",   boost::program_options::value<bool>(&floatOtherPOIs_)->default_value(floatOtherPOIs_), "POIs other than the selected ones will be kept freely floating (1) or fixed (0, default)")
        ("squareDistPoiStep","POI step size based on distance from midpoint (max-min)/2 rather than linear")
        ("skipInitialFit","Skip initial fit (save time if snapshot is loaded from previous fit)")
        ("points",  boost::program_options::value<unsigned int>(&points_)->default_value(points_), "Points to use for grid or contour scans")
        ("gridPoints",  boost::program_options::value<std::string>(&gridPoints_)->default_value(gridPoints_), "Comma separated list of points per POI for multidimensional grid scans. When set, --points is ignored.")
        ("firstPoint",  boost::program_options::value<unsigned int>(&firstPoint_)->default_value(firstPoint_), "First point to use")
        ("lastPoint",  boost::program_options::value<unsigned int>(&lastPoint_)->default_value(lastPoint_), "Last point to use")
        ("autoRange", boost::program_options::value<float>(&autoRange_)->default_value(autoRange_), "Set to any X >= 0 to do the scan in the +/- X sigma range (where the sigma is from the initial fit, so it may be fairly approximate)")
	("fixedPointPOIs",   boost::program_options::value<std::string>(&fixedPointPOIs_)->default_value(""), "Parameter space point for --algo=fixed")
        ("centeredRange", boost::program_options::value<float>(&centeredRange_)->default_value(centeredRange_), "Set to any X >= 0 to do the scan in the +/- X range centered on the nominal value")
        ("fastScan", "Do a fast scan, evaluating the likelihood without profiling it.")
        ("maxDeltaNLLForProf",  boost::program_options::value<float>(&maxDeltaNLLForProf_)->default_value(maxDeltaNLLForProf_), "Last point to use")
	("saveSpecifiedNuis",   boost::program_options::value<std::string>(&saveSpecifiedNuis_)->default_value(""), "Save specified parameters (default = none)")
	("saveSpecifiedFunc",   boost::program_options::value<std::string>(&saveSpecifiedFuncs_)->default_value(""), "Save specified function values (default = none)")
	("saveSpecifiedIndex",   boost::program_options::value<std::string>(&saveSpecifiedIndex_)->default_value(""), "Save specified indexes/discretes (default = none)")
	("saveInactivePOI",   boost::program_options::value<bool>(&saveInactivePOI_)->default_value(saveInactivePOI_), "Save inactive POIs in output (1) or not (0, default)")
	("startFromPreFit",   boost::program_options::value<bool>(&startFromPreFit_)->default_value(startFromPreFit_), "Start each point of the likelihood scan from the pre-fit values")
    ("alignEdges",   boost::program_options::value<bool>(&alignEdges_)->default_value(alignEdges_), "Align the grid points such that the endpoints of the ranges are included")
    ("setParametersForGrid", boost::program_options::value<std::string>(&setParametersForGrid_)->default_value(""), "Set the values of relevant physics model parameters. Give a comma separated list of parameter value assignments. Example: CV=1.0,CF=1.0")
	("saveFitResult",  "Save RooFitResult to multidimfit.root")
    ("out", boost::program_options::value<std::string>(&out_)->default_value(out_), "Directory to put the diagnostics output file in")
    ("robustHesse",  boost::program_options::value<bool>(&robustHesse_)->default_value(robustHesse_),  "Use a more robust calculation of the hessian/covariance matrix")
    ("robustHesseLoad",  boost::program_options::value<std::string>(&robustHesseLoad_)->default_value(robustHesseLoad_),  "Load the pre-calculated Hessian")
    ("robustHesseSave",  boost::program_options::value<std::string>(&robustHesseSave_)->default_value(robustHesseSave_),  "Save the calculated Hessian")
      ;
}

void MultiDimFit::applyOptions(const boost::program_options::variables_map &vm) 
{
    applyOptionsBase(vm);
    std::string algo = vm["algo"].as<std::string>();
    if (algo == "none") {
        algo_ = None;
    } else if (algo == "singles") {
        algo_ = Singles;
    } else if (algo == "cross") {
        algo_ = Cross;
    } else if (algo == "grid" || algo == "grid3x3" ) {
        algo_ = Grid; gridType_ = G1x1;
        if (algo == "grid3x3") gridType_ = G3x3;
    } else if (algo == "fixed") {
        algo_ = FixedPoint;
    } else if (algo == "random") {
        algo_ = RandomPoints;
    } else if (algo == "contour2d") {
        algo_ = Contour2D;
    } else if (algo == "stitch2d") {
        algo_ = Stitch2D;
    } else if (algo == "impact") {
        algo_ = Impact;
        if (vm["floatOtherPOIs"].defaulted()) floatOtherPOIs_ = true;
        if (vm["saveInactivePOI"].defaulted()) saveInactivePOI_ = true;
    } else throw std::invalid_argument(std::string("Unknown algorithm: "+algo));
    fastScan_ = (vm.count("fastScan") > 0);
    squareDistPoiStep_ = (vm.count("squareDistPoiStep") > 0);
    skipInitialFit_ = (vm.count("skipInitialFit") > 0);
    hasMaxDeltaNLLForProf_ = !vm["maxDeltaNLLForProf"].defaulted();
    loadedSnapshot_ = !vm["snapshotName"].defaulted();
    savingSnapshot_ = vm.count("saveWorkspace");
    name_ = vm["name"].as<std::string>();
    massName_ = vm["massName"].as<std::string>();
    toyName_ = vm["toyName"].as<std::string>();
    saveFitResult_ = (vm.count("saveFitResult") > 0);
}

bool MultiDimFit::runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) { 
    // one-time initialization of POI variables, TTree branches, ...
    Combine::toggleGlobalFillTree(true);

    static int isInit = false;
    if (!isInit) { initOnce(w, mc_s); isInit = true; }

    // Get PDF
    RooAbsPdf &pdf = *mc_s->GetPdf();

    // Process POI not in list
    nOtherFloatingPoi_ = 0;
    deltaNLL_ = 0;
    int nConstPoi=0;
    RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
    std::string setConstPOI;
    for (RooAbsArg *a = (RooAbsArg*) iterP.Next(); a != 0; a = (RooAbsArg*) iterP.Next()) {
        if (poiList_.contains(*a)) continue;
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv == 0) { std::cerr << "MultiDimFit: Parameter of interest " << a->GetName() << " which is not a RooRealVar will be ignored" << std::endl; continue; }
        rrv->setConstant(!floatOtherPOIs_);
	if (!floatOtherPOIs_) {
	 setConstPOI+=std::string(rrv->GetName())+", ";
	 nConstPoi++;
	}
        if (floatOtherPOIs_) nOtherFloatingPoi_++;
    }
    if (nConstPoi>0) std::cout << "Following POIs have been set constant (use --floatOtherPOIs to let them float): " << setConstPOI << std::endl;
 
    // start with a best fit
    const RooCmdArg &constrainCmdArg = withSystematics  ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooCmdArg();
    std::auto_ptr<RooFitResult> res;
    if (verbose <= 3) RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
    bool doHesse = (algo_ == Singles || algo_ == Impact) || (saveFitResult_) ;
    if ( !skipInitialFit_){
        res.reset(doFit(pdf, data, (doHesse ? poiList_ : RooArgList()), constrainCmdArg, (saveFitResult_ && !robustHesse_), 1, true, false));
        if (!res.get()) {
            std::cout << "\n " <<std::endl;
            std::cout << "\n ---------------------------" <<std::endl;
            std::cout << "\n WARNING: MultiDimFit failed" <<std::endl;
            std::cout << "\n ---------------------------" <<std::endl;
            std::cout << "\n " <<std::endl;
        }
        if (algo_ == Impact && res.get()) {
            // Set the floating parameters back to the best-fit value
            // before we write an entry into the output TTree
            w->allVars().assignValueOnly(res.get()->floatParsFinal());
        }
    } else {
        std::cout << "MultiDimFit -- Skipping initial global fit" << std::endl;
        // must still create the NLL
        nll.reset(pdf.createNLL(data, constrainCmdArg, RooFit::Extended(pdf.canBeExtended()), RooFit::Offset(true)));
    }

    //if(w->var("r")) {w->var("r")->Print();}
    if ( loadedSnapshot_ || res.get() || keepFailures_) {
        for (int i = 0, n = poi_.size(); i < n; ++i) {
            if (res.get() && doHesse ){
	    	// (res.get())->Print("v");
                RooAbsArg *rfloat = (res.get())->floatParsFinal().find(poi_[i].c_str());
                if (!rfloat) {
                    rfloat = (*res).constPars().find(poi_[i].c_str());
                }
                RooRealVar *rf = dynamic_cast<RooRealVar*>(rfloat);
                poiVals_[i] = rf->getVal();//for Singles we store the RooFitResults values
            }
            else poiVals_[i] = poiVars_[i]->getVal();
        }
        //if (algo_ != None) {
	for(unsigned int j=0; j<specifiedNuis_.size(); j++){
		specifiedVals_[j]=specifiedVars_[j]->getVal();
	}
	for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
		specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
	}
	for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
		specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
	}
	Combine::commitPoint(/*expected=*/false, /*quantile=*/-1.); // Combine will not commit a point anymore at -1 so can do it here 
	//}
    }

    if (robustHesse_) {
        RobustHesse robustHesse(*nll, verbose - 1);
        robustHesse.ProtectArgSet(*mc_s->GetParametersOfInterest());
        if (robustHesseSave_ != "") {
          robustHesse.SaveHessianToFile(robustHesseSave_);
        }
        if (robustHesseLoad_ != "") {
          robustHesse.LoadHessianFromFile(robustHesseLoad_);
        }
        robustHesse.hesse();
        if (saveFitResult_) {
            res.reset(robustHesse.GetRooFitResult(res.get()));
        }
        robustHesse.WriteOutputFile("robustHesse"+name_+".root");
    }

   
    //set snapshot for best fit
    if (savingSnapshot_) w->saveSnapshot("MultiDimFit",utils::returnAllVars(w));
    
    if (autoRange_ > 0) {
        std::cout << "Adjusting range of POIs to +/- " << autoRange_ << " standard deviations" << std::endl;
        for (int i = 0, n = poi_.size(); i < n; ++i) {
            double val = poiVars_[i]->getVal(), err = poiVars_[i]->getError(), min0 = poiVars_[i]->getMin(), max0 = poiVars_[i]->getMax();
            double min1 = std::max(min0, val - autoRange_ * err);
            double max1 = std::min(max0, val + autoRange_ * err);
            std::cout << poi_[i] << ": " << val << " +/- " << err << " [ " << min0 << " , " << max0 << " ] ==> [ " << min1 << " , " << max1 << " ]" << std::endl;
            poiVars_[i]->setRange(min1, max1);
        }
    }
    if (centeredRange_ > 0) {
        std::cout << "Adjusting range of POIs to +/- " << centeredRange_ << std::endl;
        for (int i = 0, n = poi_.size(); i < n; ++i) {
            double val = poiVars_[i]->getVal(), min0 = poiVars_[i]->getMin(), max0 = poiVars_[i]->getMax();
            double min1 = std::max(min0, val - centeredRange_);
            double max1 = std::min(max0, val + centeredRange_);
            std::cout << poi_[i] << ": " << val << " [ " << min0 << " , " << max0 << " ] ==> [ " << min1 << " , " << max1 << " ]" << std::endl;
            poiVars_[i]->setRange(min1, max1);
        }
    }

    switch(algo_) {
        case None: 
	  {
            std::cout << "\n --- MultiDimFit ---" << std::endl;
            std::cout << "best fit parameter values: "  << std::endl;
            int len = poi_[0].length();
            for (int i = 0, n = poi_.size(); i < n; ++i) {
                len = std::max<int>(len, poi_[i].length());
            }
            for (int i = 0, n = poi_.size(); i < n; ++i) {
                printf("   %*s :  %+8.3f\n", len, poi_[i].c_str(), poiVals_[i]);
            }
	  }
          if(res.get() && saveFitResult_) saveResult(*res);
          break;
        case Singles: if (res.get()) { doSingles(*res); if (saveFitResult_) {saveResult(*res);} } break;
        case Cross: doBox(*nll, cl, "box", true); break;
        case Grid: doGrid(w,*nll); break;
        case RandomPoints: doRandomPoints(w,*nll); break;
        case FixedPoint: doFixedPoint(w,*nll); break;
        case Contour2D: doContour2D(w,*nll); break;
        case Stitch2D: doStitch2D(w,*nll); break;
        case Impact: if (res.get()) doImpact(*res, *nll); break;
    }
    
    Combine::toggleGlobalFillTree(false);
    return true;
}

void MultiDimFit::initOnce(RooWorkspace *w, RooStats::ModelConfig *mc_s) {

    // Tell combine not to Fill its tree, we'll do it here;

    RooArgSet mcPoi(*mc_s->GetParametersOfInterest());
    if (poi_.empty()) {
        RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
        for (RooAbsArg *a = (RooAbsArg*) iterP.Next(); a != 0; a = (RooAbsArg*) iterP.Next()) {
            poi_.push_back(a->GetName());
        }
    }
    for (std::vector<std::string>::const_iterator it = poi_.begin(), ed = poi_.end(); it != ed; ++it) {
        RooAbsArg *a = mcPoi.find(it->c_str());
	bool isPoi=true;
        if (a == 0) { 
		a = w->arg(it->c_str());  // look for the parameter elsewhere, but remember to clear its optimizeBounds attribute 
		isPoi = false;
	}
        if (a == 0) throw std::invalid_argument(std::string("Parameter of interest ")+*it+" not in model.");
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv == 0) throw std::invalid_argument(std::string("Parameter of interest ")+*it+" not a RooRealVar.");
	if (!isPoi) {
		if (rrv->getAttribute("optimizeBounds") ) {
            rrv->setAttribute("optimizeBounds",false);
            if (rrv->hasRange("optimizeBoundRange")) rrv->setRange(rrv->getMin("optimizeBoundRange"), rrv->getMax("optimizeBoundRange"));
        }
	}
        poiVars_.push_back(rrv);
        poiVals_.push_back(rrv->getVal());
        poiList_.add(*rrv);
    }

    if(saveSpecifiedFuncs_!=""){
	    char tmp[10240] ;
	    strlcpy(tmp,saveSpecifiedFuncs_.c_str(),10240) ;
	    char* token = strtok(tmp,",") ;
	    while(token) {
		    RooAbsArg *a = w->arg(token);
		    if (a == 0) throw std::invalid_argument(std::string("function ")+token+" not in model.");
		    RooAbsReal *rrv = dynamic_cast<RooAbsReal*>(a);
		    if (rrv == 0) throw std::invalid_argument(std::string("function ")+token+" not a RooAbsReal.");
		    specifiedFuncNames_.push_back(token);
		    specifiedFunc_.push_back(rrv);
		    specifiedFuncVals_.push_back(rrv->getVal());
		    specifiedFuncList_.add(*rrv);
		    token = strtok(0,",") ; 
	    }
    }
    if(saveSpecifiedIndex_!=""){
	    char tmp[10240] ;
	    strlcpy(tmp,saveSpecifiedIndex_.c_str(),10240) ;
	    char* token = strtok(tmp,",") ;
	    while(token) {
		    RooCategory *rrv = w->cat(token);
		    if (rrv == 0) throw std::invalid_argument(std::string("function ")+token+" not a RooCategory.");
		    specifiedCatNames_.push_back(token);
		    specifiedCat_.push_back(rrv);
		    specifiedCatVals_.push_back(rrv->getIndex());
		    specifiedCatList_.add(*rrv);
		    token = strtok(0,",") ; 
	    }
    }

    if(saveSpecifiedNuis_!="" && withSystematics){
	    RooArgSet mcNuis(*mc_s->GetNuisanceParameters());
	    if(saveSpecifiedNuis_=="all"){
		    specifiedNuis_.clear();
		    RooLinkedListIter iterN = mc_s->GetNuisanceParameters()->iterator();
		    for (RooAbsArg *a = (RooAbsArg*) iterN.Next(); a != 0; a = (RooAbsArg*) iterN.Next()) {
			    if (poiList_.contains(*a)) continue;
			    specifiedNuis_.push_back(a->GetName());
		    }
	    }else{
		    char tmp[10240] ;
		    strlcpy(tmp,saveSpecifiedNuis_.c_str(),10240) ;
		    char* token = strtok(tmp,",") ;
		    while(token) {
			    const RooArgSet* group = mc_s->GetWS()->set((std::string("group_") + token).data());
			    if (group){
				    RooLinkedListIter iterN = group->iterator();
				    for (RooAbsArg *a = (RooAbsArg*) iterN.Next(); a != 0; a = (RooAbsArg*) iterN.Next()) {
					    specifiedNuis_.push_back(a->GetName());
				    }
			    }else if (!poiList_.find(token)){
				    specifiedNuis_.push_back(token);
			    }
			    token = strtok(0,",") ; 
		    }
	    }
	    for (std::vector<std::string>::const_iterator it = specifiedNuis_.begin(), ed = specifiedNuis_.end(); it != ed; ++it) {
		    RooAbsArg *a = mcNuis.find(it->c_str());
		    if (a == 0) throw std::invalid_argument(std::string("Nuisance Parameter ")+*it+" not in model.");
		    if (poiList_.contains(*a)) continue;
		    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
		    if (rrv == 0) throw std::invalid_argument(std::string("Nuisance Parameter ")+*it+" not a RooRealVar.");
		    specifiedVars_.push_back(rrv);
		    specifiedVals_.push_back(rrv->getVal());
		    specifiedList_.add(*rrv);
	    }
    }
    if(saveInactivePOI_){
	    RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
	    for (RooAbsArg *a = (RooAbsArg*) iterP.Next(); a != 0; a = (RooAbsArg*) iterP.Next()) {
		    if (poiList_.contains(*a)) continue;
		    if (specifiedList_.contains(*a)) continue;
		    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
		    specifiedNuis_.push_back(a->GetName());
		    specifiedVars_.push_back(rrv);
		    specifiedVals_.push_back(rrv->getVal());
		    specifiedList_.add(*rrv);
	    }
    }

    // then add the branches to the tree (at the end, so there are no resizes)
    for (int i = 0, n = poi_.size(); i < n; ++i) {
        Combine::addBranch(poi_[i].c_str(), &poiVals_[i], (poi_[i]+"/F").c_str()); 
    }
    for (int i = 0, n = specifiedNuis_.size(); i < n; ++i) {
	Combine::addBranch(specifiedNuis_[i].c_str(), &specifiedVals_[i], (specifiedNuis_[i]+"/F").c_str()); 
    }
    for (int i = 0, n = specifiedFuncNames_.size(); i < n; ++i) {
	Combine::addBranch(specifiedFuncNames_[i].c_str(), &specifiedFuncVals_[i], (specifiedFuncNames_[i]+"/F").c_str()); 
    }
    for (int i = 0, n = specifiedCatNames_.size(); i < n; ++i) {
	Combine::addBranch(specifiedCatNames_[i].c_str(), &specifiedCatVals_[i], (specifiedCatNames_[i]+"/I").c_str()); 
    }
    Combine::addBranch("deltaNLL", &deltaNLL_, "deltaNLL/F");
}

void MultiDimFit::doSingles(RooFitResult &res)
{
    std::cout << "\n --- MultiDimFit ---" << std::endl;
    std::cout << "best fit parameter values and profile-likelihood uncertainties: "  << std::endl;
    int len = poi_[0].length();
    for (int i = 0, n = poi_.size(); i < n; ++i) {
        len = std::max<int>(len, poi_[i].length());
    }
    for (int i = 0, n = poi_.size(); i < n; ++i) {
	RooAbsArg *rfloat = res.floatParsFinal().find(poi_[i].c_str());
	if (!rfloat) {
		rfloat = res.constPars().find(poi_[i].c_str());
	}
        RooRealVar *rf = dynamic_cast<RooRealVar*>(rfloat);
        double bestFitVal = rf->getVal();

        double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal : rf->getAsymErrorHi());
        double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal : rf->getAsymErrorLo());
        double maxError = std::max<double>(std::max<double>(hiErr, loErr), rf->getError());
        if (fabs(hiErr) < 0.001*maxError){ 
		std::cout << " Warning - No valid high-error found, will report difference to maximum of range for : " << rf->GetName() << std::endl;
		hiErr = -bestFitVal + rf->getMax();
	}
        if (fabs(loErr) < 0.001*maxError) {
		std::cout << " Warning - No valid low-error found, will report difference to minimum of range for : " << rf->GetName() << std::endl;
		loErr = +bestFitVal - rf->getMin();
	}
        
	poiVals_[i] = bestFitVal - loErr; Combine::commitPoint(true, /*quantile=*/-0.32);
        poiVals_[i] = bestFitVal + hiErr; Combine::commitPoint(true, /*quantile=*/0.32);

        double hiErr95 = +(do95_ && rf->hasRange("err95") ? rf->getMax("err95") - bestFitVal : 0);
        double loErr95 = -(do95_ && rf->hasRange("err95") ? rf->getMin("err95") - bestFitVal : 0);
        double maxError95 = std::max<double>(std::max<double>(hiErr95, loErr95), 2*rf->getError());

        if (do95_ && rf->hasRange("err95")) {
            if (fabs(hiErr95) < 0.001*maxError95){ 
		std::cout << " Warning - No valid high-error (for 95%) found, will report difference to maximum of range for : " << rf->GetName() << std::endl;
		hiErr95 = -bestFitVal + rf->getMax();
	    }
            if (fabs(loErr95) < 0.001*maxError95) {
		std::cout << " Warning - No valid low-error (for 95%) found, will report difference to minimum of range for : " << rf->GetName() << std::endl;
		loErr95 = +bestFitVal - rf->getMin();
	    }
	    poiVals_[i] = bestFitVal - loErr95; Combine::commitPoint(true, /*quantile=*/-0.05);
            poiVals_[i] = bestFitVal + hiErr95; Combine::commitPoint(true, /*quantile=*/0.05);
            //poiVals_[i] = rf->getMax("err95"); Combine::commitPoint(true, /*quantile=*/-0.05);
            //poiVals_[i] = rf->getMin("err95"); Combine::commitPoint(true, /*quantile=*/0.05);
            poiVals_[i] = bestFitVal;
            printf("   %*s :  %+8.3f   %+6.3f/%+6.3f (68%%)    %+6.3f/%+6.3f (95%%) \n", len, poi_[i].c_str(), 
                    poiVals_[i], -loErr, hiErr, -loErr95, hiErr95);
        } else {
            poiVals_[i] = bestFitVal;
            printf("   %*s :  %+8.3f   %+6.3f/%+6.3f (68%%)\n", len, poi_[i].c_str(), 
                    poiVals_[i], -loErr, hiErr);
        }
    }
}

void MultiDimFit::doImpact(RooFitResult &res, RooAbsReal &nll) {
  std::cout << "\n --- MultiDimFit ---" << std::endl;
  std::cout << "Parameter impacts: " << std::endl;

  // Save the initial parameters here to reset between NPs
  std::auto_ptr<RooArgSet> params(nll.getParameters((const RooArgSet *)0));
  RooArgSet init_snap;
  params->snapshot(init_snap);

  // Save the best-fit values of the saved parameters
  // we want to measure the impacts on
  std::vector<float> specifiedVals = specifiedVals_;
  std::vector<float> impactLo = specifiedVals_;
  std::vector<float> impactHi = specifiedVals_;

  int len = 9;
  for (int i = 0, n = poi_.size(); i < n; ++i) {
    len = std::max<int>(len, poi_[i].length());
  }
  printf("  %-*s :   %-21s", len, "Parameter", "Best-fit");
  for (int i = 0, n = specifiedNuis_.size(); i < n; ++i) {
    printf("  %-13s", specifiedNuis_[i].c_str());
  }
  printf("\n");


  for (int i = 0, n = poi_.size(); i < n; ++i) {
    RooAbsArg *rfloat = res.floatParsFinal().find(poi_[i].c_str());
    if (!rfloat) {
      rfloat = res.constPars().find(poi_[i].c_str());
    }
    RooRealVar *rf = dynamic_cast<RooRealVar *>(rfloat);
    double bestFitVal = rf->getVal();

    double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal
                                           : rf->getAsymErrorHi());
    double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal
                                           : rf->getAsymErrorLo());
      printf("  %-*s : %+8.3f  %+6.3f/%+6.3f", len, poi_[i].c_str(),
                    bestFitVal, -loErr, hiErr);
    // Reset all parameters to initial state
    *params = init_snap;
    // Then set this NP constant
    poiVars_[i]->setConstant(true);
    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    //minim.setStrategy(minimizerStrategy_);
    // Another snapshot to reset between high and low fits
    RooArgSet snap;
    params->snapshot(snap);
    std::vector<double> doVals = {bestFitVal - loErr, bestFitVal + hiErr};
    for (unsigned x = 0; x < doVals.size(); ++x) {
      *params = snap;
      poiVals_[i] = doVals[x];
      poiVars_[i]->setVal(doVals[x]);
      bool ok = minim.minimize(verbose - 1);
      if (ok) {
        for (unsigned int j = 0; j < poiVars_.size(); j++) {
          poiVals_[j] = poiVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        Combine::commitPoint(true, /*quantile=*/0.32);
      }
      for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
        if (x == 0) {
          impactLo[j] = specifiedVars_[j]->getVal() - specifiedVals[j];
        } else if (x == 1) {
          impactHi[j] = specifiedVars_[j]->getVal() - specifiedVals[j];
        }
      }
    }
    for (unsigned j = 0; j < specifiedVals.size(); ++j) {
        printf("  %+6.3f/%+6.3f", impactLo[j], impactHi[j]);
    }
    printf("\n");
  }
}


void MultiDimFit::doGrid(RooWorkspace *w, RooAbsReal &nll) 
{
    unsigned int n = poi_.size();
    //if (poi_.size() > 2) throw std::logic_error("Don't know how to do a grid with more than 2 POIs.");
    double nll0 = nll.getVal();

    if (setParametersForGrid_ != "") {
       RooArgSet allParams(w->allVars());
       allParams.add(w->allCats());
       utils::setModelParameters( setParametersForGrid_, allParams);
    }

    if (startFromPreFit_) w->loadSnapshot("clean");

    std::vector<double> p0(n), pmin(n), pmax(n);
    for (unsigned int i = 0; i < n; ++i) {
        p0[i] = poiVars_[i]->getVal();
        pmin[i] = poiVars_[i]->getMin();
        pmax[i] = poiVars_[i]->getMax();
        poiVars_[i]->setConstant(true);
    std::cout<<" POI: "<<poiVars_[i]->GetName()<<"= "<<p0[i]<<" -> ["<<pmin[i]<<","<<pmax[i]<<"]"<<std::endl;
    }


    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    if (!autoBoundsPOIs_.empty()) minim.setAutoBounds(&autoBoundsPOISet_); 
    if (!autoMaxPOIs_.empty()) minim.setAutoMax(&autoMaxPOISet_); 
    //minim.setStrategy(minimizerStrategy_);
    std::auto_ptr<RooArgSet> params(nll.getParameters((const RooArgSet *)0));
    RooArgSet snap; params->snapshot(snap);
    //snap.Print("V");

    // check if gridPoints are defined
    std::vector<unsigned int> pointsPerPoi;
    if (!gridPoints_.empty()) {
        splitGridPoints(gridPoints_, pointsPerPoi);
        if (pointsPerPoi.size() != n) {
            throw std::logic_error("Number of passed gridPoints "
                + std::to_string(pointsPerPoi.size()) + " does not match number of POIs "
                + std::to_string(n));
        }
        std::cout << "Parsed number of points per POI: ";
        for (unsigned int i = 0; i < n; i++) {
            std::cout << poi_[i] << " -> " << pointsPerPoi[i];
            if (i < n - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }

    if (n == 1) {
        unsigned int points = pointsPerPoi.size() == 0 ? points_ : pointsPerPoi[0];

        double xspacing = (pmax[0]-pmin[0]) / points;
        double xspacingOffset = 0.5;
        if (alignEdges_) {
          xspacing = (pmax[0]-pmin[0]) / (points - 1);
          if (points == 1) xspacing = 0;
          xspacingOffset = 0.0;
        }
        // can do a more intellegent spacing of points
        double xbestpoint = (p0[0] - pmin[0]) / xspacing;
        if (lastPoint_ == std::numeric_limits<unsigned int>::max()) {
          lastPoint_ = points - 1;
        }
        for (unsigned int i = 0; i < points; ++i) {
          if (i < firstPoint_) continue;
          if (i > lastPoint_)  break;
          double x = pmin[0] + (i + xspacingOffset) * xspacing;
          // If we're aligning with the edges and this is the last point,
          // set x to pmax[0] exactly
          if (alignEdges_ && i == (points - 1)) {
            x = pmax[0];
          }
          if (xbestpoint > lastPoint_) {
            int ireverse = lastPoint_ - i + firstPoint_;
            x = pmin[0] + (ireverse + xspacingOffset) * xspacing;
          }

          if (squareDistPoiStep_) {
            // distance between steps goes as ~square of distance from middle or range (could this be changed to from best fit value?)
            double phalf = (pmax[0] - pmin[0]) / 2;
            if (x < (pmin[0] + phalf)) {
              x = pmin[0] + TMath::Sqrt((x - pmin[0]) / phalf) * phalf;
            } else {
              x = pmax[0] - TMath::Sqrt((pmax[0] - x) / phalf) * phalf;
            }
          }

            //if (verbose > 1) std::cout << "Point " << i << "/" << points << " " << poiVars_[0]->GetName() << " = " << x << std::endl;
             std::cout << "Point " << i << "/" << points << " " << poiVars_[0]->GetName() << " = " << x << std::endl;
            *params = snap;
            poiVals_[0] = x;
            poiVars_[0]->setVal(x);
            // now we minimize
            nll.clearEvalErrorLog();
            deltaNLL_ = nll.getVal() - nll0;
            if (nll.numEvalErrors() > 0) {
                deltaNLL_ = 9990;
		for(unsigned int j=0; j<specifiedNuis_.size(); j++){
			specifiedVals_[j]=specifiedVars_[j]->getVal();
		}
		for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
			specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
		}
                Combine::commitPoint(true, /*quantile=*/0);
                continue;
            }
            bool ok = fastScan_ || (hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_) || utils::countFloating(*params)==0 ? 
                        true : 
                        minim.minimize(verbose-1);
            if (ok) {
                deltaNLL_ = nll.getVal() - nll0;
                double qN = 2*(deltaNLL_);
                double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
		for(unsigned int j=0; j<specifiedNuis_.size(); j++){
			specifiedVals_[j]=specifiedVars_[j]->getVal();
		}
		for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
			specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
		}
		for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
			specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
		}
                Combine::commitPoint(true, /*quantile=*/prob);
            }
        }
    } else if (n == 2) {
        RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
        CloseCoutSentry sentry(verbose < 2);

        // get number of points per axis
        unsigned int nX, nY;
        if (pointsPerPoi.size() == 0) {
            // same number of points per axis ("old" behavior)
            unsigned int sqrn = ceil(sqrt(double(points_)));
            nX = nY = sqrn;
        } else {
            // number of points different per axis
            nX = pointsPerPoi[0];
            nY = pointsPerPoi[1];
        }
        unsigned int nTotal = nX * nY;

        // determine grid variables
        double deltaX, deltaY, spacingOffsetX, spacingOffsetY;
        if (nX == 1) {
            deltaX = 0;
            spacingOffsetX = 0;
        } else if (alignEdges_) {
            deltaX = (pmax[0] - pmin[0]) / (nX - 1);
            spacingOffsetX = 0;
        } else {
            deltaX = (pmax[0] - pmin[0]) / nX;
            spacingOffsetX = 0.5;
        }
        if (nY == 1) {
            deltaY = 0;
            spacingOffsetY = 0;
        } else if (alignEdges_) {
            deltaY = (pmax[1] - pmin[1]) / (nY - 1);
            spacingOffsetY = 0;
        } else {
            deltaY = (pmax[1] - pmin[1]) / nY;
            spacingOffsetY = 0.5;
        }
        unsigned int ipoint = 0, nprint = ceil(0.005 * nTotal);

        // loop through the grid
        for (unsigned int i = 0; i < nX; ++i) {
            for (unsigned int j = 0; j < nY; ++j, ++ipoint) {
                if (ipoint < firstPoint_) continue;
                if (ipoint > lastPoint_)  break;
                *params = snap;
                double x =  pmin[0] + (i + spacingOffsetX) * deltaX;
                double y =  pmin[1] + (j + spacingOffsetY) * deltaY;
                if (verbose && (ipoint % nprint == 0)) {
                         fprintf(sentry.trueStdOut(), "Point %d/%d, (i,j) = (%d,%d), %s = %f, %s = %f\n",
                                        ipoint,nTotal, i,j, poiVars_[0]->GetName(), x, poiVars_[1]->GetName(), y);
                }
                poiVals_[0] = x;
                poiVals_[1] = y;
                poiVars_[0]->setVal(x);
                poiVars_[1]->setVal(y);
                nll.clearEvalErrorLog(); nll.getVal();
                if (nll.numEvalErrors() > 0) { 
			for(unsigned int j=0; j<specifiedNuis_.size(); j++){
				specifiedVals_[j]=specifiedVars_[j]->getVal();
			}
			for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
				specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
			}
			for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
				specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
			}
                    deltaNLL_ = 9999; Combine::commitPoint(true, /*quantile=*/0); 
                    if (gridType_ == G3x3) {
                        for (int i2 = -1; i2 <= +1; ++i2) {
                            for (int j2 = -1; j2 <= +1; ++j2) {
                                if (i2 == 0 && j2 == 0) continue;
                                poiVals_[0] = x + 0.33333333*i2*deltaX;
                                poiVals_[1] = y + 0.33333333*j2*deltaY;
				for(unsigned int j=0; j<specifiedNuis_.size(); j++){
					specifiedVals_[j]=specifiedVars_[j]->getVal();
				}
				for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
					specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
				}
				for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
					specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
				}
                                deltaNLL_ = 9999; Combine::commitPoint(true, /*quantile=*/0); 
                            }
                        }
                    }
                    continue;
                }
                // now we minimize
                bool skipme = hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_;
                bool ok = fastScan_ || skipme ? true :  minim.minimize(verbose-1);
                if (ok) {
                    deltaNLL_ = nll.getVal() - nll0;
                    double qN = 2*(deltaNLL_);
                    double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
		    for(unsigned int j=0; j<specifiedNuis_.size(); j++){
			    specifiedVals_[j]=specifiedVars_[j]->getVal();
		    }
		    for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
			    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
		    }
		    for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
			    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
		    }
                    Combine::commitPoint(true, /*quantile=*/prob);
                }
                if (gridType_ == G3x3) {
                    bool forceProfile = !fastScan_ && std::min(fabs(deltaNLL_ - 1.15), fabs(deltaNLL_ - 2.995)) < 0.5;
                    utils::CheapValueSnapshot center(*params);
                    double x0 = x, y0 = y;
                    for (int i2 = -1; i2 <= +1; ++i2) {
                        for (int j2 = -1; j2 <= +1; ++j2) {
                            if (i2 == 0 && j2 == 0) continue;
                            center.writeTo(*params);
                            x = x0 + 0.33333333*i2*deltaX;
                            y = y0 + 0.33333333*j2*deltaY;
                            poiVals_[0] = x; poiVars_[0]->setVal(x);
                            poiVals_[1] = y; poiVars_[1]->setVal(y);
                            nll.clearEvalErrorLog(); nll.getVal();
                            if (nll.numEvalErrors() > 0) { 
				    for(unsigned int j=0; j<specifiedNuis_.size(); j++){
					    specifiedVals_[j]=specifiedVars_[j]->getVal();
				    }
				    for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
					    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
				    }
				    for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
					    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
				    }
                                deltaNLL_ = 9999; Combine::commitPoint(true, /*quantile=*/0); 
                                continue;
                            }
                            deltaNLL_ = nll.getVal() - nll0;
                            if (forceProfile || (!fastScan_ && std::min(fabs(deltaNLL_ - 1.15), fabs(deltaNLL_ - 2.995)) < 0.5)) {
                                minim.minimize(verbose-1);
                                deltaNLL_ = nll.getVal() - nll0;
                            }
                            double qN = 2*(deltaNLL_);
                            double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
			    for(unsigned int j=0; j<specifiedNuis_.size(); j++){
				    specifiedVals_[j]=specifiedVars_[j]->getVal();
			    }
			    for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
				    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
			    }
			    for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
				    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
			    }
                            Combine::commitPoint(true, /*quantile=*/prob);
                        }
                    }
                }
            }
        }

    } else { // Use utils routine if n > 2
        RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
        CloseCoutSentry sentry(verbose < 2);

        // get number of points per axis
        std::vector<int> axis_points;
        if (pointsPerPoi.size() == 0) {
            // same number of points per axis ("old" behavior)
            unsigned int rootn = ceil(TMath::Power(double(points_),double(1./n)));
            axis_points.resize(n, (int)rootn);
        } else {
            for (auto p : pointsPerPoi) {
                axis_points.push_back(p);
            }
        }
        unsigned int nTotal = 1;
        for (auto p : axis_points) nTotal *= p;
        unsigned int ipoint = 0, nprint = ceil(0.005*nTotal);

        // Create permutations
        std::vector<std::vector<int> > permutations = utils::generateCombinations(axis_points);

        // Step through points
        std::vector<std::vector<int> >::iterator perm_it = permutations.begin();
        int npermutations = permutations.size();
        for (;perm_it!=permutations.end(); perm_it++) {
            if (ipoint < firstPoint_) {
                ipoint++;
                continue;
            }
            if (ipoint > lastPoint_) break;
            *params = snap;

            if (verbose && (ipoint % nprint == 0)) {
                fprintf(sentry.trueStdOut(), "Point %d/%d, ", ipoint,npermutations);
            }
            for (unsigned int poi_i=0;poi_i<n;poi_i++) {
                int ip = (*perm_it)[poi_i];
                double deltaXi = (pmax[poi_i]-pmin[poi_i])/axis_points[poi_i];
                double spacingOffset = 0.5;
                if (alignEdges_) {
                    deltaXi = (pmax[poi_i] - pmin[poi_i]) / (axis_points[poi_i] - 1);
                    if (axis_points[poi_i] == 1) {
                        deltaXi = 0.;
                    }
                    spacingOffset = 0.0;
                }
                double xi = pmin[poi_i] + deltaXi * (ip + spacingOffset);
                poiVals_[poi_i] = xi; poiVars_[poi_i]->setVal(xi);
                if (verbose && (ipoint % nprint == 0)) {
                    fprintf(sentry.trueStdOut(), " %s = %f ", poiVars_[poi_i]->GetName(), xi);
                }
            }
            if (verbose && (ipoint % nprint == 0)) fprintf(sentry.trueStdOut(), "\n");

            nll.clearEvalErrorLog(); nll.getVal();
            if (nll.numEvalErrors() > 0) {
                for (unsigned int j=0; j<specifiedNuis_.size(); j++) {
                    specifiedVals_[j]=specifiedVars_[j]->getVal();
                }
                for (unsigned int j=0; j<specifiedFuncNames_.size(); j++) {
                    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
                }
                for (unsigned int j=0; j<specifiedCatNames_.size(); j++) {
                    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
                }
                deltaNLL_ = 9999; Combine::commitPoint(true, /*quantile=*/0);
                ipoint++;
                continue;
            }

            // now we minimize
            bool skipme = hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_;
            bool ok = fastScan_ || skipme ? true :  minim.minimize(verbose-1);
            if (ok) {
                deltaNLL_ = nll.getVal() - nll0;
                double qN = 2*(deltaNLL_);
                double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
                for (unsigned int j=0; j<specifiedNuis_.size(); j++) {
                    specifiedVals_[j]=specifiedVars_[j]->getVal();
                }
                for (unsigned int j=0; j<specifiedFuncNames_.size(); j++) {
                    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
                }
                for (unsigned int j=0; j<specifiedCatNames_.size(); j++) {
                    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
                }
                Combine::commitPoint(true, /*quantile=*/prob);
            }
            ipoint++;
        }
    }
}

void MultiDimFit::doRandomPoints(RooWorkspace *w, RooAbsReal &nll) 
{
    double nll0 = nll.getVal();
    if (startFromPreFit_) w->loadSnapshot("clean");
    for (unsigned int i = 0, n = poi_.size(); i < n; ++i) {
        poiVars_[i]->setConstant(true);
    }

    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    if (!autoBoundsPOIs_.empty()) minim.setAutoBounds(&autoBoundsPOISet_); 
    if (!autoMaxPOIs_.empty()) minim.setAutoMax(&autoMaxPOISet_); 
    //minim.setStrategy(minimizerStrategy_);
    unsigned int n = poi_.size();
    for (unsigned int j = 0; j < points_; ++j) {
        for (unsigned int i = 0; i < n; ++i) {
            poiVars_[i]->randomize();
            poiVals_[i] = poiVars_[i]->getVal(); 
        }
        // now we minimize
        {   
            CloseCoutSentry sentry(verbose < 3);    
            bool ok = minim.minimize(verbose-1);
            if (ok) {
                double qN = 2*(nll.getVal() - nll0);
                double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
		for(unsigned int j=0; j<specifiedNuis_.size(); j++){
			specifiedVals_[j]=specifiedVars_[j]->getVal();
		}
		for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
			specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
		}
		for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
			specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
		}
                Combine::commitPoint(true, /*quantile=*/prob);
            }
        } 
    }
}
void MultiDimFit::doFixedPoint(RooWorkspace *w, RooAbsReal &nll) 
{
    double nll0 = nll.getVal();
    if (startFromPreFit_) w->loadSnapshot("clean");
    for (unsigned int i = 0, n = poi_.size(); i < n; ++i) {
        poiVars_[i]->setConstant(true);
    }

    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    if (!autoBoundsPOIs_.empty()) minim.setAutoBounds(&autoBoundsPOISet_); 
    if (!autoMaxPOIs_.empty()) minim.setAutoMax(&autoMaxPOISet_); 
    //minim.setStrategy(minimizerStrategy_);
    unsigned int n = poi_.size();

    //for (unsigned int i = 0; i < n; ++i) {
    //        std::cout<<" Before setting fixed point "<<poiVars_[i]->GetName()<<"= "<<poiVals_[i]<<std::endl;
    //}
    if (fixedPointPOIs_ != "") {
	    utils::setModelParameters( fixedPointPOIs_, w->allVars());
    } else if (setPhysicsModelParameterExpression_ != "") {
            std::cout << " --fixedPointPOIs option not used, so will use the argument of --setParameters instead" << std::endl;
	    utils::setModelParameters( setPhysicsModelParameterExpression_, w->allVars());
    }   

    for (unsigned int i = 0; i < n; ++i) {
	    poiVars_[i] -> setConstant(true);
	    poiVals_[i] = poiVars_[i]->getVal(); 
            std::cout<<" Evaluating fixed point with "<<poiVars_[i]->GetName()<<"= "<<poiVals_[i]<<std::endl;
    }
    // now we minimize
    {   
	    CloseCoutSentry sentry(verbose < 3);    
	    bool ok = minim.minimize(verbose-1);
	    if (ok) {
		    nll0Value_ = nll0;
		    nllValue_ = nll.getVal();
		    deltaNLL_ = nll.getVal() - nll0;
		    double qN = 2*(nll.getVal() - nll0);
		    double prob = ROOT::Math::chisquared_cdf_c(qN, n+nOtherFloatingPoi_);
		    for(unsigned int j=0; j<specifiedNuis_.size(); j++){
			    specifiedVals_[j]=specifiedVars_[j]->getVal();
		    }
		    for(unsigned int j=0; j<specifiedFuncNames_.size(); j++){
			    specifiedFuncVals_[j]=specifiedFunc_[j]->getVal();
		    }
		    for(unsigned int j=0; j<specifiedCatNames_.size(); j++){
			    specifiedCatVals_[j]=specifiedCat_[j]->getIndex();
		    }
		    Combine::commitPoint(true, /*quantile=*/prob);
    //for (unsigned int i = 0; i < n; ++i) {
    //        std::cout<<" after the fit "<<poiVars_[i]->GetName()<<"= "<<poiVars_[i]->getVal()<<std::endl;
    //}
        }
    } 
}

void MultiDimFit::doContour2D(RooWorkspace *, RooAbsReal &nll) 
{
    if (poi_.size() != 2) throw std::logic_error("Contour2D works only in 2 dimensions");
    RooRealVar *xv = poiVars_[0]; double x0 = poiVals_[0]; float &x = poiVals_[0];
    RooRealVar *yv = poiVars_[1]; double y0 = poiVals_[1]; float &y = poiVals_[1];

    double threshold = nll.getVal() + 0.5*ROOT::Math::chisquared_quantile_c(1-cl,2+nOtherFloatingPoi_);
    if (verbose>0) std::cout << "Best fit point is for " << xv->GetName() << ", "  << yv->GetName() << " =  " << x0 << ", " << y0 << std::endl;

    // make a box
    doBox(nll, cl, "box");
    double xMin = xv->getMin("box"), xMax = xv->getMax("box");
    double yMin = yv->getMin("box"), yMax = yv->getMax("box");

    verbose--; // reduce verbosity to avoid messages from findCrossing
    // ===== Get relative min/max of x for several fixed y values =====
    yv->setConstant(true);
    for (unsigned int j = 0; j <= points_; ++j) {
        if (j < firstPoint_) continue;
        if (j > lastPoint_)  break;
        // take points uniformly spaced in polar angle in the case of a perfect circle
        double yc = 0.5*(yMax + yMin), yr = 0.5*(yMax - yMin);
        yv->setVal( yc + yr * std::cos(j*M_PI/double(points_)) );
        // ===== Get the best fit x (could also do without profiling??) =====
        xv->setConstant(false);  xv->setVal(x0);
        CascadeMinimizer minimXI(nll, CascadeMinimizer::Unconstrained, xv);
        if (!autoBoundsPOIs_.empty()) minimXI.setAutoBounds(&autoBoundsPOISet_); 
        if (!autoMaxPOIs_.empty()) minimXI.setAutoMax(&autoMaxPOISet_); 
        //minimXI.setStrategy(minimizerStrategy_);
        {
            CloseCoutSentry sentry(verbose < 3);    
            minimXI.minimize(verbose-1);
        }
        double xc = xv->getVal(); xv->setConstant(true);
        if (verbose>-1) std::cout << "Best fit " << xv->GetName() << " for  " << yv->GetName() << " = " << yv->getVal() << " is at " << xc << std::endl;
        // ===== Then get the range =====
        CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
        if (!autoBoundsPOIs_.empty()) minim.setAutoBounds(&autoBoundsPOISet_); 
        if (!autoMaxPOIs_.empty()) minim.setAutoMax(&autoMaxPOISet_); 
        double xup = findCrossing(minim, nll, *xv, threshold, xc, xMax);
        if (!std::isnan(xup)) { 
            x = xup; y = yv->getVal(); Combine::commitPoint(true, /*quantile=*/1-cl);
            if (verbose>-1) std::cout << "Minimum of " << xv->GetName() << " at " << cl << " CL for " << yv->GetName() << " = " << y << " is " << x << std::endl;
        }
        
        double xdn = findCrossing(minim, nll, *xv, threshold, xc, xMin);
        if (!std::isnan(xdn)) { 
            x = xdn; y = yv->getVal(); Combine::commitPoint(true, /*quantile=*/1-cl);
            if (verbose>-1) std::cout << "Maximum of " << xv->GetName() << " at " << cl << " CL for " << yv->GetName() << " = " << y << " is " << x << std::endl;
        }
    }

    verbose++; // restore verbosity
}

void MultiDimFit::doStitch2D(RooWorkspace *, RooAbsReal &nll)
{
    if (poi_.size() != 2) throw std::logic_error("Contour2D works only in 2 dimensions");
    //RooRealVar *xv = poiVars_[0]; double x0 = poiVals_[0]; float &x = poiVals_[0];
    //RooRealVar *yv = poiVars_[1]; double y0 = poiVals_[1]; float &y = poiVals_[1];

    //double threshold = nll.getVal() + 0.5*ROOT::Math::chisquared_quantile_c(1-cl,2+nOtherFloatingPoi_);
    //if (verbose>0) std::cout << "Best fit point is for " << xv->GetName() << ", "  << yv->GetName() << " =  " << x0 << ", " << y0 << std::endl;

    // make a box
    //doBox(nll, cl, "box");
    //double xMin = xv->getMin("box"), xMax = xv->getMax("box");
    //double yMin = yv->getMin("box"), yMax = yv->getMax("box");

//    verbose--; // reduce verbosity to avoid messages from findCrossing
//    // ===== Get relative min/max of x for several fixed y values =====
//    yv->setConstant(true);
//    for (unsigned int j = 0; j <= points_; ++j) {
//        if (j < firstPoint_) continue;
//        if (j > lastPoint_)  break;
//        // take points uniformly spaced in polar angle in the case of a perfect circle
//        double yc = 0.5*(yMax + yMin), yr = 0.5*(yMax - yMin);
//        yv->setVal( yc + yr * std::cos(j*M_PI/double(points_)) );
//        // ===== Get the best fit x (could also do without profiling??) =====
//        xv->setConstant(false);  xv->setVal(x0);
//        CascadeMinimizer minimXI(nll, CascadeMinimizer::Unconstrained, xv);
//        minimXI.setStrategy(minimizerStrategy_);
//        {
//            CloseCoutSentry sentry(verbose < 3);
//            minimXI.minimize(verbose-1);
//        }
//        double xc = xv->getVal(); xv->setConstant(true);
//        if (verbose>-1) std::cout << "Best fit " << xv->GetName() << " for  " << yv->GetName() << " = " << yv->getVal() << " is at " << xc << std::endl;
//        // ===== Then get the range =====
//        CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
//        double xup = findCrossing(minim, nll, *xv, threshold, xc, xMax);
//        if (!std::isnan(xup)) {
//            x = xup; y = yv->getVal(); Combine::commitPoint(true, /*quantile=*/1-cl);
//            if (verbose>-1) std::cout << "Minimum of " << xv->GetName() << " at " << cl << " CL for " << yv->GetName() << " = " << y << " is " << x << std::endl;
//        }
//
//        double xdn = findCrossing(minim, nll, *xv, threshold, xc, xMin);
//        if (!std::isnan(xdn)) {
//            x = xdn; y = yv->getVal(); Combine::commitPoint(true, /*quantile=*/1-cl);
//            if (verbose>-1) std::cout << "Maximum of " << xv->GetName() << " at " << cl << " CL for " << yv->GetName() << " = " << y << " is " << x << std::endl;
//        }
//    }
//
//    verbose++; // restore verbosity
}


void MultiDimFit::doBox(RooAbsReal &nll, double cl, const char *name, bool commitPoints)  {
    unsigned int n = poi_.size();
    double nll0 = nll.getVal(), threshold = nll0 + 0.5*ROOT::Math::chisquared_quantile_c(1-cl,n+nOtherFloatingPoi_);

    std::vector<double> p0(n);
    for (unsigned int i = 0; i < n; ++i) {
        p0[i] = poiVars_[i]->getVal();
        poiVars_[i]->setConstant(false);
    }

    verbose--; // reduce verbosity due to findCrossing
    for (unsigned int i = 0; i < n; ++i) {
        RooRealVar *xv = poiVars_[i];
        xv->setConstant(true);
        CascadeMinimizer minimX(nll, CascadeMinimizer::Constrained);
        if (!autoBoundsPOIs_.empty()) minimX.setAutoBounds(&autoBoundsPOISet_); 
        if (!autoMaxPOIs_.empty()) minimX.setAutoMax(&autoMaxPOISet_); 
        //minimX.setStrategy(minimizerStrategy_);

        for (unsigned int j = 0; j < n; ++j) poiVars_[j]->setVal(p0[j]);
        double xMin = findCrossing(minimX, nll, *xv, threshold, p0[i], xv->getMin()); 
        if (!std::isnan(xMin)) { 
            if (verbose > -1) std::cout << "Minimum of " << xv->GetName() << " at " << cl << " CL for all others floating is " << xMin << std::endl;
            for (unsigned int j = 0; j < n; ++j) poiVals_[j] = poiVars_[j]->getVal();
            if (commitPoints) Combine::commitPoint(true, /*quantile=*/1-cl);
        } else {
            xMin = xv->getMin();
            for (unsigned int j = 0; j < n; ++j) poiVals_[j] = poiVars_[j]->getVal();
            double prob = ROOT::Math::chisquared_cdf_c(2*(nll.getVal() - nll0), n+nOtherFloatingPoi_);
            if (commitPoints) Combine::commitPoint(true, /*quantile=*/prob);
            if (verbose > -1) std::cout << "Minimum of " << xv->GetName() << " at " << cl << " CL for all others floating is " << xMin << " (on the boundary, p-val " << prob << ")" << std::endl;
        }
        
        for (unsigned int j = 0; j < n; ++j) poiVars_[j]->setVal(p0[j]);
        double xMax = findCrossing(minimX, nll, *xv, threshold, p0[i], xv->getMax()); 
        if (!std::isnan(xMax)) { 
            if (verbose > -1) std::cout << "Maximum of " << xv->GetName() << " at " << cl << " CL for all others floating is " << xMax << std::endl;
            for (unsigned int j = 0; j < n; ++j) poiVals_[j] = poiVars_[j]->getVal();
            if (commitPoints) Combine::commitPoint(true, /*quantile=*/1-cl);
        } else {
            xMax = xv->getMax();
            double prob = ROOT::Math::chisquared_cdf_c(2*(nll.getVal() - nll0), n+nOtherFloatingPoi_);
            for (unsigned int j = 0; j < n; ++j) poiVals_[j] = poiVars_[j]->getVal();
            if (commitPoints) Combine::commitPoint(true, /*quantile=*/prob);
            if (verbose > -1) std::cout << "Maximum of " << xv->GetName() << " at " << cl << " CL for all others floating is " << xMax << " (on the boundary, p-val " << prob << ")" << std::endl;
        }

        xv->setRange(name, xMin, xMax);
        xv->setConstant(false);
    }
    verbose++; // restore verbosity 
}

void MultiDimFit::saveResult(RooFitResult &res) {
    if (verbose>2) res.Print();
    if (out_ == "none") return;
    const bool longName = runtimedef::get(std::string("longName"));
    std::string mdname(out_+"/multidimfit"+name_);
    if (longName)
        mdname += "."+massName_+toyName_+"root";
    else
        mdname += ".root";
    fitOut.reset(TFile::Open(mdname.c_str(), "RECREATE"));
    fitOut->WriteTObject(&res,"fit_mdf");
    fitOut->cd();
    fitOut.release()->Close();
}

void MultiDimFit::splitGridPoints(const std::string& s, std::vector<unsigned int>& points) const {
    // split by comma
    std::vector<std::string> strPoints;
    boost::split(strPoints, s, boost::is_any_of(","));

    // convert to int and add
    for (const auto strPoint : strPoints) {
        points.push_back(std::stoul(strPoint));
    }
}
