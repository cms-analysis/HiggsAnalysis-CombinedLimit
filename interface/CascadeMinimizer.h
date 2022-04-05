#ifndef HiggsAnalysis_CombinedLimit_CascadeMinimizer_h
#define HiggsAnalysis_CombinedLimit_CascadeMinimizer_h

class RooAbsReal;
class RooArgSet;
class RooRealVar;
#include <RooArgSet.h>
#include <RooListProxy.h>
#include <RooSetProxy.h>
//#include "HiggsAnalysis/CombinedLimit/interface/RooMinimizerOpt.h"
#include "RooMinimizer.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>


class CascadeMinimizer {
    public:
        enum Mode { Constrained, Unconstrained };
        CascadeMinimizer(RooAbsReal &nll, Mode mode, RooRealVar *poi=0) ;
        // do a new minimization, assuming the initial state is random
        bool minimize(int verbose=0, bool cascade=true);
        // run minos
        bool minos(const RooArgSet &, int verbose = 0 );
        // run hesse
        bool hesse(int verbose = 0 );
        // do a new minimization, assuming a plausible initial state
        bool improve(int verbose=0, bool cascade=true, bool forceResetMinimizer=false);
        // declare nuisance parameters for pre-fit
        void setNuisanceParameters(const RooArgSet *nuis) { nuisances_ = nuis; }
        RooMinimizer & minimizer() { return *minimizer_; }
        RooFitResult *save() { return minimizer().save(); }
        void  setStrategy(int strategy) { strategy_ = strategy; }
        void  setErrorLevel(float errorLevel) { minimizer_->setErrorLevel(errorLevel); }
        static void  initOptions() ;
        static void  applyOptions(const boost::program_options::variables_map &vm) ;
        static const boost::program_options::options_description & options() { return options_; }
        void trivialMinimize(const RooAbsReal &nll, RooRealVar &r, int points=100) const ;
        //void collectIrrelevantNuisances(RooAbsCollection &irrelevant) const ;
	bool freezeDiscParams(const bool);
        void setAutoBounds(const RooArgSet *pois) ;
        void setAutoMax(const RooArgSet *pois) ; 
	double tolerance() {return defaultMinimizerTolerance_;};
	std::string algo() {return defaultMinimizerAlgo_;};
    private:
        RooAbsReal & nll_;
        std::auto_ptr<RooMinimizer> minimizer_;
        Mode         mode_;
        static int          strategy_;
        RooRealVar * poi_; 
        const RooArgSet *nuisances_;
        /// automatically enlarge bounds for POIs if they're within 10% from the boundary
        bool autoBounds_;
        const RooArgSet *poisForAutoBounds_, *poisForAutoMax_;

        bool improveOnce(int verbose, bool noHesse=false);
        bool autoBoundsOk(int verbose) ;

	bool multipleMinimize(const RooArgSet &,bool &,double &,int,bool,int
		,std::vector<std::vector<bool> > & );
       
        bool iterativeMinimize(double &,int,bool); 

        void remakeMinimizer() ;

        /// options configured from command line
        static boost::program_options::options_description options_;
        /// compact information about an algorithm
        struct Algo { 
            Algo() : type(), algo(), tolerance(), strategy(-1) {}
            Algo(const std::string &tystr, const std::string &str, float tol=-1.f, int strategy=-1) :type(tystr), algo(str), tolerance(tol), strategy(strategy) {}
            std::string type; std::string algo; float tolerance; int strategy;
            static float default_tolerance() { return 0.1; }
            static int   default_strategy() { return -1; }
        };
        /// list of algorithms to run if the default one fails
        static std::vector<Algo> fallbacks_;
        /// do a pre-scan
        static bool preScan_;
        /// do a pre-fit (with larger tolerance)
        static double approxPreFitTolerance_;
        /// do a pre-fit (with larger tolerance)
        static int approxPreFitStrategy_;
        /// do a pre-fit (w/o nuisances)
        static int preFit_;
        /// do first a fit of only the POI
        static bool poiOnlyFit_;
        /// do first a minimization of each nuisance individually 
        static bool singleNuisFit_;
        /// do first a minimization of each nuisance individually 
        static float nuisancePruningThreshold_;
        /// do first a fit of only the POI
        static bool setZeroPoint_;
        /// don't do old fallback using robustMinimize 
        static bool oldFallback_;
        /// call Hesse before or after the minimization 
        static bool firstHesse_, lastHesse_;
        /// storage level for minuit2 (toggles storing of intermediate covariances)
        static int minuit2StorageLevel_;

	static double discreteMinTol_;

	static std::string defaultMinimizerType_;
	static std::string defaultMinimizerAlgo_;
	static double 	   defaultMinimizerTolerance_;
	static double 	   defaultMinimizerPrecision_;
	//static int 	   defaultMinimizerStrategy_;

    	static bool runShortCombinations; 
        //static void setDefaultIntegrator(RooCategory &cat, const std::string & val) ;

	static std::map<std::string,std::vector<std::string> > const minimizerAlgoMap_;
	static bool checkAlgoInType(std::string, std::string );
};

// Singleton Class inside!
class CascadeMinimizerGlobalConfigs{

	private:
	  CascadeMinimizerGlobalConfigs(){};
	  CascadeMinimizerGlobalConfigs(const CascadeMinimizerGlobalConfigs&) {};
	  CascadeMinimizerGlobalConfigs& operator=(const CascadeMinimizerGlobalConfigs&);

	public:

	  //RooCategory* x;
	  RooListProxy pdfCategories; 
	  RooListProxy nuisanceParameters; 
	  RooListProxy allFloatingParameters; 
	  RooListProxy parametersOfInterest; 
	  RooListProxy allRooMultiPdfParams;
	  RooListProxy allRooMultiPdfs;

	  static CascadeMinimizerGlobalConfigs& O(){

		static CascadeMinimizerGlobalConfigs singleton;
		return singleton;
	  }

};

#endif
