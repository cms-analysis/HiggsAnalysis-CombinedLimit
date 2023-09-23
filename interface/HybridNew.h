#ifndef HiggsAnalysis_CombinedLimit_HybridNew_h
#define HiggsAnalysis_CombinedLimit_HybridNew_h
/** \class HybridNew
 *
 * Module to compute limits by tossing toys (CLs, CLs+b, Feldman-Cousins), and related significances
 *
 * \author Luca Lista (INFN), Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "LimitAlgo.h"
#include <algorithm> 
#include <RooStats/ModelConfig.h>
#include <RooStats/HybridCalculator.h>
#include <RooStats/ToyMCSampler.h>
#include <TF1.h>

class RooRealVar;
class TGraphErrors;
class TDirectory;

class HybridNew : public LimitAlgo {
public:
  HybridNew() ; 
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;
  virtual void applyDefaultOptions() ; 

  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual bool runLimit(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual bool runSignificance(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual bool runSinglePoint(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual bool runTestStatistics(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual const std::string & name() const {
    static const std::string name("HybridNew");
    return name;
  }
  
  // I don't understand why if it's declared private it gcc does not like it
  enum WorkingMode { MakeLimit, MakeSignificance, MakePValues, MakeTestStatistics, MakeSignificanceTestStatistics };
private:
  static WorkingMode workingMode_;
  static unsigned int nToys_;
  static double clsAccuracy_, rAbsAccuracy_, rRelAccuracy_, interpAccuracy_;
  static bool CLs_;
  static std::string rule_, testStat_;
  static bool reportPVal_;
  static bool genNuisances_, genGlobalObs_, fitNuisances_;
  static std::string         rValue_;   // name=value,name=value,... or just value for a 1D model
  static RooArgSet           rValues_;  // values of the parameters
  static unsigned int iterations_;
  static bool saveHybridResult_, readHybridResults_; 
  static std::string gridFile_;
  static bool expectedFromGrid_, clsQuantiles_; 
  static float quantileForExpectedFromGrid_;
  static bool fullBToys_; 
  static bool fullGrid_; 
  static bool saveGrid_; 
  static bool noUpdateGrid_; 
  static unsigned int nCpu_, fork_;
  static bool importanceSamplingNull_, importanceSamplingAlt_;
  static std::string algo_;
  static std::string mode_;
  static std::string plot_;
//  static std::string minimizerAlgo_;
//  static float       minimizerTolerance_;

  static bool optimizeProductPdf_;
  static bool optimizeTestStatistics_;
  static bool newToyMCSampler_;
  static bool rMinSet_;
  static bool rMaxSet_;
  float mass_;
  static std::string	scaleAndConfidenceSelection_; 
  static float maxProbability_;
  static float confidenceToleranceForToyScaling_;
  static float adaptiveToys_;

  static double EPS;
  // graph, used to compute the limit, not just for plotting!
  std::unique_ptr<TGraphErrors> limitPlot_;
 
  // performance counter: remember how many toys have been thrown
  unsigned int perf_totalToysRun_;

  //----- extra variables used for cross-checking the implementation of frequentist toy tossing in RooStats
  // mutable RooAbsData *realData_;
  // std::unique_ptr<RooAbsCollection>  snapGlobalObs_;

  struct Setup {
    RooStats::ModelConfig modelConfig, modelConfig_bonly;
    std::unique_ptr<RooStats::TestStatistic> qvar;
    std::unique_ptr<RooStats::ToyMCSampler>  toymcsampler;
    std::unique_ptr<RooStats::ProofConfig> pc;
    RooArgSet cleanupList;
  };

  void validateOptions() ;

  // make sure our rValues_ is contains all pois in the model, and does not contain anything else
  void setupPOI(RooStats::ModelConfig *mc_s) ;

  /// generic version
  std::pair<double,double> eval(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, const RooAbsCollection & rVals, bool adaptive=false, double clsTarget=-1) ;
  /// specific version used in unidimensional searches (calls the generic one)
  std::pair<double,double> eval(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double rVal, bool adaptive=false, double clsTarget=-1)  ;

  /// generic version
  std::unique_ptr<RooStats::HybridCalculator> create(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, const RooAbsCollection & rVals, Setup &setup);
  /// specific version used in unidimensional searches (calls the generic one)
  std::unique_ptr<RooStats::HybridCalculator> create(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double rVal, Setup &setup);

  std::pair<double,double> eval(RooStats::HybridCalculator &hc, const RooAbsCollection & rVals, bool adaptive=false, double clsTarget=-1) ;

  // Compute CLs or Pmu from a HypoTestResult. In the case of the LHCFC test statistics, do the proper transformation to LHC or Profile (which needs to know the POI)
  std::pair<double,double> eval(const RooStats::HypoTestResult &hcres, const RooAbsCollection & rVals) ;
  std::pair<double,double> eval(const RooStats::HypoTestResult &hcres, double rVal) ;

  void applyExpectedQuantile(RooStats::HypoTestResult &hcres);
  void applyClsQuantile(RooStats::HypoTestResult &hcres);
  void applySignalQuantile(RooStats::HypoTestResult &hcres);
  RooStats::HypoTestResult *evalGeneric(RooStats::HybridCalculator &hc, bool forceNoFork=false);
  RooStats::HypoTestResult *evalWithFork(RooStats::HybridCalculator &hc);
  // RooStats::HypoTestResult *evalFrequentist(RooStats::HybridCalculator &hc);  // cross-check implementation, 
  RooStats::HypoTestResult *readToysFromFile(const RooAbsCollection & rVals);

  std::map<double, RooStats::HypoTestResult *> grid_;

  void clearGrid(); 
  void readGrid(TDirectory *directory, double rMin, double rMax); 
  void updateGridData(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, bool smart, double clsTarget); 
  void updateGridDataFC(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, bool smart, double clsTarget); 
  std::pair<double,double> updateGridPoint(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, std::map<double, RooStats::HypoTestResult *>::iterator point);
  std::pair<double,double> interpolateAndUncert(TGraphErrors *gr, double clsTarget);
  std::vector<std::pair<double, double> > findIntervalsFromSplines(TGraphErrors *limitPlot_,double clsTarget);


  void useGrid();

  bool doFC_;
  
};

#endif
