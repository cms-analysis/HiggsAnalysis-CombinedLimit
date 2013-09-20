#ifndef HiggsAnalysis_CombinedLimit_Asymptotic_h
#define HiggsAnalysis_CombinedLimit_Asymptotic_h
/** \class Asymptotic
 *
 * CLs asymptotic limits 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "../interface/LimitAlgo.h"
#include "../interface/utils.h"
#include <memory>
class RooRealVar;
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooFitResult.h>

#include "TFile.h"
#include "TTree.h"

class Asymptotic : public LimitAlgo {
public:
  Asymptotic() ; 
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;
  virtual void applyDefaultOptions() ; 

  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual bool runLimit(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  std::vector<std::pair<float,float> > runLimitExpected(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) ;

  float findExpectedLimitFromCrossing(RooAbsReal &nll, RooRealVar *r, double rMin, double rMax, double nll0, double quantile) ; 

  virtual const std::string& name() const { static std::string name_ = "Asymptotic"; return name_; }
private:
  static double rAbsAccuracy_, rRelAccuracy_;
  static std::string what_;
  static std::string gridFileName_;
  static bool qtilde_; 
  static bool picky_; 
  static bool noFitAsimov_; 
  static bool useGrid_; 
  static bool newExpected_; 
  static std::string minosAlgo_;
  static std::string minimizerAlgo_;
  static float       minimizerTolerance_;
  static int         minimizerStrategy_;

  static double rValue_;

  static bool   strictBounds_;

  bool    hasFloatParams_;
  bool    hasDiscreteParams_;
  mutable std::auto_ptr<RooArgSet>  params_;
  mutable std::auto_ptr<RooAbsReal> nllD_, nllA_; 
  //mutable std::auto_ptr<RooFitResult> fitFreeD_, fitFreeA_;
  //mutable std::auto_ptr<RooFitResult> fitFixD_,  fitFixA_;
  utils::CheapValueSnapshot fitFreeD_, fitFreeA_, fitFixD_,  fitFixA_;

  mutable double                      minNllD_,  minNllA_;
  mutable RooArgSet snapGlobalObsData, snapGlobalObsAsimov;

  float calculateLimitFromGrid(RooRealVar *, double, double);

  RooAbsData *asimovDataset(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data);
  double getCLs(RooRealVar &r, double rVal, bool getAlsoExpected=false, double *limit=0, double *limitErr=0);
  
  TFile *gridFile_;
  TTree *limitsTree_;
  double readCL_;
  double readMU_;

};

#endif
