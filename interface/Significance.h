#ifndef HiggsAnalysis_CombinedLimit_Significance_h
#define HiggsAnalysis_CombinedLimit_Significance_h
/** \class Significance
 *
 * Class for computing limits and significances from naive Profile Likelihood asymptotics (i.e. no Asimov dataset) 
 *
 * \author Giovanni Petrucciani (UCSD) 
 *
 *
 */
#include "LimitAlgo.h"

class RooAbsPdf; class RooRealVar; class RooAbsData; class RooArgSet;

class Significance : public LimitAlgo {
public:
  Significance() ;
  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual const std::string & name() const {
    static const std::string name("Significance");
    return name;
  }
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;

  /// Setup Minimizer configuration on creation, reset the previous one on destruction.
  class MinimizerSentry {
     public:
        MinimizerSentry(const std::string &algo, double tolerance);
        ~MinimizerSentry();
     private:
        std::string minimizerTypeBackup, minimizerAlgoBackup;
        double minimizerTollBackup;
  };

protected:
  //static std::string minimizerAlgo_, 
  static std::string minimizerAlgoForBF_;
  //static float       minimizerTolerance_, 
  static float minimizerToleranceForBF_;

  static bool useMinos_, bruteForce_;
  static std::string bfAlgo_;
  static int  points_;

  // ----- options for handling cases where the likelihood fit misbihaves ------
  /// compute the limit N times
  static int         tries_;
  /// trying up to M times from different points
  static int         maxTries_;
  /// maximum relative deviation of the different points from the median to accept 
  static float       maxRelDeviation_;
  /// Ignore up to this fraction of results if they're too far from the median
  static float       maxOutlierFraction_;
  /// Stop trying after finding N outliers
  static int         maxOutliers_;
  /// Try first a plain fit
  static bool        preFit_;

  /// Report p-value instead of significance
  static bool reportPVal_;
  static bool uncapped_;

  static float signalForSignificance_;
  static float mass_;

  static std::string plot_;

  bool runSignificance(RooWorkspace *w, RooStats::ModelConfig *mc, RooAbsData &data, double &limit, double &limitErr);
  bool runLimit(RooWorkspace *w, RooStats::ModelConfig *mc, RooAbsData &data, double &limit, double &limitErr);

  double upperLimitWithMinos(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance, double cl) const ;
  std::pair<double,double> upperLimitBruteForce(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance, double cl) const ;
  double significanceBruteForce(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance) const ;
  double significanceFromScan(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance, int npoints) const ;
};

#endif
