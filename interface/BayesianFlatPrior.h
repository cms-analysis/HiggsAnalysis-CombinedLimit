#ifndef HiggsAnalysis_CombinedLimit_BayesianFlatPrior_h
#define HiggsAnalysis_CombinedLimit_BayesianFlatPrior_h
/** \class BayesianFlatPrior
 *
 * abstract interface for physics objects
 *
 * \author Luca Lista (INFN), from initial implementation by Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "LimitAlgo.h"

class BayesianFlatPrior : public LimitAlgo {
public:
  BayesianFlatPrior() ;
  bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) override;
  const std::string & name() const override {
    static const std::string name("BayesianSimple");
    return name;
  }
private:
  static int maxDim_;
};

#endif
