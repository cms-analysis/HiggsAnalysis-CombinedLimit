#ifndef HiggsAnalysis_CombinedLimit_GenerateOnly_h
#define HiggsAnalysis_CombinedLimit_GenerateOnly_h
/** \class GenerateOnly
 *
 * Class for generation of toy samples without any actual limit computation
 *
 * \author Luca Lista (INFN)
 *
 *
 */
#include "LimitAlgo.h"

class RooAbsPdf; class RooRealVar; class RooAbsData; class RooArgSet;

class GenerateOnly : public LimitAlgo {
public:
  GenerateOnly() ;
  bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) override;
  const std::string & name() const override {
    static const std::string name("GenerateOnly");
    return name;
  }
  void applyOptions(const boost::program_options::variables_map &vm) override ;
};

#endif
