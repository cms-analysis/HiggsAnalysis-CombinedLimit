#ifndef HiggsAnalysis_CombinedLimit_ChannelCompatibilityCheck_h
#define HiggsAnalysis_CombinedLimit_ChannelCompatibilityCheck_h
/** \class ChannelCompatibilityCheck
 *
 * Do a ML fit of the data with background and signal+background hypothesis and print out diagnostics plots 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "FitterAlgoBase.h"

class ChannelCompatibilityCheck : public FitterAlgoBase {
public:
  ChannelCompatibilityCheck() ;
  const std::string & name() const override {
    static const std::string name("ChannelCompatibilityCheck");
    return name;
  }
  void applyOptions(const boost::program_options::variables_map &vm) override ;

protected:
  std::string nameForLabel(const char *label) ;

  static float mu_;
  static bool  fixedMu_;

  static bool runMinos_;
  static bool saveFitResult_;

  static std::vector<std::string> groups_;
  static std::map<TString,std::pair<double,double>> groupRanges_;

  bool runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) override;
};


#endif
