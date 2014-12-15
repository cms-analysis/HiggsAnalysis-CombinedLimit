#ifndef ICHiggsTauTau_CombineTools_MorphFunctions_h
#define ICHiggsTauTau_CombineTools_MorphFunctions_h
#include <algorithm>
#include <vector>
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "CombineTools/interface/CombineHarvester.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "HiggsAnalysis/CombinedLimit/interface/RooMorphingPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/AsymPow.h"
#pragma GCC diagnostic pop

namespace ch {

inline TH1F AsTH1F(TH1 const* hist) {
  TH1F res;
  TH1F const* test_f = dynamic_cast<TH1F const*>(hist);
  TH1D const* test_d = dynamic_cast<TH1D const*>(hist);
  if (test_f) {
    test_f->Copy(res);
  } else if (test_d) {
    test_d->Copy(res);
  } else {
    throw std::runtime_error("TH1 is not a TH1F or a TH1D");
  }
  return res;
}

void BuildRooMorphing(RooWorkspace& ws, CombineHarvester& cb,
                      RooAbsReal& mass_var, bool verbose,
                      std::string norm_postfix = "_norm");

void BuildRooMorphing(RooWorkspace& ws, CombineHarvester& cb,
                      std::string const& bin, std::string const& process,
                      RooAbsReal& mass_var, std::string const& mass_min,
                      std::string const& mass_max,
                      bool allow_morph, bool verbose);
}
#endif
