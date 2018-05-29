#ifndef HiggsAnalysis_CombinedLimit_RobustHesse_h
#define HiggsAnalysis_CombinedLimit_RobustHesse_h

#include <vector>
#include <string>
#include <unordered_map>
// #include <TGraphAsymmErrors.h>
// #include <TString.h>
// #include <RooHistError.h>
// #include <RooFitResult.h>
// #include <TH1.h>
#include "RooAbsReal.h"
#include "RooRealVar.h"

class RobustHesse {
 public:
  RobustHesse(RooAbsReal &nll);

  void SaveHessianToFile(std::string const& filename);
  void LoadHessianFromFile(std::string const& filename);

  int hesse();

 private:
  int factorial(int n) {
    if (n > 1)
      return n * factorial(n - 1);
    else
      return 1;
  }

  void initialize();

  double deltaNLL();
  double deltaNLL(std::vector<unsigned> const& indices, std::vector<double> const& vals);

  int setParameterStencil(unsigned i);


  std::pair<int, double> findBound(unsigned i, double x, double initialDelta, double initialMult, double scaleMult, double threshold, double hardBound, unsigned maxIters);

  double improveWithBisect(unsigned i, double x, double max, double target, double targetLo, double targetHi, unsigned maxIters);

  std::vector<double> getFdCoeffs(unsigned n, std::vector<double> const& stencil);
  RooAbsReal * nll_;
  double nll0_;
  std::vector<RooRealVar *> v_;
  std::vector<double> nominal_;
  std::vector<std::vector<double>> stencils_;
  std::vector<std::vector<double>> d1coeffs_;
  std::vector<std::vector<double>> d2coeffs_;
  std::vector<double> rescales_;
  std::vector<std::vector<double>> sampled_;


  // Parameters that control the behaviour
  double targetNllForStencils_;
  double minNllForStencils_;
  double maxNllForStencils_;

  bool doRescale_;

  std::string saveFile_;
  std::string loadFile_;

  std::map<std::pair<unsigned, double>, double> nllcache_;
};

#endif
