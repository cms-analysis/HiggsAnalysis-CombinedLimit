#ifndef HiggsAnalysis_CombinedLimit_RobustHesse_h
#define HiggsAnalysis_CombinedLimit_RobustHesse_h

#include <vector>
#include <string>
#include <unordered_map>
// #include <TGraphAsymmErrors.h>
// #include <TString.h>
// #include <RooHistError.h>
#include "RooFitResult.h"
// #include <TH1.h>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMatrixDSymEigen.h"

class RooFitResultBuilder : public RooFitResult {
 public:
  RooFitResultBuilder() : RooFitResult(){
    this->RooFitResult::setInitParList(RooArgList());
    this->RooFitResult::setConstParList(RooArgList());
  };

  RooFitResultBuilder(RooFitResult const& other) : RooFitResult(other) {}

  void setFinalParList(RooArgList const& pars) { this->RooFitResult::setFinalParList(pars); }

  void setCovarianceMatrix(TMatrixDSym & matrix) { this->RooFitResult::setCovarianceMatrix(matrix);}

  RooFitResult Get() { return *this; }
};

class RobustHesse {
 public:
  RobustHesse(RooAbsReal &nll, unsigned verbose = 0);

  void SaveHessianToFile(std::string const& filename);
  void LoadHessianFromFile(std::string const& filename);

  void ProtectArgSet(RooArgSet const& set);
  void ProtectVars(std::vector<std::string> const& names);

  int hesse();

  void WriteOutputFile(std::string const& outputFileName) const;
  RooFitResult * GetRooFitResult(RooFitResult const* current) const;

 private:
  int factorial(int n) {
    if (n > 1)
      return n * factorial(n - 1);
    else
      return 1;
  }

  struct Var {
    RooRealVar * v;
    double nominal;
    std::vector<double> stencil;
    std::vector<double> d1coeffs;
    std::vector<double> d2coeffs;
    double rescale;
  };

  void initialize();

  double deltaNLL();
  double deltaNLL(std::vector<unsigned> const& indices, std::vector<double> const& vals);

  int setParameterStencil(unsigned i);


  std::pair<int, double> findBound(unsigned i, double x, double initialDelta, double initialMult, double scaleMult, double threshold, double hardBound, unsigned maxIters);

  double improveWithBisect(unsigned i, double x, double max, double target, double targetLo, double targetHi, unsigned maxIters);

  std::vector<double> getFdCoeffs(unsigned n, std::vector<double> const& stencil);

  void RemoveFromHessian(std::vector<unsigned> const& ids);

  void ReplaceVars(std::vector<Var> newVars) {
    cVars_ = newVars;
    nllcache_.clear();
    nllEvals_ = 0;
    nllEvalsCached_ = 0;
  }



  RooAbsReal * nll_;
  double nll0_;

  std::vector<Var> cVars_;
  std::vector<Var> invalidStencilVars_;
  std::vector<Var> removedFromHessianVars_;

  // Parameters that control the behaviour
  double targetNllForStencils_;
  double minNllForStencils_;
  double maxNllForStencils_;
  unsigned maxRemovalsFromHessian_;

  bool doRescale_;

  std::string saveFile_;
  std::string loadFile_;

  std::unique_ptr<TMatrixDSym> hessian_;
  std::unique_ptr<TMatrixDSym> covariance_;

  std::map<std::pair<unsigned, double>, double> nllcache_;
  unsigned nllEvals_;
  unsigned nllEvalsCached_;

  std::set<std::string> proctected_;

  int verbosity_;
};

#endif
