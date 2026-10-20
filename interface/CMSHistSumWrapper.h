#ifndef CMSHistSumWrapper_h
#define CMSHistSumWrapper_h
#include <ostream>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "FastTemplate_Old.h"
#include "SimpleCacheSentry.h"
#include "CMSHistSum.h"

// Per-process view of a CMSHistSum: evaluates the morphed template of process
// idx including its share of the bin-wise (autoMCStats) shifts. Plays the role
// that CMSHistFuncWrapper has for CMSHistErrorPropagator. Either stored in the
// workspace by text2workspace or created on the fly by CMSHistSum::wrapperList().
class CMSHistSumWrapper : public RooAbsReal {
 public:
  CMSHistSumWrapper();

  CMSHistSumWrapper(const char* name, const char* title, RooAbsReal& x,
                    CMSHistSum& sum, unsigned idx);

  CMSHistSumWrapper(CMSHistSumWrapper const& other, const char* name = 0);

  TObject* clone(const char* newname) const override {
    return new CMSHistSumWrapper(*this, newname);
  }
  ~CMSHistSumWrapper() override {}

  Double_t evaluate() const override;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const override;

  // x_ is also reached via sum_, which would otherwise make
  // RooRealIntegral fall back to numeric integration over x
  Bool_t forceAnalyticalInt(const RooAbsArg& dep) const override {
    return dep.namePtr() == x_.arg().namePtr();
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const override;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const override;

  void updateCache() const;

  inline FastHisto const& cache() const { return cache_; }

  inline unsigned idx() const { return idx_; }

 protected:
  RooRealProxy x_;
  RooRealProxy sum_;

  mutable FastHisto cache_;
  unsigned idx_;
  mutable SimpleCacheSentry sentry_; //!

 private:
  ClassDefOverride(CMSHistSumWrapper, 1)
  mutable CMSHistSum const* psum_; //!
  mutable bool initialized_; //! not to be serialized

  void initialize() const;
};

#endif
