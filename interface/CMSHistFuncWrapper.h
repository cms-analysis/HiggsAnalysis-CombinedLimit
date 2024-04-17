#ifndef CMSHistFuncWrapper_h
#define CMSHistFuncWrapper_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooArgProxy.h"
#include "Rtypes.h"
#include "FastTemplate_Old.h"
#include "SimpleCacheSentry.h"
#include "CMSHistErrorPropagator.h"
#include "CMSHistFunc.h"
#include "CMSHistV.h"

class CMSHistFuncWrapper : public RooAbsReal {
 public:
  CMSHistFuncWrapper();

  CMSHistFuncWrapper(const char* name, const char* title, RooRealVar& x,
              CMSHistFunc & func, CMSHistErrorPropagator & err, unsigned idx);

  CMSHistFuncWrapper(CMSHistFuncWrapper const& other, const char* name = 0);

  TObject* clone(const char* newname) const override {
    return new CMSHistFuncWrapper(*this, newname);
  }
  ~CMSHistFuncWrapper() override {}

  Double_t evaluate() const override;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const override;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const override;

  void updateCache() const;

  inline FastHisto const& cache() const { return cache_; }

  friend class CMSHistV<CMSHistFuncWrapper>;

 protected:
  RooRealProxy x_;
  RooRealProxy func_;
  RooArgProxy err_;

  mutable FastHisto cache_;
  unsigned idx_;
  mutable SimpleCacheSentry sentry_; //!

 private:
  ClassDefOverride(CMSHistFuncWrapper, 1)
  mutable CMSHistFunc const* pfunc_;
  mutable CMSHistErrorPropagator * perr_;
  mutable bool initialized_; //! not to be serialized

  void initialize() const;
};

#endif
