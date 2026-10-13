#include "../interface/CMSHistSumWrapper.h"
#include <ostream>
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"

CMSHistSumWrapper::CMSHistSumWrapper()
    : psum_(nullptr), initialized_(false) {
  idx_ = 0;
}

CMSHistSumWrapper::CMSHistSumWrapper(const char* name, const char* title, RooAbsReal& x,
                                     CMSHistSum& sum, unsigned idx)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      sum_("sum", "", this, sum),
      idx_(idx),
      sentry_(TString(name) + "_sentry", ""),
      psum_(nullptr),
      initialized_(false) {
  cache_ = sum.cache();
  cache_.Clear();
}

CMSHistSumWrapper::CMSHistSumWrapper(CMSHistSumWrapper const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      sum_("sum", this, other.sum_),
      cache_(other.cache_),
      idx_(other.idx_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      psum_(nullptr),
      initialized_(false) {
}

void CMSHistSumWrapper::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  psum_ = dynamic_cast<CMSHistSum const*>(&(sum_.arg()));
  auto sentry_args = psum_->getSentryArgs();
  for (RooAbsArg *arg : *sentry_args) {
    sentry_.addArg(*arg);
  }
  sentry_.setValueDirty();
  initialized_ = true;
}

void CMSHistSumWrapper::updateCache() const {
  initialize();

  // The CMSHistSum sentries will send a dirty flag whenever we need to update the cache
  if (!sentry_.good()) {
    psum_->fillProcessShape(idx_, cache_);
  }
  sentry_.reset();
}

Double_t CMSHistSumWrapper::evaluate() const {
  updateCache();
  return cache_.GetAt(x_);
}

void CMSHistSumWrapper::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  cache_.Dump();
  std::cout << ">> Sentry: " << sentry_.good() << "\n";
  sentry_.Print("v");
}

Int_t CMSHistSumWrapper::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistSumWrapper::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  switch (code) {
    case 1: {
      updateCache();
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}
