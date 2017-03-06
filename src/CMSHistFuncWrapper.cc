#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFuncWrapper.h"
#include <vector>
#include <ostream>
#include "RooRealProxy.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"

CMSHistFuncWrapper::CMSHistFuncWrapper()
    : pfunc_(nullptr), perr_(nullptr), v(0), initialized_(false) {
  v = 0;
  idx_ = 0;
}

CMSHistFuncWrapper::CMSHistFuncWrapper(const char* name, const char* title, RooRealVar& x,
              CMSHistFunc & func, CMSHistErrorPropagator & err, unsigned idx)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      func_("func", "", this, func),
      err_("err", "", this, err, true, false),
      idx_(idx),
      sentry_(TString(name) + "_sentry", ""),
      pfunc_(nullptr),
      perr_(nullptr),
      v(0),
      initialized_(false) {
  cache_ = func.getCacheHisto();
}

CMSHistFuncWrapper::CMSHistFuncWrapper(CMSHistFuncWrapper const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      func_("func", this, other.func_),
      err_("err", this, other.err_),
      cache_(other.cache_),
      idx_(other.idx_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.sentry_.GetName()), ""),
      pfunc_(nullptr),
      perr_(nullptr),
      v(other.v),
      initialized_(false) {
}

void CMSHistFuncWrapper::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  pfunc_ = dynamic_cast<CMSHistFunc const*>(&(func_.arg()));
  perr_ = dynamic_cast<CMSHistErrorPropagator *>(err_.absArg());
  // sentry_.addArg(*err_.absArg());
  auto sentry_args = perr_->getSentryArgs();
  RooFIter iter = sentry_args->fwdIterator() ;
  RooAbsArg* arg;
  while((arg = iter.next())) {
    sentry_.addArg(*arg);
  }
  // sentry_args->
  // sentry_.addArg(*perr_->getSentryArgs()->first());
  sentry_.setValueDirty();
  initialized_ = true;
}

void CMSHistFuncWrapper::updateCache() const {
  // FNLOGC(std::cout, v) << "CMSHistFuncWrapper::updateCache()\n";
  // FNLOGC(std::cout, v) << pfunc_ << "\t" << perr_ << "\n";

  initialize();

  // FNLOGC(std::cout, v) << "Sentry: " << sentry_.good() << "\n";

  // The ErrorPropagator will send a dirty flag whenever we need to update the cache
  if (!sentry_.good()) {
    perr_->applyErrorShifts(idx_, pfunc_->getCacheHisto(), cache_);
    cache_.CropUnderflows();
    // if (v > 1) {
    //   FNLOG(std::cout) << "Updated cache from CMSHistFunc:\n";
    //   pfunc_->getCacheHisto().Dump();
    //   FNLOG(std::cout) << "After shifts and cropping:\n";
    //   cache_.Dump();
    // }
  }
  sentry_.reset();
}

Double_t CMSHistFuncWrapper::evaluate() const {
  // LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  updateCache();
  return cache_.GetAt(x_);
  // Just for testing - bypass everything
  // initialize();
  // return pfunc_->evaluate();
}

void CMSHistFuncWrapper::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  cache_.Dump();
  std::cout << ">> Sentry: " << sentry_.good() << "\n";
  sentry_.Print("v");
}

Int_t CMSHistFuncWrapper::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistFuncWrapper::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache();
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}


CMSHistFuncWrapperV::CMSHistFuncWrapperV(const CMSHistFuncWrapper& hpdf, const RooAbsData& data,
                           bool includeZeroWeights)
    : hpdf_(hpdf), begin_(0), end_(0) {
  hpdf.updateCache();
  std::vector<int> bins;
  RooArgSet obs(hpdf.x_.arg());
  const RooRealVar& x = static_cast<const RooRealVar&>(*obs.first());
  bool aligned = true;
  for (int i = 0, n = data.numEntries(); i < n; ++i) {
    obs = *data.get(i);
    if (data.weight() == 0 && !includeZeroWeights) continue;
    int idx = hpdf.cache_.FindBin(x.getVal());
    if (!bins.empty() && idx != bins.back() + 1) aligned = false;
    bins.push_back(idx);
  }
  if (bins.empty()) {
    // nothing to do.
  } else if (aligned) {
    begin_ = bins.front();
    end_ = bins.back() + 1;
    // std::cout << "Created CMSHistFuncWrapperV from " << hpdf.GetName() << ",
    // aligned, " << (end_-begin_) << " bins." << std::endl;
  } else {
    nbins_ = bins.size();
    bins_.swap(bins);
    blocks_.clear();
    int start = bins_[0], istart = 0;
    for (int i = 1, n = bins_.size(); i < n; ++i) {
      if (bins_[i] != bins_[i - 1] + 1) {
        blocks_.push_back(Block(istart, start, bins_[i - 1] + 1));
        start = bins_[i];
        istart = i;
      }
    }
    blocks_.push_back(Block(istart, start, bins_.back() + 1));
    if (blocks_.size() < 4 * bins_.size()) {
      // std::cout << "Created CMSHistFuncWrapperV from " << hpdf.GetName() << ",
      // block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      bins_.clear();
    } else {
      // std::cout << "Created CMSHistFuncWrapperV from " << hpdf.GetName() << ",
      // non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      blocks_.clear();
    }
  }
}

void CMSHistFuncWrapperV::fill(std::vector<Double_t>& out) const {
  hpdf_.updateCache();
  if (begin_ != end_) {
    out.resize(end_ - begin_);
    std::copy(&hpdf_.cache_.GetBinContent(begin_),
              &hpdf_.cache_.GetBinContent(end_), out.begin());
  } else if (!blocks_.empty()) {
    out.resize(nbins_);
    for (auto b : blocks_)
      std::copy(&hpdf_.cache_.GetBinContent(b.begin),
                &hpdf_.cache_.GetBinContent(b.end), out.begin() + b.index);
  } else {
    out.resize(bins_.size());
    for (int i = 0, n = bins_.size(); i < n; ++i) {
      out[i] = hpdf_.cache_.GetBinContent(bins_[i]);
    }
  }
}
