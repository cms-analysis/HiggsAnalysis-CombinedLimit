#include "../interface/CMSHistFuncWrapper.h"
#include <vector>
#include <ostream>
#include "RooRealProxy.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"

#define HFVERBOSE 0

CMSHistFuncWrapper::CMSHistFuncWrapper()
    : pfunc_(nullptr), perr_(nullptr), initialized_(false) {
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
      initialized_(false) {
  cache_ = func.cache();
}

CMSHistFuncWrapper::CMSHistFuncWrapper(CMSHistFuncWrapper const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      func_("func", this, other.func_),
      err_("err", this, other.err_),
      cache_(other.cache_),
      idx_(other.idx_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      pfunc_(nullptr),
      perr_(nullptr),
      initialized_(false) {
}

void CMSHistFuncWrapper::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  pfunc_ = dynamic_cast<CMSHistFunc const*>(&(func_.arg()));
  perr_ = dynamic_cast<CMSHistErrorPropagator *>(err_.absArg());
  auto sentry_args = perr_->getSentryArgs();
  for (RooAbsArg *arg : *sentry_args) {
    sentry_.addArg(*arg);
  }
  sentry_.setValueDirty();
  initialized_ = true;
}

void CMSHistFuncWrapper::updateCache() const {
  initialize();

  // The ErrorPropagator will send a dirty flag whenever we need to update the cache
  if (!sentry_.good()) {
    perr_->applyErrorShifts(idx_, pfunc_->cache(), cache_);
    // cache_.CropUnderflows();
#if HFVERBOSE > 2
    std::cout << "Updated cache from CMSHistFunc:\n";
    pfunc_->cache().Dump();
    std::cout << "After shifts and cropping:\n";
    cache_.Dump();
#endif
  }
  sentry_.reset();
}

Double_t CMSHistFuncWrapper::evaluate() const {
  updateCache();
  return cache_.GetAt(x_);

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
  switch (code) {
    case 1: {
      updateCache();
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitHS3/JSONIO.h"
#include "static_execute.h"

namespace {
using RooFit::Detail::JSONNode;

// Wrappers are plotting-only (not in the likelihood) but always exported by ROOT: write a frozen "step" snapshot, as a declining exporter would leave a broken stub.
class CMSHistFuncWrapperStreamer : public RooFit::JSONIO::Exporter {
public:
   // func_/err_ have no exporters of their own and are already inlined in the snapshot.
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "step";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto const *wrapper = dynamic_cast<const CMSHistFuncWrapper *>(p);
     if (!wrapper) return false;

     wrapper->updateCache();
     auto const &c = wrapper->cache();

     // exportHisto() takes the bin count from the variable, which need not match the cache: use a throwaway clone.
     auto const *xOrig = dynamic_cast<RooRealVar const *>(&wrapper->getXVar());
     if (!xOrig) return false;
     RooRealVar xTmp(xOrig->GetName(), xOrig->GetTitle(), xOrig->getMin(), xOrig->getMax());
     xTmp.setBins(static_cast<int>(c.fullsize()));
     RooArgSet obs(xTmp);
     tool->exportHisto(obs, c.fullsize(), c.Values().data(), elem["data"].set_map());
     // The snapshot reimports as a RooHistFunc of the real observable, which ROOT exports as a
     // dependant; export it here too, so the first export already matches.
     tool->queueExport(*xOrig);
     elem["type"] << key();
     return true;
   }
};

STATIC_EXECUTE([]() {
  using namespace RooFit::JSONIO;
  registerExporter<CMSHistFuncWrapperStreamer>(CMSHistFuncWrapper::Class(), false);
});
}  // namespace

#undef HFVERBOSE
