#include "HiggsAnalysis/CombinedLimit/interface/RooMorphingPdf2.h"
#include <stdexcept>
#include <vector>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"

CMSHistFunc::CMSHistFunc() {}

CMSHistFunc::CMSHistFunc(const char* name, const char* title, RooRealVar& x,
                         TH1 const& hist)
    : RooAbsReal(name, title), x_("x", "", this, x), cache_(hist) {
  //
}

CMSHistFunc::CMSHistFunc(CMSHistFunc const& other, const char* name)
    : RooAbsReal(other, name), x_("x", this, other.x_), cache_(other.cache_) {}

Double_t CMSHistFunc::evaluate() const {
  return cache_.GetAt(x_);
}
