#include "../interface/CMSHistBarlowBeestonGradientTerm.h"

#include "../interface/CMSHistHelpers.h"

#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RtypesCore.h>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

static void gradientTermsOf(RooAbsReal const& hist, std::vector<double>& slopes, std::vector<double>& ratios) {
  cmsHistVisit(hist, [&](auto const& h) { h.barlowBeestonGradientTerms(slopes, ratios); });
}

CMSHistBarlowBeestonGradientTerm::CMSHistBarlowBeestonGradientTerm(const char* name,
                                                                   const char* title,
                                                                   RooAbsReal& x,
                                                                   RooAbsReal const& hist)
    : RooAbsReal(name, title),
      m_x(/*theName=*/"x", /*desc=*/"observable", this, x),
      m_hist(&hist),
      m_slopes(/*inName=*/"slopes", /*desc=*/"slopes", this),
      m_ratios(/*inName=*/"ratios", /*desc=*/"reference ratios", this) {
  auto const& binTypes = cmsHistVisit(hist, [](auto const& h) -> auto const& { return h.binTypes(); });
  RooArgSet owned;
  for (unsigned j = 0; j < binTypes.size(); ++j) {
    if (binTypes[j][0] != 1)
      continue;
    m_bins.push_back(j);
    std::string const bin = std::to_string(j);
    auto slope = std::make_unique<RooRealVar>((std::string(name) + "_slope" + bin).c_str(), /*title=*/"", /*value=*/0.);
    auto ratio = std::make_unique<RooRealVar>((std::string(name) + "_ratio" + bin).c_str(), /*title=*/"", /*value=*/0.);
    slope->setConstant(/*value=*/true);
    ratio->setConstant(/*value=*/true);
    m_slopes.add(*slope);
    m_ratios.add(*ratio);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 26, 0)
    owned.addOwned(std::move(slope));
    owned.addOwned(std::move(ratio));
  }
  addOwnedComponents(std::move(owned));
#else
    owned.addOwned(*slope.release());
    owned.addOwned(*ratio.release());
  }
  // Older ROOT has no overloads that transfer the ownership out of `owned`.
  addOwnedComponents(owned);
  owned.releaseOwnership();
#endif
}

CMSHistBarlowBeestonGradientTerm::CMSHistBarlowBeestonGradientTerm(CMSHistBarlowBeestonGradientTerm const& other,
                                                                   const char* name)
    : RooAbsReal(other, name),
      m_x(/*theName=*/"x", this, other.m_x),
      m_hist(other.m_hist),
      m_slopes(/*inName=*/"slopes", this, other.m_slopes),
      m_ratios(/*inName=*/"ratios", this, other.m_ratios),
      m_bins(other.m_bins) {}

void CMSHistBarlowBeestonGradientTerm::setFromHistogram() {
  std::vector<double> slopes;
  std::vector<double> ratios;
  gradientTermsOf(*m_hist, slopes, ratios);
  for (std::size_t k = 0; k < m_bins.size(); ++k) {
    dynamic_cast<RooRealVar&>(m_slopes[static_cast<Int_t>(k)]).setVal(slopes[k]);
    dynamic_cast<RooRealVar&>(m_ratios[static_cast<Int_t>(k)]).setVal(ratios[k]);
  }
}

void CMSHistBarlowBeestonGradientTerm::clearSlopes() {
  for (RooAbsArg* slope : m_slopes) {
    dynamic_cast<RooRealVar&>(*slope).setVal(/*value=*/0.);
  }
}

Int_t CMSHistBarlowBeestonGradientTerm::getAnalyticalIntegral(RooArgSet& allVars,
                                                              RooArgSet& analVars,
                                                              const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, m_x))
    return 1;
  return 0;
}

double CMSHistBarlowBeestonGradientTerm::analyticalIntegral(Int_t /*code*/, const char* /*rangeName*/) const {
  std::vector<double> slopes;
  std::vector<double> ratios;
  gradientTermsOf(*m_hist, slopes, ratios);
  double sum = 0.;
  for (std::size_t k = 0; k < m_bins.size(); ++k) {
    sum += dynamic_cast<RooAbsReal const&>(m_slopes[static_cast<Int_t>(k)]).getVal() *
           (ratios[k] - dynamic_cast<RooAbsReal const&>(m_ratios[static_cast<Int_t>(k)]).getVal());
  }
  return sum;
}
