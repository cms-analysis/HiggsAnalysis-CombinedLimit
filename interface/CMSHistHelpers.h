#ifndef HiggsAnalysis_CombinedLimit_CMSHistHelpers_h
#define HiggsAnalysis_CombinedLimit_CMSHistHelpers_h

#include "CMSHistErrorPropagator.h"
#include "CMSHistSum.h"

#include <stdexcept>
#include <string>

/// Whether a node is the histogram node of a template channel: a CMSHistSum or
/// a CMSHistErrorPropagator. The two classes share the interface of the
/// Barlow-Beeston minimization, but no base class, see cmsHistVisit().
inline bool isCMSHist(RooAbsArg const &arg) {
  return dynamic_cast<CMSHistSum const *>(&arg) || dynamic_cast<CMSHistErrorPropagator const *>(&arg);
}

/// Calls f with the histogram node cast to its class.
template <class F>
decltype(auto) cmsHistVisit(RooAbsReal const &hist, F &&f) {
  if (auto const *sum = dynamic_cast<CMSHistSum const *>(&hist))
    return std::forward<F>(f)(*sum);
  if (auto const *prop = dynamic_cast<CMSHistErrorPropagator const *>(&hist))
    return std::forward<F>(f)(*prop);
  throw std::invalid_argument(std::string(hist.GetName()) + " is neither a CMSHistSum nor a CMSHistErrorPropagator");
}

#endif
