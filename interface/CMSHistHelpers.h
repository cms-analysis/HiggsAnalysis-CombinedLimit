#ifndef HiggsAnalysis_CombinedLimit_CMSHistHelpers_h
#define HiggsAnalysis_CombinedLimit_CMSHistHelpers_h

#include "CMSHistErrorPropagator.h"
#include "CMSHistSum.h"

#include <stdexcept>
#include <string>
#include <vector>

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


/// Ingredients for the gradient of the analytically profiled likelihood of a
/// histogram node (CMSHistSum::barlowBeestonGradientTerms() and the same for
/// CMSHistErrorPropagator), one entry per bin with a single Gaussian parameter
/// for the total, in bin order: the ratio of the summed yield to the total
/// error (without the bin parameters), and the derivative of the profiled
/// likelihood with respect to that ratio. The derivative is zero unless the
/// parameter is profiled analytically and its solution sits at the boundary
/// of zero expected yield. There the likelihood has a kink and the derivative
/// of the profiled likelihood is not the partial derivative at the solution:
/// the missing term is the slope times the derivative of the ratio.
template <class ValSum, class TotErr, class BarlowBeeston>
void cmsHistBarlowBeestonGradientTerms(std::vector<std::vector<unsigned>> const &binTypes,
                                       ValSum const &valsum,
                                       TotErr const &toterr,
                                       BarlowBeeston const &bb,
                                       std::vector<double> &slopes,
                                       std::vector<double> &ratios) {
  slopes.clear();
  ratios.clear();
  unsigned k = 0;  // position in bb.use, which lists the profiled bins in increasing order
  for (unsigned j = 0; j < binTypes.size(); ++j) {
    if (binTypes[j][0] != 1)
      continue;
    ratios.push_back(toterr[j] > 0. ? valsum[j] / toterr[j] : 0.);
    double slope = 0.;
    if (bb.init && k < bb.use.size() && bb.use[k] == j) {
      // The solution is at the boundary if the shifted yield is cropped
      // (compare updateCache()). The solution is then minus the ratio, the
      // Poisson term is flat, and only the Gaussian constraint (unit width)
      // contributes to the derivative: (res - gobs) * d(res)/d(ratio).
      if (valsum[j] + (toterr[j] * bb.res[k]) < 1e-9)
        slope = bb.gobs[k] - bb.res[k];
      ++k;
    }
    slopes.push_back(slope);
  }
}

#endif
