#ifndef HiggsAnalysis_CombinedLimit_CMSHistHelpers_h
#define HiggsAnalysis_CombinedLimit_CMSHistHelpers_h

#include "CMSHistErrorPropagator.h"
#include "CMSHistSum.h"
#include "FastTemplate_Old.h"

#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 32, 0)
#include "RooFit/EvalContext.h"
#endif
#include "RooAbsArg.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
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

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 32, 0)
/// Fills the output of a batched RooFit evaluation (RooAbsReal::doEval) of one
/// of the CMSHist classes from its histogram cache: for each entry, the content
/// of the bin that the observable falls into. The observable is a scalar when
/// it is not read from a dataset, e.g. when the function is evaluated by a
/// numeric integral. The doEval() overrides serve the RooFit Evaluator ("cpu"
/// backend and the value of the "codegen" backend before ROOT 6.40); the
/// fallback RooAbsReal::doEval() would loop over all servers and all bins and
/// call evaluate() per bin, which costs more than the update of the cache.
inline void cmsHistFillBatch(RooFit::EvalContext &ctx, RooRealProxy const &x, FastHisto const &hist) {
  auto xVals = ctx.at(x);
  auto output = ctx.output();
  for (std::size_t i = 0; i < output.size(); ++i) {
    output[i] = hist.GetAt(xVals.size() == 1 ? xVals[0] : xVals[i]);
  }
}
#endif

#endif
