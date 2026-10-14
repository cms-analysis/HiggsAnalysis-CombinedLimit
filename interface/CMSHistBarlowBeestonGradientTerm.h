#ifndef HiggsAnalysis_CombinedLimit_CMSHistBarlowBeestonGradientTerm_h
#define HiggsAnalysis_CombinedLimit_CMSHistBarlowBeestonGradientTerm_h

#include <RooAbsReal.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooListProxy.h>
#include <RooRealProxy.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TObject.h>

#include <vector>

/// Correction of the analytic gradient of a likelihood whose Barlow-Beeston
/// parameters are profiled in closed form. Where the closed-form solution sits
/// at the boundary of zero expected yield, the derivative of the profiled
/// likelihood picks up an extra term: the slope of the constraint at the
/// solution times the derivative of the yield-to-error ratio of the bin, see
/// cmsHistBarlowBeestonGradientTerms().
///
/// This function has the value zero and the analytical integral
/// sum_j slope_j * (ratio_j - ratio0_j) over the bins with one Gaussian bin
/// parameter. Added with a unit coefficient to the extended RooRealSumPdf of a
/// channel, it adds exactly that sum to the extended likelihood. Before a
/// gradient call, cacheutils::ChannelSumNLL sets the slopes and the reference
/// ratios `ratio0` to the current ratios, so the sum vanishes but its
/// derivative in the generated code is the missing term; at all other times
/// the slopes are zero. Slopes and reference ratios are constant RooRealVars,
/// so that the generated code reads them at run time.
///
/// The histogram is deliberately not a server: a branch that depends on the
/// observable would make RooRealIntegral integrate numerically instead of
/// calling analyticalIntegral(). The generated integral refers to the arrays
/// of the histogram by name.
class CMSHistBarlowBeestonGradientTerm : public RooAbsReal {
 public:
  CMSHistBarlowBeestonGradientTerm() = default;

  /// `hist` is a CMSHistSum or a CMSHistErrorPropagator with observable `x`.
  CMSHistBarlowBeestonGradientTerm(const char* name, const char* title, RooAbsReal& x, RooAbsReal const& hist);

  CMSHistBarlowBeestonGradientTerm(CMSHistBarlowBeestonGradientTerm const& other, const char* name = nullptr);

  ~CMSHistBarlowBeestonGradientTerm() override = default;
  CMSHistBarlowBeestonGradientTerm& operator=(CMSHistBarlowBeestonGradientTerm const&) = delete;
  CMSHistBarlowBeestonGradientTerm(CMSHistBarlowBeestonGradientTerm&&) = delete;
  CMSHistBarlowBeestonGradientTerm& operator=(CMSHistBarlowBeestonGradientTerm&&) = delete;

  // The clone() contract of RooFit is a raw owning pointer.
  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  TObject* clone(const char* newname) const override { return new CMSHistBarlowBeestonGradientTerm(*this, newname); }

  RooAbsReal const& histogram() const { return *m_hist; }
  /// The slopes and reference ratios, one per bin with a single Gaussian bin
  /// parameter of the histogram, in bin order.
  RooArgList const& slopes() const { return m_slopes; }
  RooArgList const& ratios() const { return m_ratios; }
  /// The indices of these bins in the histogram.
  std::vector<unsigned> const& bins() const { return m_bins; }

  /// Take the slopes and the reference ratios from the current state of the
  /// histogram.
  void setFromHistogram();
  /// Set all slopes to zero.
  void clearSlopes();

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = nullptr) const override;
  double analyticalIntegral(Int_t code, const char* rangeName = nullptr) const override;

 protected:
  double evaluate() const override { return 0.; }

 private:
  RooRealProxy m_x;
  RooAbsReal const* m_hist = nullptr;  //! not a server, see above
  RooListProxy m_slopes;
  RooListProxy m_ratios;
  std::vector<unsigned> m_bins;

  ClassDefOverride(CMSHistBarlowBeestonGradientTerm, 1)
};

#endif
