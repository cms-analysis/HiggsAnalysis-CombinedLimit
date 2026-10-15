#ifndef HiggsAnalysis_CombinedLimit_ChannelSumNLL_h
#define HiggsAnalysis_CombinedLimit_ChannelSumNLL_h

// ROOT_VERSION: ROOT/RVersion.hxx exists since ROOT 6.30, older versions only have RVersion.h.
#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif

// The RooFit evaluation backends that this class builds on exist since ROOT 6.30.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#include "SimNLLBase.h"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooListProxy.h>
#include <RooRealVar.h>
#include <RooSetProxy.h>
#include <RooSimultaneous.h>

class CMSHistBarlowBeestonGradientTerm;
class RooMultiPdf;
class RooPoisson;

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace cacheutils {

  /// Likelihood of a RooSimultaneous: one RooFit likelihood per channel plus one
  /// for the constraint terms, summed here.
  ///
  /// Each term comes from RooAbsPdf::createNLL() with the requested RooFit
  /// backend and owns its evaluator, so a masked channel can be left out of the
  /// sum, of getParameters() and of the evaluation, which a single evaluator
  /// can't do. The constraint terms are factorized out of the channels. The
  /// simple Gaussian and Poisson ones (what text2workspace.py writes) are summed
  /// here in closed form, see addFastConstraint(); the others are evaluated by
  /// a never-masked likelihood: a RooUniform of a dummy observable with a
  /// one-entry dataset (exactly zero) and the constraint pdfs as external
  /// constraints. Like CachingSimNLL, the RooMultiPdf penalty is added per
  /// channel, since createNLL() only applies the correction of the top-level
  /// pdf. With "codegen" the analytic gradients of the terms are summed;
  /// hasGradient() requires all terms to have one.
  ///
  /// A RooFit::Evaluator stores a data token in every non-RooRealVar node, so a
  /// RooCategory shared by several channels (a RooMultiPdf index) can't be in
  /// two evaluators. Such channels get a cloned pdf with private category copies,
  /// synchronized with the originals before every evaluation; getParameters()
  /// reports the originals. The per-channel datasets from splitting the input
  /// are owned here, because the terms keep their data by pointer.
  ///
  /// The channel masks of a RooSimultaneousOpt (text2workspace.py
  /// --channel-masks) are honored like in CachingSimNLL: a channel whose mask
  /// is non-zero is left out of the sum, but not of getParameters().
  ///
  /// setAnalyticBarlowBeeston() works like in CachingSimNLL: the bin-by-bin
  /// parameters of the histogram nodes of the original channel pdfs are set
  /// constant and their closed-form best fit is pushed into them before every
  /// evaluation, which the RooFit likelihoods see like any parameter change,
  /// because they share the parameters with the originals.
  class ChannelSumNLL : public SimNLLBase {
  public:
    /// \param simPdf     The simultaneous pdf, usually a RooSimultaneousOpt
    ///                   whose channel pdfs still contain the constraint terms.
    /// \param data       The combined dataset, containing the index category.
    /// \param constrain  Parameters to constrain (the nuisance parameters), or
    ///                   nullptr for RooFit's automatic selection.
    /// \param evalBackend  Name of the RooFit evaluation backend.
    ChannelSumNLL(const char *name,
                  const char *title,
                  RooSimultaneous &simPdf,
                  RooAbsData &data,
                  const RooArgSet *constrain,
                  std::string evalBackend);
    ChannelSumNLL(const ChannelSumNLL &other, const char *name = nullptr);
    ~ChannelSumNLL() override = default;
    ChannelSumNLL &operator=(const ChannelSumNLL &) = delete;
    ChannelSumNLL(ChannelSumNLL &&) = delete;
    ChannelSumNLL &operator=(ChannelSumNLL &&) = delete;
    // The clone() contract of RooFit is a raw owning pointer.
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    ChannelSumNLL *clone(const char *name = nullptr) const override { return new ChannelSumNLL(*this, name); }

    bool isDerived() const override { return true; }
    double defaultErrorLevel() const override { return 0.5; }
    /// Always evaluate: the terms track their own dirty state, and the private
    /// category copies are synchronized in evaluate().
    double getValV(const RooArgSet *) const override { return evaluate(); }

    /// Replace the data of all channels. The new dataset must have the same
    /// structure as the original one (same columns and index category).
    bool setData(RooAbsData &data, bool cloneData = true) override;
    bool getParameters(const RooArgSet *observables,
                       RooArgSet &outputSet,
                       bool stripDisconnected = true) const override;
    void applyWeightSquared(bool flag) override;

    bool hasGradient() const override { return m_hasGradient; }
    /// Gradient in the order of getParameters(nullptr, set, /*stripDisconnected=*/false),
    /// as RooMinimizer expects it. Masked channels don't contribute.
    void gradient(double *out) const override;

    void setZeroPoint() override;
    void clearZeroPoint() override;
    /// There is no constant zero point in this class, nothing to do.
    void clearConstantZeroPoint() override {}
    void setMaskNonDiscreteChannels(bool mask) override;
    void setAnalyticBarlowBeeston(bool flag) override;

  protected:
    double evaluate() const override;

  private:
    struct CategorySync {
      RooCategory *original = nullptr;
      RooCategory *clone = nullptr;
    };
    // Constraint terms evaluated in closed form, see addFastConstraint().
    struct GaussianConstraint {
      RooRealVar const *x = nullptr;
      RooRealVar const *mean = nullptr;
      RooAbsReal const *sigma = nullptr;
    };
    struct PoissonConstraint {
      RooPoisson const *pdf = nullptr;
      RooAbsReal const *count = nullptr;
      RooRealVar const *mean = nullptr;
      bool normalized = false;  // RooFit would normalize it over the mean
      // The normalization integral only depends on the count and on the range
      // of the mean, cached here.
      mutable bool normValid = false;
      mutable double normCount = 0.;
      mutable double normMin = 0.;
      mutable double normMax = 0.;
      mutable double logNorm = 0.;
    };

    void setup(RooAbsData &data);
    std::vector<std::unique_ptr<RooAbsData>> splitData(RooAbsData &data) const;
    /// The pdf to build the channel likelihood from: the channel pdf itself, or
    /// a clone with private copies of the categories it depends on.
    RooAbsPdf &channelPdfForNLL(RooAbsPdf &pdf, std::vector<CategorySync> &syncs);
    void syncCategories() const;
    void replaceClonesByOriginals(RooArgSet &params) const;
    /// For the "codegen" backend: the channel pdf with a gradient term per
    /// histogram node, see CMSHistBarlowBeestonGradientTerm. Returns the input
    /// pdf if there is no histogram node.
    RooAbsPdf &addBarlowBeestonGradientTerms(RooAbsPdf &pdf, std::vector<RooAbsReal const *> const &hists);
    void updateGradientMaps() const;
    /// Take over a constraint pdf that is evaluated in closed form, with the
    /// same value as its RooFit likelihood would have: a SimpleGaussianConstraint
    /// of a parameter with a variable as mean, or a RooPoisson of a constant
    /// count with a parameter as mean. Returns false for anything else.
    bool addFastConstraint(RooAbsPdf const &pdf);
    double fastConstraintsValue() const;
    static double poissonConstraintValue(PoissonConstraint const &c);
    void fastConstraintsGradient(double *out) const;
    /// Masked by the discrete profiling scan (setMaskNonDiscreteChannels()) or
    /// by the channel mask of the RooSimultaneousOpt.
    bool isMasked(std::size_t i) const {
      return i < m_nChannels && ((!m_internalMasks.empty() && !m_internalMasks[i]) ||
                                 (m_channelMaskVars[i] != nullptr && m_channelMaskVars[i]->getVal() != 0.));
    }
    /// Hand the observed counts of a channel to its histogram nodes, which the
    /// analytic Barlow-Beeston minimization needs.
    void setHistogramData(std::size_t iChannel, RooAbsData const &channelData) const;
    /// The Barlow-Beeston step before an evaluation, for the active channels.
    void runAnalyticBarlowBeeston() const;

    RooSimultaneous *m_simPdf = nullptr;  // not owned
    RooAbsData *m_data = nullptr;         // not owned, kept for cloning
    bool m_hasConstrain = false;
    RooArgSet m_constrain;
    std::string m_evalBackend;

    std::unique_ptr<RooAbsPdf> m_factorizedPdfOwner;  // set if the factorization created a new pdf
    RooSimultaneous *m_factorizedPdf = nullptr;       // channel pdfs without constraint terms
    RooArgList m_constraints;

    std::vector<std::string> m_channelNames;
    std::vector<RooAbsReal const *> m_channelMaskVars;  // the channel masks of a RooSimultaneousOpt, or nullptr
    std::vector<std::unique_ptr<RooAbsPdf>> m_privateChannelPdfs;
    std::vector<std::unique_ptr<RooCategory>> m_categoryClones;
    std::vector<std::vector<CategorySync>> m_categorySyncs;
    RooSetProxy m_originalCategories;  // the original categories as servers, for graph traversals
    std::vector<std::unique_ptr<RooAbsData>> m_datasets;
    std::unique_ptr<RooRealVar> m_constraintDummyObs;
    std::unique_ptr<RooAbsPdf> m_constraintPdf;
    std::unique_ptr<RooDataSet> m_constraintData;

    // The likelihood terms: one per channel, then the one of the constraint
    // terms if there are any. Masks and penalties apply to the channels only.
    std::vector<std::unique_ptr<RooAbsReal>> m_terms;
    std::size_t m_nChannels = 0;
    std::vector<RooArgSet> m_channelParams;                            // all parameters of each channel likelihood
    std::vector<std::vector<RooMultiPdf const *>> m_channelMultiPdfs;  // for the penalty terms
    // The histogram nodes (CMSHistSum or CMSHistErrorPropagator) of the
    // original channel pdfs, for the analytic Barlow-Beeston minimization.
    std::vector<std::vector<RooAbsReal const *>> m_channelHists;
    bool m_analyticBarlowBeeston = false;
    // With an analytic gradient, every histogram node gets a gradient term
    // (see CMSHistBarlowBeestonGradientTerm) inside the channel pdf, fed from
    // the original node before each gradient call.
    std::vector<std::vector<CMSHistBarlowBeestonGradientTerm *>> m_channelGradientTerms;  // owned by their pdfs
    RooArgSet m_gradientTermVars;  // their inputs, hidden from getParameters()
    std::vector<double> m_zeroPoints;

    // The constraint terms evaluated in closed form and their zero point.
    std::vector<RooAbsPdf const *> m_fastConstraintPdfs;
    std::vector<GaussianConstraint> m_gaussianConstraints;
    std::vector<PoissonConstraint> m_poissonConstraints;
    double m_fastConstraintZeroPoint = 0;

    // All terms as servers, for dirty-state propagation and graph traversal.
    RooListProxy m_components;

    bool m_hasGradient = false;

    // Channel masking, see setMaskNonDiscreteChannels(). Empty if no mask is
    // set, otherwise true for the channels that are kept.
    std::vector<bool> m_internalMasks;
    RooArgSet m_activeParameters;
    double m_maskingOffset = 0;      // keeps the value continuous when channels are masked
    double m_maskingOffsetZero = 0;  // and the corresponding zero point

    // Maps from the parameter index of each term to the index in our own
    // parameter list, -1 if the parameter is not in our list (masked).
    mutable bool m_gradientMapsDirty = true;
    mutable std::size_t m_nGradientParams = 0;
    mutable std::vector<std::vector<int>> m_gradientMaps;
    // Positions of the closed-form constraint variables in our parameter list,
    // -1 if absent: x, mean and sigma of each Gaussian, the mean of each Poisson.
    mutable std::vector<std::array<int, 3>> m_gaussianGradientIndices;
    mutable std::vector<int> m_poissonGradientIndices;
    mutable std::vector<double> m_gradientBuffer;
  };

}  // namespace cacheutils

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#endif
