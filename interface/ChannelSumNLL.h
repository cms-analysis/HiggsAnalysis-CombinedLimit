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

class RooMultiPdf;

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
  /// can't do. The constraint terms are factorized out of the channels and
  /// evaluated by a never-masked likelihood: a RooUniform of a dummy observable
  /// with a one-entry dataset (exactly zero) and the constraint pdfs as external
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

  protected:
    double evaluate() const override;

  private:
    struct CategorySync {
      RooCategory *original = nullptr;
      RooCategory *clone = nullptr;
    };

    void setup(RooAbsData &data);
    std::vector<std::unique_ptr<RooAbsData>> splitData(RooAbsData &data) const;
    /// The pdf to build the channel likelihood from: the channel pdf itself, or
    /// a clone with private copies of the categories it depends on.
    RooAbsPdf &channelPdfForNLL(RooAbsPdf &pdf, std::vector<CategorySync> &syncs);
    void syncCategories() const;
    void replaceClonesByOriginals(RooArgSet &params) const;
    void updateGradientMaps() const;
    bool isMasked(std::size_t i) const {
      return i < m_nChannels && !m_internalMasks.empty() && !m_internalMasks[i];
    }

    RooSimultaneous *m_simPdf = nullptr;  // not owned
    RooAbsData *m_data = nullptr;         // not owned, kept for cloning
    bool m_hasConstrain = false;
    RooArgSet m_constrain;
    std::string m_evalBackend;

    std::unique_ptr<RooAbsPdf> m_factorizedPdfOwner;  // set if the factorization created a new pdf
    RooSimultaneous *m_factorizedPdf = nullptr;       // channel pdfs without constraint terms
    RooArgList m_constraints;

    std::vector<std::string> m_channelNames;
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
    std::vector<double> m_zeroPoints;

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
    mutable std::vector<double> m_gradientBuffer;
  };

}  // namespace cacheutils

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#endif
