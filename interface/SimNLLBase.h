#ifndef HiggsAnalysis_CombinedLimit_SimNLLBase_h
#define HiggsAnalysis_CombinedLimit_SimNLLBase_h

#include <RooAbsReal.h>

namespace cacheutils {

  /// Interface of CachingSimNLL and ChannelSumNLL for the CascadeMinimizer:
  /// the zero-point subtraction that keeps the values seen by Minuit close to
  /// zero, and the channel masking of the discrete profiling scans.
  class SimNLLBase : public RooAbsReal {
  public:
    SimNLLBase() = default;
    SimNLLBase(const char *name, const char *title) : RooAbsReal(name, title) {}
    SimNLLBase(const SimNLLBase &other, const char *name = nullptr) : RooAbsReal(other, name) {}
    ~SimNLLBase() override = default;
    // Like RooAbsArg, copies are made with the copy constructor and clone() only.
    SimNLLBase &operator=(const SimNLLBase &) = delete;
    SimNLLBase(SimNLLBase &&) = delete;
    SimNLLBase &operator=(SimNLLBase &&) = delete;

    /// Subtract the current value of every term from the likelihood, so that
    /// the minimizer works with numbers close to zero.
    virtual void setZeroPoint() = 0;
    /// Undo setZeroPoint().
    virtual void clearZeroPoint() = 0;
    void updateZeroPoint() {
      clearZeroPoint();
      setZeroPoint();
    }
    /// Remove the constant offset applied at construction, so that the value is
    /// the full NLL (needed to compare likelihoods of different models, e.g. in
    /// the goodness-of-fit tests).
    virtual void clearConstantZeroPoint() = 0;

    /// Keep only the channels with a floating RooCategory (the ones in a discrete
    /// profiling scan) in the sum and in getParameters(). An offset keeps the
    /// value continuous, so masked and unmasked values compare. `mask = false`
    /// restores all channels.
    virtual void setMaskNonDiscreteChannels(bool mask) = 0;
  };

}  // namespace cacheutils

#endif
