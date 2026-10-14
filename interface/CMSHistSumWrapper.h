#ifndef CMSHistSumWrapper_h  // NOLINT(llvm-header-guard)
#define CMSHistSumWrapper_h

#include <ostream>

#include "CMSHistSum.h"
#include "FastTemplate_Old.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "RtypesCore.h"
#include "SimpleCacheSentry.h"
#include "TString.h"

// Per-process view of a CMSHistSum: evaluates the morphed template of process
// idx including its share of the bin-wise (autoMCStats) shifts. Plays the role
// that CMSHistFuncWrapper has for CMSHistErrorPropagator. Either stored in the
// workspace by text2workspace or created on the fly by CMSHistSum::wrapperList().
class CMSHistSumWrapper : public RooAbsReal {  // NOLINT(cppcoreguidelines-special-member-functions)
public:
  CMSHistSumWrapper() = default;

  CMSHistSumWrapper(const char* name, const char* title, RooAbsReal& x, CMSHistSum& sum, unsigned idx);

  CMSHistSumWrapper(CMSHistSumWrapper const& other, const char* name = nullptr);

  TObject* clone(const char* newname) const override {
    return new CMSHistSumWrapper(*this, newname);  // NOLINT(cppcoreguidelines-owning-memory)
  }
  ~CMSHistSumWrapper() override = default;

  Double_t evaluate() const override;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose, TString indent) const override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = nullptr) const override;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = nullptr) const override;

  void updateCache() const;

  FastHisto const& cache() const { return m_cache; }

  unsigned idx() const { return m_idx; }

private:
  RooRealProxy m_x;
  RooRealProxy m_sum;

  mutable FastHisto m_cache;
  unsigned m_idx = 0;
  mutable SimpleCacheSentry m_sentry;          //!
  mutable CMSHistSum const* m_psum = nullptr;  //!
  mutable bool m_initialized = false;          //! not to be serialized

  void initialize() const;

  ClassDefOverride(CMSHistSumWrapper, 1)  // NOLINT
};

#endif
