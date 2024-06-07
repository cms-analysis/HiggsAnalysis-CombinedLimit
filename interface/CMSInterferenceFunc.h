#ifndef CMSInterferenceFunc_h
#define CMSInterferenceFunc_h
#include "RooListProxy.h"
#include "SimpleCacheSentry.h"

#include "CMSExternalMorph.h"

class _InterferenceEval;

class CMSInterferenceFunc : public CMSExternalMorph {
  public:
    CMSInterferenceFunc();
    CMSInterferenceFunc(CMSInterferenceFunc const& other, const char* name = 0);
    /*
     * For a coefficients list of length n and edges array of length b+1,
     * the binscaling nested vector should have b entries (for b bins) with
     * each being of length n*(n+1)/2, corresponding to the lower triangular
     * elements of the scaling matrix, i.e. (m_00, m_10, m_11, m_20, m_21, m_22, ...)
     */
    CMSInterferenceFunc(
        const char* name,
        const char* title,
        RooRealVar& x,
        const std::vector<double>& edges,
        const RooArgList& coefficients,
        const std::vector<std::vector<double>> binscaling
        );
    ~CMSInterferenceFunc() override;

    TObject* clone(const char* newname) const override {
      return new CMSInterferenceFunc(*this, newname);
    };

    void printMultiline(
        std::ostream& os, Int_t contents, Bool_t verbose, TString indent
        ) const override;

    bool hasChanged() const override { return !sentry_.good(); };
    const std::vector<double>& batchGetBinValues() const override;

  protected:
    RooListProxy coefficients_;
    std::vector<std::vector<double>> binscaling_;

    mutable SimpleCacheSentry sentry_; //!
    mutable std::unique_ptr<_InterferenceEval> evaluator_; //!

  private:
    void initialize() const;
    void updateCache() const;

    ClassDefOverride(CMSInterferenceFunc, 1)
};

#endif // CMSInterferenceFunc_h
