#ifndef CMSInterferenceFunc_h
#define CMSInterferenceFunc_h
#include <vector>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "SimpleCacheSentry.h"

class _InterferenceEval;

class CMSInterferenceFunc : public RooAbsReal {
  public:
    CMSInterferenceFunc();
    CMSInterferenceFunc(CMSInterferenceFunc const& other, const char* name = 0);
    /*
     * For a coefficients list of length n and edges array of length b+1,
     * the binscaling nested vector should have b entries (for b bins) with
     * each being of length n*(n+1)/2, corresponding to the upper triangular
     * elements of the scaling matrix, i.e. (m_00, m_01, m_02, ..., m_11, m_12, ...)
     */
    CMSInterferenceFunc(
        const char* name,
        const char* title,
        RooRealVar& x,
        const RooArgList& coefficients,
        const std::vector<double>& edges,
        const std::vector<std::vector<double>> binscaling
        );
    virtual ~CMSInterferenceFunc();

    virtual TObject* clone(const char* newname) const override {
      return new CMSInterferenceFunc(*this, newname);
    };

    void printMultiline(
        std::ostream& os, Int_t contents, Bool_t verbose, TString indent
        ) const override;

    // batch accessor for CMSHistFunc / CMSHistSum
    bool hasChanged() const { return !sentry_.good(); };
    const std::vector<double>& batchGetBinValues() const;

  protected:
    RooRealProxy x_;
    RooListProxy coefficients_;
    std::vector<double> edges_;
    std::vector<std::vector<double>> binscaling_;

    mutable SimpleCacheSentry sentry_; //!
    mutable std::unique_ptr<_InterferenceEval> evaluator_; //!

    double evaluate() const override;

  private:
    void initialize() const;
    void updateCache() const;

    ClassDefOverride(CMSInterferenceFunc, 1)
};

#endif // CMSInterferenceFunc_h
