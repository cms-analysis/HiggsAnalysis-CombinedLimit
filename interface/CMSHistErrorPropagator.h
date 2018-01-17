#ifndef CMSHistErrorPropagator_h
#define CMSHistErrorPropagator_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate_Old.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistV.h"

class CMSHistErrorPropagator : public RooAbsReal {
private:
  struct BarlowBeeston {
    bool init = false;
    std::vector<unsigned> use;
    std::vector<double> dat;
    std::vector<double> valsum;
    std::vector<double> toterr;
    std::vector<double> err;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> tmp;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> res;
    std::vector<double> gobs;
    std::set<RooAbsArg*> dirty_prop;
    std::vector<RooRealVar*> push_res;
  };
public:
  CMSHistErrorPropagator();

  CMSHistErrorPropagator(const char* name, const char* title, RooRealVar& x,
                         RooArgList const& funcs, RooArgList const& coeffs);

  CMSHistErrorPropagator(CMSHistErrorPropagator const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistErrorPropagator(*this, newname);
  }

  virtual ~CMSHistErrorPropagator() {;}

  void applyErrorShifts(unsigned idx, FastHisto const& nominal, FastHisto & result);

  Double_t evaluate() const;

  RooArgList * setupBinPars(double poissonThreshold);

  std::unique_ptr<RooArgSet> getSentryArgs() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  void setData(RooAbsData const& data) const;

  void setAnalyticBarlowBeeston(bool flag) const;

  inline FastHisto const& cache() const { return cache_; }

  RooArgList wrapperList() const;
  RooArgList const& coefList() const { return coeffs_; }
  RooArgList const& funcList() const { return funcs_; }

  friend class CMSHistV<CMSHistErrorPropagator>;

 protected:
  RooRealProxy x_;
  RooListProxy funcs_;
  RooListProxy coeffs_;
  RooListProxy binpars_;
  mutable std::vector<CMSHistFunc const*> vfuncs_; //!
  mutable std::vector<RooAbsReal const*> vcoeffs_; //!
  mutable std::vector<std::vector<RooAbsReal *>> vbinpars_; //!
  std::vector<std::vector<unsigned>> bintypes_;

  mutable std::vector<double> coeffvals_; //!
  mutable FastHisto valsum_; //!
  mutable FastHisto cache_; //!
  mutable std::vector<double> err2sum_; //!
  mutable std::vector<double> toterr_; //!
  mutable std::vector<std::vector<double>> binmods_; //!
  mutable std::vector<std::vector<double>> scaledbinmods_; //!
  mutable SimpleCacheSentry sentry_; //!
  mutable SimpleCacheSentry binsentry_; //!
  mutable std::vector<double> data_; //!

  mutable BarlowBeeston bb_; //!

  mutable bool initialized_; //! not to be serialized

  mutable int last_eval_; //! not to be serialized

  mutable bool analytic_bb_; //! not to be serialized

  void initialize() const;
  void updateCache(int eval = 1) const;

  void runBarlowBeeston() const;


 private:
  ClassDef(CMSHistErrorPropagator,1)
};

#endif
