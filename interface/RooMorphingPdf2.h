#ifndef RooMorphingPdf2_h
#define RooMorphingPdf2_h
#include <map>
#include <ostream>
#include <vector>
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logging.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooArgProxy.h"
#include "Rtypes.h"
#include "TH1F.h"

class CMSHistFunc : public RooAbsReal {
 private:
  struct GlobalCache {
    std::vector<double> bedgesn;
    bool single_point = true;
    unsigned p1 = 0;
    unsigned p2 = 0;
  };

  struct Cache {
    FastTemplate cdf;
    double integral = 0.;

    // The interpolation points between this cache and another one
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> y;

    FastTemplate sum;
    FastTemplate diff;

    FastTemplate step1;
    FastTemplate step2;

    bool cdf_set = false;
    bool interp_set = false;
    // bool step1_set = false;
    // bool step2_set = false;
  };

 public:
  CMSHistFunc();

  CMSHistFunc(const char* name, const char* title, RooRealVar& x,
              TH1 const& hist);

  CMSHistFunc(CMSHistFunc const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistFunc(*this, newname);
  }
  virtual ~CMSHistFunc() {}

  void addHorizontalMorph(RooAbsReal& hvar, TVectorD hpoints);

  void setVerticalMorphs(RooArgList const& vvars);

  void prepareStorage();

  void setShape(unsigned hindex, unsigned hpoint, unsigned vindex,
                unsigned vpoint, TH1 const& hist);

  Double_t evaluate() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  inline void setMorphStrategy(unsigned strategy) {
    morph_strategy_ = strategy;
  };
  inline void setEvalVerbose(unsigned val) { veval = val; };

  inline FastTemplate const& getBinErrors() const { return binerrors_; }
  inline FastHisto const& getCacheHisto() const { return cache_; }

  /*

  – RooAbsArg::setVerboseEval(Int_t level) • Level 0 – No messages
   Level 1 – Print one-line message each time a normalization integral is
  recalculated
   Level 2 – Print one-line message each time a PDF is recalculated
   Level 3 – Provide details of convolution integral recalculations
  */
 protected:
  RooRealProxy x_;
  RooListProxy vmorphs_;
  RooListProxy hmorphs_;
  std::vector<std::vector<double>> hpoints_;
  mutable SimpleCacheSentry vmorph_sentry_;  //! not to be serialized
  mutable SimpleCacheSentry hmorph_sentry_;  //! not to be serialized

  mutable FastHisto cache_;
  FastTemplate binerrors_;
  std::vector<FastHisto> storage_;

  mutable GlobalCache global_;         //! not to be serialized
  mutable std::vector<Cache> mcache_;  //! not to be serialized

  unsigned morph_strategy_;
  int veval;

 private:
  unsigned getIdx(unsigned hindex, unsigned hpoint, unsigned vindex,
                  unsigned vpoint) const;

  inline double smoothStepFunc(double x) const;

  void setCdf(Cache& c, FastHisto const& h) const;

  void prepareInterpCache(Cache& c1, Cache const& c2) const;

  FastTemplate cdfMorph(unsigned idx, double par1, double par2,
                        double parinterp) const;

  ClassDef(CMSHistFunc, 1)
};

class CMSHistErrorPropagator : public RooAbsReal {
public:
  CMSHistErrorPropagator();
  CMSHistErrorPropagator(const char* name, const char* title, RooArgSet const& obs,
                         RooArgList const& funcs, RooArgList const& coeffs,
                         RooArgList const& binpars);

  CMSHistErrorPropagator(CMSHistErrorPropagator const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistErrorPropagator(*this, newname);
  }

  virtual ~CMSHistErrorPropagator() {;}


  // std::vector<double> getErrorMultipliers(unsigned idx);
  // std::vector<double> const& getErrorAdditions(unsigned idx);
  void applyErrorShifts(unsigned idx, FastHisto const& nominal, FastHisto & result);

  Double_t evaluate() const { return 0.; }


 protected:
  RooListProxy funcs_;
  RooListProxy coeffs_;
  RooListProxy binpars_;
  std::vector<CMSHistFunc const*> vfuncs_; //!
  std::vector<RooAbsReal const*> vcoeffs_; //!
  std::vector<RooAbsReal const*> vbinpars_; //!

  std::vector<double> coeffvals_; //!
  std::vector<std::vector<double>> valvec_; //!
  std::vector<std::vector<double>> err2vec_; //!
  std::vector<double> valsum_; //!
  std::vector<double> err2sum_; //!
  std::vector<double> toterr_; //!
  std::vector<std::vector<double>> binmods_; //!
  std::vector<std::vector<double>> scaledbinmods_; //!
  mutable SimpleCacheSentry sentry_;
  mutable SimpleCacheSentry binsentry_;

  int v;

  void fillSumAndErr();


 private:
  ClassDef(CMSHistErrorPropagator,1)
};

class CMSHistFuncWrapper : public RooAbsReal {
 public:
  CMSHistFuncWrapper();

  CMSHistFuncWrapper(const char* name, const char* title, RooRealVar& x,
              CMSHistFunc & func, CMSHistErrorPropagator & err, unsigned idx);

  CMSHistFuncWrapper(CMSHistFuncWrapper const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistFuncWrapper(*this, newname);
  }
  virtual ~CMSHistFuncWrapper() {}

  Double_t evaluate() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  inline void setEvalVerbose(unsigned val) { veval = val; };

 protected:
  RooRealProxy x_;
  RooRealProxy func_;
  RooArgProxy err_;

  mutable FastHisto cache_;
  int veval;
  unsigned idx_;
  mutable SimpleCacheSentry sentry_;

 private:
  ClassDef(CMSHistFuncWrapper, 1)
  mutable CMSHistFunc const* pfunc_;
  mutable CMSHistErrorPropagator * perr_;
  int v;
};

#endif
