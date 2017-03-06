#ifndef CMSHistFunc_h
#define CMSHistFunc_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooArgProxy.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logging.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"

class CMSHistFuncV;

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

  void updateCache() const;

  std::unique_ptr<RooArgSet> getSentryArgs() const;

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

  friend class CMSHistFuncV;

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
  mutable bool initialized_; //! not to be serialized

 private:
  void initialize() const;

  unsigned getIdx(unsigned hindex, unsigned hpoint, unsigned vindex,
                  unsigned vpoint) const;

  inline double smoothStepFunc(double x) const;

  void setCdf(Cache& c, FastHisto const& h) const;

  void prepareInterpCache(Cache& c1, Cache const& c2) const;

  FastTemplate cdfMorph(unsigned idx, double par1, double par2,
                        double parinterp) const;

  ClassDef(CMSHistFunc, 1)
};

class CMSHistFuncV {
 public:
  CMSHistFuncV(const CMSHistFunc &,
                              const RooAbsData& data,
                              bool includeZeroWeights = false);
  void fill(std::vector<Double_t>& out) const;

 private:
  const CMSHistFunc& hpdf_;
  int begin_, end_, nbins_;
  struct Block {
    int index, begin, end;
    Block(int i, int begin_, int end_) : index(i), begin(begin_), end(end_) {}
  };
  std::vector<Block> blocks_;
  std::vector<int> bins_;
};

#endif
