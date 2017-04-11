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
#include "RooAbsData.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TMatrix.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistV.h"
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logging.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"

class CMSHistFunc : public RooAbsReal {
 private:
  struct GlobalCache {
    bool single_point = true;
    unsigned p1 = 0;
    unsigned p2 = 0;

    TMatrixD m;
    std::vector<double> c_scale;  // coeffs for scaling
    std::vector<double> c_sum;    // coeffs for summing
    std::vector<double> means;
    std::vector<double> sigmas;
    std::vector<double> slopes;
    std::vector<double> offsets;
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

    // For moment morphing
    double mean;
    double sigma;

    bool meansig_set = false;
  };

 public:
  enum HorizontalType { Closest, Integral, Moment };

  enum MomentSetting {
    Linear,
    NonLinear,
    NonLinearPosFractions,
    NonLinearLinFractions,
    SineLinear
  };

  CMSHistFunc();

  CMSHistFunc(const char* name, const char* title, RooRealVar& x,
              TH1 const& hist, bool divideByWidth = true);

  CMSHistFunc(CMSHistFunc const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistFunc(*this, newname);
  }
  virtual ~CMSHistFunc() {}

  void addHorizontalMorph(RooAbsReal& hvar, TVectorD hpoints);

  void setVerticalMorphs(RooArgList const& vvars);

  void prepareStorage();

  void setActiveBins(unsigned bins);

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
    resetCaches();
  };

  inline void setHorizontalType(CMSHistFunc::HorizontalType htype) {
    htype_ = htype;
    resetCaches();
  };

  inline void setMomentType(CMSHistFunc::MomentSetting mtype) {
    mtype_ = mtype;
    resetCaches();
  };

  inline void setEvalVerbose(unsigned val) { veval = val; };

  inline FastTemplate const& getBinErrors() const { return binerrors_; }
  inline FastHisto const& getCacheHisto() const { return cache_; }

  friend class CMSHistV<CMSHistFunc>;

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
  std::vector<FastTemplate> storage_;

  mutable GlobalCache global_;  //! not to be serialized
  mutable std::vector<Cache> mcache_;  //! not to be serialized

  unsigned morph_strategy_;
  int veval;
  mutable bool initialized_; //! not to be serialized

  HorizontalType htype_;
  MomentSetting mtype_;

  bool divide_by_width_;

 private:
  void initialize() const;
  void setGlobalCache() const;

  void resetCaches();

  unsigned getIdx(unsigned hindex, unsigned hpoint, unsigned vindex,
                  unsigned vpoint) const;

  inline double smoothStepFunc(double x) const;

  void setCdf(Cache& c, FastTemplate const& h) const;
  void setMeanSig(Cache& c, FastTemplate const& h) const;
  void updateMomentFractions(double m) const;

  void prepareInterpCache(Cache& c1, Cache const& c2) const;

  FastTemplate cdfMorph(unsigned idx, double par1, double par2,
                        double parinterp) const;

  double integrateTemplate(FastTemplate const& t) const;

  ClassDef(CMSHistFunc, 1)
};


#endif
