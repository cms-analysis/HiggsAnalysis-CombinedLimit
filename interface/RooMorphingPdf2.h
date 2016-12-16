#ifndef RooMorphingPdf2_h
#define RooMorphingPdf2_h
#include <map>
#include <vector>
#include <ostream>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"
#include "Rtypes.h"
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logging.h"

struct CMSVMorphPersistCache {
  FastHisto sum;
  FastHisto diff;
};

struct HMorphCache {
  Int_t nbn;
  Double_t xminn;
  Double_t xmaxn;
  std::vector<Double_t> sigdis1;
  std::vector<Double_t> sigdis2;
  std::vector<Double_t> sigdisn;
  std::vector<Double_t> xdisn;
  std::vector<Double_t> sigdisf;
};

struct HMorphNewGlobalCache {
  std::vector<double> bedgesn;
};

struct HMorphNewCDFCache {
  double integral;
  FastTemplate cdf;
  // The interpolation points between this cache and another one
  std::vector<double> x1;
  std::vector<double> x2;
  std::vector<double> y;
};


class CMSHistFunc : public RooAbsReal {
 public:
  CMSHistFunc();

  CMSHistFunc(const char* name, const char* title, RooRealVar & x,
              TH1 const& hist);

  CMSHistFunc(CMSHistFunc const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistFunc(*this, newname);
  }
  virtual ~CMSHistFunc() {}

  Double_t evaluate() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  void addHorizontalMorph(RooAbsReal & hvar, TVectorD hpoints);
  void setVerticalMorphs(RooArgList const& vvars);

  void prepareStorage();

  void setShape(unsigned hindex, unsigned hpoint, unsigned vindex, unsigned vpoint, TH1 const& hist);


/*

– RooAbsArg::setVerboseEval(Int_t level) • Level 0 – No messages
 Level 1 – Print one-line message each time a normalization integral is recalculated
 Level 2 – Print one-line message each time a PDF is recalculated
 Level 3 – Provide details of convolution integral recalculations
*/
 protected:
  RooRealProxy x_;
  RooListProxy vmorphs_;
  RooListProxy hmorphs_;
  mutable FastHisto cache_;
  std::vector<FastHisto> storage_;
  mutable std::vector<CMSVMorphPersistCache> vmorph_cache_; // !not to be serialized
  unsigned n_h_;
  unsigned n_v_;
  mutable HMorphCache mc_; //! not to be serialized
  mutable HMorphNewGlobalCache hmorph_global_cache_; //! not to be serialized
  mutable std::vector<HMorphNewCDFCache> hmorph_cache_; // !not to be serialized

 private:
  ClassDef(CMSHistFunc, 1)

  unsigned getIdx(unsigned hindex, unsigned hpoint, unsigned vindex, unsigned vpoint) const;

  inline double smoothStepFunc(double x) const {
    double _smoothRegion = 1.0;
    if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }

  mutable SimpleCacheSentry vmorph_sentry_; //! not to be serialized
  mutable SimpleCacheSentry hmorph_sentry_; //! not to be serialized

  FastTemplate morph(FastTemplate const& hist1, FastTemplate const& hist2,
                     double par1, double par2, double parinterp) const;

  FastTemplate morph2(FastTemplate const& hist1, FastTemplate const& hist2,
                     double par1, double par2, double parinterp) const;

  void prepareInterpCaches();
  void prepareInterpCache(HMorphNewCDFCache &c1, HMorphNewCDFCache const&c2);
  std::vector<std::vector<double>> hpoints_;
};

// class RooMorphingPdf : public RooAbsPdf {
//  protected:

//   // Store morphing parameters and allocate arrays
//   // so we don't have to do it every time in ::morph
//   struct MorphCache {
//     Int_t nbn;
//     Double_t xminn;
//     Double_t xmaxn;
//     std::vector<Double_t> sigdis1;
//     std::vector<Double_t> sigdis2;
//     std::vector<Double_t> sigdisn;
//     std::vector<Double_t> xdisn;
//     std::vector<Double_t> sigdisf;
//   };

//   mutable MorphCache mc_; //! not to be serialized


//   RooRealProxy x_;                // The x-axis variable
//   RooRealProxy mh_;               // The mass variable
//   RooListProxy pdfs_;             // pdfs
//   std::vector<double> masses_;    // mass points
//   mutable SimpleCacheSentry sentry_; //! not to be serialized
//   bool can_morph_;                // Allowed to do horizontal morphing

//   TArrayI rebin_;                 // Rebinning scheme

//   TAxis target_axis_;             // Target axis
//   TAxis morph_axis_;              // Morphing axis


//   typedef std::map<double, FastVerticalInterpHistPdf2*> MassMap;
//   typedef MassMap::const_iterator MassMapIter;
//   mutable MassMap hmap_;          //! not to be serialized
//   mutable bool init_;             //! not to be serialized
//   mutable FastHisto cache_;       //! not to be serialized

//   // Store some transient info about the current point
//   // Only need to refresh this when the mh parameter changes
//   mutable bool single_point_;  //! not to be serialized
//   mutable FastVerticalInterpHistPdf2 const* p1_;  //! not to be serialized
//   mutable FastVerticalInterpHistPdf2 const* p2_;  //! not to be serialized
//   mutable double mh_lo_;  //! not to be serialized
//   mutable double mh_hi_;  //! not to be serialized

//   void SetAxisInfo();
//   void Init() const;

//   FastTemplate morph(FastTemplate const& hist1, FastTemplate const& hist2,
//                      double par1, double par2, double parinterp) const;

//  public:
//   // Default constructor
//   RooMorphingPdf();
//   // Standard constructor
//   RooMorphingPdf(const char* name, const char* title, RooRealVar& x,
//                  RooAbsReal& mh, RooArgList const& pdfs,
//                  std::vector<double> const& masses, bool const& can_morph,
//                  TAxis const& target_axis, TAxis const& morph_axis);
//   // Copy constructor
//   RooMorphingPdf(const RooMorphingPdf& other, const char* name = 0);
//   // Destructor
//   virtual ~RooMorphingPdf() {}
//   // Clone method
//   virtual TObject* clone(const char* newname) const {
//     return new RooMorphingPdf(*this, newname);
//   }
//   // void AddPoint(double point, FastVerticalInterpHistPdf & hist);

//   Bool_t selfNormalized() const { return kTRUE; }

//   virtual Double_t evaluate() const;

//  public:
//   ClassDef(RooMorphingPdf, 1);
// };

#endif
