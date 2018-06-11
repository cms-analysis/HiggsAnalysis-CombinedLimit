#ifndef RooMorphingPdf_h
#define RooMorphingPdf_h
#include <map>
#include <vector>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"
#include "Rtypes.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"

class RooMorphingPdf : public RooAbsPdf {
 protected:

  // Store morphing parameters and allocate arrays
  // so we don't have to do it every time in ::morph
  struct MorphCache {
    Int_t nbn;
    Double_t xminn;
    Double_t xmaxn;
    std::vector<Double_t> sigdis1;
    std::vector<Double_t> sigdis2;
    std::vector<Double_t> sigdisn;
    std::vector<Double_t> xdisn;
    std::vector<Double_t> sigdisf;
  };

  mutable MorphCache mc_; //! not to be serialized


  RooRealProxy x_;                // The x-axis variable
  RooRealProxy mh_;               // The mass variable
  RooListProxy pdfs_;             // pdfs
  std::vector<double> masses_;    // mass points
  mutable SimpleCacheSentry sentry_; //! not to be serialized
  bool can_morph_;                // Allowed to do horizontal morphing

  TArrayI rebin_;                 // Rebinning scheme

  TAxis target_axis_;             // Target axis
  TAxis morph_axis_;              // Morphing axis


  typedef std::map<double, FastVerticalInterpHistPdf2*> MassMap;
  typedef MassMap::const_iterator MassMapIter;
  mutable MassMap hmap_;          //! not to be serialized
  mutable bool init_;             //! not to be serialized
  mutable FastHisto cache_;       //! not to be serialized

  // Store some transient info about the current point
  // Only need to refresh this when the mh parameter changes
  mutable bool single_point_;  //! not to be serialized
  mutable FastVerticalInterpHistPdf2 const* p1_;  //! not to be serialized
  mutable FastVerticalInterpHistPdf2 const* p2_;  //! not to be serialized
  mutable double mh_lo_;  //! not to be serialized
  mutable double mh_hi_;  //! not to be serialized

  void SetAxisInfo();
  void Init() const;

  FastTemplate morph(FastTemplate const& hist1, FastTemplate const& hist2,
                     double par1, double par2, double parinterp) const;

 public:
  // Default constructor
  RooMorphingPdf();
  // Standard constructor
  RooMorphingPdf(const char* name, const char* title, RooRealVar& x,
                 RooAbsReal& mh, RooArgList const& pdfs,
                 std::vector<double> const& masses, bool const& can_morph,
                 TAxis const& target_axis, TAxis const& morph_axis);
  // Copy constructor
  RooMorphingPdf(const RooMorphingPdf& other, const char* name = 0);
  // Destructor
  virtual ~RooMorphingPdf() {}
  // Clone method
  virtual TObject* clone(const char* newname) const {
    return new RooMorphingPdf(*this, newname);
  }
  // void AddPoint(double point, FastVerticalInterpHistPdf & hist);

  Bool_t selfNormalized() const { return kTRUE; }

  virtual Double_t evaluate() const;

 public:
  ClassDef(RooMorphingPdf, 1);
};

#endif
