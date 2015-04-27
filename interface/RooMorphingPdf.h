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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "../interface/VerticalInterpHistPdf.h"
#pragma GCC diagnostic pop

class RooMorphingPdf : public RooAbsPdf {
 protected:
  RooRealProxy x_;                // The x-axis variable
  RooRealProxy mh_;               // The mass variable
  RooListProxy pdfs_;             // pdfs
  std::vector<double> masses_;    // mass points
  mutable double current_mh_;     // The last-used value of the mass
  bool can_morph_;                // Allowed to do horizontal morphing

  TArrayI rebin_;                 // Rebinning scheme

  TAxis target_axis_;             // Target axis
  TAxis morph_axis_;              // Morphing axis


  typedef std::map<double, FastVerticalInterpHistPdf2*> MassMap;
  typedef MassMap::const_iterator MassMapIter;
  mutable MassMap hmap_;          //! not to be serialized

  mutable bool init_;             //! not to be serialized
  mutable FastHisto cache_;       //! not to be serialized


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
