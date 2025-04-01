#ifndef ROO_REAL_SUM_PDF_SAFE
#define ROO_REAL_SUM_PDF_SAFE

/** Vertical interpolation between multiple histograms (or pdfs).
    Based on RooRealSumPdf */

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooObjCacheManager.h"
#include "RtypesCore.h"

#include "CombineCodegenImpl.h"

class VerticalInterpPdf : public RooAbsPdf {
public:

  VerticalInterpPdf() = default;
  VerticalInterpPdf(const char *name, const char *title, const RooArgList& funcList, const RooArgList& coefList, Double_t quadraticRegion=0., Int_t quadraticAlgo=0) ;
  VerticalInterpPdf(const VerticalInterpPdf& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new VerticalInterpPdf(*this,newname) ; }

  Double_t evaluate() const override ;
  Bool_t checkObservables(const RooArgSet* nset) const override ;	
  COMBINE_DECLARE_TRANSLATE

  Bool_t forceAnalyticalInt(const RooAbsArg&) const override { return kTRUE ; }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& numVars, const RooArgSet* normSet, const char* rangeName=0) const override ;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const override ;
  COMBINE_DECLARE_ANALYTICAL_INTEGRAL

  const RooArgList& funcList() const { return _funcList ; }
  const RooArgList& funcIntListFromCache() const;
  const RooArgList& coefList() const { return _coefList ; }

  const Double_t quadraticRegion() const { return _quadraticRegion; }
  const Int_t quadraticAlgo() const { return _quadraticAlgo; }

  Double_t pdfFloorVal() const { return _pdfFloorVal; }
  Double_t integralFloorVal() const { return _integralFloorVal; }

  void setFloorVals(Double_t const& pdf_val, Double_t const& integral_val);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,34,06)
  std::unique_ptr<RooAbsArg> compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext & ctx) const override;
#endif

protected:
  
  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem()  {} ;
    RooArgList containedArgs(Action) override { RooArgList ret(_funcIntList) ; ret.add(_funcNormList) ; return ret ; }
    RooArgList _funcIntList ;
    RooArgList _funcNormList ;
  } ;
  mutable RooObjCacheManager _normIntMgr ; // The integration cache manager

  RooListProxy _funcList ;   //  List of component FUNCs
  RooListProxy _coefList ;  //  List of coefficients
  Double_t     _quadraticRegion = 0;
  Int_t        _quadraticAlgo;

  Double_t _pdfFloorVal = 1e-15; // PDF floor should be customizable, default is 1e-15
  Double_t _integralFloorVal = 1e-10; // PDF integral floor should be customizable, default is 1e-10

  Double_t interpolate(Double_t coeff, Double_t central, RooAbsReal *fUp, RooAbsReal *fDown) const ; 

  bool isConditionalProdPdf(RooAbsReal *pdf) const;
  RooAbsReal* makeConditionalProdPdfIntegral(RooAbsPdf* pdf, RooArgSet const& analVars) const;

  ClassDefOverride(VerticalInterpPdf,3) // PDF constructed from a sum of (non-pdf) functions
};

COMBINE_DECLARE_CODEGEN_IMPL(VerticalInterpPdf);
COMBINE_DECLARE_CODEGEN_INTEGRAL_IMPL(VerticalInterpPdf);

#endif
