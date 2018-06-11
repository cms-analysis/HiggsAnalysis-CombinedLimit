#ifndef VVHZZ4L_ROOSPINZEROPDF_1D_FAST
#define VVHZZ4L_ROOSPINZEROPDF_1D_FAST

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"


class VVHZZ4L_RooSpinZeroPdf_1D_fast : public RooAbsPdf{
protected:
  RooRealProxy fai1;

  RooListProxy obsList; //  List of observables
  RooListProxy coefList; //  List of pdf components

public:
  VVHZZ4L_RooSpinZeroPdf_1D_fast();
  VVHZZ4L_RooSpinZeroPdf_1D_fast(
    const char *name, const char *title,
    RooAbsReal& in_fai1,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );

  VVHZZ4L_RooSpinZeroPdf_1D_fast(const VVHZZ4L_RooSpinZeroPdf_1D_fast& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new VVHZZ4L_RooSpinZeroPdf_1D_fast(*this, newname); }
  inline virtual ~VVHZZ4L_RooSpinZeroPdf_1D_fast(){}

  Float_t interpolateFcn(Int_t code, const char* rangeName=0) const;
  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:
  ClassDef(VVHZZ4L_RooSpinZeroPdf_1D_fast, 1)

};
 
#endif
