#ifndef VBFHZZ4L_ROOSPINZEROPDF_FAST
#define VBFHZZ4L_ROOSPINZEROPDF_FAST

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"


class VBFHZZ4L_RooSpinZeroPdf_fast : public RooAbsPdf{
protected:
  RooRealProxy a1;
  RooRealProxy ai1;

  RooListProxy obsList; //  List of observables
  RooListProxy coefList; //  List of pdf components

public:
  VBFHZZ4L_RooSpinZeroPdf_fast();
  VBFHZZ4L_RooSpinZeroPdf_fast(
    const char *name, const char *title,
    RooAbsReal& in_a1,
    RooAbsReal& in_ai1,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );

  VBFHZZ4L_RooSpinZeroPdf_fast(const VBFHZZ4L_RooSpinZeroPdf_fast& other, const char* name=0);
  TObject* clone(const char* newname) const override { return new VBFHZZ4L_RooSpinZeroPdf_fast(*this, newname); }
  inline ~VBFHZZ4L_RooSpinZeroPdf_fast() override{}

  Float_t interpolateFcn(Int_t code, const char* rangeName=0) const;
  Double_t evaluate() const override;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

protected:
  ClassDefOverride(VBFHZZ4L_RooSpinZeroPdf_fast, 1)

};
 
#endif
