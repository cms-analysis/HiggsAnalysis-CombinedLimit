#ifndef HZZ4L_ROOSPINZEROPDF_2D_FAST
#define HZZ4L_ROOSPINZEROPDF_2D_FAST

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"


class HZZ4L_RooSpinZeroPdf_2D_fast : public RooAbsPdf{
protected:
  RooRealProxy fai1;
  RooRealProxy fai2;
  RooRealProxy phi1;
  RooRealProxy phi2;

  RooListProxy obsList; //  List of observables
  RooListProxy coefList; //  List of pdf components

public:
  HZZ4L_RooSpinZeroPdf_2D_fast();
  HZZ4L_RooSpinZeroPdf_2D_fast(
    const char *name, const char *title,
    RooAbsReal& in_fai1,
    RooAbsReal& in_fai2,
    RooAbsReal& in_phi1,
    RooAbsReal& in_phi2,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );

  HZZ4L_RooSpinZeroPdf_2D_fast(const HZZ4L_RooSpinZeroPdf_2D_fast& other, const char* name=0);
  TObject* clone(const char* newname) const override { return new HZZ4L_RooSpinZeroPdf_2D_fast(*this, newname); }
  inline ~HZZ4L_RooSpinZeroPdf_2D_fast() override{}

  Float_t interpolateFcn(Int_t code, const char* rangeName=0) const;
  Double_t evaluate() const override;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

protected:
  ClassDefOverride(HZZ4L_RooSpinZeroPdf_2D_fast, 1)

};
 
#endif
