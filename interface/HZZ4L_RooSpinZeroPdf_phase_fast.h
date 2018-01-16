#ifndef HZZ4L_ROOSPINZEROPDF_PHASE_FAST
#define HZZ4L_ROOSPINZEROPDF_PHASE_FAST

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"


class HZZ4L_RooSpinZeroPdf_phase_fast : public RooAbsPdf{
protected:
  RooRealProxy fai1;
  RooRealProxy phi1;

  RooListProxy obsList; //  List of observables
  RooListProxy coefList; //  List of pdf components

public:
  HZZ4L_RooSpinZeroPdf_phase_fast();
  HZZ4L_RooSpinZeroPdf_phase_fast(
    const char *name, const char *title,
    RooAbsReal& in_fai1,
    RooAbsReal& in_phi1,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );

  HZZ4L_RooSpinZeroPdf_phase_fast(const HZZ4L_RooSpinZeroPdf_phase_fast& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new HZZ4L_RooSpinZeroPdf_phase_fast(*this, newname); }
  inline virtual ~HZZ4L_RooSpinZeroPdf_phase_fast(){}

  Float_t interpolateFcn(Int_t code, const char* rangeName=0) const;
  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:
  ClassDef(HZZ4L_RooSpinZeroPdf_phase_fast, 1)

};
 
#endif
