#ifndef VVHZZ4L_ROOSPINZEROPDF_1D_FAST_BASE
#define VVHZZ4L_ROOSPINZEROPDF_1D_FAST_BASE

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooListProxy.h"


class VVHZZ4L_RooSpinZeroPdf_1D_fast_base : public RooAbsPdf{
protected:
  RooListProxy obsList; //  List of observables
  RooListProxy coefList; //  List of pdf components

  virtual double a1Val() const = 0;
  virtual double ai1Val() const = 0;
  virtual bool isAnomalousCouplingValueValid() const = 0;

public:
  VVHZZ4L_RooSpinZeroPdf_1D_fast_base();
  VVHZZ4L_RooSpinZeroPdf_1D_fast_base(
    const char *name, const char *title,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );

  VVHZZ4L_RooSpinZeroPdf_1D_fast_base(const VVHZZ4L_RooSpinZeroPdf_1D_fast_base& other, const char* name=0);
  virtual TObject* clone(const char* newname) const = 0;
  inline virtual ~VVHZZ4L_RooSpinZeroPdf_1D_fast_base(){}

  Float_t interpolateFcn(Int_t code, const char* rangeName=0) const;
  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
};
 
#endif
