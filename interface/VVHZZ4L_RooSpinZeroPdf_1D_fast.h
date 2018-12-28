#ifndef VVHZZ4L_ROOSPINZEROPDF_1D_FAST
#define VVHZZ4L_ROOSPINZEROPDF_1D_FAST

#include "VVHZZ4L_RooSpinZeroPdf_1D_fast_base.h"

#include "RooRealProxy.h"

class VVHZZ4L_RooSpinZeroPdf_1D_fast : public VVHZZ4L_RooSpinZeroPdf_1D_fast_base{
protected:
  RooRealProxy fai1;
  virtual double a1Val() const;
  virtual double ai1Val() const;
  virtual bool isAnomalousCouplingValueValid() const;
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

protected:
  ClassDef(VVHZZ4L_RooSpinZeroPdf_1D_fast, 1)
};

#endif
