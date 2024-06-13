#ifndef HZZ4L_ROOSPINZEROPDF_1D_FAST
#define HZZ4L_ROOSPINZEROPDF_1D_FAST

#include "HZZ4L_RooSpinZeroPdf_1D_fast_base.h"

#include "RooRealProxy.h"

class HZZ4L_RooSpinZeroPdf_1D_fast : public HZZ4L_RooSpinZeroPdf_1D_fast_base{
protected:
  RooRealProxy fai1;
  virtual double a1Val() const;
  virtual double ai1Val() const;
  virtual bool isAnomalousCouplingValueValid() const;
public:
  HZZ4L_RooSpinZeroPdf_1D_fast();
  HZZ4L_RooSpinZeroPdf_1D_fast(
    const char *name, const char *title,
    RooAbsReal& in_fai1,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );
  HZZ4L_RooSpinZeroPdf_1D_fast(const HZZ4L_RooSpinZeroPdf_1D_fast& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new HZZ4L_RooSpinZeroPdf_1D_fast(*this, newname); }

protected:
  ClassDef(HZZ4L_RooSpinZeroPdf_1D_fast, 1)
};

#endif
