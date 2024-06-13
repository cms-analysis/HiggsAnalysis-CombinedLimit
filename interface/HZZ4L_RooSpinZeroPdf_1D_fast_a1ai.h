#ifndef HZZ4L_ROOSPINZEROPDF_1D_FAST_A1AI
#define HZZ4L_ROOSPINZEROPDF_1D_FAST_A1AI

#include "HZZ4L_RooSpinZeroPdf_1D_fast_base.h"

#include "RooRealProxy.h"

class HZZ4L_RooSpinZeroPdf_1D_fast_a1ai : public HZZ4L_RooSpinZeroPdf_1D_fast_base{
protected:
  RooRealProxy a1, ai1;
  virtual double a1Val() const;
  virtual double ai1Val() const;
  virtual bool isAnomalousCouplingValueValid() const;
public:
  HZZ4L_RooSpinZeroPdf_1D_fast_a1ai();
  HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(
    const char *name, const char *title,
    RooAbsReal& in_a1,
    RooAbsReal& in_ai1,
    const RooArgList& inObsList,
    const RooArgList& inCoefList
    );
  HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(const HZZ4L_RooSpinZeroPdf_1D_fast_a1ai& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(*this, newname); }

protected:
  ClassDef(HZZ4L_RooSpinZeroPdf_1D_fast_a1ai, 1)
};

#endif
