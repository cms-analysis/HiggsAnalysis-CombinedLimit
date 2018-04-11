#ifndef ROONCSPLINE_1D_FAST
#define ROONCSPLINE_1D_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooNCSplineCore.h"


class RooNCSpline_1D_fast : public RooNCSplineCore{
protected:
  BoundaryCondition const bcBeginX;
  BoundaryCondition const bcEndX;

  std::vector<T> FcnList; // List of function values

  std::vector<T> kappaX;
  std::vector<std::vector<T>> coefficients;

public:
  RooNCSpline_1D_fast();
  RooNCSpline_1D_fast(
    const char* name,
    const char* title
    );
  RooNCSpline_1D_fast(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    const std::vector<T>& inXList,
    const std::vector<T>& inFcnList,
    RooNCSplineCore::BoundaryCondition const bcBeginX_=RooNCSplineCore::bcNaturalSpline,
    RooNCSplineCore::BoundaryCondition const bcEndX_=RooNCSplineCore::bcNaturalSpline,
    Bool_t inUseFloor=true,
    T inFloorEval=0,
    T inFloorInt=0
    );
  RooNCSpline_1D_fast(const RooNCSpline_1D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new RooNCSpline_1D_fast(*this, newname); }
	inline virtual ~RooNCSpline_1D_fast(){}

  void setRangeValidity(const T valmin, const T valmax, const Int_t whichDirection=0);

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

protected:
  virtual void emptyFcnList(){ std::vector<T> tmp; FcnList.swap(tmp); }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const;
  virtual void getKappas(std::vector<T>& kappas, const Int_t whichDirection);

  Bool_t testRangeValidity(const T& val, const Int_t whichDirection=0)const;
  void cropValueForRange(T& val, const Int_t whichDirection=0)const;

  virtual T interpolateFcn(Int_t code, const char* rangeName=0)const;

  virtual Double_t evaluate()const;


  ClassDef(RooNCSpline_1D_fast, 2)

};
 
#endif
