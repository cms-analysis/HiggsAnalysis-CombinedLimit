#include "../interface/AsymQuad.h"

#include <cmath>
#include <cassert>
#include <cstdio>

AsymQuad::AsymQuad() :
RooAbsReal(),
_funcList("funcList", "List of functions", this),
_coefList("coefList", "List of coefficients", this),
smoothRegion_(0),
smoothAlgo_(0)
{
  _funcIter  = _funcList.createIterator();
  _coefIter  = _coefList.createIterator();
}

AsymQuad::AsymQuad(const char *name, const char *title, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t smoothRegion, Int_t smoothAlgo) :
RooAbsReal(name, title),
_funcList("funcList", "List of functions", this),
_coefList("coefList", "List of coefficients", this),
smoothRegion_(smoothRegion),
smoothAlgo_(smoothAlgo)
{
  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "AsymQuad::AsymQuad(" << GetName()
      << ") number of functions and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl;
    assert(0);
  }

  TIterator* funcIter = inFuncList.createIterator();
  RooAbsArg* func;
  while ((func = (RooAbsArg*)funcIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: AsymQuad::AsymQuad(" << GetName() << ") function  " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _funcList.add(*func);
  }
  delete funcIter;

  TIterator* coefIter = inCoefList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR: AsymQuad::AsymQuad(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*coef);
  }
  delete coefIter;

  _funcIter  = _funcList.createIterator();
  _coefIter = _coefList.createIterator();
}

AsymQuad::AsymQuad(const AsymQuad& other, const char* name):
RooAbsReal(other, name),
_funcList("!funcList", this, other._funcList),
_coefList("!coefList", this, other._coefList),
smoothRegion_(other.smoothRegion_),
smoothAlgo_(other.smoothAlgo_)
{
  _funcIter  = _funcList.createIterator();
  _coefIter = _coefList.createIterator();
}

AsymQuad::~AsymQuad() {
  delete _funcIter;
  delete _coefIter;
}

Double_t AsymQuad::evaluate() const {
  Double_t result(0);

  _funcIter->Reset();
  _coefIter->Reset();
  RooAbsReal* coef;
  RooAbsReal* func = (RooAbsReal*)_funcIter->Next();

  Double_t central = func->getVal();
  result = central;

  while ((coef=(RooAbsReal*)_coefIter->Next())) {
    Double_t coefVal = coef->getVal();
    RooAbsReal* funcUp = (RooAbsReal*)_funcIter->Next();
    RooAbsReal* funcDn = (RooAbsReal*)_funcIter->Next();
    result += interpolate(coefVal, central, funcUp->getVal(), funcDn->getVal());
  }

  return result;
}

Double_t AsymQuad::interpolate(Double_t theta_, Double_t valueCenter_, Double_t valueHigh_, Double_t valueLow_) const {
  if (smoothAlgo_<0) return 0;
  else{
    if (fabs(theta_)>=smoothRegion_) return theta_ * (theta_ > 0 ? valueHigh_ - valueCenter_ : valueCenter_ - valueLow_);
    if (smoothAlgo_ == 0) {
      // Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries
      Double_t c_up  = +theta_ * (smoothRegion_ + theta_) / (2 * smoothRegion_);
      Double_t c_dn  = -theta_ * (smoothRegion_ - theta_) / (2 * smoothRegion_);
      Double_t c_cen = -theta_ * theta_ / smoothRegion_;
      return c_up * valueHigh_ + c_dn * valueLow_ + c_cen * valueCenter_;
    }
    else if (smoothAlgo_ == 1){
      // Quadratic interpolation that is everywhere differentiable but not null at zero
      Double_t c_up  = (smoothRegion_ + theta_) * (smoothRegion_ + theta_) / (4 * smoothRegion_);
      Double_t c_dn  = (smoothRegion_ - theta_) * (smoothRegion_ - theta_) / (4 * smoothRegion_);
      Double_t c_cen = -c_up - c_dn;
      return c_up * valueHigh_ + c_dn * valueLow_ + c_cen * valueCenter_;
    }
    else/* if (smoothAlgo_ == 2)*/{
      // Quadratic interpolation that is everywhere differentiable and null at zero
      Double_t cnorm = theta_/smoothRegion_;
      Double_t cnorm2 = pow(cnorm, 2);
      Double_t hi = valueHigh_ - valueCenter_;
      Double_t lo = valueLow_ - valueCenter_;
      Double_t sum = hi+lo;
      Double_t diff = hi-lo;
      Double_t a = theta_/2.; // cnorm*smoothRegion_
      Double_t b = 0.125 * cnorm * (cnorm2 * (3.*cnorm2 - 10.) + 15.);
      Double_t result = a*(diff + b*sum);
      return result;
    }
  }
}

ClassImp(AsymQuad)
