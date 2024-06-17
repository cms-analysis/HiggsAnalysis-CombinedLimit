#include "../interface/HMuMuRooPdfs.h"
//#include "HMuMuRooPdfs.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "RooMath.h"
#include "Math/SpecFuncMathMore.h"
#include "TError.h"
#include <cmath>
#include <complex>
#include <iostream>

/////////////////////
ClassImp(RooModZPdf)

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a):
RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b{"b", "b", this},
  c{"c", "c", this},
  m{"m", "m", this},
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b("b", "b", this, _b),
  c("c", "c", this, _c),
  m{"m", "m", this},
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _c):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b{"b", "b", this},
  c("c", "c", this, _c),
  m{"m", "m", this},
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, const RooArgList& _coef):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b{"b", "b", this},
  c{"c", "c", this},
  m{"m", "m", this},
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
  for (RooAbsArg* coef : _coef) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << "RooBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl ;
      R__ASSERT(0) ;
    }
    bernCoef.add(*coef);
  }
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b("b", "b", this, _b),
  c("c", "c", this, _c),
  m{"m", "m", this},
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
  for (RooAbsArg* coef : _coef) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << "RooBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl ;
      R__ASSERT(0) ;
    }
    bernCoef.add(*coef);
  }
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b("b", "b", this, _b),
  c("c", "c", this, _c),
  m("m", "m", this, _m),
  w{"w", "w", this},
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w):
  RooAbsPdf(name, title),
  x("x", "x", this, _x),
  a("a", "a", this, _a),
  b("b", "b", this, _b),
  c("c", "c", this, _c),
  m("m", "m", this, _m),
  w("w", "w", this, _w),
  bernCoef("coefficients", "List of Bernstein coefficients", this)
{
}

RooModZPdf::RooModZPdf(const RooModZPdf& other, const char* name):
  RooAbsPdf(other, name),
  x("x", this, other.x),
  a("a", this, other.a),
  b("b", this, other.b),
  c("c", this, other.c),
  m("m", this, other.m),
  w("w", this, other.w),
  bernCoef("coefficients", this, other.bernCoef)
{
}

double RooModZPdf::evaluate() const {

  const double xmin = x.min();
  const double xmax = x.max();

  // default values
  double zm = 91.2;
  double zw = 2.5;
  double bv = 0.;
  double cv = 2.;

  if(&m.arg()) zm = m;
  if(&w.arg()) zw = w;
  if(&b.arg()) bv = b;
  if(&c.arg()) cv = c;

  double val = 1.0;
  val *= exp(a*x+bv*x*x);
  val /= (pow(x*x-zm*zm,cv) + pow(zw*zm,cv));
  
  Int_t degree = bernCoef.getSize();
  if (degree <= 0) return val;
  
  double xv = (x - xmin) / (xmax - xmin);
  double bernval = 1.;
  double coefsum = 0.;
  double coef    = 0.;

  for (Int_t i = 0; i < degree; i++) {
    coef = static_cast<RooAbsReal &>(bernCoef[i]).getVal();
    coefsum -= coef;
    bernval += (degree+1.) * TMath::Binomial(degree, i) * pow(xv, degree-i) * pow(1.-xv, i) * coef;
  }
  bernval += coefsum * pow(1.-xv, degree);  
  return val * bernval;
  
}

/////////
ClassImp(RooExpPdf)

RooExpPdf::RooExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, const bool & _offset):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  m {"m" , "m" , this},
  offset(_offset)
{
}

RooExpPdf::RooExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m):
  RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  m ("m" , "m" , this, _m ),
  offset(true)
{
}

RooExpPdf::RooExpPdf(const RooExpPdf& other, const char* name):
  RooAbsPdf(other, name),
  x ("x" , this, other.x ),
  a1("a1", this, other.a1),
  m ("m" , this, other.m ),
  offset(other.offset)
{
}

Double_t RooExpPdf::evaluate() const {

  const double xmin = x.min();
  const double xmax = x.max();

  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  // Normalization coefficients for each term of the PDF --> this ensures that the overall PDF is normalized correctly
  double retval1 = 0.;
  if(a1 == 0) retval1 = 1/(xmax-xmin);
  else        retval1 = a1/(exp(a1*(xmax-zm))-exp(a1*(xmin-zm)));
  
  return retval1*exp(a1*(x-zm));
  
}

Int_t RooExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const {
  if (matchArgs(allVars, analVars, x)) return 1;
  else return 0;
}

Double_t RooExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

  // Assert if code is not one
  R__ASSERT(code == 1);

  // take the integration range
  const Double_t xmin = x.min(rangeName);
  const Double_t xmax = x.max(rangeName);

  // if integral over full range return 1.
  if (xmin == x.min() && xmax == x.max()) return 1.;

  // correction when PDF integrated in a smaller range
  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  // scale the integral by the relative fractions
  Double_t integral_part1 = 0.;
  if (a1 != 0.) {
    integral_part1 *= (exp(a1*(xmax-zm)) - exp(a1*(xmin-zm)));
    integral_part1 /= (exp(a1*(x.max()-zm)) - exp(a1*(x.min()-zm)));
  }
  else {
    integral_part1 *= (xmax - xmin);
    integral_part1 /= (x.max()-x.min());
  }
  
  return integral_part1;
}

/////////
ClassImp(RooSumTwoExpPdf)

RooSumTwoExpPdf::RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, const bool & _offset):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  a2("a2", "a2", this, _a2),
  f ("f" , "f" , this, _f ),
  m {"m" , "m" , this},
  offset(_offset)
{
}

RooSumTwoExpPdf::RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, RooAbsReal& _m):
  RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  a2("a2", "a2", this, _a2),
  f ("f" , "f" , this, _f ),
  m ("m" , "m" , this, _m ),
  offset(true)
{
}

RooSumTwoExpPdf::RooSumTwoExpPdf(const RooSumTwoExpPdf& other, const char* name):
  RooAbsPdf(other, name),
  x ("x" , this, other.x ),
  a1("a1", this, other.a1),
  a2("a2", this, other.a2),
  f ("f" , this, other.f ),
  m ("m" , this, other.m ),
  offset(other.offset)
{
}

Double_t RooSumTwoExpPdf::evaluate() const {

  const double xmin = x.min();
  const double xmax = x.max();

  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  // Normalization coefficients for each term of the PDF --> this ensures that the overall PDF is normalized correctly
  double retval1 = 0.;
  double retval2 = 0.;
  if(a1 == 0) retval1 = 1/(xmax-xmin);
  else        retval1 = a1/(exp(a1*(xmax-zm))-exp(a1*(xmin-zm)));
  if(a2 == 0) retval2 = 1/(xmax-xmin);
  else        retval2 = a2/(exp(a2*xmax)-exp(a2*xmin));

  return retval1*f*exp(a1*(x-zm))+retval2*(1.-f)*exp(a2*x);     
  
}

Int_t RooSumTwoExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const {
  if (matchArgs(allVars, analVars, x)) return 1;
  else return 0;
}

Double_t RooSumTwoExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

  // Assert if code is not one
  R__ASSERT(code == 1);
  
  // take the integration range
  const Double_t xmin = x.min(rangeName);
  const Double_t xmax = x.max(rangeName);

  // if integral over full range return 1.
  if (xmin == x.min() && xmax == x.max()) return 1.;

  // correction when PDF integrated in a smaller range
  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  // scale the integral by the relative fractions
  Double_t integral_part1 = f;
  if (a1 != 0.) {
    integral_part1 *= (exp(a1*(xmax-zm)) - exp(a1*(xmin-zm)));
    integral_part1 /= (exp(a1*(x.max()-zm)) - exp(a1*(x.min()-zm)));
  }
  else {
    integral_part1 *= (xmax - xmin);
    integral_part1 /= (x.max()-x.min());
  }
  
  // scale the integral by the relative fractions
  Double_t integral_part2 = (1. - f);
  if (a2 != 0.) {
    integral_part2 *= (exp(a2*xmax) - exp(a2*xmin));
    integral_part2 /= (exp(a2*x.max()) - exp(a2*x.min()));
  }
  else {
    integral_part2 *= (xmax - xmin);
    integral_part2 /= (x.max()-x.min());
  }

  return integral_part1 + integral_part2;
}


//////////
ClassImp(RooPowerLawPdf)

RooPowerLawPdf::RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, const bool & _offset):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  m{"m", "m", this},
  offset(_offset)
{}

RooPowerLawPdf::RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  m ("m", "m", this, _m),
  offset(true)
{}

RooPowerLawPdf::RooPowerLawPdf(const RooPowerLawPdf& other, const char* name):
  RooAbsPdf(other, name),
  x ("x" , this, other.x ),
  a1("a1", this, other.a1),
  m ("m", this, other.m),
  offset(other.offset)
{
}

Double_t RooPowerLawPdf::evaluate() const {

  const Double_t xmin = x.min();
  const Double_t xmax = x.max();

  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  double retval = 0.;
  if(a1 == 0.)       retval = 1/(xmax-xmin);
  else if (a1 == -1) retval = 1/(log(xmax-zm)-log(xmin-zm));
  else               retval = (a1+1)/(pow(xmax-zm,a1+1)-pow(xmin-zm,a1+1));

  double val = pow(x-zm,a1)*retval;  
  return val;
}


Int_t RooPowerLawPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const {
  
  if (matchArgs(allVars, analVars, x)) return 1;
  else return 0;
}

Double_t RooPowerLawPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

  // Assert if code is not one
  R__ASSERT(code == 1);

  // take the integration range
  const Double_t xmin = x.min(rangeName);
  const Double_t xmax = x.max(rangeName);

  // if integral over full range return 1.
  if (xmin == x.min() && xmax == x.max()) return 1.;

  // correction when PDF integrated in a smaller range
  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  double integral = 0.;
  if(a1 == 0.)       integral = (xmax-xmin)/(x.max()-x.min());
  else if (a1 == -1) integral = (log(xmax-zm)-log(xmin-zm))/(log(x.max()-zm)-log(x.min()-zm));
  else               integral = (pow(xmax-zm,a1+1)-pow(xmin-zm,a1+1))/(pow(x.max()-zm,a1+1)-pow(x.min()-zm,a1+1));

  return integral;
}



///////////
ClassImp(RooSumTwoPowerLawPdf)

RooSumTwoPowerLawPdf::RooSumTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, const bool & _offset):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  a2("a2", "a2", this, _a2),
  f("f", "f", this, _f),
  m{"m", "m", this},
  offset(_offset)
{
}

 RooSumTwoPowerLawPdf::RooSumTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, RooAbsReal& _m):
RooAbsPdf(name, title),
  x ("x" , "x" , this, _x),
  a1("a1", "a1", this, _a1),
  a2("a2", "a2", this, _a2),
  f("f", "f", this, _f),
  m("m", "m", this, _m),
  offset(true)
{
}

RooSumTwoPowerLawPdf::RooSumTwoPowerLawPdf(const RooSumTwoPowerLawPdf& other, const char* name):
  RooAbsPdf(other, name),
  x ("x" , this, other.x ),
  a1("a1", this, other.a1),
  a2("a2", this, other.a2),
  f("f", this, other.f),
  m("m", this, other.m),
  offset(other.offset)
{
}

Double_t RooSumTwoPowerLawPdf::evaluate() const {

  const Double_t xmin = x.min();
  const Double_t xmax = x.max();

  double zm = 91.2;
  if(&m.arg()) zm = m;
  if(not offset) zm = 0;

  double retval1 = 0.;
  if(a1 == 0)       retval1 = 1/(xmax-xmin);
  else if(a1 == -1) retval1 = 1/(log(xmax-zm)-log(xmin-zm));
  else              retval1 = (a1+1)/(pow(xmax-zm,a1+1)-pow(xmin-zm,a1+1));

  double retval2 = 0.;
  if(a2 == 0)       retval2 = 1/(xmax-xmin);
  else if(a2 == -1) retval2 = 1/(log(xmax)-log(xmin));
  else              retval2 = (a2+1)/(pow(xmax,a2+1)-pow(xmin,a2+1));

  return f*retval1*pow(x-zm,a1)+(1-f)*retval2*pow(x,a2);
  
}

Int_t RooSumTwoPowerLawPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

    if (matchArgs(allVars, analVars, x)) return 1;
    else return 0;
}

Double_t RooSumTwoPowerLawPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

    R__ASSERT(code == 1);

    const Double_t xmin = x.min(rangeName);
    const Double_t xmax = x.max(rangeName);

    if (xmin == x.min() && xmax == x.max()) return 1.;

    double zm = 91.2;
    if(&m.arg()) zm = m;
    if(not offset) zm = 0;
  
    double integral_part1 = 0.;
    if(a1 == 0)       integral_part1 = (xmax-xmin)/(x.max()-x.min());
    else if(a1 == -1) integral_part1 = (log(xmax-zm)-log(xmin-zm))/(log(x.max()-zm)-log(x.min()-zm));
    else              integral_part1 = (pow(xmax-zm,a1+1)-pow(xmin-zm,a1+1))/(pow(x.max()-zm,a1+1)-pow(x.min()-zm,a1+1));

    double integral_part2 = 0.;
    if(a2 == 0)       integral_part2 = (xmax-xmin)/(x.max()-x.min());
    else if(a2 == -1) integral_part2 = (log(xmax-zm)-log(xmin-zm))/(log(x.max()-zm)-log(x.min()-zm));
    else              integral_part2 = (pow(xmax-zm,a2+1)-pow(xmin-zm,a2+1))/(pow(x.max()-zm,a2+1)-pow(x.min()-zm,a2+1));
    
    return integral_part1 + integral_part2;

}

