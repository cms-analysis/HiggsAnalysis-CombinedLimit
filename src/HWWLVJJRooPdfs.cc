#include "Riostream.h" 

#include "HiggsAnalysis/CombinedLimit/interface/HWWLVJJRooPdfs.h" 
#include "RooAbsReal.h" 
//#include "RooAbsCategory.h" 
#include "RooMath.h"
#include <math.h> 
#include "TMath.h" 
#include "TFile.h"

ClassImp(RooChebyshevPDF)

RooChebyshevPDF::RooChebyshevPDF(const char *name, const char *title, 
				 RooAbsReal& var, RooArgList& coefList) :
  RooAbsPdf(name,title),
  x("x", "x", this, var),
  coefs("coefs", "coefs", this) {

  TIterator *cx = coefList.createIterator();
  RooAbsReal *coef;
  while ((coef = (RooAbsReal*)cx->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cerr << "Coefficient " << coef->GetName() << " is not good." << std::endl;
      assert(0);
    }
    coefs.add(*coef);
  }
  delete cx;
}

RooChebyshevPDF::RooChebyshevPDF(const RooChebyshevPDF& other, 
				 const char *name) :
  RooAbsPdf(other, name),
  x("x", this, other.x),
  coefs("coefs", this, other.coefs) {

}

RooChebyshevPDF::~RooChebyshevPDF() {
}

Double_t RooChebyshevPDF::ChebyshevP(Int_t order, Double_t v) {
  if (order == 0) return 1.0;
  if (order == 1) return v;
  if (order == 2) return 2.0*v*v - 1.0;
  if (order == 3) return 4.0*v*v*v - 3.0*v;
  if (order == 4) return 8.0*v*v*v*v - 8.0*v*v + 1.0;
  if (order == 5) return 16.*v*v*v*v*v - 20.*v*v*v + 5.*v;
  if (order == 6) return 32.*v*v*v*v*v*v - 48.*v*v*v*v + 18.*v*v - 1.0;
  if (order == 7) return 64.*v*v*v*v*v*v*v - 112.*v*v*v*v*v + 56.*v*v*v - 7.*v;
  if (order > 7) return 2.0*v*ChebyshevP(order-1,v)-ChebyshevP(order-2,v);
  assert(order > -1);
  return 0.0;
}

Double_t RooChebyshevPDF::evaluate() const {
  RooAbsReal *coef;
  Int_t tord = 1;
  Double_t val = 1.0;
  Double_t v = 2.0/(x.max()-x.min())*(x - (x.min() + (x.max()-x.min())/2.0));
  for (tord = 1; tord <= coefs.getSize(); ++tord) {
    coef = dynamic_cast<RooAbsReal*>(coefs.at(tord-1));
    val += coef->getVal()*ChebyshevP(tord,v);
  }
  return val;
}

Int_t RooChebyshevPDF::getAnalyticalIntegral(RooArgSet& allVars, 
					  RooArgSet& analVars,
					  const char */*rangeName*/) const {
  if (matchArgs(allVars, analVars, x)) return 1;
  return 0;
}

Double_t RooChebyshevPDF::analyticalIntegral(Int_t code, 
					     const char *rangeName) const {
  assert(code);
  RooAbsReal *coef;
  Int_t tord = 1;
  Double_t dv = (x.max()-x.min())/2.0;
  Double_t vmax = 2.0/(x.max()-x.min())*(x.max(rangeName) - 
					 (x.min() + (x.max()-x.min())/2.0));
  Double_t vmin = 2.0/(x.max()-x.min())*(x.min(rangeName) - 
					 (x.min() + (x.max()-x.min())/2.0));
  Double_t val = vmax-vmin;
  for (tord = 1; tord <= coefs.getSize(); ++tord) {
    coef = dynamic_cast<RooAbsReal*>(coefs.at(tord-1));
    if (tord == 1) {
      val += coef->getVal()*vmax*vmax/2.0;
      val -= coef->getVal()*vmin*vmin/2.0;
    } else {
      val += coef->getVal()*(tord*ChebyshevP(tord+1, vmax)/(tord*tord-1) - 
			     vmax*ChebyshevP(tord, vmax)/(tord-1));
      val -= coef->getVal()*(tord*ChebyshevP(tord+1, vmin)/(tord*tord-1) - 
			     vmin*ChebyshevP(tord, vmin)/(tord-1));
    }
  }
  return dv*val;
}

ClassImp(RooErfPdf) 

RooErfPdf::RooErfPdf(const char *name, const char *title, 
		     RooAbsReal& _x,
		     RooAbsReal& _turnOn,
		     RooAbsReal& _width,
		     int _onOff) : RooAbsPdf(name,title), 
  x("x","x",this,_x),
  turnOn("turnOn","turnOn",this,_turnOn),
  width("width","width",this,_width),
  onOff(_onOff)
{ 
  if (_onOff < 0)
    onOff = -1;
  else
    onOff = 1;
} 


RooErfPdf::RooErfPdf(const RooErfPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  turnOn("turnOn",this,other.turnOn),
  width("width",this,other.width),
  onOff(other.onOff)
{ 
} 



 Double_t RooErfPdf::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return (1.+onOff*TMath::Erf((x-turnOn)/width))*0.5 ; 
 }



 Int_t RooErfPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
 { 
   // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
   // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
   // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
   // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
   // EXPRESSION MULTIPLE TIMES

   if (matchArgs(allVars,analVars,x)) return 1 ; 
   return 0 ; 
 } 



 Double_t RooErfPdf::analyticalIntegral(Int_t code, const char* rangeName) const  
 { 
   // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
   // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
   // BOUNDARIES FOR EACH OBSERVABLE x

   if (code==1) {
     
     return  0.5*(x.max(rangeName)-x.min(rangeName) + 
		  onOff*(indefErfIntegral(x.max(rangeName)) - 
			 indefErfIntegral(x.min(rangeName))));
   } 
   return 0 ; 
 } 


double RooErfPdf::indefErfIntegral(double val) const {
  static double const rootpi = TMath::Sqrt(TMath::Pi());
  return (val-turnOn)*TMath::Erf((val-turnOn)/width) + 
    width/rootpi*TMath::Exp(-(val-turnOn)*(val-turnOn)/width/width);
}

void RooErfPdf::printMultiline(std::ostream& os, Int_t contents, 
			       Bool_t verbose, TString indent) const {
  RooAbsPdf::printMultiline(os, contents,verbose,indent);
  os << indent << "--- RooErfPdf --" << '\n';
  os << indent << "onOff: " << onOff << std::endl;
}

ClassImp(RooPowerExpPdf) 

RooPowerExpPdf::RooPowerExpPdf(const char *name, const char *title, 
			       RooAbsReal& _x,
			       RooAbsReal& _c,
			       RooAbsReal& _power) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  c("c","c",this,_c),
  power("power","power",this,_power)
{ 
} 


RooPowerExpPdf::RooPowerExpPdf(const RooPowerExpPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  c("c",this,other.c),
  power("power",this,other.power)
{ 
} 



Double_t RooPowerExpPdf::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  return TMath::Power(x,power)*TMath::Exp(c*x) ; 
} 



Int_t RooPowerExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const  
{ 
  // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
  // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
  // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
  // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
  // EXPRESSION MULTIPLE TIMES
  
  if ((-1*c*x.min(rangeName)> 0) && (matchArgs(allVars,analVars,x)))
    return 1 ; 
  return 0 ; 
} 

double ExpIntegralE(double p, double x) {
    
  // std::cout << "ExpIntegralE(" << p << "," << x << ")\n";
  double sum(0.), term(10.);
  int k(0);
  do {
    term = TMath::Exp(k*TMath::Log(x)-TMath::LnGamma(2-p+k));
    sum += term;
    ++k;
    //std::cout << "  term " << k << ": " << term << '\n';
  } while ((TMath::Abs(term/sum) > 1e-9)&&(k<1000));
    
  double ret(TMath::Gamma(1-p)*(TMath::Power(x, p-1) - TMath::Exp(-x)*sum));
  //std::cout << "returns: " << ret << '\n';
  return ret;
}

double incomplete_gamma(double s, double x) {
  double ret(TMath::Gamma(s, x)*TMath::Gamma(s));
  // std::cout << "incomplete_gamma(" << s << ',' << x << ") = " << ret << '\n';
  return ret;
}

Double_t RooPowerExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const  
{ 
  // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
  // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
  // BOUNDARIES FOR EACH OBSERVABLE x

  if (code==1) { 
    double maxTerm = -1*TMath::Power(x.max(rangeName), power+1) * 
      ExpIntegralE(-1*power, -1*c*x.max(rangeName));
    double minTerm = -1*TMath::Power(x.min(rangeName), power+1) * 
      ExpIntegralE(-1*power, -1*c*x.min(rangeName));
    // double maxTerm = -1*incomplete_gamma(power+1, -1*c*x.max(rangeName)) *
    //   TMath::Power(x.max(rangeName), power+1) * 
    //   TMath::Power((-1*c*x.max(rangeName)), -1*power -1);
    // double minTerm = -1*incomplete_gamma(power+1, -1*c*x.min(rangeName)) *
    //   TMath::Power(x.min(rangeName), power+1) * 
    //   TMath::Power((-1*c*x.min(rangeName)), -1*power -1);
    // std::cout << "c " << c << " power " << power << '\n';
    // std::cout << "integral: " << maxTerm << " - " << minTerm << " = " << 
    //   (maxTerm-minTerm) << '\n';
    return (maxTerm-minTerm) ; 
  } 
  return 0 ; 
}

ClassImp(RooTH1DPdf) 
  
RooTH1DPdf::RooTH1DPdf(const char *name, const char *title, 
		       RooAbsReal& _x,
		       TH1D& _hist, bool _interpolate) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  hist(_hist),
  interpolate(_interpolate)
{ 
  hist.SetDirectory(0);
} 


RooTH1DPdf::RooTH1DPdf(const RooTH1DPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  hist(other.hist),
  interpolate(other.interpolate)
{ 
  hist.SetDirectory(0);
} 



Double_t RooTH1DPdf::evaluate() const 
{ 
  if (interpolate)
    return TMath::Max(hist.Interpolate(x), 0.);
  else
    return TMath::Max(hist.GetBinContent(hist.FindBin(x)), 0.);
}

Int_t RooTH1DPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
{ 
  // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
  // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
  // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
  // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
  // EXPRESSION MULTIPLE TIMES
  
  if (matchArgs(allVars,analVars,x)) return 1 ; 
  return 0 ; 
} 



Double_t RooTH1DPdf::analyticalIntegral(Int_t code, const char* rangeName) const  
{ 
  // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
  // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
  // BOUNDARIES FOR EACH OBSERVABLE x

  if (code==1) {
    
    return hist.Integral(hist.FindBin(x.min(rangeName)), 
			 hist.FindBin(x.max(rangeName)), "width");
  } 
  return 0 ; 
} 

ClassImp(RooPowerLaw)


RooPowerLaw::RooPowerLaw(const char * name, const char * title, RooAbsReal& _x,
			 RooAbsReal& _power) : 
  RooAbsPdf(name, title),
  x("x", "dependent", this, _x), 
  power("power", "power", this, _power) {
}

RooPowerLaw::RooPowerLaw(const RooPowerLaw& other, const char * name) :
  RooAbsPdf(other, name), 
  x("x", this, other.x), 
  power("power", this, other.power) {
}

Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& anVars, 
					 const char * rangeName) const {
  if ((power == -1.0) && (x.min(rangeName) < 0.))
    return 0;
  if (matchArgs(allVars, anVars, x)) return 1 ;
  return 0;
}

Double_t RooPowerLaw::analyticalIntegral(Int_t code, 
					 const char * rangeName) const {
  switch(code) {
  case 1:
    Double_t ret(0.);
    if (power == -1.) {
      ret = TMath::Log(x.max(rangeName)) - TMath::Log(x.min(rangeName));
    } else {
      ret = (TMath::Power(x.max(rangeName), power + 1) - 
	     TMath::Power(x.min(rangeName), power + 1))/(power + 1);
    }
    return ret;
    break;
  }

  assert(0);
  return 0.;
}

Double_t RooPowerLaw::evaluate() const {
  return TMath::Power(x, power);
}

ClassImp(RooPowerFunction) 

 RooPowerFunction::RooPowerFunction(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _r) :
   RooAbsReal(name,title), 
   x("x","x",this,_x),
   r("r","r",this,_r)
 { 
 } 


 RooPowerFunction::RooPowerFunction(const RooPowerFunction& other, const char* name) :  
   RooAbsReal(other,name), 
   x("x",this,other.x),
   r("r",this,other.r)
 { 
 } 



 Double_t RooPowerFunction::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return pow(x,r) ; 
 } 



 Int_t RooPowerFunction::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
 { 
   // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
   // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
   // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
   // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
   // EXPRESSION MULTIPLE TIMES

   // if (matchArgs(allVars,analVars,x)) return 1 ; 
   return 0 ; 
 } 



 Double_t RooPowerFunction::analyticalIntegral(Int_t code, const char* rangeName) const  
 { 
   // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
   // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
   // BOUNDARIES FOR EACH OBSERVABLE x

   // assert(code==1) ; 
   // return (x.max(rangeName)-x.min(rangeName)) ; 
   return 0 ; 
 } 

ClassImp(RooExpPoly) 


RooExpPoly::RooExpPoly(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooArgList& _coefs) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  coefs("coefs","coefs",this)
{
  TIterator *cx = _coefs.createIterator();
  RooAbsReal *coef;
  while ((coef = (RooAbsReal*)cx->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cerr << "Coefficient " << coef->GetName() << " is not good." << std::endl;
      assert(0);
    }
    coefs.add(*coef);
  }
  delete cx;
}


RooExpPoly::RooExpPoly(const RooExpPoly& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  coefs("coefs",this,other.coefs)
{ 
} 



Double_t RooExpPoly::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  Double_t exponent(0.);
  Int_t order(0);
  RooAbsReal *coef;
  for (order = 1; order <= coefs.getSize(); ++order) {
    coef = dynamic_cast<RooAbsReal *>(coefs.at(order-1));
    exponent += coef->getVal()*TMath::Power(x, order);
  }
  return TMath::Exp(exponent) ; 
} 

Int_t RooExpPoly::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
					const char */*rangeName*/) const {
  if ( (coefs.getSize()<2) && (matchArgs(allVars, analVars, x)) )
    return 1;
  if ( (coefs.getSize()==2) && 
       // (dynamic_cast<RooAbsReal *>(coefs.at(1))->getVal() <= 0) && 
       (matchArgs(allVars, analVars, x)) )
    return 1;

  return 0;
}

Double_t RooExpPoly::analyticalIntegral(Int_t code, 
					const char *rangeName) const {
  assert(code);
  double coef1(0.), coef2(0.);
  static Double_t const rootpi(TMath::Sqrt(TMath::Pi()));
  
  Double_t val(x.max(rangeName) - x.min(rangeName));
  if (coefs.getSize() > 0) {
    coef1 = dynamic_cast<RooAbsReal *>(coefs.at(0))->getVal();
    if (coef1 != 0.)
      val = (TMath::Exp(coef1*x.max(rangeName)) -
	     TMath::Exp(coef1*x.min(rangeName)))/coef1;
    if (coefs.getSize() > 1) {
      coef2 = dynamic_cast<RooAbsReal *>(coefs.at(1))->getVal();
      double absrootc2(TMath::Sqrt(TMath::Abs(coef2)));
      if (coef2 != 0) {
	// std::cout << "coef1: " << coef1 << " coef2: " << coef2 << '\n';
	RooComplex c1(coef1, 0); RooComplex c2(coef2, 0);
	RooComplex rootc2((coef2 > 0) ? RooComplex(TMath::Sqrt(coef2),0) :
			  RooComplex(0,TMath::Sqrt(-1*coef2)));
	RooComplex zmax = c1*0.5/rootc2 + rootc2*x.max(rangeName);
	RooComplex zmin = c1*0.5/rootc2 + rootc2*x.min(rangeName);
	double erfimax((coef2 > 0) ? erfi(zmax).re() : erfi(zmax).im());
	double erfimin((coef2 > 0) ? erfi(zmin).re() : erfi(zmin).im());
	// if (coef2 < 0) {
	//   erfimax = TMath::Erf(absrootc2*x.max(rangeName) - 
	// 		       coef1*0.5/absrootc2);
	//   erfimin = TMath::Erf(absrootc2*x.min(rangeName) - 
	// 		       coef1*0.5/absrootc2);
	//   // std::cout << "erfimax: " << erfimax << '\n'
	//   // 	    << "erfimin: " << erfimin << '\n';
	// }
	double zval_max(erfimax*0.5*rootpi*
			TMath::Exp(-coef1*coef1/4./coef2)/absrootc2);
	double zval_min(erfimin*0.5*rootpi*
			TMath::Exp(-coef1*coef1/4./coef2)/absrootc2);
	// std::cout << "c1: "; c1.Print();
	// std::cout << "c2: "; c2.Print();
	// std::cout << "rootc2: "; rootc2.Print();
	// std::cout << "zmax: "; zmax.Print();
	// std::cout << "zval_max: " << zval_max << '\n';
	// std::cout << "zmin: "; zmin.Print();
	// std::cout << "zval_min: " << zval_min << '\n';
	val = zval_max - zval_min;
      }
    }
  }
  return val;
}

RooComplex RooExpPoly::erfi(RooComplex z) {
  static RooComplex myi(0,1);
  static RooComplex myone(1,0);

  RooComplex zsq(z*z);
  RooComplex ret(1,0);
  // std::cout << "z: "; z.Print();
  if (z.im() > 5.)
    ret =  RooComplex(0., 1.);
  else if (z.im() < -5.)
    ret =  RooComplex(0.,-1.);
  else {
    RooComplex w(RooMath::ComplexErrFunc(z));
    // std::cout << "w(z): "; w.Print();
    // std::cout << "zsq: "; zsq.Print();
    // std::cout << "zsq.exp(): "; zsq.exp().Print();
    ret = myi*(myone - ((zsq.exp())*(w)));
  }

  // std::cout << "erfi: "; ret.Print();
  return ret;
}
