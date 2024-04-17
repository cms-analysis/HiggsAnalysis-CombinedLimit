// -*- mode: C++ -*-

#ifndef ROOWWLNUJJPDFS_H
#define ROOWWLNUJJPDFS_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
//#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
//#include "RooAbsCategory.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include <complex> 

class RooChebyshevPDF : public RooAbsPdf {
public:

  RooChebyshevPDF() {} ;
  RooChebyshevPDF(const char *name, const char *title, RooAbsReal& var,
                  RooArgList& coefList);
  RooChebyshevPDF(const RooChebyshevPDF& other, const char *name=0);
  TObject* clone(const char *newname) const override {
    return new RooChebyshevPDF(*this, newname); }

  ~RooChebyshevPDF() override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, 
                              const char *rangeName = 0) const override;
  Double_t analyticalIntegral(Int_t code, const char *rangeName = 0) const override;

  static Double_t ChebyshevP(Int_t order, Double_t v);

protected:

  RooRealProxy x;
  RooListProxy coefs;

  Double_t evaluate() const override;

private:

  ClassDefOverride(RooChebyshevPDF,1) //Chebyshev polynomial implementation.
};


/*
  PDF for an erf((x-turnOn)/width)
 */
 
class RooErfPdf : public RooAbsPdf {
public:
  RooErfPdf() {} ; 
  RooErfPdf(const char *name, const char *title,
	    RooAbsReal& _x,
	    RooAbsReal& _turnOn,
	    RooAbsReal& _width,
	    int _onOff = 1);
  RooErfPdf(const RooErfPdf& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooErfPdf(*this,newname); }
  inline ~RooErfPdf() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;
  void	printMultiline(std::ostream& os, Int_t contents, Bool_t verbose = kFALSE, TString indent = "") const override ;

protected:

  double indefErfIntegral(double val) const ;

  RooRealProxy x ;
  RooRealProxy turnOn ;
  RooRealProxy width ;
  int onOff;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooErfPdf,1) // Your description goes here...
};
 
/*
  PDF for an x^(power) * exp(c*x)
 */
 
class RooPowerExpPdf : public RooAbsPdf {
public:
  RooPowerExpPdf() {} ; 
  RooPowerExpPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _c,
	      RooAbsReal& _power);
  RooPowerExpPdf(const RooPowerExpPdf& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooPowerExpPdf(*this,newname); }
  inline ~RooPowerExpPdf() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy power ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooPowerExpPdf,1) // Your description goes here...
};
 

/*
  PDF using a TH1D for the density function.  It is quick and dirty not
  with a lot of precision.
 */


class RooTH1DPdf : public RooAbsPdf {
public:
  RooTH1DPdf() {} ; 
  RooTH1DPdf(const char *name, const char *title,
	     RooAbsReal& _x,
	     TH1D& _hist, bool _interpolate = false);
  RooTH1DPdf(const RooTH1DPdf& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooTH1DPdf(*this,newname); }
  inline ~RooTH1DPdf() override { }

  TH1D & getHist() { return hist ; }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x ;
  mutable TH1D hist ;
  bool interpolate ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooTH1DPdf,1) // Your description goes here...
};
 
class RooPowerFunction : public RooAbsReal {
public:
  RooPowerFunction() {} ; 
  RooPowerFunction(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _r);
  RooPowerFunction(const RooPowerFunction& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooPowerFunction(*this,newname); }
  inline ~RooPowerFunction() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x ;
  RooRealProxy r ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooPowerFunction,1) // Your description goes here...
};

/*
  PDF for x^(power)
 */

class RooPowerLaw : public RooAbsPdf {
public:
  RooPowerLaw() { } ;
  RooPowerLaw(const char * name, const char * title, RooAbsReal& _x,
	      RooAbsReal& _power) ;
  RooPowerLaw(const RooPowerLaw& other, const char * name = 0) ;

  TObject * clone(const char * newname) const override { return new RooPowerLaw(*this, newname); }

  inline ~RooPowerLaw() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& anVars, 
  			      const char * rangeName = 0) const override ;
  Double_t analyticalIntegral(Int_t code, const char * rangeName = 0) const override ;

protected:
  RooRealProxy x;
  RooRealProxy power;

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooPowerLaw, 1) // Power law PDF
};
 
class RooExpPoly : public RooAbsPdf {
public:
  RooExpPoly() {} ; 
  RooExpPoly(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooArgList& _coefs);
  RooExpPoly(const RooExpPoly& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooExpPoly(*this,newname); }
  inline ~RooExpPoly() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, 
			      const char *rangeName = 0) const override;
  Double_t analyticalIntegral(Int_t code, const char *rangeName = 0) const override;

  //Emuation of deprecated RooComplex (that used Double_t)
  typedef std::complex<Double_t> DoubleComplex_t;

  static DoubleComplex_t erfi(DoubleComplex_t xval);
protected:

  RooRealProxy x ;
  RooListProxy coefs ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooExpPoly,1) // Your description goes here...
};

#endif
