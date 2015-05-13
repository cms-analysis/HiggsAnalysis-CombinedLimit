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
#include "RooComplex.h"
 
class RooChebyshevPDF : public RooAbsPdf {
public:

  RooChebyshevPDF() {} ;
  RooChebyshevPDF(const char *name, const char *title, RooAbsReal& var,
                  RooArgList& coefList);
  RooChebyshevPDF(const RooChebyshevPDF& other, const char *name=0);
  virtual TObject* clone(const char *newname) const {
    return new RooChebyshevPDF(*this, newname); }

  virtual ~RooChebyshevPDF();

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, 
                              const char *rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char *rangeName = 0) const;

  static Double_t ChebyshevP(Int_t order, Double_t v);

protected:

  RooRealProxy x;
  RooListProxy coefs;

  Double_t evaluate() const;

private:

  ClassDef(RooChebyshevPDF,1) //Chebyshev polynomial implementation.
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
  virtual TObject* clone(const char* newname) const { return new RooErfPdf(*this,newname); }
  inline virtual ~RooErfPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  virtual void	printMultiline(std::ostream& os, Int_t contents, Bool_t verbose = kFALSE, TString indent = "") const ;

protected:

  double indefErfIntegral(double val) const ;

  RooRealProxy x ;
  RooRealProxy turnOn ;
  RooRealProxy width ;
  int onOff;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooErfPdf,1) // Your description goes here...
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
  virtual TObject* clone(const char* newname) const { return new RooPowerExpPdf(*this,newname); }
  inline virtual ~RooPowerExpPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy power ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooPowerExpPdf,1) // Your description goes here...
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
  virtual TObject* clone(const char* newname) const { return new RooTH1DPdf(*this,newname); }
  inline virtual ~RooTH1DPdf() { }

  TH1D & getHist() { return hist ; }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  mutable TH1D hist ;
  bool interpolate ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooTH1DPdf,1) // Your description goes here...
};
 
class RooPowerFunction : public RooAbsReal {
public:
  RooPowerFunction() {} ; 
  RooPowerFunction(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _r);
  RooPowerFunction(const RooPowerFunction& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPowerFunction(*this,newname); }
  inline virtual ~RooPowerFunction() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy r ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooPowerFunction,1) // Your description goes here...
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

  virtual TObject * clone(const char * newname) const { return new RooPowerLaw(*this, newname); }

  inline virtual ~RooPowerLaw() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& anVars, 
  			      const char * rangeName = 0) const ;
  Double_t analyticalIntegral(Int_t code, const char * rangeName = 0) const ;

protected:
  RooRealProxy x;
  RooRealProxy power;

  Double_t evaluate() const;

private:
  ClassDef(RooPowerLaw, 1) // Power law PDF
};
 
class RooExpPoly : public RooAbsPdf {
public:
  RooExpPoly() {} ; 
  RooExpPoly(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooArgList& _coefs);
  RooExpPoly(const RooExpPoly& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpPoly(*this,newname); }
  inline virtual ~RooExpPoly() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, 
			      const char *rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char *rangeName = 0) const;

  static RooComplex erfi(RooComplex xval);
protected:

  RooRealProxy x ;
  RooListProxy coefs ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooExpPoly,1) // Your description goes here...
};

#endif
