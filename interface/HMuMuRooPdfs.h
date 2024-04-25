#ifndef HMUMUROOPDFS
#define HMUMUROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooListProxy.h"
#include "RooDataSet.h"

// Modified Breit-Wigner: BWZ, BWZRedux, BWZ x Bern, BWZRedux x Bern
class RooModZPdf : public RooAbsPdf {

 public:
  RooModZPdf() {};  
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a); // BWZ with mean and width fixed to Z-boson
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _c); // BWZRedux without second exp param
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c); // BWZRedux with mean and width fixed to Z-boson
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m); // BWZRedux with BW mean as parameter, width fixed
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w); // BWZRedux with BW mean and width as parameters
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, const RooArgList& _coef); // BWZ x Bernstein with mean and width fixed to Z-boson
  RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef); // BWZRedux x Bernstein with mean and width fixed to Z-boson 
  RooModZPdf(const RooModZPdf& other, const char* name=0) ;

  virtual TObject* clone(const char* newname) const {
    return new RooModZPdf(*this,newname);
  }
  
  inline virtual ~RooModZPdf() {}
  
 protected:
  RooRealProxy x;
  RooRealProxy a;
  RooRealProxy b;
  RooRealProxy c;
  RooRealProxy m;
  RooRealProxy w;
  
  RooListProxy bernCoef;  
  Double_t evaluate() const ;
  
 private:
  ClassDef(RooModZPdf,1)
};


/////////
class RooExpPdf : public RooAbsPdf {

 public:
  RooExpPdf() {};  
  RooExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, const bool & _offset = true);   // standard exponential centered at Z-boson mass
  RooExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m);  // exponential with free paramaters for center x0
  RooExpPdf(const RooExpPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const {
    return new RooExpPdf(*this, newname);
  }
  
  inline virtual ~RooExpPdf() {
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  
 protected:
  RooRealProxy x;
  RooRealProxy a1;
  RooRealProxy m;
  bool offset;
  Double_t evaluate() const ;
  
 private:  
  ClassDef(RooExpPdf,1)

};


///////////
class RooSumTwoExpPdf : public RooAbsPdf {

 public:
  RooSumTwoExpPdf() {};  
  RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, const bool & _offset = true);
  RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, RooAbsReal& _m);  
  RooSumTwoExpPdf(const RooSumTwoExpPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const {
    return new RooSumTwoExpPdf(*this, newname);
  }
  
  inline virtual ~RooSumTwoExpPdf() {
  }
    
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

 protected:
  RooRealProxy x;
  RooRealProxy a1;
  RooRealProxy a2;
  RooRealProxy f;
  RooRealProxy m;
  bool offset;
  Double_t evaluate() const ;
  
 private:  
  ClassDef(RooSumTwoExpPdf,1)
    
};

/////////
class RooPowerLawPdf : public RooAbsPdf {

 public:
  RooPowerLawPdf() {};  
  RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, const bool & _offset = true);  
  RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m);  
  RooPowerLawPdf(const RooPowerLawPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const {
    return new RooPowerLawPdf(*this, newname);
  }
  
  inline virtual ~RooPowerLawPdf() {
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  
 protected:
  RooRealProxy x;
  RooRealProxy a1;
  RooRealProxy m;
  bool offset;
  Double_t evaluate() const ;

 private:  
  ClassDef(RooPowerLawPdf,1)

};

////////////
class RooSumTwoPowerLawPdf : public RooAbsPdf {

 public:
  RooSumTwoPowerLawPdf() {};  
  RooSumTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, const bool & _offset = true);
  RooSumTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f, RooAbsReal& _m);
  RooSumTwoPowerLawPdf(const RooSumTwoPowerLawPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const {
    return new RooSumTwoPowerLawPdf(*this, newname);
  }
  
  inline virtual ~RooSumTwoPowerLawPdf() {
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
    
 protected:
  RooRealProxy x;
  RooRealProxy a1;
  RooRealProxy a2;
  RooRealProxy f;
  RooRealProxy m;
  bool offset;
  Double_t evaluate() const ;
  
 private:  
  ClassDef(RooSumTwoPowerLawPdf,1)

};


#endif
