#ifndef HZZ2L2QROOPDFS
#define HZZ2L2QROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooCB : public RooAbsPdf {
 public:
  RooCB();
  RooCB(const char *name, const char *title,
        RooAbsReal& _x,
        RooAbsReal& _mean,
        RooAbsReal& _width,
        RooAbsReal& _alpha,
        RooAbsReal& _n,
        RooAbsReal& _theta
	);
  RooCB(const RooCB& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooCB(*this,newname); }
  inline ~RooCB() override { }

 protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha;
  RooRealProxy n;
  RooRealProxy theta;

  Double_t evaluate() const override ;

 private:

  ClassDefOverride(RooCB,1)
    };

 
class RooDoubleCB : public RooAbsPdf {
public:
  RooDoubleCB();
  RooDoubleCB(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mean,
	      RooAbsReal& _width,
	      RooAbsReal& _alpha1,
	      RooAbsReal& _n1,
	      RooAbsReal& _alpha2,
	      RooAbsReal& _n2
	   );
  RooDoubleCB(const char *name, const char *title,
        RooAbsReal& _x,
        RooAbsReal& _xp,
        RooAbsReal& _mean,
        RooAbsReal& _width,
        RooAbsReal& _alpha1,
        RooAbsReal& _n1,
        RooAbsReal& _alpha2,
        RooAbsReal& _n2
  );
  RooDoubleCB(const RooDoubleCB& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooDoubleCB(*this,newname); }
  inline ~RooDoubleCB() override { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x;
  RooRealProxy xp;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooDoubleCB,2)
};
 
class RooFermi : public RooAbsPdf {
public:
  RooFermi();
  RooFermi(const char *name, const char *title,
 	    RooAbsReal& _x,
            RooAbsReal& _cutOff,
	   RooAbsReal& _beta
	   );
  RooFermi(const RooFermi& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooFermi(*this,newname); }
  inline ~RooFermi() override { }

protected:

  RooRealProxy x ;
  RooRealProxy cutOff ;
  RooRealProxy beta ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooFermi,1) 
};
  
class RooRelBW : public RooAbsPdf {
public:
  RooRelBW();
  RooRelBW(const char *name, const char *title,
	   RooAbsReal& _x,
	   RooAbsReal& _mean,
	   RooAbsReal& _width,
	   RooAbsReal& _n
	   );
  RooRelBW(const RooRelBW& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooRelBW(*this,newname); }
  inline ~RooRelBW() override { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy n ;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooRelBW,1)
};
 

class Triangle : public RooAbsPdf {
public:
  Triangle();
  Triangle(const char *name, const char *title,                
	   RooAbsReal& _m,
	   RooAbsReal& _start,
	   RooAbsReal& _turn,
	   RooAbsReal& _stop
	   );	
  
  Triangle(const Triangle& other, const char* name = 0);
  TObject* clone(const char* newname) const override { 
    return new Triangle(*this,newname); }

  inline ~Triangle() override { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy m;
  RooRealProxy start;
  RooRealProxy turn;
  RooRealProxy stop;
  
  Double_t evaluate() const override;

private:
  
  ClassDefOverride(Triangle,1)
};


class RooLevelledExp : public RooAbsPdf {
 public:
  RooLevelledExp();
  RooLevelledExp(const char *name, const char *title,
		 RooAbsReal& _x,
		 RooAbsReal& _sigma,
		 RooAbsReal& _alpha,
		 RooAbsReal& _m,
		 //RooAbsReal& _k,
		 RooAbsReal& _theta
		);

  RooLevelledExp(const RooLevelledExp& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooLevelledExp(*this,newname); }
  inline ~RooLevelledExp() override { }

 protected:

  RooRealProxy x ;
  RooRealProxy sigma;
  RooRealProxy alpha;
  RooRealProxy m;
  // RooRealProxy k;
  RooRealProxy theta;
  

  Double_t evaluate() const override ;

 private:

  ClassDefOverride(RooLevelledExp,1)
    };


#endif
