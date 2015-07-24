#ifndef ROOPARAMETRICHIST
#define ROOPARAMETRICHIST

#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooAbsReal.h"
//#include "RooAbsData.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1.h"
#include "TString.h"
  
class RooParametricHist : public RooAbsPdf {
public:
  
  RooParametricHist() {} ;
  RooParametricHist (const char *name, const char *title, RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape);
  
  RooParametricHist(const RooParametricHist& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooParametricHist(*this,newname); }
  inline virtual ~RooParametricHist (){};
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t, const char* rangeName=0) const ;   
  
  //RooAddition & getYieldVar(){return sum;};
  
protected:
  
  RooRealProxy x;
  //RooAddition sum;
  RooListProxy pars;
  mutable int N_bins; 
  mutable std::vector<double> bins; 
  mutable std::vector<double> widths; 

  void initializeBins(const TH1&) const;
  //void initializeNorm();
  
  Double_t evaluate() const ;
  double getFullSum() const ;
 
  mutable double cval;
  void update_cval(double r){cval=r;};
private:
   ClassDef(RooParametricHist, 1) 
};

#endif
