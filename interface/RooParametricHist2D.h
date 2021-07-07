#ifndef ROOPARAMETRICHIST2D
#define ROOPARAMETRICHIST2D

#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooAbsReal.h"
//#include "RooAbsData.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH2D.h"
  
class RooParametricHist2D : public RooAbsPdf {
public:
  
  RooParametricHist2D() {} ;
  RooParametricHist2D (const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _y, RooArgList& _pars, const TH2 &_shape);
  
  RooParametricHist2D(const RooParametricHist2D& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooParametricHist2D(*this,newname); }
  inline virtual ~RooParametricHist2D (){};
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t, const char* rangeName=0) const ;   
  
  //RooAddition & getYieldVar(){return sum;};
  
protected:
  
  RooRealProxy x;
  RooRealProxy y;
  //RooAddition sum;
  RooListProxy pars;
  mutable int N_bins; 
  mutable int N_bins_x;
  mutable int N_bins_y;
  mutable std::vector<double> bins_x; 
  mutable std::vector<double> bins_y;
  mutable std::vector<double> widths_x;
  mutable std::vector<double> widths_y;

  void initializeBins(const TH2&) const;
  //void initializeNorm();
  
  Double_t evaluate() const ;
  double getFullSum() const ;
 
  mutable double cval;
  void update_cval(double r){cval=r;};
private:
   ClassDef(RooParametricHist2D, 1) 
};

#endif
