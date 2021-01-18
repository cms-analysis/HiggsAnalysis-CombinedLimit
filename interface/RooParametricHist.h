#ifndef ROOPARAMETRICHIST
#define ROOPARAMETRICHIST

#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooRealVar.h"
//#include "RooAbsData.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1.h"
#include "TString.h"
  
class RooParametricHist : public RooAbsPdf {
public:
  
  RooParametricHist() {} ;
  RooParametricHist (const char *name, const char *title, RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape);
  RooParametricHist (const char *name, const char *title, RooArgList& _pars, const RooParametricHist& other);
  
  RooParametricHist(const RooParametricHist& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooParametricHist(*this,newname); }
  inline virtual ~RooParametricHist (){};

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t, const char* rangeName=0) const ;   

  RooAbsArg  & getBinVar(const int i) const ;
  RooArgList & getAllBinVars() const ;
  
  RooRealVar & getObs() const { return (RooRealVar&)x; };
  const std::vector<double>  getBins()   const { return bins;   };
  const std::vector<double>  getWidths() const { return widths; };
 
  const double quickSum() const {return getFullSum() ;}
  //RooAddition & getYieldVar(){return sum;};
  
  void addMorphs(RooArgList _morphPdfs, RooArgList _coeffs)

protected:
  
  // return a smooth function that is equal to +/-1 for |x| >= smoothRegion_ and it's null in zero
  RooRealProxy x;
  //RooAddition sum;
  RooListProxy pars;
  mutable int N_bins; 
  mutable std::vector<double> bins; 
  mutable std::vector<double> widths; 

  // access morphing parameters for systematic functions 
  const double _smoothRegion;
  mutable std::vector<double> diffs; 
  mutable std::vector<double> sums; 
  

  void initializeBins(const TH1&) const;
  //void initializeNorm();
 
  RooListProxy _coefList; 

  Double_t evaluate() const ;
  double getFullSum() const ;
 
  mutable double cval;
  void update_cval(double r){cval=r;};
  
  inline double smoothStepFunc(double x) const { 
    if (fabs(x) >= _smoothRegion) return x > 0. ? +1. : -1.;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15.);
  }

private:
   ClassDef(RooParametricHist, 1) 
};

#endif
