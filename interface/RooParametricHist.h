#ifndef ROOPARAMETRICHIST
#define ROOPARAMETRICHIST

#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
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
  TObject* clone(const char* newname) const override { return new RooParametricHist(*this,newname); }
  inline ~RooParametricHist () override{};

  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t, const char* rangeName=0) const override ;

  RooAbsArg  & getBinVar(const int i) const ;
  RooArgList & getAllBinVars() const ;

  RooRealVar & getObs() const { return (RooRealVar&)x; };
  RooAbsReal& observable() const { return const_cast<RooAbsReal&>(static_cast<const RooAbsReal&>(x.arg())); }
  const std::vector<double>  getBins()   const { return bins;   };
  const int getNBins() const { return N_bins; };
  const std::vector<double>  getWidths() const { return widths; };

  const RooArgList& getPars() const { return pars; };
  const RooArgList& getCoeffList() const { return _coeffList; };

  const double quickSum() const;
  //RooAddition & getYieldVar(){return sum;};

  // how can we pass this version? is there a Collection object for RooDataHists?
  //void addMorphs(RooArgList &_morphPdfsUp, RooArgList &_morphPdfsDown, RooArgList &_coeffs, double smoothRegion);
  void addMorphs(RooDataHist&, RooDataHist&, RooRealVar&, double );
  Double_t evaluate() const override;

  // Accessors for evaluation data
  double getX() const { return x; }
  double getSmoothRegion() const { return _smoothRegion; }
  bool hasMorphs() const { return _hasMorphs; }

  double getParVal(int bin_i) const;

  // Utility functions for data extraction
  const std::vector<double>& getParVals() const;
  const std::vector<double>& getCoeffs() const;
  void getFlattenedMorphs(std::vector<double>& diffs_flat, std::vector<double>& sums_flat) const;

protected:

  // return a smooth function that is equal to +/-1 for |x| >= smoothRegion and it's null in zero
  RooRealProxy x;
  //RooAddition sum;
  RooListProxy pars;
  RooListProxy _coeffList;
  mutable int N_bins;
  mutable std::vector<double> bins;
  mutable std::vector<double> widths;

  // access morphing parameters for systematic functions
  mutable double _smoothRegion;
  mutable bool   _hasMorphs;
  mutable std::vector<std::vector <double> > _diffs;
  mutable std::vector<std::vector<double> > _sums;

  mutable std::vector<double> pars_vals_;   //! Don't serialize me
  mutable std::vector<double> coeffs_;      //! Don't serialize me
  mutable std::vector<double> diffs_flat_;  //! Don't serialize me
  mutable std::vector<double> sums_flat_;   //! Don't serialize me

  void initializeBins(const TH1&) const;
  //void initializeNorm();

  mutable double _cval;
  void update_cval(double r){_cval=r;};

  inline double smoothStepFunc(double x) const {
    if (fabs(x) >= _smoothRegion) return x > 0. ? +1. : -1.;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15.);
  }

private:
   ClassDefOverride(RooParametricHist, 2)
};

#endif
