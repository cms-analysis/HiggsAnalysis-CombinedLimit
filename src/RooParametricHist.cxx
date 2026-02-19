/*****************************************************************************
 *****************************************************************************/


#include "Riostream.h"

#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "../interface/RooParametricHist.h"
#include "../interface/CombineMathFuncs.h"

#include <math.h>
#include "TMath.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooFit.h"

#include "TFile.h"

ClassImp(RooParametricHist)

RooParametricHist::RooParametricHist(const char *name,
						 const char *title,
						 RooAbsReal& _x,
						 RooArgList& _pars,
						 const TH1 &_shape  // only need this to initialize bins
						 ) :
  RooAbsPdf(name,title),
  x("observable","observable",this,_x),
  pars("pars","pars",this),
  _coeffList("_coeffList","_coeffList",this)
  //SM_shape("SM_shape","SM_shape",this,_SM_shape),
{
  pars.add(_pars);
  if ( pars.getSize() != _shape.GetNbinsX() ){
	std::cerr << " Warning, number of parameters not equal to number of bins in shape histogram! " << std::endl;
	assert(0);
  }
  initializeBins(_shape);
//  initializeNorm();
  _cval = -1;
  _smoothRegion = 0.;
  _hasMorphs = false;
}

//_____________________________________________________________________________
RooParametricHist::RooParametricHist(const RooParametricHist& other, const char* name) :
 RooAbsPdf(other, name),x("observable",this,other.x),pars("pars",this,other.pars),_coeffList("_coeffList",this,other._coeffList)
{

  N_bins = other.N_bins;
  _smoothRegion=other._smoothRegion;
  _hasMorphs=other._hasMorphs;
  _cval = other._cval;

  for(int i=0; i<=N_bins; i++) {
     bins.push_back(other.bins[i]);
     if (i<N_bins) {
      widths.push_back(other.widths[i]);
      if (other._hasMorphs){
        std::vector<double> su;
        std::vector<double> di;
        for (int j=0; j<other._coeffList.getSize();j++){
          su.push_back(other._sums[i][j]);
          di.push_back(other._diffs[i][j]);
        }
        _sums.push_back(su);
        _diffs.push_back(di);
      }
     }
  }
}
void RooParametricHist::initializeBins(const TH1 &shape) const {

  N_bins = shape.GetNbinsX();
  for(int i=1; i<=N_bins+1; ++i) {
     bins.push_back(shape.GetBinLowEdge(i));
     if (i<=N_bins) widths.push_back(shape.GetBinWidth(i));
  }
}

RooAbsArg & RooParametricHist::getBinVar(const int i) const {
  if (i > N_bins ) std::cerr  << " Error in RooParametricHist::getBinBar -- Asked for bin " << i << " which is more than N_bins-1 -> " << N_bins << std::endl;
  return *pars.at(i);
}

RooArgList & RooParametricHist::getAllBinVars() const {
  return (RooArgList&)pars;
}

const double RooParametricHist::quickSum() const {
  getParVals();
  getCoeffs();
  getFlattenedMorphs(diffs_flat_, sums_flat_);

  return RooFit::Detail::MathFuncs::parametricHistFullSum(pars_vals_.data(),
                                                          N_bins,
                                                          _hasMorphs,
                                                          _coeffList.getSize(),
                                                          coeffs_.data(),
                                                          _hasMorphs ? diffs_flat_.data() : nullptr,
                                                          _hasMorphs ? sums_flat_.data() : nullptr,
                                                          _smoothRegion);
}

Int_t RooParametricHist::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet & analVars, const char*) const  {
  if (matchArgs(allVars,analVars,x)){
    return 1;
  }
  return 0;
}

Double_t RooParametricHist::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code == 1);
  getParVals();
  getCoeffs();
  getFlattenedMorphs(diffs_flat_, sums_flat_);
  const double* bins_ptr = bins.data();
  const double* widths_ptr = widths.data();
  return RooFit::Detail::MathFuncs::parametricHistIntegral(pars_vals_.data(),
                                                           bins_ptr,
                                                           N_bins,
                                                           coeffs_.data(),
                                                           _coeffList.getSize(),
                                                           _hasMorphs ? diffs_flat_.data() : nullptr,
                                                           _hasMorphs ? sums_flat_.data() : nullptr,
                                                           widths_ptr,
                                                           _smoothRegion,
                                                           rangeName,
                                                           x.min(rangeName),
                                                           x.max(rangeName));
}

void RooParametricHist::addMorphs(RooDataHist &hpdfU, RooDataHist &hpdfD, RooRealVar &cVar, double smoothRegion){
 
  if (!_hasMorphs){
    for (int i=0;i<N_bins;i++){
      std::vector<double> difs;
      std::vector<double> sums;
      _diffs.push_back(difs);
      _sums.push_back(sums);
    }
  }
  for (int i=0;i<N_bins;i++){
    double f0 = static_cast<RooAbsReal*>(pars.at(i))->getVal();
   
    hpdfU.get(i); hpdfD.get(i);
    double dh = (hpdfU.weight()-f0);
    double dl = (hpdfD.weight()-f0);
    _diffs[i].push_back(dh-dl);
    _sums[i].push_back(dh+dl);
  }
  _coeffList.add(cVar);
  _hasMorphs = true;
  _smoothRegion = smoothRegion;
}

Double_t RooParametricHist::evaluate() const {
  // Find which bin we are in and lookup the parameter value
  double xVal = getX();
  int bin_i = RooFit::Detail::MathFuncs::parametricHistFindBin(N_bins, bins, xVal);
  if (bin_i < 0)
    return 0.0;  // Out of range
  pars_vals_.resize(N_bins);
  pars_vals_[bin_i] = getParVal(bin_i);
  if (bins.empty() || widths.empty())
    return 0.0;
  const double* bins_ptr = bins.data();
  const double* widths_ptr = widths.data();
  int nMorphs = _coeffList.getSize();
  getFlattenedMorphs(diffs_flat_, sums_flat_);

  if (_hasMorphs) {
    getCoeffs();
  }

  return RooFit::Detail::MathFuncs::parametricHistEvaluate(bin_i,
                                                           pars_vals_.data(),
                                                           bins_ptr,
                                                           N_bins,
                                                           _hasMorphs ? coeffs_.data() : nullptr,
                                                           nMorphs,
                                                           _hasMorphs ? diffs_flat_.data() : nullptr,
                                                           _hasMorphs ? sums_flat_.data() : nullptr,
                                                           widths_ptr,
                                                           _smoothRegion);
}

double RooParametricHist::getParVal(int bin_i) const { return static_cast<RooAbsReal*>(pars.at(bin_i))->getVal(); }

const std::vector<double>& RooParametricHist::getParVals() const {
  pars_vals_.clear();
  pars_vals_.reserve(pars.getSize());
  for (int i = 0; i < pars.getSize(); ++i) {
    pars_vals_.push_back(static_cast<RooAbsReal*>(pars.at(i))->getVal());
  }
  return pars_vals_;
}

const std::vector<double>& RooParametricHist::getCoeffs() const {
  coeffs_.clear();
  coeffs_.reserve(_coeffList.getSize());
  for (int i = 0; i < _coeffList.getSize(); ++i) {
    coeffs_.push_back(static_cast<RooRealVar*>(_coeffList.at(i))->getVal());
  }
  return coeffs_;
}

void RooParametricHist::getFlattenedMorphs(std::vector<double>& diffs_flat, std::vector<double>& sums_flat) const {
  if (!_hasMorphs)
    return;
  // _diffs/_sums are immutable after construction, so only flatten once
  if (!diffs_flat.empty())
    return;

  int nMorphs = _coeffList.getSize();
  diffs_flat.reserve(N_bins * nMorphs);
  sums_flat.reserve(N_bins * nMorphs);

  for (int i = 0; i < N_bins; ++i) {
    for (int j = 0; j < nMorphs; ++j) {
      diffs_flat.push_back(_diffs[i][j]);
      sums_flat.push_back(_sums[i][j]);
    }
  }
}
