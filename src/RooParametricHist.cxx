/*****************************************************************************
 *****************************************************************************/


#include "Riostream.h"

#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "../interface/RooParametricHist.h"

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
  pars("pars","pars",this)
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
 RooAbsPdf(other, name),x("observable",this,other.x),pars("_pars",this,RooListProxy()),_coeffList("_coeffList",this,RooListProxy())
{

  N_bins = other.N_bins;
  _smoothRegion=other._smoothRegion;
  _hasMorphs=other._hasMorphs;
  _cval = other._cval;

  pars.add(other.pars);
  _coeffList.add(other._coeffList);

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

double RooParametricHist::getFullSum() const {
    double sum=0;
    for (int i = 0; i < pars.getSize(); ++i) {
	  double thisVal = static_cast<RooAbsReal&>(pars[i]).getVal();
	  if (_hasMorphs) thisVal*=evaluateMorphFunction(i);
	  sum+=thisVal;
	  i++;
    }
    return sum;
}

Int_t RooParametricHist::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet & analVars, const char*) const  {
  if (matchArgs(allVars,analVars,x)){
    return 1;
  }
  return 0;
}

Double_t RooParametricHist::analyticalIntegral(Int_t code, const char* rangeName) const
{
 assert(code==1) ;

 // Case without range is trivial: p.d.f is by construction normalized
 if (!rangeName) {
   //return 1;//getFullSum() ;
   return getFullSum();
 }
 // Case with ranges, calculate integral explicitly
 double xmin = x.min(rangeName) ;
 double xmax = x.max(rangeName) ;
 double sum=0 ;
 int i ;
 for (i=1 ; i<=N_bins ; i++) {
   double binVal = (static_cast<RooAbsReal*>(pars.at(i-1))->getVal())/widths[i-1];
   if (_hasMorphs) binVal*= evaluateMorphFunction(i-1);
   if (bins[i-1]>=xmin && bins[i]<=xmax) {
      // Bin fully in the integration domain
      sum += (bins[i]-bins[i-1])*binVal ;
   } else if (bins[i-1]<xmin && bins[i]>xmax) {
      // Domain is fully contained in this bin
      sum += (xmax-xmin)*binVal ;
      // Exit here, this is the last bin to be processed by construction
      return sum/getFullSum() ;
   } else if (bins[i-1]<xmin && bins[i]<=xmax && bins[i]>xmin) {
      // Lower domain boundary is in bin
      sum +=  (bins[i]-xmin)*binVal ;
   } else if (bins[i-1]>=xmin && bins[i]>xmax && bins[i-1]<xmax) {
      sum +=  (xmax-bins[i-1])*binVal ;
      // Upper domain boundary is in bin
      // Exit here, this is the last bin to be processed by construction
      return sum ;
   }
 }
 return sum;
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
  smoothRegion = _smoothRegion;
}

double RooParametricHist::evaluateMorphFunction(int j) const
{
    double scale=1.0;
    if (!_hasMorphs) return scale;
    
    int ndim = _coeffList.getSize();
    double f0 = static_cast<RooAbsReal*>(pars.at(j))->getVal();
    // apply all morphs one by one to the bin
    // almost certaintly a faster way to do this in a vectorized way ....
    for (int i = 0; i < ndim; ++i) {
        double x = (dynamic_cast<RooRealVar*>(_coeffList.at(i)))->getVal();
        double a = 0.5*x, b = smoothStepFunc(x);
	scale *= 1+(1./f0) * a*(_diffs[j][i] + b*_sums[j][i]);
	//std::cout << " at coeff " << (dynamic_cast<RooRealVar*>(_coeffList.at(i)))->GetName() << " = " << x << std::endl;
	//std::cout << " .... scale is now " << scale << std::endl;
    }
    return scale;
}
double RooParametricHist::evaluatePartial() const
{
  auto it = std::upper_bound(std::begin(bins), std::end(bins), x);
  if ( it == std::begin(bins) ) {
    // underflow
    return 0;
  }
  else if ( it == std::end(bins) ) {
    // overflow
    return 0;
  }
  size_t bin_i = std::distance(std::begin(bins), it) - 1;
  RooAbsReal *retVar = (RooAbsReal*)pars.at(bin_i);

  double ret = retVar->getVal();
  ret /= widths[bin_i];
  return ret;
}

double RooParametricHist::evaluateFull() const
{
  int bin_i;
  if (x < bins[0]) 	   			return 0;    // should set to 0 instead?
  else if (x >= bins[N_bins])   return 0;

  else {
    for(bin_i=0; bin_i<N_bins; bin_i++) {   // faster way to loop through ?
      if (x>=bins[bin_i] && x < bins[bin_i+1] ) break;
    }
  }
  double mVar = evaluateMorphFunction(bin_i);
  RooAbsReal *retVar = (RooAbsReal*)pars.at(bin_i);

  double ret = retVar->getVal()*mVar;
  ret /= widths[bin_i];
  return ret;
}

Double_t RooParametricHist::evaluate() const
{
  double ret = _hasMorphs ? evaluateFull() : evaluatePartial() ;
  _cval=ret;
  return ret > 0 ? ret : 0;
}
