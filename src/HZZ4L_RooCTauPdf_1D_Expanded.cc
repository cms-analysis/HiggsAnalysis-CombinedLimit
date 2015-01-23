#include "Riostream.h" 
#include "../interface/HZZ4L_RooCTauPdf_1D_Expanded.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;

ClassImp(HZZ4L_RooCTauPdf_1D_Expanded) 

HZZ4L_RooCTauPdf_1D_Expanded::HZZ4L_RooCTauPdf_1D_Expanded(
	  const char *name,
	  const char *title,
	  RooAbsReal& _kd,
	  RooAbsReal& _ctau,
    const RooArgList& inFuncList,
    const RooArgList& inCoefList,

	  double _ctau_min,
	  double _ctau_max,

    Double_t smoothRegion_ , Int_t smoothAlgo_
	  ):
RooAbsPdf(name, title),
kd("kd", "kd", this, _kd),
ctau("ctau", "ctau", this, _ctau),
_funcList("funcList", "List of funcficients", this),
_coefList("coefList", "List of coefficients", this),
ctau_min(_ctau_min),
ctau_max(_ctau_max),
smoothRegion(smoothRegion_),
smoothAlgo(smoothAlgo_)
{
  nCoef=0;

  TIterator* coefIter = inCoefList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR::HZZ4L_RooCTauPdf_1D_Expanded(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    _coefList.add(*coef);
    nCoef++;
  }
  delete coefIter;

  if (nCoef>2){
    coutE(InputArguments) << "ERROR::HZZ4L_RooCTauPdf_1D_Expanded(" << GetName() << ") number coefficients " << nCoef << " is not supported." << endl;
    assert(0);
  }
//  else cout << "Number of coefficients: " << nCoef << endl;

  TIterator* funcIter = inFuncList.createIterator();
	RooAbsArg* func;
	while ((func = (RooAbsArg*)funcIter->Next())) {
		if (!dynamic_cast<RooAbsReal*>(func)) {
			coutE(InputArguments) << "ERROR::HZZ4L_RooCTauPdf_1D_Expanded(" << GetName() << ") funcfions " << func->GetName() << " is not of type RooAbsReal" << endl;
			assert(0);
		}
		_funcList.add(*func);
	}
	delete funcIter;

  int integralArraySize = TMath::Min(_funcList.getSize(), 5*101);
	nbins_ctau = _funcList.getSize()/(2*nCoef+1);
//  cout << "Number of ctau bins: " << nbins_ctau << endl;
  for (int mp=0; mp<integralArraySize; mp++) Integral_T[mp] = dynamic_cast<const RooHistFunc*>(_funcList.at(mp))->analyticalIntegral(1000);
	for(int mp=nbins_ctau;mp<(5*101);mp++) Integral_T[mp] = Integral_T[nbins_ctau-1];

//	cout << "Finishing constructor, obtained the following integrals:" << endl;
//	for(int mp=0;mp<5*101;mp++) cout << "Integral[" << mp << "]: " << Integral_T[mp] << endl;
}

HZZ4L_RooCTauPdf_1D_Expanded::HZZ4L_RooCTauPdf_1D_Expanded(
	const HZZ4L_RooCTauPdf_1D_Expanded& other,
	const char* name
	):
RooAbsPdf(other, name),
kd("kd", this, other.kd),
ctau("ctau", this, other.ctau),
_funcList("funcList", this, other._funcList),
_coefList("coefList", this, other._coefList),
ctau_min(other.ctau_min),
ctau_max(other.ctau_max),
smoothRegion(other.smoothRegion),
smoothAlgo(other.smoothAlgo)
{
  nCoef = _coefList.getSize();
  nbins_ctau = _funcList.getSize()/(2*nCoef+1);
  for (int mp=0; mp<(5*101); mp++) Integral_T[mp] = (other.Integral_T)[mp];
}

int HZZ4L_RooCTauPdf_1D_Expanded::findNeighborBins() const{
	double gridsize = (ctau_max-ctau_min)/(nbins_ctau-1);
	int bincode = (int) (ctau-ctau_min)/gridsize;
//  cout << "Bin code " << bincode << "@ctau=" << ctau << endl;
	return bincode;
}
Double_t HZZ4L_RooCTauPdf_1D_Expanded::interpolateBin() const{
	int bincode = findNeighborBins();
	double gridsize = (ctau_max-ctau_min)/(nbins_ctau-1);
	Double_t result = 1.0e-100;
	int lowbin,highbin;
	if (bincode < 0){
		lowbin = 0; highbin = 1;
	}
	else if (bincode >= (nbins_ctau - 1)){
		lowbin = nbins_ctau - 2; highbin = nbins_ctau - 1;
	}
	else{
		lowbin = bincode; highbin = lowbin + 1;
	}

	Double_t v1_nominal = dynamic_cast<const RooHistFunc*>(_funcList.at(lowbin))->getVal();
	Double_t v2_nominal = dynamic_cast<const RooHistFunc*>(_funcList.at(highbin))->getVal();
  Double_t v1 = v1_nominal;
  Double_t v2 = v2_nominal;

  TIterator* coefIter = _coefList.createIterator();
  RooAbsReal* coef;
  int numCoef=0;
  while ((coef=(RooAbsReal*) coefIter->Next())) {
    Double_t coefVal = coef->getVal();
    numCoef++;
    Double_t varUp = dynamic_cast<const RooHistFunc*>(_funcList.at(lowbin + ((2*numCoef-1)*nbins_ctau)))->getVal();
    Double_t varDn = dynamic_cast<const RooHistFunc*>(_funcList.at(lowbin + (2*numCoef*nbins_ctau)))->getVal();
    v1 += interpolateVariation(coefVal, v1_nominal, varUp, varDn);
//    cout << "v1 up: " << varUp << ", down: " << varDn << ". theta: " << coefVal << " so " << v1_nominal << "->" << v1 << endl;
    varUp = dynamic_cast<const RooHistFunc*>(_funcList.at(highbin + ((2*numCoef-1)*nbins_ctau)))->getVal();
    varDn = dynamic_cast<const RooHistFunc*>(_funcList.at(highbin + (2*numCoef*nbins_ctau)))->getVal();
    v2 += interpolateVariation(coefVal, v2_nominal, varUp, varDn);
//    cout << "v2 up: " << varUp << ", down: " << varDn << ". theta: " << coefVal << " so " << v2_nominal << "->" << v2 << endl;
  }
  delete coefIter;

	Double_t fdistance = ( ( (ctau-ctau_min)/gridsize ) - ((double) lowbin) ) / (highbin-lowbin);
//  cout << "ctau=" << ctau << " => fraction=" << fdistance << endl;

	result = (1.-fdistance)*v1 + fdistance*v2;
	return result;
}
Double_t HZZ4L_RooCTauPdf_1D_Expanded::interpolateIntegral() const{
	int bincode = findNeighborBins();
	double gridsize = (ctau_max-ctau_min)/(nbins_ctau-1);
	Double_t result = 1.0e-100;
	int lowbin,highbin;
	if (bincode < 0){
		lowbin = 0; highbin = 1;
	}
	else if (bincode >= (nbins_ctau - 1)){
		lowbin = nbins_ctau - 2; highbin = nbins_ctau - 1;
	}
	else{
		lowbin = bincode; highbin = lowbin + 1;
	}

	Double_t v1_nominal = Integral_T[lowbin];
	Double_t v2_nominal = Integral_T[highbin];

  Double_t v1 = v1_nominal;
  Double_t v2 = v2_nominal;

  TIterator* coefIter = _coefList.createIterator();
  RooAbsReal* coef;
  int numCoef=0;
  while ((coef=(RooAbsReal*)coefIter->Next())) {
    Double_t coefVal = coef->getVal();
    numCoef++;
    Double_t varUp = Integral_T[lowbin + ((2*numCoef-1)*nbins_ctau)];
    Double_t varDn = Integral_T[lowbin + (2*numCoef*nbins_ctau)];
    v1 += interpolateVariation(coefVal, v1_nominal, varUp, varDn);
    varUp = Integral_T[highbin + ((2*numCoef-1)*nbins_ctau)];
    varDn = Integral_T[highbin + (2*numCoef*nbins_ctau)];
    v2 += interpolateVariation(coefVal, v2_nominal, varUp, varDn);
  }
  delete coefIter;


	Double_t fdistance = ( ( (ctau-ctau_min)/gridsize ) - ((double) lowbin) ) / (highbin-lowbin);

	result = (1.-fdistance)*v1 + fdistance*v2;
	return result;
}
Double_t HZZ4L_RooCTauPdf_1D_Expanded::interpolateVariation(Double_t theta_, Double_t valueCenter_, Double_t valueHigh_, Double_t valueLow_) const {
  if (smoothAlgo<0) return 0;
  else{
    if (fabs(theta_)>=smoothRegion) return theta_ * (theta_ > 0 ? valueHigh_ - valueCenter_ : valueCenter_ - valueLow_);

    Double_t c_up  = 0;
    Double_t c_dn  = 0;
    Double_t c_cen = 0;
    Double_t addVal = 0;

    if (smoothAlgo != 1) {
      // Quadratic interpolation null at zero and continuous at boundaries but not smooth at boundaries
      c_up  = +theta_ * (smoothRegion + theta_) / (2 * smoothRegion);
      c_dn  = -theta_ * (smoothRegion - theta_) / (2 * smoothRegion);
      c_cen = -theta_ * theta_ / smoothRegion;
    }
    else{
      // Quadratic interpolation that is everywhere differentiable but not null at zero
      c_up  = (smoothRegion + theta_) * (smoothRegion + theta_) / (4 * smoothRegion);
      c_dn  = (smoothRegion - theta_) * (smoothRegion - theta_) / (4 * smoothRegion);
      c_cen = -c_up - c_dn;
    }
//    cout << theta_ << ": " << addVal << endl;
    addVal = c_up * valueHigh_ + c_dn * valueLow_ + c_cen * valueCenter_;
    return addVal;
  }
}

Double_t HZZ4L_RooCTauPdf_1D_Expanded::evaluate() const
{
	double value = interpolateBin();
	if (value <= 0) return 1.0e-15;
	return value;
}
Int_t HZZ4L_RooCTauPdf_1D_Expanded::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg()))) return 1;
  return 0 ;
}
Double_t HZZ4L_RooCTauPdf_1D_Expanded::analyticalIntegral(Int_t code, const char* rangeName) const
{
	switch (code)
	{
	case 1:
		{
			double integral = interpolateIntegral();
			if( integral <= 0 ) integral = 1.0e-10;
			return integral;
		}
	default:
		cerr << "getAnalyticalIntegral failed, so analytical integration did not complete!" << endl;
		assert(0);
	}
}

