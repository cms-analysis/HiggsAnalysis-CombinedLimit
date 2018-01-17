#include "Riostream.h" 
#include "../interface/HZZ4L_RooCTauPdf_1D_Expanded.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include <vector>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;
using namespace RooFit;
using namespace std;

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

    Double_t smoothRegion_, Int_t smoothAlgo_, Int_t nLinearVariations_
	  ):
RooAbsPdf(name, title),
kd("kd", "kd", this, _kd),
ctau("ctau", "ctau", this, _ctau),
_funcList("funcList", "List of funcficients", this),
_coefList("coefList", "List of coefficients", this),
ctau_min(_ctau_min),
ctau_max(_ctau_max),
smoothRegion(smoothRegion_),
smoothAlgo(smoothAlgo_),
nLinearVariations(nLinearVariations_)
{
  nCoef=0;
  nFuncs=0;

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
			coutE(InputArguments) << "ERROR::HZZ4L_RooCTauPdf_1D_Expanded(" << GetName() << ") function " << func->GetName() << " is not of type RooAbsReal" << endl;
			assert(0);
		}
		_funcList.add(*func);
    nFuncs++;
  }
	delete funcIter;

  if (nLinearVariations>nCoef){
    coutE(InputArguments) << "ERROR::HZZ4L_RooCTauPdf_1D_Expanded(" << GetName() << ") number of linear variations (" << nLinearVariations << ") is more than the number of coefficients (" << nCoef << ")!!!" << endl;
    assert(0);
  }

  int nTrueSyst = nCoef - nLinearVariations;
  int integralArraySize = TMath::Min(nFuncs, 9*101);
  nbins_ctau = nFuncs/((2*nTrueSyst+1)*(2*nLinearVariations+1));

//  cout << "Number of ctau bins: " << nbins_ctau << endl;
  for (int mp=0; mp<integralArraySize; mp++) Integral_T[mp] = dynamic_cast<const RooHistFunc*>(_funcList.at(mp))->analyticalIntegral(1000);
	for(int mp=nbins_ctau;mp<(9*101);mp++) Integral_T[mp] = Integral_T[nbins_ctau-1];

//	cout << "Finishing constructor, obtained the following integrals:" << endl;
//	for(int mp=0;mp<9*101;mp++) cout << "Integral[" << mp << "]: " << Integral_T[mp] << endl;
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
smoothAlgo(other.smoothAlgo),
nLinearVariations(other.nLinearVariations)
{
  nCoef = _coefList.getSize();
  nFuncs = _funcList.getSize();
  int nTrueSyst = nCoef - nLinearVariations;
  nbins_ctau = nFuncs/((2*nTrueSyst+1)*(2*nLinearVariations+1));
  for (int mp=0; mp<(9*101); mp++) Integral_T[mp] = (other.Integral_T)[mp];
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

  int nTrueSyst = nCoef - nLinearVariations;
  std::vector<Double_t> v1_array;
  std::vector<Double_t> v2_array;
  for (int ss=0; ss<(2*nLinearVariations+1); ss++){
    int iArray[2] ={ lowbin, highbin };
    int addIndex = (2*nTrueSyst+1)*nbins_ctau*ss;
    for (int iarr=0; iarr<2; iarr++) iArray[iarr] = iArray[iarr]+addIndex;

    Double_t v1Nom = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[0]))->getVal();
    Double_t v2Nom = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[1]))->getVal();
    Double_t v1Inst = v1Nom;
    Double_t v2Inst = v2Nom;

    for (int iSyst=1; iSyst<=nTrueSyst; iSyst++){
      RooAbsReal* coef = dynamic_cast<RooAbsReal*>(_coefList.at(iSyst-1));
      Double_t coefVal = coef->getVal();

      int iSystUp = (2*iSyst-1)*nbins_ctau;
      int iSystDn = 2*iSyst*nbins_ctau;

      Double_t v1_Up = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[0] + iSystUp))->getVal();
      Double_t v1_Dn = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[0] + iSystDn))->getVal();
      Double_t v2_Up = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[1] + iSystUp))->getVal();
      Double_t v2_Dn = dynamic_cast<const RooHistFunc*>(_funcList.at(iArray[1] + iSystDn))->getVal();

      v1Inst += interpolateVariation(coefVal, v1Nom, v1_Up, v1_Dn);
      v2Inst += interpolateVariation(coefVal, v2Nom, v2_Up, v2_Dn);
    }
    v1_array.push_back(v1Inst);
    v2_array.push_back(v2Inst);
  }

  Double_t v1_nominal=v1_array[0];
  Double_t v2_nominal=v2_array[0];
  Double_t v1=v1_nominal;
  Double_t v2=v2_nominal;
  for (int ss=0; ss<nLinearVariations; ss++){
    RooAbsReal* coef = dynamic_cast<RooAbsReal*>( _coefList.at( nTrueSyst + ss ) );
    Double_t coefVal = coef->getVal();

    int iLinear_Up = 2*ss+1;
    int iLinear_Dn = 2*ss+2;

    Double_t v1_variation_up = v1_array[iLinear_Up];
    Double_t v2_variation_up = v2_array[iLinear_Up];
    Double_t v1_variation_dn = v1_array[iLinear_Dn];
    Double_t v2_variation_dn = v2_array[iLinear_Dn];

    v1 += interpolateVariation(coefVal, v1_nominal, v1_variation_up, v1_variation_dn);
    v2 += interpolateVariation(coefVal, v2_nominal, v2_variation_up, v2_variation_dn);
  }
  v1_array.clear();
  v2_array.clear();

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

  int nTrueSyst = nCoef - nLinearVariations;
  std::vector<Double_t> v1_array;
  std::vector<Double_t> v2_array;
  for (int ss=0; ss<(2*nLinearVariations+1); ss++){
    int iArray[2] ={ lowbin, highbin };
    int addIndex = (2*nTrueSyst+1)*nbins_ctau*ss;
    for (int iarr=0; iarr<2; iarr++) iArray[iarr] = iArray[iarr]+addIndex;

    Double_t v1Nom = Integral_T[ iArray[0] ];
    Double_t v2Nom = Integral_T[ iArray[1] ];
    Double_t v1Inst = v1Nom;
    Double_t v2Inst = v2Nom;

    for (int iSyst=1; iSyst<=nTrueSyst; iSyst++){
      RooAbsReal* coef = dynamic_cast<RooAbsReal*>(_coefList.at(iSyst-1));
      Double_t coefVal = coef->getVal();

      int iSystUp = (2*iSyst-1)*nbins_ctau;
      int iSystDn = 2*iSyst*nbins_ctau;

//      cout << "Production array: " << ss << "\tSystematic: " << coef->GetName() << " (" << iSyst << ", " << coefVal << ")" << endl;

      Double_t v1_Up = Integral_T[ iArray[0] + iSystUp ];
      Double_t v1_Dn = Integral_T[ iArray[0] + iSystDn ];
      Double_t v2_Up = Integral_T[ iArray[1] + iSystUp ];
      Double_t v2_Dn = Integral_T[ iArray[1] + iSystDn ];

      v1Inst += interpolateVariation(coefVal, v1Nom, v1_Up, v1_Dn);
      v2Inst += interpolateVariation(coefVal, v2Nom, v2_Up, v2_Dn);
    }
    v1_array.push_back(v1Inst);
    v2_array.push_back(v2Inst);
  }

  Double_t v1_nominal=v1_array[0];
  Double_t v2_nominal=v2_array[0];
  Double_t v1=v1_nominal;
  Double_t v2=v2_nominal;
  for (int ss=0; ss<nLinearVariations; ss++){
    RooAbsReal* coef = dynamic_cast<RooAbsReal*>(_coefList.at(nTrueSyst + ss));
    Double_t coefVal = coef->getVal();

    int iLinear_Up = 2*ss+1;
    int iLinear_Dn = 2*ss+2;

//    cout << "Production: " << ss << " -> " << coef->GetName() << " (" << iLinear_Up << ", " << iLinear_Dn << ", " << coefVal << ")" << endl;

    Double_t v1_variation_up = v1_array[iLinear_Up];
    Double_t v2_variation_up = v2_array[iLinear_Up];
    Double_t v1_variation_dn = v1_array[iLinear_Dn];
    Double_t v2_variation_dn = v2_array[iLinear_Dn];

    v1 += interpolateVariation(coefVal, v1_nominal, v1_variation_up, v1_variation_dn);
    v2 += interpolateVariation(coefVal, v2_nominal, v2_variation_up, v2_variation_dn);
  }
  v1_array.clear();
  v2_array.clear();

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

