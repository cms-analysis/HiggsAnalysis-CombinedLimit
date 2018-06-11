#include "Riostream.h" 
#include "../interface/HZZ4L_RooCTauPdf_1D.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;
using namespace RooFit;
using namespace std;

ClassImp(HZZ4L_RooCTauPdf_1D) 

HZZ4L_RooCTauPdf_1D::HZZ4L_RooCTauPdf_1D(
	  const char *name,
	  const char *title,
	  RooAbsReal& _kd,
	  RooAbsReal& _ctau,
	  const RooArgList& inCoefList,

	  double _ctau_min,
	  double _ctau_max
	  ):
RooAbsPdf(name, title),
kd("kd", "kd", this, _kd),
ctau("ctau", "ctau", this, _ctau),
_coefList("coefList", "List of funcficients", this),
ctau_min(_ctau_min),
ctau_max(_ctau_max)
{
	TIterator* coefIter = inCoefList.createIterator();
	RooAbsArg* func;
	while ((func = (RooAbsArg*)coefIter->Next())) {
		if (!dynamic_cast<RooAbsReal*>(func)) {
			coutE(InputArguments) << "ERROR: :HZZ4L_RooCTauPdf_1D(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
			assert(0);
		}
		_coefList.add(*func);
	}
	delete coefIter;

	nbins_ctau = _coefList.getSize();
	for(int mp=0;mp<(nbins_ctau>101 ? 101 : nbins_ctau);mp++) Integral_T[mp] = dynamic_cast<const RooHistFunc*>(_coefList.at(mp))->analyticalIntegral(1000);
	for(int mp=nbins_ctau;mp<101;mp++) Integral_T[mp] = Integral_T[nbins_ctau-1];

//	cout << "Finishing constructor, obtained the following integrals:" << endl;
//	for(int mp=0;mp<101;mp++) cout << mp << '\t' << Integral_T[mp] << endl;
}

HZZ4L_RooCTauPdf_1D::HZZ4L_RooCTauPdf_1D(
	const HZZ4L_RooCTauPdf_1D& other,
	const char* name
	):
RooAbsPdf(other, name),
kd("kd", this, other.kd),
ctau("ctau", this, other.ctau),
_coefList("coefList", this, other._coefList),
ctau_min(other.ctau_min),
ctau_max(other.ctau_max)
{
	nbins_ctau = _coefList.getSize();
	for(int mp=0;mp<101;mp++) Integral_T[mp] = (other.Integral_T)[mp];
}

int HZZ4L_RooCTauPdf_1D::findNeighborBins() const{
	double gridsize = (ctau_max-ctau_min)/(nbins_ctau-1);
	int bincode = (int) (ctau-ctau_min)/gridsize;
	return bincode;
}
Double_t HZZ4L_RooCTauPdf_1D::interpolateBin() const{
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

	Double_t v1 = dynamic_cast<const RooHistFunc*>(_coefList.at(lowbin))->getVal();
	Double_t v2 = dynamic_cast<const RooHistFunc*>(_coefList.at(highbin))->getVal();

	Double_t fdistance = ( ( (ctau-ctau_min)/gridsize ) - ((double) lowbin) ) / (highbin-lowbin);

	result = (1.-fdistance)*v1 + fdistance*v2;
	return result;
}
Double_t HZZ4L_RooCTauPdf_1D::interpolateIntegral() const{
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

	Double_t v1 = Integral_T[lowbin];
	Double_t v2 = Integral_T[highbin];

	Double_t fdistance = ( ( (ctau-ctau_min)/gridsize ) - ((double) lowbin) ) / (highbin-lowbin);

	result = (1.-fdistance)*v1 + fdistance*v2;
	return result;
}


Double_t HZZ4L_RooCTauPdf_1D::evaluate() const
{
	double value = interpolateBin();
	if (value <= 0) return 1.0e-15;
	return value;
}
Int_t HZZ4L_RooCTauPdf_1D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg()))) return 1;
  return 0 ;
}
Double_t HZZ4L_RooCTauPdf_1D::analyticalIntegral(Int_t code, const char* rangeName) const
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

