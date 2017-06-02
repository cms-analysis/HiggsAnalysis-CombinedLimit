//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "HiggsAnalysis/CombinedLimit/interface/RooParametricShapeBinPdf.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"
#include "RooTFnBinding.h"

using namespace std;
using namespace RooFit;

ClassImp(RooParametricShapeBinPdf)
//---------------------------------------------------------------------------
/*
RooParametricShapeBinPdf::RooParametricShapeBinPdf(const char *name, const char *title, const char *formula, 
			       RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape ) : RooAbsPdf(name, title), 
  x("x", "x Observable", this, _x),
  pars("pars","pars",this),
  mypdf("mypdf","mypdf",this),
  xBins(0),
  xMax(0),
  xMin(0),
  relTol(1E-12),
  absTol(1E-12),
  nPars(0)
{
  memset(&xArray, 0, sizeof(xArray));
  TIterator *varIter=_pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*)varIter->Next()) ){
	pars.add(*fVar);
  }
  setTH1Binning(_shape);
  myfunc = new TF1("myfunc",formula,xMin,xMax);
  mypdf.setArg((RooAbsPdf&)*bindFunction(myfunc,(RooAbsReal&)x.arg(),pars));
  //RooListProxy obs;
  //obs.add(x.arg());
  //mypdf.setArg((RooAbsPdf&)* new RooTFnBinding("mypdf","mypdf",myfunc,obs,pars));
  nPars = pars.getSize();
}
*/
//---------------------------------------------------------------------------
RooParametricShapeBinPdf::RooParametricShapeBinPdf(const char *name, const char *title, RooAbsReal& _pdf, 
			       RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape ) : RooAbsPdf(name, title), 
  x("x", "x Observable", this, _x),
  pars("pars","pars",this),
  mypdf("mypdf","mypdf", this, _pdf),
  xBins(0),
  xMax(0),
  xMin(0),
  relTol(1E-12),
  absTol(1E-12),
  nPars(0)
{
  memset(&xArray, 0, sizeof(xArray));
  TIterator *varIter=_pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*)varIter->Next()) ){
	pars.add(*fVar);
  }
  setTH1Binning(_shape);
  RooListProxy obs;
  obs.add(x.arg());
  myfunc = _pdf.asTF(obs,pars);
  nPars = pars.getSize();
}
//---------------------------------------------------------------------------
RooParametricShapeBinPdf::RooParametricShapeBinPdf(const RooParametricShapeBinPdf& other, const char* name) : RooAbsPdf(other, name), 
   x("x", this, other.x),
   pars("pars",this,RooListProxy()),
   mypdf("mypdf",this,other.mypdf),
   xBins(other.xBins),
   xMax(other.xMax),
   xMin(other.xMin),
   relTol(other.relTol),
   absTol(other.absTol),
   nPars(other.nPars)
{
  //memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
  
  TIterator *varIter=other.pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*) varIter->Next()) ){
	pars.add(*fVar);
  }
  myfunc = new TF1(*(other.myfunc));
}
//---------------------------------------------------------------------------
void RooParametricShapeBinPdf::setTH1Binning(const TH1 &_Hnominal){
  xBins = _Hnominal.GetXaxis()->GetNbins();
  xMin = _Hnominal.GetXaxis()->GetBinLowEdge(1);
  xMax = _Hnominal.GetXaxis()->GetBinUpEdge(xBins);
  memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal.GetXaxis()->GetBinLowEdge(i+1);
  }
}
//---------------------------------------------------------------------------
void RooParametricShapeBinPdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooParametricShapeBinPdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
/// Return the parameteric p.d.f
RooAbsPdf* RooParametricShapeBinPdf::getPdf() const {
  return mypdf ? ((RooAbsPdf*)mypdf.absArg()) : 0 ;
}
//---------------------------------------------------------------------------
Double_t RooParametricShapeBinPdf::evaluate() const
{
  Double_t integral = 0.0;
  Int_t iBin;
  for(iBin=0; iBin<xBins; iBin++) {  
    if (x>=xArray[iBin] && x < xArray[iBin+1] ) break;
  }
  
  if(iBin < 0 || iBin >= xBins) {
    //cout << "in bin " << iBin << " which is outside of range" << endl;
    return 0.0;
  }
  
  Double_t xLow = xArray[iBin];
  Double_t xHigh = xArray[iBin+1];
    
  // define the function to be integrated numerically  
  ROOT::Math::WrappedTF1 func(*myfunc);
  double *params = myfunc->GetParameters();  
  TIterator *varIter=pars.createIterator(); 
  RooAbsReal *fVar;
  int iPar = 0;
  while ( (fVar = (RooAbsReal*) varIter->Next()) ){
    params[iPar] = fVar->getVal();
    iPar+=1;
  }
  myfunc->SetParameters(params);
  func.SetParameters(params);

  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
  ig.SetFunction(func,false);
  
  integral = ig.Integral(xLow,xHigh) / (xHigh-xLow); //return integral as a density 
  //Double_t total_integral = ig.Integral(xMin,xMax);

  if (integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooParametricShapeBinPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooParametricShapeBinPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t xRangeMin = x.min(rangeName); Double_t xRangeMax = x.max(rangeName);

   Double_t integral = 0.0;
   
   // define the function to be integrated numerically  
   ROOT::Math::WrappedTF1 func(*myfunc);
   double *params = myfunc->GetParameters();
   TIterator *varIter=pars.createIterator(); 
   RooAbsReal *fVar;
   int iPar = 0;
   while ( (fVar = (RooAbsReal*) varIter->Next()) ){
     params[iPar] = fVar->getVal();
     iPar+=1;
   }
   myfunc->SetParameters(params);
   func.SetParameters(params);
   
   ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
   ig.SetFunction(func,false);
    
   if (code==1 && xRangeMin<=xMin && xRangeMax>=xMax){
     integral = ig.Integral(xMin,xMax);
     return integral;
   }
   else if(code==1) {     
     //Double_t integral = ig.Integral(xRangeMin,xRangeMax);
     //return integral;
     for (Int_t iBin=0; iBin<xBins; iBin++){       
       Double_t xLow = xArray[iBin];
       Double_t xHigh = xArray[iBin+1];
       Double_t partial_integral = ig.Integral(xLow,xHigh);
       if (xLow>=xRangeMin && xHigh<=xRangeMax) {
     	 // Bin fully in the integration domain
     	 integral += partial_integral;
       } else if (xLow<xRangeMin && xHigh>xRangeMax) {
     	 // Domain is fully contained in this bin
     	 integral += (xRangeMax-xRangeMin)*partial_integral/(xHigh-xLow);
     	 // Exit here, this is the last bin to be processed by construction
     	 return integral;
       } else if (xLow<xRangeMin && xHigh<=xRangeMax && xHigh>xRangeMin) {
     	 // Lower domain boundary is in bin
     	 integral += (xHigh-xRangeMin)*partial_integral/(xHigh-xLow);
       } else if (xLow>=xRangeMin && xHigh>xRangeMax && xLow<xRangeMax) {
     	 // Upper domain boundary is in bin
     	 integral +=  (xRangeMax-xLow)*partial_integral/(xHigh-xLow);
     	 // Exit here, this is the last bin to be processed by construction
     	 return integral;
       }
     }
     return integral;
   } else {
     cout << "WARNING IN RooParametricShapeBinPdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (integral>0.0) {
     return integral;
   } else return 1.0;
}
// //---------------------------------------------------------------------------

