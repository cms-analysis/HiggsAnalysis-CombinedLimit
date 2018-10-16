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
#include "RooRealProxy.h"
#include "RooListProxy.h"

using namespace std;
using namespace RooFit;

ClassImp(RooParametricShapeBinPdf)
//---------------------------------------------------------------------------

RooParametricShapeBinPdf::RooParametricShapeBinPdf(const char *name, const char *title, RooAbsReal& _pdf, 
			       RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape ) : RooAbsPdf(name, title), 
  x("x", "x Observable", this, _x),
  pars("pars","pars",this),
  mypdf("mypdf","mypdf", this, _pdf),
  myintegrals("myintegrals","myintegrals",this),
  xBins(0),
  xMax(0),
  xMin(0)
{
  memset(&xArray, 0, sizeof(xArray));
  TIterator *varIter=_pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*)varIter->Next()) ){
    pars.add(*fVar);
  }
  setTH1Binning(_shape);
  RooAbsReal* myintegral;
  RooListProxy obs;
  obs.add(x.arg());
  for (Int_t iBin=0; iBin<xBins; iBin++){
    std::string rangeName  = Form("%s_%s_range_bin%d", GetName(), x.GetName(), iBin);
    if (!x.arg().hasRange(rangeName.c_str())) {
      RooRealVar x_rrv = dynamic_cast<const RooRealVar &>(x.arg());
      Double_t xLow = xArray[iBin];
      Double_t xHigh = xArray[iBin+1];
      x_rrv.setRange(rangeName.c_str(),xLow,xHigh);
    } 
    myintegral = getPdf()->createIntegral(obs,Range(rangeName.c_str()));
    myintegrals.add(*myintegral);
  }
  //last one is full integral (no range)
  myintegral = getPdf()->createIntegral(obs);
  myintegrals.add(*myintegral);
}
//---------------------------------------------------------------------------
RooParametricShapeBinPdf::RooParametricShapeBinPdf(const RooParametricShapeBinPdf& other, const char* name) : RooAbsPdf(other, name), 
   x("x", this, other.x),
   pars("pars",this,RooListProxy()),
   mypdf("mypdf",this,other.mypdf),
   myintegrals("myintegrals",this,RooListProxy()),
   xBins(other.xBins),
   xMax(other.xMax),
   xMin(other.xMin)
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

  TIterator *intIter=other.myintegrals.createIterator(); 
  RooAbsReal *fInt;
  while ( (fInt = (RooAbsReal*) intIter->Next()) ){
    myintegrals.add(*fInt);
  }
  
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
/// Return the parameteric p.d.f
RooAbsPdf* RooParametricShapeBinPdf::getPdf() const {
  return mypdf ? ((RooAbsPdf*)mypdf.absArg()) : 0 ;
}
//---------------------------------------------------------------------------
/// Return the bin-by-bin integrals
RooAbsReal* RooParametricShapeBinPdf::getIntegral(int index) const {
  return myintegrals.at(index) ? ((RooAbsReal*)myintegrals.at(index)) : 0;
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

  // check again if x variable has the right range already defined 
  // needed when combining multiple workspaces, and taking variable x from only one of them!
  std::string rangeName  = Form("%s_%s_range_bin%d", GetName(), x.GetName(), iBin);
  if (!x.arg().hasRange(rangeName.c_str())) {
    RooRealVar x_rrv = dynamic_cast<const RooRealVar &>(x.arg());
    Double_t xLow = xArray[iBin];
    Double_t xHigh = xArray[iBin+1];
    x_rrv.setRange(rangeName.c_str(),xLow,xHigh);
  } 
  integral = getIntegral(iBin)->getVal() / (xHigh-xLow);
  
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
  
  RooListProxy obs;
  obs.add(x.arg());
  
  if (code==1 && xRangeMin<=xMin && xRangeMax>=xMax){
    integral = getIntegral(xBins)->getVal();
    return integral;
  }
  else if(code==1) {     
    RooAbsReal* myintegral = getPdf()->createIntegral(obs,Range(rangeName));
    integral = myintegral->getVal();
    return integral;
   } else {
    cout << "WARNING IN RooParametricShapeBinPdf: integration code is not correct" << endl;
    cout << "                           what are you integrating on?" << endl;
    return 1.0;
  }
}
// //---------------------------------------------------------------------------

