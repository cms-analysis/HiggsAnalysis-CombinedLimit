/***************************************************************************** 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h" 

#include <math.h> 
#include "TMath.h" 
#include "RooFormulaVar.h"
#include "RooAbsReal.h"
#include "RooFit.h"

#include "TFile.h"

//using namespace RooFit ;

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
  TIterator *varIter=_pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*)varIter->Next()) ){
	pars.add(*fVar);
  }
  if ( pars.getSize() != _shape.GetNbinsX() ){
	std::cout << " Warning, number of parameters not equal to number of bins in shape histogram! " << std::endl;
  }
  initializeBins(_shape);
//  initializeNorm();
  cval = -1; 
} 

//_____________________________________________________________________________
RooParametricHist::RooParametricHist(const RooParametricHist& other, const char* name) :
 RooAbsPdf(other, name),x("observable",this,other.x),pars("_pars",this,RooListProxy())
{
  
  N_bins = other.N_bins;
  //sum    = other.sum;
  TIterator *varIter=other.pars.createIterator(); 
  RooAbsReal *fVar;
  while ( (fVar = (RooAbsReal*) varIter->Next()) ){
	pars.add(*fVar);
  }

  for(int i=0; i<=N_bins; ++i) {
     bins.push_back(other.bins[i]);
     if (i<N_bins) widths.push_back(other.widths[i]);
  }

  cval = other.cval; 

}
void RooParametricHist::initializeBins(const TH1 &shape) const {

  ///std::vector<double> bins;
  //std::vector<double> widths;
  N_bins = shape.GetNbinsX();
  for(int i=1; i<=N_bins+1; ++i) {
     bins.push_back(shape.GetBinLowEdge(i));
     if (i<=N_bins) widths.push_back(shape.GetBinWidth(i));
  }
}

double RooParametricHist::getFullSum() const {
    double sum=0;
    TIterator *varIter=pars.createIterator(); 
    RooAbsReal *fVar;
    while ( (fVar = (RooAbsReal*) varIter->Next()) ){
      sum+=fVar->getVal();	
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


Double_t RooParametricHist::evaluate() const 
{ 
  int bin_i;
  if (x < bins[0]) 	   	return 0;    // should set to 0 instead?  
  else if (x >= bins[N_bins])   return 0;

  else {
    for(bin_i=0; bin_i<N_bins; bin_i++) {   // faster way to loop through ?
      if (x>=bins[bin_i] && x < bins[bin_i+1] ) break;
    }
  }
  RooAbsReal *retVar = (RooAbsReal*)pars.at(bin_i);
  double ret = retVar->getVal() / widths[bin_i];
  cval=ret;
  return ret; 
}

