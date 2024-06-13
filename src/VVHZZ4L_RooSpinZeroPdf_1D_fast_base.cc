#include <HiggsAnalysis/CombinedLimit/interface/VVHZZ4L_RooSpinZeroPdf_1D_fast_base.h>
#include <HiggsAnalysis/CombinedLimit/interface/FastTemplateFunc.h>
#include <cmath>
#include <cassert>
#include "TMath.h"
#include "Riostream.h"

using namespace TMath;
using namespace std;
using namespace RooFit;


VVHZZ4L_RooSpinZeroPdf_1D_fast_base::VVHZZ4L_RooSpinZeroPdf_1D_fast_base() :
RooAbsPdf(),
obsList("obsList", "List of pdf observables", this),
coefList("coefList", "List of pdf components", this)
{}

VVHZZ4L_RooSpinZeroPdf_1D_fast_base::VVHZZ4L_RooSpinZeroPdf_1D_fast_base(
  const char *name, const char *title,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
  ) :
  RooAbsPdf(name, title),
  obsList("obsList", "List of pdf observables", this),
  coefList("coefList", "List of pdf components", this)
{
  TIterator* coefIter = inObsList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "VVHZZ4L_RooSpinZeroPdf_1D_fast_base(" << GetName() << ") observable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    obsList.add(*coef);
  }
  delete coefIter;

  coefIter = inCoefList.createIterator();
  coef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (dynamic_cast<FastTemplateFunc_f*>(coef)==0 && dynamic_cast<FastTemplateFunc_d*>(coef)==0){
      coutE(InputArguments) << "VVHZZ4L_RooSpinZeroPdf_1D_fast_base(" << GetName() << ") component " << coef->GetName() << " is not of type FastTemplateFunc_f" << endl;
      assert(0);
    }
    coefList.add(*coef);
  }
  delete coefIter;
}


VVHZZ4L_RooSpinZeroPdf_1D_fast_base::VVHZZ4L_RooSpinZeroPdf_1D_fast_base(
  const VVHZZ4L_RooSpinZeroPdf_1D_fast_base& other, const char* name
  ) :
  RooAbsPdf(other, name),
  obsList("obsList", this, other.obsList),
  coefList("coefList", this, other.coefList)
{}


Float_t VVHZZ4L_RooSpinZeroPdf_1D_fast_base::interpolateFcn(Int_t code, const char* rangeName) const{
  if (!isAnomalousCouplingValueValid()) return 0;

  DefaultAccumulator<Float_t> value = 0;

  vector<Float_t> coefs; coefs.reserve(5);
  coefs.push_back((Float_t)pow(a1Val(), 4)); // a1**4
  coefs.push_back((Float_t)pow(a1Val(), 3)*ai1Val()); // a1**3 x ai1
  coefs.push_back((Float_t)pow(a1Val()*ai1Val(), 2)); // a1**2 x ai1**2
  coefs.push_back((Float_t)a1Val()*pow(ai1Val(), 3)); // a1 x ai1**3
  coefs.push_back((Float_t)pow(ai1Val(), 4)); // ai1**4

  if (coefList.getSize() != (Int_t)coefs.size()){
    cerr << "VVHZZ4L_RooSpinZeroPdf_1D_fast_base::interpolateFcn: coefList.getSize()=" << coefList.getSize() << " != coefs.size()=" << coefs.size() << endl;
    assert(0);
  }

  if (code==0){
    for (int ic=0; ic<coefList.getSize(); ic++) value += (Float_t)((dynamic_cast<const RooAbsReal*>(coefList.at(ic))->getVal())*coefs.at(ic));
  }
  else{
    for (int ic=0; ic<coefList.getSize(); ic++) value += (Float_t)((dynamic_cast<const RooAbsReal*>(coefList.at(ic))->analyticalIntegral(code, rangeName))*coefs.at(ic));
  }

  Float_t result = value.sum();
  return result;
}
Double_t VVHZZ4L_RooSpinZeroPdf_1D_fast_base::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-100;
  return value;
}
Int_t VVHZZ4L_RooSpinZeroPdf_1D_fast_base::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  Int_t code = 0;
  if (coefList.getSize()>0){
    if (dynamic_cast<const FastTemplateFunc_f*>(coefList.at(0))!=0) code = dynamic_cast<const FastTemplateFunc_f*>(coefList.at(0))->getAnalyticalIntegral(allVars, analVars, rangeName);
    else if (dynamic_cast<const FastTemplateFunc_d*>(coefList.at(0))!=0) code = dynamic_cast<const FastTemplateFunc_d*>(coefList.at(0))->getAnalyticalIntegral(allVars, analVars, rangeName);
  }
  return code;
}
Double_t VVHZZ4L_RooSpinZeroPdf_1D_fast_base::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) return 1e-100;
  return value;
}


ClassImp(VVHZZ4L_RooSpinZeroPdf_1D_fast_base)

