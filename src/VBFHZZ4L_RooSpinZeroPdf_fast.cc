#include <HiggsAnalysis/CombinedLimit/interface/VBFHZZ4L_RooSpinZeroPdf_fast.h>
#include <HiggsAnalysis/CombinedLimit/interface/FastTemplateFunc.h>
#include <cmath>
#include <cassert>
#include "TMath.h"
#include "Riostream.h"

using namespace TMath;
using namespace std;
using namespace RooFit;


VBFHZZ4L_RooSpinZeroPdf_fast::VBFHZZ4L_RooSpinZeroPdf_fast() :
RooAbsPdf(),
a1("a1", "a1", this),
ai1("ai1", "ai1", this),
obsList("obsList", "List of pdf observables", this),
coefList("coefList", "List of pdf components", this)
{}

VBFHZZ4L_RooSpinZeroPdf_fast::VBFHZZ4L_RooSpinZeroPdf_fast(
  const char *name, const char *title,
  RooAbsReal& in_a1,
  RooAbsReal& in_ai1,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
  ) :
  RooAbsPdf(name, title),
  a1("a1", "a1", this, in_a1),
  ai1("ai1", "ai1", this, in_ai1),
  obsList("obsList", "List of pdf observables", this),
  coefList("coefList", "List of pdf components", this)
{
  TIterator* coefIter = inObsList.createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "VBFHZZ4L_RooSpinZeroPdf_fast(" << GetName() << ") observable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    obsList.add(*coef);
  }
  delete coefIter;

  coefIter = inCoefList.createIterator();
  coef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (dynamic_cast<FastTemplateFunc_f*>(coef)==0 && dynamic_cast<FastTemplateFunc_d*>(coef)==0){
      coutE(InputArguments) << "VBFHZZ4L_RooSpinZeroPdf_fast(" << GetName() << ") component " << coef->GetName() << " is not of type FastTemplateFunc_f" << endl;
      assert(0);
    }
    coefList.add(*coef);
  }
  delete coefIter;
}


VBFHZZ4L_RooSpinZeroPdf_fast::VBFHZZ4L_RooSpinZeroPdf_fast(
  const VBFHZZ4L_RooSpinZeroPdf_fast& other, const char* name
  ) :
  RooAbsPdf(other, name),
  a1("a1", this, other.a1),
  ai1("ai1", this, other.ai1),
  obsList("obsList", this, other.obsList),
  coefList("coefList", this, other.coefList)
{}


Float_t VBFHZZ4L_RooSpinZeroPdf_fast::interpolateFcn(Int_t code, const char* rangeName) const{
  DefaultAccumulator<Float_t> value = 0;

  vector<Float_t> coefs; coefs.reserve(5);
  coefs.push_back((Float_t)pow(a1, 4)); // a1**4
  coefs.push_back((Float_t)pow(a1, 3)*ai1); // a1**3 x ai1
  coefs.push_back((Float_t)pow(a1*ai1, 2)); // a1**2 x ai1**2
  coefs.push_back((Float_t)a1*pow(ai1, 3)); // a1 x ai1**3
  coefs.push_back((Float_t)pow(ai1, 4)); // ai1**4

  if (coefList.getSize() != (Int_t)coefs.size()){
    cerr << "VBFHZZ4L_RooSpinZeroPdf_fast::interpolateFcn: coefList.getSize()=" << coefList.getSize() << " != coefs.size()=" << coefs.size() << endl;
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
Double_t VBFHZZ4L_RooSpinZeroPdf_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-100;
  return value;
}
Int_t VBFHZZ4L_RooSpinZeroPdf_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  Int_t code = 0;
  if (coefList.getSize()>0){
    if (dynamic_cast<const FastTemplateFunc_f*>(coefList.at(0))!=0) code = dynamic_cast<const FastTemplateFunc_f*>(coefList.at(0))->getAnalyticalIntegral(allVars, analVars, rangeName);
    else if (dynamic_cast<const FastTemplateFunc_d*>(coefList.at(0))!=0) code = dynamic_cast<const FastTemplateFunc_d*>(coefList.at(0))->getAnalyticalIntegral(allVars, analVars, rangeName);
  }
  return code;
}
Double_t VBFHZZ4L_RooSpinZeroPdf_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) return 1e-100;
  return value;
}


ClassImp(VBFHZZ4L_RooSpinZeroPdf_fast)

