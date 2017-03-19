#include <HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_2D_fast.h>
#include <HiggsAnalysis/CombinedLimit/interface/FastTemplateFunc.h>
#include <cmath>
#include <cassert>
#include "TMath.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "Riostream.h"

using namespace TMath;
using namespace std;
using namespace RooFit;


HZZ4L_RooSpinZeroPdf_2D_fast::HZZ4L_RooSpinZeroPdf_2D_fast() :
RooAbsPdf(),
fai1("fai1", "fai1", this),
fai2("fai2", "fai2", this),
phi1("phi1", "phi1", this),
phi2("phi2", "phi2", this),
obsList("obsList", "List of pdf observables", this),
coefList("coefList", "List of pdf components", this),
fullIntCode(1)
{}

HZZ4L_RooSpinZeroPdf_2D_fast::HZZ4L_RooSpinZeroPdf_2D_fast(
  const char *name, const char *title,
  RooAbsReal& in_fai1,
  RooAbsReal& in_fai2,
  RooAbsReal& in_phi1,
  RooAbsReal& in_phi2,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
  ) :
  RooAbsPdf(name, title),
  fai1("fai1", "fai1", this, in_fai1),
  fai2("fai2", "fai2", this, in_fai2),
  phi1("phi1", "phi1", this, in_phi1),
  phi2("phi2", "phi2", this, in_phi2),
  obsList("obsList", "List of pdf observables", this),
  coefList("coefList", "List of pdf components", this),
  fullIntCode(1)
{
  const Int_t code_prime[3]={ 2, 3, 5 };
  TIterator* coefIter = inObsList.createIterator();
  RooAbsArg* coef;
  RooArgSet intSet;
  unsigned int icoef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<RooAbsReal*>(coef)){
      coutE(InputArguments) << "HZZ4L_RooSpinZeroPdf_2D_fast(" << GetName() << ") observable " << coef->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    obsList.add(*coef);
    intSet.add(*coef);
    if (icoef<3) fullIntCode*=code_prime[icoef];
    else{
      fullIntCode=1;
      if (icoef==3) cout << "HZZ4L_RooSpinZeroPdf_2D_fast(" << GetName() << ") number of observables is grater than 3. Full integration will switch to analytical." << endl;
    }
    icoef++;
  }
  delete coefIter;

  coefIter = inCoefList.createIterator();
  coef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    if (!dynamic_cast<FastTemplateFunc_f*>(coef)){
      coutE(InputArguments) << "HZZ4L_RooSpinZeroPdf_2D_fast(" << GetName() << ") component " << coef->GetName() << " is not of type FastTemplateFunc_f" << endl;
      assert(0);
    }
    coefList.add(*coef);
    fullIntegral.push_back((Float_t)dynamic_cast<const FastTemplateFunc_f*>(coef)->analyticalIntegral(fullIntCode));
  }
  delete coefIter;
}


HZZ4L_RooSpinZeroPdf_2D_fast::HZZ4L_RooSpinZeroPdf_2D_fast(
  const HZZ4L_RooSpinZeroPdf_2D_fast& other, const char* name
  ) :
  RooAbsPdf(other, name),
  fai1("fai1", this, other.fai1),
  fai2("fai2", this, other.fai2),
  phi1("phi1", this, other.phi1),
  phi2("phi2", this, other.phi2),
  obsList("obsList", this, other.obsList),
  coefList("coefList", this, other.coefList),
  fullIntegral(other.fullIntegral),
  fullIntCode(other.fullIntCode)
{}


Float_t HZZ4L_RooSpinZeroPdf_2D_fast::interpolateFcn(Int_t code, const char* /*rangeName*/) const{
  Float_t fa1 = (1.-fabs(fai1) - fabs(fai2)); if (fa1<0.) return 0;

  DefaultAccumulator<Float_t> value = 0;
  Float_t sgn_fai1 = (fai1<0. ? 1. : -1.);
  Float_t sgn_fai2 = (fai2<0. ? 1. : -1.);
  Float_t absfai1 = fabs(fai1);
  Float_t absfai2 = fabs(fai2);
  vector<Float_t> coefs; coefs.reserve(9);
  coefs.push_back((Float_t)fa1);
  coefs.push_back((Float_t)absfai1);
  coefs.push_back((Float_t)absfai2);
  coefs.push_back((Float_t)sgn_fai1*sqrt(fa1*absfai1)*cos(phi1));
  coefs.push_back((Float_t)sgn_fai1*sqrt(fa1*absfai2)*cos(phi2));
  coefs.push_back((Float_t)sgn_fai1*sgn_fai2*sqrt(absfai1*absfai2)*cos(phi2-phi1));
  coefs.push_back((Float_t)sgn_fai1*sqrt(fa1*absfai1)*sin(phi1));
  coefs.push_back((Float_t)sgn_fai1*sqrt(fa1*absfai2)*sin(phi2));
  coefs.push_back((Float_t)sgn_fai1*sgn_fai2*sqrt(absfai1*absfai2)*sin(phi2-phi1));
  if (code==0){
    if (coefList.getSize() != (Int_t)coefs.size()){
      cerr << "HZZ4L_RooSpinZeroPdf_2D_fast::interpolateFcn: coefList.getSize() (=" << coefList.getSize() << ")!=coefs.size() (=" << coefs.size() << ")" << endl;
      assert(0);
    }
    for (int ic=0; ic<coefList.getSize(); ic++) value += (Float_t)((dynamic_cast<const RooAbsReal*>(coefList.at(ic))->getVal())*coefs.at(ic));
  }
  else if (code==fullIntCode){
    if (fullIntegral.size() != coefs.size()){
      cerr << "HZZ4L_RooSpinZeroPdf_2D_fast::interpolateFcn: fullIntegral.size() (=" << fullIntegral.size() << ")!=coefs.size() (=" << coefs.size() << ")" << endl;
      assert(0);
    }
    for (unsigned int ic=0; ic<fullIntegral.size(); ic++) value += (Float_t)(fullIntegral.at(ic)*coefs.at(ic));
  }

  Float_t result = value.sum();
  return result;
}
Double_t HZZ4L_RooSpinZeroPdf_2D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (value<=0.) return 1e-15;
  return value;
}
Int_t HZZ4L_RooSpinZeroPdf_2D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = 1;
  const Int_t code_prime[3]={ 2, 3, 5 };
  for (int ic=0; ic<obsList.getSize(); ic++){
    if (ic>=3){ code=0; break; }
    if (dynamic_cast<RooRealVar*>(obsList.at(ic))!=0){
      if (matchArgs(allVars, analVars, RooArgSet(*(obsList.at(ic))))) code*=code_prime[ic];
    }
  }
  if (code==1) code=0;
  return code;
}
Double_t HZZ4L_RooSpinZeroPdf_2D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (value<=0.) return 1e-10;
  return value;
}


ClassImp(HZZ4L_RooSpinZeroPdf_2D_fast)

