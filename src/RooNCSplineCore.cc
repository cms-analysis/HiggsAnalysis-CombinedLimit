#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineCore.h" 
#include <cmath>
#include "TMath.h"
#include "TIterator.h"
#include "Riostream.h"

using namespace TMath;
using namespace RooFit;
using namespace std;


ClassImp(RooNCSplineCore)

RooNCSplineCore::RooNCSplineCore() :
RooAbsReal(),
verbosity(RooNCSplineCore::kSilent),
useFloor(true), floorEval(0), floorInt(0),
rangeXmin(1), rangeXmax(-1),
theXVar("theXVar", "theXVar", this),
leafDepsList("leafDepsList", "leafDepsList", this)
{}

RooNCSplineCore::RooNCSplineCore(
  const char* name,
  const char* title
  ) :
  RooAbsReal(name, title),
  verbosity(RooNCSplineCore::kSilent),
  useFloor(true), floorEval(0), floorInt(0),
  rangeXmin(1), rangeXmax(-1),
  theXVar("theXVar", "theXVar", this),
  leafDepsList("leafDepsList", "leafDepsList", this)
{}

RooNCSplineCore::RooNCSplineCore(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<T>& inXList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooAbsReal(name, title),
  verbosity(RooNCSplineCore::kSilent),
  useFloor(inUseFloor), floorEval(inFloorEval), floorInt(inFloorInt),
  rangeXmin(1), rangeXmax(-1),
  theXVar("theXVar", "theXVar", this, inXVar),
  leafDepsList("leafDepsList", "leafDepsList", this),
  XList(inXList)
{}

RooNCSplineCore::RooNCSplineCore(
  const RooNCSplineCore& other,
  const char* name
  ) :
  RooAbsReal(other, name),
  verbosity(other.verbosity),
  useFloor(other.useFloor), floorEval(other.floorEval), floorInt(other.floorInt),
  rangeXmin(other.rangeXmin), rangeXmax(other.rangeXmax),
  theXVar("theXVar", this, other.theXVar),
  leafDepsList("leafDepsList", this, other.leafDepsList),
  XList(other.XList)
{}

void RooNCSplineCore::setVerbosity(VerbosityLevel flag){ verbosity = flag; }
void RooNCSplineCore::setEvalFloor(RooNCSplineCore::T val){ floorEval = val; }
void RooNCSplineCore::setIntFloor(RooNCSplineCore::T val){ floorInt = val; }
void RooNCSplineCore::doFloor(Bool_t flag){ useFloor = flag; }

void RooNCSplineCore::getBArray(const std::vector<RooNCSplineCore::T>& kappas, const vector<RooNCSplineCore::T>& fcnList, std::vector<RooNCSplineCore::T>& BArray)const{
  BArray.clear();
  int npoints=kappas.size();
  if (npoints!=(int)fcnList.size()){
    coutE(InputArguments) << "RooNCSplineCore::getBArray: Dim(kappas)=" << npoints << " != Dim(fcnList)=" << fcnList.size() << endl;
    assert(0);
  }
  if (npoints>1){
    BArray.push_back(3.*(fcnList.at(1)-fcnList.at(0)));
    for (int j=1; j<npoints-1; j++){
      RooNCSplineCore::T val_j = fcnList.at(j);
      RooNCSplineCore::T val_jpo = fcnList.at(j+1);
      RooNCSplineCore::T val_jmo = fcnList.at(j-1);
      RooNCSplineCore::T kappa_j = kappas.at(j);
      RooNCSplineCore::T kappa_jmo = kappas.at(j-1);
      RooNCSplineCore::T rsq = pow(kappa_j/kappa_jmo, 2);
      RooNCSplineCore::T Bval = val_jpo*rsq + val_j*(1.-rsq) - val_jmo;
      Bval *= 3.;
      BArray.push_back(Bval);
    }
    BArray.push_back(3.*(fcnList.at(npoints-1)-fcnList.at(npoints-2)));
  }
  else if (npoints==1) BArray.push_back(0);
}
void RooNCSplineCore::getAArray(const vector<RooNCSplineCore::T>& kappas, vector<vector<RooNCSplineCore::T>>& AArray)const{
  AArray.clear();
  Int_t npoints = kappas.size();
  for (int i=0; i<npoints; i++){
    vector<RooNCSplineCore::T> Ai(npoints, 0.);
    if (npoints==1) Ai[0]=1;
    else if (i==0){ Ai[0]=2; Ai[1]=kappas.at(1)/kappas.at(0); }
    else if (i==npoints-1){ Ai[npoints-2]=1; Ai[npoints-1]=2.*kappas.at(npoints-1)/kappas.at(npoints-2); }
    else{
      RooNCSplineCore::T kappa_j = kappas.at(i);
      RooNCSplineCore::T kappa_jmo = kappas.at(i-1);
      RooNCSplineCore::T kappa_jpo = kappas.at(i+1);

      Ai[i-1]=1;
      Ai[i]=2.*kappa_j/kappa_jmo*(1.+kappa_j/kappa_jmo);
      Ai[i+1]=kappa_j*kappa_jpo/pow(kappa_jmo, 2);
    }
    AArray.push_back(Ai);
  }
}

vector<RooNCSplineCore::T> RooNCSplineCore::getCoefficients(const TVector_t& S, const vector<RooNCSplineCore::T>& kappas, const vector<RooNCSplineCore::T>& fcnList, const Int_t& bin)const{
  DefaultAccumulator<RooNCSplineCore::T> A, B, C, D;
  vector<RooNCSplineCore::T> res;

  const int fcnsize = fcnList.size();
  if (fcnsize>bin){
    A=fcnList.at(bin);
    B=S[bin];
    if (fcnsize>(bin+1)){
      DefaultAccumulator<RooNCSplineCore::T> dFcn = fcnList.at(bin+1);
      dFcn -= A;

      C += RooNCSplineCore::T(3.)*dFcn;
      C -= 2.*B;
      C -= S[bin+1]*kappas.at(bin+1)/kappas.at(bin);

      D += RooNCSplineCore::T(-2.)*dFcn;
      D += B;
      D += S[bin+1]*kappas.at(bin+1)/kappas.at(bin);
    }
  }

  res.push_back(A);
  res.push_back(B);
  res.push_back(C);
  res.push_back(D);
  return res;
}
vector<vector<RooNCSplineCore::T>> RooNCSplineCore::getCoefficientsAlongDirection(const std::vector<RooNCSplineCore::T>& kappas, const TMatrix_t& Ainv, const vector<RooNCSplineCore::T>& fcnList, const Int_t pickBin)const{
  vector<RooNCSplineCore::T> BArray;
  getBArray(kappas, fcnList, BArray);

  Int_t npoints = BArray.size();
  TVector_t Btrans(npoints);
  for (int i=0; i<npoints; i++) Btrans[i]=BArray.at(i);
  TVector_t Strans = Ainv*Btrans;

  vector<vector<RooNCSplineCore::T>> coefs;
  for (Int_t bin=0; bin<(npoints>1 ? npoints-1 : 1); bin++){
    if (pickBin>=0 && bin!=pickBin) continue;
    vector<RooNCSplineCore::T> coef = getCoefficients(Strans, kappas, fcnList, bin);
    coefs.push_back(coef);
  }
  return coefs;
}

RooNCSplineCore::T RooNCSplineCore::evalSplineSegment(const std::vector<RooNCSplineCore::T>& coefs, const RooNCSplineCore::T& kappa, const RooNCSplineCore::T& tup, const RooNCSplineCore::T& tdn, Bool_t doIntegrate)const{
  DefaultAccumulator<RooNCSplineCore::T> res;
  for (unsigned int ic=0; ic<coefs.size(); ic++){
    if (doIntegrate) res += coefs.at(ic)*(pow(tup, (int)(ic+1))-pow(tdn, (int)(ic+1)))/((RooNCSplineCore::T)(ic+1));
    else res += coefs.at(ic)*pow(tup, (int)ic);
  }
  if (doIntegrate) res /= kappa;
  return res;
}

void RooNCSplineCore::getLeafDependents(RooRealProxy& proxy, RooArgSet& set){
  RooArgSet deps;
  proxy.absArg()->leafNodeServerList(&deps, 0, true);
  set.add(deps);
}
void RooNCSplineCore::addLeafDependents(RooArgSet& set){
  TIterator* iter = set.createIterator();
  RooAbsArg* absarg;
  while ((absarg = (RooAbsArg*)iter->Next())){ if (dynamic_cast<RooRealVar*>(absarg)) leafDepsList.add(*absarg); }
  delete iter;
}
