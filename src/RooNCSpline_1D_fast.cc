#include "HiggsAnalysis/CombinedLimit/interface/RooNCSpline_1D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;


ClassImp(RooNCSpline_1D_fast)

RooNCSpline_1D_fast::RooNCSpline_1D_fast() :
  RooNCSplineCore(),
  bcBeginX(RooNCSplineCore::bcNaturalSpline), bcEndX(RooNCSplineCore::bcNaturalSpline)
{}

RooNCSpline_1D_fast::RooNCSpline_1D_fast(
  const char* name,
  const char* title
  ) :
  RooNCSplineCore(name, title),
  bcBeginX(RooNCSplineCore::bcNaturalSpline), bcEndX(RooNCSplineCore::bcNaturalSpline)
{}

RooNCSpline_1D_fast::RooNCSpline_1D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<T>& inXList,
  const std::vector<T>& inFcnList,
  RooNCSplineCore::BoundaryCondition const bcBeginX_,
  RooNCSplineCore::BoundaryCondition const bcEndX_,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooNCSplineCore(name, title, inXVar, inXList, inUseFloor, inFloorEval, inFloorInt),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  FcnList(inFcnList)
{
  if (npointsX()>1){
    int npoints;

    vector<vector<RooNCSplineCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA, bcBeginX, bcEndX);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    Double_t det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSpline_1D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    coefficients = getCoefficientsAlongDirection(kappaX, xAinv, FcnList, bcBeginX, bcEndX, -1);
  }
  else assert(0);

  RooArgSet leafset;
  getLeafDependents(theXVar, leafset);
  addLeafDependents(leafset);

  emptyFcnList();
}

RooNCSpline_1D_fast::RooNCSpline_1D_fast(
  const RooNCSpline_1D_fast& other,
  const char* name
  ) :
  RooNCSplineCore(other, name),
  bcBeginX(other.bcBeginX), bcEndX(other.bcEndX),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  coefficients(other.coefficients)
{}

RooNCSplineCore::T RooNCSpline_1D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<RooNCSplineCore::T> res;

  if (verbosity==RooNCSplineCore::kVerbose) cout << "RooNCSpline_1D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl;

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1;
  RooNCSplineCore::T tx=0, txmin=0, txmax=0;
  if (code==0 || code%2!=0){ // Case to just compute the value at x
    if (!testRangeValidity(theXVar)) return 0;
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Case to integrate along x
    RooNCSplineCore::T coordmin = theXVar.min(rangeName); cropValueForRange(coordmin);
    RooNCSplineCore::T coordmax = theXVar.max(rangeName); cropValueForRange(coordmax);
    xbinmin = getWhichBin(coordmin, 0);
    txmin = getTVar(kappaX, coordmin, xbinmin, 0);
    xbinmax = getWhichBin(coordmax, 0);
    txmax = getTVar(kappaX, coordmax, xbinmax, 0);
  }

  int nxbins = (int)coefficients.size();
  for (int ix=0; ix<nxbins; ix++){
    if (
      (xbin>=0 && ix!=xbin)
      ||
      (xbinmin>=0 && xbinmax>=xbinmin && !(xbinmin<=ix && ix<=xbinmax))
      ) continue;

    RooNCSplineCore::T txlow=0, txhigh=1;
    if (code>0 && code%2==0){
      if (ix==xbinmin) txlow=txmin;
      if (ix==xbinmax) txhigh=txmax;
    }
    else txhigh=tx;

    // Get the x coefficients at bin ix and evaluate value of spline at x
    res += evalSplineSegment(coefficients.at(ix), kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void RooNCSpline_1D_fast::getKappas(vector<RooNCSplineCore::T>& kappas, const Int_t /*whichDirection*/){
  kappas.clear();
  RooNCSplineCore::T kappa=1;

  Int_t npoints;
  vector<RooNCSplineCore::T> const* coord;
  npoints=npointsX();
  coord=&XList;

  for (Int_t j=0; j<npoints-1; j++){
    RooNCSplineCore::T val_j = coord->at(j);
    RooNCSplineCore::T val_jpo = coord->at(j+1);
    RooNCSplineCore::T val_diff = (val_jpo-val_j);
    if (fabs(val_diff)>RooNCSplineCore::T(0)) kappa = 1./val_diff;
    else kappa = 0;
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t RooNCSpline_1D_fast::getWhichBin(const RooNCSplineCore::T& val, const Int_t /*whichDirection*/)const{
  Int_t bin=-1;
  RooNCSplineCore::T valj, valjpo;
  Int_t npoints;
  vector<RooNCSplineCore::T> const* coord;
  coord=&XList;
  npoints=npointsX();

  if (npoints<=1) bin=0;
  else{
    valjpo = coord->at(0);
    for (Int_t j=0; j<npoints-1; j++){
      valj = coord->at(j);
      valjpo = coord->at(j+1);
      if (val<valjpo && val>=valj){ bin=j; break; }
    }
    if (bin==-1 && val>=valjpo) bin=npoints-2;
    else if (bin==-1) bin=0;
  }

  return bin;
}
RooNCSplineCore::T RooNCSpline_1D_fast::getTVar(const vector<RooNCSplineCore::T>& kappas, const RooNCSplineCore::T& val, const Int_t& bin, const Int_t /*whichDirection*/)const{
  const RooNCSplineCore::T& K=kappas.at(bin);
  return (val-XList.at(bin))*K;
}

Double_t RooNCSpline_1D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=RooNCSplineCore::kError) coutE(Eval) << "RooNCSpline_1D_fast ERROR::RooNCSpline_1D_fast(" << GetName() << ") evaluation returned " << value << " at x = " << theXVar << endl;
    value = floorEval;
  }
  if (verbosity==RooNCSplineCore::kVerbose){
    cout << "RooNCSpline_1D_fast(" << GetName() << ")::evaluate = " << value << " at x = " << theXVar << endl;
    RooArgSet Xdeps; theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);
    TIterator* iter = Xdeps.createIterator();
    RooAbsArg* var;
    while ((var = (RooAbsArg*)iter->Next())){
      cout << var->GetName() << " value = " << dynamic_cast<RooAbsReal*>(var)->getVal() << endl;
    }
    delete iter;
    cout << endl;
  }
  return value;
}
Int_t RooNCSpline_1D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (_forceNumInt) return 0;
  Int_t code=0;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code=2;
  }
  return code;
}
Double_t RooNCSpline_1D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=RooNCSplineCore::kError) coutE(Integration) << "RooNCSpline_1D_fast ERROR::RooNCSpline_1D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==RooNCSplineCore::kVerbose){ cout << "RooNCSpline_1D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

Bool_t RooNCSpline_1D_fast::testRangeValidity(const T& val, const Int_t /*whichDirection*/) const{
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  return (*(range[0])>*(range[1]) || (val>=*(range[0]) && val<=*(range[1])));
}
void RooNCSpline_1D_fast::setRangeValidity(const T valmin, const T valmax, const Int_t /*whichDirection*/){
  T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  *(range[0])=valmin;
  *(range[1])=valmax;
}
void RooNCSpline_1D_fast::cropValueForRange(T& val, const Int_t /*whichDirection*/)const{
  if (testRangeValidity(val)) return;
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  if (val<*(range[0])) val = *(range[0]);
  if (val>*(range[1])) val = *(range[1]);
}
