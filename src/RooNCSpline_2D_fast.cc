#include "HiggsAnalysis/CombinedLimit/interface/RooNCSpline_2D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;


ClassImp(RooNCSpline_2D_fast)

RooNCSpline_2D_fast::RooNCSpline_2D_fast() :
  RooNCSplineCore(),
  rangeYmin(1), rangeYmax(-1),
  bcBeginX(RooNCSplineCore::bcNaturalSpline), bcEndX(RooNCSplineCore::bcNaturalSpline),
  bcBeginY(RooNCSplineCore::bcNaturalSpline), bcEndY(RooNCSplineCore::bcNaturalSpline),
  theYVar("theYVar", "theYVar", this)
{}

RooNCSpline_2D_fast::RooNCSpline_2D_fast(
  const char* name,
  const char* title
  ) :
  RooNCSplineCore(name, title),
  rangeYmin(1), rangeYmax(-1),
  bcBeginX(RooNCSplineCore::bcNaturalSpline), bcEndX(RooNCSplineCore::bcNaturalSpline),
  bcBeginY(RooNCSplineCore::bcNaturalSpline), bcEndY(RooNCSplineCore::bcNaturalSpline),
  theYVar("theYVar", "theYVar", this)
{}

RooNCSpline_2D_fast::RooNCSpline_2D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  const std::vector<T>& inXList,
  const std::vector<T>& inYList,
  const std::vector<std::vector<T>>& inFcnList,
  RooNCSplineCore::BoundaryCondition const bcBeginX_,
  RooNCSplineCore::BoundaryCondition const bcEndX_,
  RooNCSplineCore::BoundaryCondition const bcBeginY_,
  RooNCSplineCore::BoundaryCondition const bcEndY_,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooNCSplineCore(name, title, inXVar, inXList, inUseFloor, inFloorEval, inFloorInt),
  rangeYmin(1), rangeYmax(-1),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  bcBeginY(bcBeginY_), bcEndY(bcEndY_),
  theYVar("theYVar", "theYVar", this, inYVar),
  YList(inYList),
  FcnList(inFcnList)
{
  if (npointsX()>1 && npointsY()>1){
    // Prepare A and kappa arrays for x and y coordinates
    int npoints;
    Double_t det;

    vector<vector<RooNCSplineCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA, bcBeginX, bcEndX);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSpline_2D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<RooNCSplineCore::T>> yA; getKappas(kappaY, 1); getAArray(kappaY, yA, bcBeginY, bcEndY);
    npoints=kappaY.size();
    TMatrix_t yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrix_t yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "RooNCSpline_2D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
    }

    // Get the grid of coefficients
    vector<vector<vector<RooNCSplineCore::T>>> coefsAlongY; // [Ax(y),Bx(y),Cx(y),Dx(y)][xbin][ybin]
    int npoldim=0;
    int nxbins=0;
    for (unsigned int j=0; j<npointsY(); j++){
      vector<vector<RooNCSplineCore::T>> xcoefsAtYj = getCoefficientsPerY(kappaX, xAinv, j, bcBeginX, bcEndX, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j
      if (j==0){
        nxbins=xcoefsAtYj.size();
        npoldim=xcoefsAtYj.at(0).size();
        for (int ipow=0; ipow<npoldim; ipow++){
          vector<vector<RooNCSplineCore::T>> dum_xycoefarray;
          for (int ix=0; ix<nxbins; ix++){
            vector<RooNCSplineCore::T> dum_ycoefarray;
            dum_xycoefarray.push_back(dum_ycoefarray);
          }
          coefsAlongY.push_back(dum_xycoefarray);
        }
      }
      if (nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()){
        coutE(InputArguments) << "RooNCSpline_2D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()!" << endl;
        assert(0);
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int ipow=0; ipow<npoldim; ipow++) coefsAlongY.at(ipow).at(ix).push_back(xcoefsAtYj.at(ix).at(ipow));
      }
    }

    for (int ix=0; ix<nxbins; ix++){
      // Get the x coefficients interpolated across y
      vector<vector<vector<RooNCSplineCore::T>>> xCoefs;
      for (int ic=0; ic<npoldim; ic++){
        vector<vector<RooNCSplineCore::T>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ic).at(ix), bcBeginY, bcEndY, -1); // [iy][A,B,C,D]
        xCoefs.push_back(yCoefs);
      }
      coefficients.push_back(xCoefs);
    }
  }
  else assert(0);

  RooArgSet leafset;
  getLeafDependents(theXVar, leafset);
  getLeafDependents(theYVar, leafset);
  addLeafDependents(leafset);

  emptyFcnList();
}

RooNCSpline_2D_fast::RooNCSpline_2D_fast(
  const RooNCSpline_2D_fast& other,
  const char* name
  ) :
  RooNCSplineCore(other, name),
  rangeYmin(other.rangeYmin), rangeYmax(other.rangeYmax),
  bcBeginX(other.bcBeginX), bcEndX(other.bcEndX),
  bcBeginY(other.bcBeginY), bcEndY(other.bcEndY),
  theYVar("theYVar", this, other.theYVar),
  YList(other.YList),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  kappaY(other.kappaY),
  coefficients(other.coefficients)
{}


RooNCSplineCore::T RooNCSpline_2D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<RooNCSplineCore::T> res;

  if (verbosity==RooNCSplineCore::kVerbose){ cout << "RooNCSpline_2D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl; }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1, ybin=-1, ybinmin=-1, ybinmax=-1;
  RooNCSplineCore::T tx=0, txmin=0, txmax=0, ty=0, tymin=0, tymax=0;
  if (code==0 || code%2!=0){ // Case to just compute the value at x
    if (!testRangeValidity(theXVar, 0)) return 0;
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Case to integrate along x
    RooNCSplineCore::T coordmin = theXVar.min(rangeName); cropValueForRange(coordmin, 0);
    RooNCSplineCore::T coordmax = theXVar.max(rangeName); cropValueForRange(coordmax, 0);
    xbinmin = getWhichBin(coordmin, 0);
    txmin = getTVar(kappaX, coordmin, xbinmin, 0);
    xbinmax = getWhichBin(coordmax, 0);
    txmax = getTVar(kappaX, coordmax, xbinmax, 0);
  }
  if (code==0 || code%3!=0){ // Case to just compute the value at y
    if (!testRangeValidity(theYVar, 1)) return 0;
    ybin = getWhichBin(theYVar, 1);
    ty = getTVar(kappaY, theYVar, ybin, 1);
  }
  else{ // Case to integrate along y
    RooNCSplineCore::T coordmin = theYVar.min(rangeName); cropValueForRange(coordmin, 1);
    RooNCSplineCore::T coordmax = theYVar.max(rangeName); cropValueForRange(coordmax, 1);
    ybinmin = getWhichBin(coordmin, 1);
    tymin = getTVar(kappaY, coordmin, ybinmin, 1);
    ybinmax = getWhichBin(coordmax, 1);
    tymax = getTVar(kappaY, coordmax, ybinmax, 1);
  }

  for (int ix=0; ix<(int)coefficients.size(); ix++){
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

    if (verbosity==RooNCSplineCore::kVerbose){
      if (code==0 || code%2!=0) cout << "Evaluating tx=" << txhigh << " in bin " << ix << endl;
      else cout << "Evaluating tx[" << txlow << ", " << txhigh << "] in bin " << ix << endl;
    }

    // Get the x coefficients interpolated across y
    vector<RooNCSplineCore::T> xCoefs;
    for (int ic=0; ic<(int)coefficients.at(ix).size(); ic++){
      const vector<vector<RooNCSplineCore::T>>& yCoefs = coefficients.at(ix).at(ic);

      if (verbosity==RooNCSplineCore::kVerbose) cout << "\tCoefficient " << ic << ":\n";

      DefaultAccumulator<RooNCSplineCore::T> theCoef;
      for (int iy=0; iy<(int)yCoefs.size(); iy++){
        if (
          (ybin>=0 && iy!=ybin)
          ||
          (ybinmin>=0 && ybinmax>=ybinmin && !(ybinmin<=iy && iy<=ybinmax))
          ) continue;

        RooNCSplineCore::T tylow=0, tyhigh=1;
        if (code>0 && code%3==0){
          if (iy==ybinmin) tylow=tymin;
          if (iy==ybinmax) tyhigh=tymax;
        }
        else tyhigh=ty;

        if (verbosity==RooNCSplineCore::kVerbose){
          if (code==0 || code%3!=0) cout << "\tEvaluating ty=" << tyhigh << " in bin " << iy << endl;
          else cout << "\tEvaluating ty[" << tylow << ", " << tyhigh << "] in bin " << iy << endl;
        }

        theCoef += evalSplineSegment(yCoefs.at(iy), kappaY.at(iy), tyhigh, tylow, (code>0 && code%3==0));
      }

      //if (code==0) cout << "\tCoefficient is " << theCoef << endl;

      xCoefs.push_back(theCoef);
    }

    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void RooNCSpline_2D_fast::getKappas(vector<RooNCSplineCore::T>& kappas, const Int_t whichDirection){
  kappas.clear();
  RooNCSplineCore::T kappa=1;

  Int_t npoints;
  vector<RooNCSplineCore::T> const* coord;
  if (whichDirection==0){
    npoints=npointsX();
    coord=&XList;
  }
  else{
    npoints=npointsY();
    coord=&YList;
  }

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
Int_t RooNCSpline_2D_fast::getWhichBin(const RooNCSplineCore::T& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  RooNCSplineCore::T valj, valjpo;
  Int_t npoints;
  vector<RooNCSplineCore::T> const* coord;
  if (whichDirection==0){
    coord=&XList;
    npoints=npointsX();
  }
  else{
    coord=&YList;
    npoints=npointsY();
  }

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
RooNCSplineCore::T RooNCSpline_2D_fast::getTVar(const vector<RooNCSplineCore::T>& kappas, const RooNCSplineCore::T& val, const Int_t& bin, const Int_t whichDirection)const{
  const RooNCSplineCore::T& K=kappas.at(bin);
  vector<RooNCSplineCore::T> const* coord;
  if (whichDirection==0) coord=&XList;
  else coord=&YList;
  return (val-coord->at(bin))*K;
}

vector<vector<RooNCSplineCore::T>> RooNCSpline_2D_fast::getCoefficientsPerY(const std::vector<RooNCSplineCore::T>& kappaX, const TMatrix_t& xAinv, const Int_t& ybin, RooNCSplineCore::BoundaryCondition const& bcBegin, RooNCSplineCore::BoundaryCondition const& bcEnd, const Int_t xbin)const{
  vector<RooNCSplineCore::T> fcnList;
  for (unsigned int bin=0; bin<npointsX(); bin++){ fcnList.push_back(FcnList.at(ybin).at(bin)); }
  vector<vector<RooNCSplineCore::T>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, bcBegin, bcEnd, xbin);
  return coefs;
}

Double_t RooNCSpline_2D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=RooNCSplineCore::kError) coutE(Eval) << "RooNCSpline_2D_fast ERROR::RooNCSpline_2D_fast(" << GetName() << ") evaluation returned " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl;
    value = floorEval;
  }
  if (verbosity==RooNCSplineCore::kVerbose){ cout << "RooNCSpline_2D_fast(" << GetName() << ")::evaluate = " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl; }
  return value;
}
Int_t RooNCSpline_2D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (_forceNumInt) return 0;

  Int_t code=1;

  RooArgSet Xdeps, Ydeps;
  RooRealVar* rrv_x = dynamic_cast<RooRealVar*>(theXVar.absArg());
  RooRealVar* rrv_y = dynamic_cast<RooRealVar*>(theYVar.absArg());
  if (rrv_x==0) theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);
  if (rrv_y==0) theYVar.absArg()->leafNodeServerList(&Ydeps, 0, true);

  if (rrv_x!=0){
    if (Ydeps.find(*rrv_x)==0 || rrv_y!=0){
      if (matchArgs(allVars, analVars, theXVar)) code*=2;
    }
  }
  if (rrv_y!=0){
    if (Xdeps.find(*rrv_y)==0 || rrv_x!=0){
      if (matchArgs(allVars, analVars, theYVar)) code*=3;
    }
  }

  if (code==1) code=0;
  return code;
}
Double_t RooNCSpline_2D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=RooNCSplineCore::kError) coutE(Integration) << "RooNCSpline_2D_fast ERROR::RooNCSpline_2D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==RooNCSplineCore::kVerbose){ cout << "RooNCSpline_2D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

Bool_t RooNCSpline_2D_fast::testRangeValidity(const T& val, const Int_t whichDirection) const{
  const T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  return (*(range[0])>*(range[1]) || (val>=*(range[0]) && val<=*(range[1])));
}
void RooNCSpline_2D_fast::setRangeValidity(const T valmin, const T valmax, const Int_t whichDirection){
  T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  *(range[0])=valmin;
  *(range[1])=valmax;
}
void RooNCSpline_2D_fast::cropValueForRange(T& val, const Int_t whichDirection)const{
  if (testRangeValidity(val, whichDirection)) return;
  const T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  if (val<*(range[0])) val = *(range[0]);
  if (val>*(range[1])) val = *(range[1]);
}
