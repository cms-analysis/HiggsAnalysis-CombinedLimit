#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_3D.h"

using namespace std;


RooNCSplineFactory_3D::RooNCSplineFactory_3D(
  RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_,
  RooNCSplineCore::BoundaryCondition const bcBeginX_,
  RooNCSplineCore::BoundaryCondition const bcEndX_,
  RooNCSplineCore::BoundaryCondition const bcBeginY_,
  RooNCSplineCore::BoundaryCondition const bcEndY_,
  RooNCSplineCore::BoundaryCondition const bcBeginZ_,
  RooNCSplineCore::BoundaryCondition const bcEndZ_
) :
  appendName(appendName_),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  bcBeginY(bcBeginY_), bcEndY(bcEndY_),
  bcBeginZ(bcBeginZ_), bcEndZ(bcEndZ_),
  XVar(&XVar_), YVar(&YVar_), ZVar(&ZVar_),
  fcn(0),
  PDF(0)
{}
RooNCSplineFactory_3D::~RooNCSplineFactory_3D(){
  destroyPDF();
}

void RooNCSplineFactory_3D::addUnique(std::vector<RooNCSplineCore::T>& list, RooNCSplineCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<splineQuadruplet_t> RooNCSplineFactory_3D::getPoints(
  const std::vector<RooNCSplineCore::T>& XList,
  const std::vector<RooNCSplineCore::T>& YList,
  const std::vector<RooNCSplineCore::T>& ZList,
  const std::vector<RooNCSplineCore::T>& FcnList
  ){
  const unsigned int nX = XList.size();
  const unsigned int nY = YList.size();
  const unsigned int nZ = ZList.size();
  const unsigned int n = FcnList.size();
  if (nX*nY*nZ!=n){
    cerr << "RooNCSplineFactory_3D::getPoints: nX=" << nX << " x nY=" << nY << " x nZ=" << nZ << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<splineQuadruplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    RooNCSplineCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      RooNCSplineCore::T yval = YList.at(iy);
      for (unsigned int iz=0; iz<nZ; iz++){
        RooNCSplineCore::T zval = ZList.at(iz);
        unsigned int ip = nZ*(nY*ix + iy) + iz;
        pList.push_back(splineQuadruplet_t(xval, yval, zval, FcnList.at(ip)));
      }
    }
  }
  return pList;
}

void RooNCSplineFactory_3D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void RooNCSplineFactory_3D::initPDF(const std::vector<splineQuadruplet_t>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  vector<RooNCSplineCore::T> XList;
  vector<RooNCSplineCore::T> YList;
  vector<RooNCSplineCore::T> ZList;
  vector<vector<vector<RooNCSplineCore::T>>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
    addUnique(ZList, (pList.at(ip))[2]);
  }
  FcnList.reserve(ZList.size());
  for (unsigned int iz=0; iz<ZList.size(); iz++){
    vector<vector<RooNCSplineCore::T>> dumz;
    dumz.reserve(YList.size());
    for (unsigned int iy=0; iy<YList.size(); iy++){
      vector<RooNCSplineCore::T> dumy;
      dumy.reserve(XList.size());
      for (unsigned int ix=0; ix<XList.size(); ix++){
        unsigned int ip = ZList.size()*(YList.size()*ix + iy) + iz;
        dumy.push_back((pList.at(ip))[3]); // Do not use unique here
      }
      dumz.push_back(dumy);
    }
    FcnList.push_back(dumz);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new RooNCSpline_3D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar, *ZVar,
    XList, YList, ZList, FcnList,
    bcBeginX, bcEndX,
    bcBeginY, bcEndY,
    bcBeginZ, bcEndZ
  );

  name.Prepend("PDF_"); title=name;
  PDF = new RooFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void RooNCSplineFactory_3D::setPoints(TTree* tree){
  vector<splineQuadruplet_t> pList;
  RooNCSplineCore::T x, y, z, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Z", &z);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(splineQuadruplet_t(x, y, z, fcn)); }
  setPoints(pList);
}

void RooNCSplineFactory_3D::setEndConditions(
  RooNCSplineCore::BoundaryCondition const bcBegin,
  RooNCSplineCore::BoundaryCondition const bcEnd,
  const unsigned int direction
){
  switch (direction){
  case 0:
    bcBeginX=bcBegin;
    bcEndX=bcEnd;
    break;
  case 1:
    bcBeginY=bcBegin;
    bcEndY=bcEnd;
    break;
  case 2:
    bcBeginZ=bcBegin;
    bcEndZ=bcEnd;
    break;
  default:
    // Do nothing
    break;
  }
}
