#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_2D.h"

using namespace std;


RooNCSplineFactory_2D::RooNCSplineFactory_2D(
  RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_,
  RooNCSplineCore::BoundaryCondition const bcBeginX_,
  RooNCSplineCore::BoundaryCondition const bcEndX_,
  RooNCSplineCore::BoundaryCondition const bcBeginY_,
  RooNCSplineCore::BoundaryCondition const bcEndY_
) :
  appendName(appendName_),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  bcBeginY(bcBeginY_), bcEndY(bcEndY_),
  XVar(&XVar_), YVar(&YVar_),
  fcn(0),
  PDF(0)
{}
RooNCSplineFactory_2D::~RooNCSplineFactory_2D(){
  destroyPDF();
}

void RooNCSplineFactory_2D::addUnique(std::vector<RooNCSplineCore::T>& list, RooNCSplineCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<splineTriplet_t> RooNCSplineFactory_2D::getPoints(
  const std::vector<RooNCSplineCore::T>& XList,
  const std::vector<RooNCSplineCore::T>& YList,
  const std::vector<RooNCSplineCore::T>& FcnList
  ){
  const unsigned int nX = XList.size();
  const unsigned int nY = YList.size();
  const unsigned int n = FcnList.size();
  if (nX*nY!=n){
    cerr << "RooNCSplineFactory_2D::getPoints: nX=" << nX << " x nY=" << nY << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<splineTriplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    RooNCSplineCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      unsigned int ip = nY*ix + iy;
      RooNCSplineCore::T yval = YList.at(iy);
      pList.push_back(splineTriplet_t(xval, yval, FcnList.at(ip)));
    }
  }
  return pList;
}

void RooNCSplineFactory_2D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void RooNCSplineFactory_2D::initPDF(const std::vector<splineTriplet_t>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  vector<RooNCSplineCore::T> XList;
  vector<RooNCSplineCore::T> YList;
  vector<vector<RooNCSplineCore::T>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
  }
  FcnList.reserve(YList.size());
  for (unsigned int iy=0; iy<YList.size(); iy++){
    vector<RooNCSplineCore::T> dum; dum.reserve(XList.size());
    for (unsigned int ix=0; ix<XList.size(); ix++){
      unsigned int ip = YList.size()*ix + iy;
      dum.push_back((pList.at(ip))[2]); // Do not use unique here
    }
    FcnList.push_back(dum);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new RooNCSpline_2D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar,
    XList, YList, FcnList,
    bcBeginX, bcEndX,
    bcBeginY, bcEndY
    );

  name.Prepend("PDF_"); title=name;
  PDF = new RooFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void RooNCSplineFactory_2D::setPoints(TTree* tree){
  vector<splineTriplet_t> pList;
  RooNCSplineCore::T x, y, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(splineTriplet_t(x, y, fcn)); }
  setPoints(pList);
}

void RooNCSplineFactory_2D::setEndConditions(
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
  default:
    // Do nothing
    break;
  }
}
