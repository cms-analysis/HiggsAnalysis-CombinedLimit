#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_1D.h"
#include <cassert>

using namespace std;


RooNCSplineFactory_1D::RooNCSplineFactory_1D(
  RooAbsReal& splineVar_, TString appendName_,
  RooNCSplineCore::BoundaryCondition const bcBeginX_,
  RooNCSplineCore::BoundaryCondition const bcEndX_
) :
  appendName(appendName_),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  splineVar(&splineVar_),
  fcn(0),
  PDF(0)
{}
RooNCSplineFactory_1D::~RooNCSplineFactory_1D(){
  destroyPDF();
}
void RooNCSplineFactory_1D::setPoints(TTree* tree){
  vector<pair<RooNCSplineCore::T, RooNCSplineCore::T>> pList;
  RooNCSplineCore::T x, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(pair<RooNCSplineCore::T, RooNCSplineCore::T>(x, fcn)); }
  setPoints(pList);
}
void RooNCSplineFactory_1D::setPoints(TGraph* tg){
  vector<pair<RooNCSplineCore::T, RooNCSplineCore::T>> pList;
  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++) pList.push_back(pair<RooNCSplineCore::T, RooNCSplineCore::T>(xx[ip], yy[ip]));
  setPoints(pList);
}
const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>> RooNCSplineFactory_1D::getPoints(const std::vector<RooNCSplineCore::T>& XList, const std::vector<RooNCSplineCore::T>& FcnList){
  const unsigned int nX = XList.size();
  const unsigned int n = FcnList.size();
  if (nX!=n){
    cerr << "RooNCSplineFactory_1D::getPoints: nX=" << nX << " != nFcn=" << n << endl;
    assert(0);
  }
  std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>> pList; pList.reserve(n);
  for (unsigned int ip=0; ip<n; ip++) pList.push_back(pair<RooNCSplineCore::T, RooNCSplineCore::T>(XList.at(ip), FcnList.at(ip)));
  return pList;
}

void RooNCSplineFactory_1D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void RooNCSplineFactory_1D::initPDF(const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  std::vector<RooNCSplineCore::T> XList;
  std::vector<RooNCSplineCore::T> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    XList.push_back(pList.at(ip).first);
    FcnList.push_back(pList.at(ip).second);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new RooNCSpline_1D_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList, FcnList,
    bcBeginX, bcEndX
    );

  name.Prepend("PDF_"); title=name;
  PDF = new RooFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void RooNCSplineFactory_1D::setEndConditions(
  RooNCSplineCore::BoundaryCondition const bcBegin,
  RooNCSplineCore::BoundaryCondition const bcEnd,
  const unsigned int /*direction*/
){
  bcBeginX=bcBegin;
  bcEndX=bcEnd;
}
