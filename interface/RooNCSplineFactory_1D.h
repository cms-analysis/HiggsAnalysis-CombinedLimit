#ifndef ROONCSPLINEFACTORY_1D
#define ROONCSPLINEFACTORY_1D

#include <vector>
#include <utility>
#include <algorithm>
#include "TGraph.h"
#include "TTree.h"
#include "RooNCSpline_1D_fast.h"
#include "RooFuncPdf.h"

class RooNCSplineFactory_1D{
protected:
  TString appendName;

  RooNCSplineCore::BoundaryCondition bcBeginX;
  RooNCSplineCore::BoundaryCondition bcEndX;

  RooAbsReal* splineVar;
  RooNCSpline_1D_fast* fcn;
  RooFuncPdf* PDF;

  const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>> getPoints(const std::vector<RooNCSplineCore::T>& XList, const std::vector<RooNCSplineCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>>& pList);

public:
  RooNCSplineFactory_1D(
    RooAbsReal& splineVar_, TString appendName_="",
    RooNCSplineCore::BoundaryCondition const bcBeginX_=RooNCSplineCore::bcNaturalSpline,
    RooNCSplineCore::BoundaryCondition const bcEndX_=RooNCSplineCore::bcNaturalSpline
  );
  ~RooNCSplineFactory_1D();

  RooNCSpline_1D_fast* getFunc(){ return fcn; }
  RooFuncPdf* getPDF(){ return PDF; }

  void setEndConditions(
    RooNCSplineCore::BoundaryCondition const bcBegin,
    RooNCSplineCore::BoundaryCondition const bcEnd,
    const unsigned int direction=0
  );

  void setPoints(TTree* tree);
  void setPoints(TGraph* tg);
  void setPoints(const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& FcnList){
    std::vector<RooNCSplineCore::T> transXList;
    std::vector<RooNCSplineCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((RooNCSplineCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((RooNCSplineCore::T)FcnList.at(ip));
    const std::vector<std::pair<RooNCSplineCore::T, RooNCSplineCore::T>> pList = getPoints(transXList, transFcnList);
    initPDF(pList);
  }

};

template void RooNCSplineFactory_1D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& FcnList);
template void RooNCSplineFactory_1D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

#endif



