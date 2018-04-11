#ifndef ROONCSPLINEFACTORY_2D
#define ROONCSPLINEFACTORY_2D

#include <vector>
#include <utility>
#include <algorithm>
#include "TTree.h"
#include "RooNCSpline_2D_fast.h"
#include "RooFuncPdf.h"


namespace NumUtils{
  template<typename T> struct triplet{
    T value[3];
    triplet(T i1, T i2, T i3){
      value[0]=i1;
      value[1]=i2;
      value[2]=i3;

    }
    triplet(T i1){
      value[0]=i1;
      value[1]=i1;
      value[2]=i1;

    }
    triplet(){}
    T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
    const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
  };

  typedef triplet<int> intTriplet_t;
  typedef triplet<float> floatTriplet_t;
  typedef triplet<double> doubleTriplet_t;
}

typedef NumUtils::triplet<RooNCSplineCore::T> splineTriplet_t;


class RooNCSplineFactory_2D{
protected:
  TString appendName;

  RooNCSplineCore::BoundaryCondition bcBeginX;
  RooNCSplineCore::BoundaryCondition bcEndX;
  RooNCSplineCore::BoundaryCondition bcBeginY;
  RooNCSplineCore::BoundaryCondition bcEndY;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooNCSpline_2D_fast* fcn;
  RooFuncPdf* PDF;

  const std::vector<splineTriplet_t> getPoints(const std::vector<RooNCSplineCore::T>& XList, const std::vector<RooNCSplineCore::T>& YList, const std::vector<RooNCSplineCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<splineTriplet_t>& pList);

  void addUnique(std::vector<RooNCSplineCore::T>& list, RooNCSplineCore::T val);

public:
  RooNCSplineFactory_2D(
    RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_="",
    RooNCSplineCore::BoundaryCondition const bcBeginX_=RooNCSplineCore::bcNaturalSpline,
    RooNCSplineCore::BoundaryCondition const bcEndX_=RooNCSplineCore::bcNaturalSpline,
    RooNCSplineCore::BoundaryCondition const bcBeginY_=RooNCSplineCore::bcNaturalSpline,
    RooNCSplineCore::BoundaryCondition const bcEndY_=RooNCSplineCore::bcNaturalSpline
  );
  ~RooNCSplineFactory_2D();

  RooNCSpline_2D_fast* getFunc(){ return fcn; }
  RooFuncPdf* getPDF(){ return PDF; }

  void setEndConditions(
    RooNCSplineCore::BoundaryCondition const bcBegin,
    RooNCSplineCore::BoundaryCondition const bcEnd,
    const unsigned int direction
  );

  void setPoints(TTree* tree);
  void setPoints(const std::vector<splineTriplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& FcnList){
    std::vector<RooNCSplineCore::T> transXList;
    std::vector<RooNCSplineCore::T> transYList;
    std::vector<RooNCSplineCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((RooNCSplineCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((RooNCSplineCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((RooNCSplineCore::T)FcnList.at(ip));
    const std::vector<splineTriplet_t> pList = getPoints(transXList, transYList, transFcnList);
    setPoints(pList);
  }

};

template void RooNCSplineFactory_2D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& FcnList);
template void RooNCSplineFactory_2D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& FcnList);

#endif



