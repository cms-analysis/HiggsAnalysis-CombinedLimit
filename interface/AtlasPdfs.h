// Adapted from ATLAS code to run with combine
// Ported by Andrea Carlo Marini 
// Date Fri Oct  1 15:19:36 CEST 2021
// --- Importing file RooBSplineBases
// @(#)root/roostats:$Id:  armbrust $
// Author: Aaron Armbruster
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ROOBSPLINEBASES
#define ROOSTATS_ROOBSPLINEBASES

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include <sstream>


class RooRealVar;
class RooArgList ;

namespace RooStats{
namespace HistFactory{

  class RooBSplineBases : public RooAbsReal {
  public:

    RooBSplineBases() ;
    RooBSplineBases(const char* name, const char* title, int order, std::vector<double>& tValues,
		    RooAbsReal& t, int nrClose=0);

    RooBSplineBases(const char *name, const char *title);
    RooBSplineBases(const RooBSplineBases&, const char*);

    virtual TObject* clone(const char* newname) const { return new RooBSplineBases(*this, newname); }
    virtual ~RooBSplineBases() ;

/*     Double_t getCurvature() const; */

    int getOrder() const {return _n;}
    //void setWeights(const RooArgList& weights);
    Double_t getBasisVal(int n, int i, bool rebuild=true) const;
    std::vector<double> getTValues() const {return _tValues;}
    std::vector<double> getTAry() {return _t_ary;}

  protected:

    void buildTAry() const;

    //RooListProxy _controlPoints;
    std::vector<double> _tValues;
    //RooListProxy _t_ary;
    int _m;
/*     mutable double* _t_ary; //[_m] */
    RooRealProxy _t;
    int _n;
    int _nrClose;
    //int _nPlusOne;
    //mutable double** _bin; //[_nPlusOne][_m]
    mutable std::vector<double> _t_ary;
    mutable std::vector<std::vector<double> > _bin;

    Double_t evaluate() const;

    ClassDef(RooStats::HistFactory::RooBSplineBases,1) // Uniform B-Spline
  };
}
}

#endif
// --- Importing file RooBSpline
// @(#)root/roostats:$Id:  armbrust $
// Author: Aaron Armbruster
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ROOBSPLINE
#define ROOSTATS_ROOBSPLINE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooSetProxy.h"


#include "RooObjCacheManager.h"
#include "RooNumIntConfig.h"


#include <sstream>


class RooRealVar;
class RooArgList ;

namespace RooStats{
namespace HistFactory{

  class RooBSpline : public RooAbsReal {
  public:

    RooBSpline() ;
    RooBSpline(const char* name, const char* title,
	       const RooArgList& controlPoints, RooBSplineBases& bases, const RooArgSet& vars);
    //RooBSpline(const char* name, const char* title, int order, vector<double>& tValues,
    //       const RooArgList& controlPoints, RooAbsReal& t, const RooArgSet& vars, int nrClose=0);

    RooBSpline(const char *name, const char *title);
    RooBSpline(const RooBSpline&, const char*);

    virtual TObject* clone(const char* newname) const { return new RooBSpline(*this, newname); }
    virtual ~RooBSpline() ;

/*     Double_t getCurvature() const; */

//    RooBSplinePenalty* getRealPenalty(int k, RooRealVar* obs, RooRealVar* beta, const char* name = "") const;


    void setWeights(const RooArgList& weights);

    Bool_t setBinIntegrator(RooArgSet& allVars) ;
    Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet,const char* rangeName=0) const ;
    Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0) const ;

    const RooArgList& getControlPoints() const {return _controlPoints;}

    RooBSplineBases* getBases() const {return (RooBSplineBases*)&_bases.arg();}
    int getOrder() const {return _n;}

  protected:

    RooListProxy _controlPoints;
    //RooListProxy _t_ary;
    int _m;
/*     double* _t_ary; //[_m] */
/*     RooRealProxy _t; */
    int _n;
    RooListProxy _weights;
    RooRealProxy _bases;
    RooSetProxy _vars;


    // Cache the integrals   
    class CacheElem : public RooAbsCacheElement {
    public:
      virtual ~CacheElem();
      // Payload
      RooArgList _I ;
      virtual RooArgList containedArgs(Action) ;
    };
    mutable RooObjCacheManager _cacheMgr ; // The cache manager


    Double_t evaluate() const;

    ClassDef(RooStats::HistFactory::RooBSpline,2) // Uniform B-Spline
  };
}
}

#endif
// --- Importing file RooParamKeysPdf
/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooParamKeysPdf.h,v 1.10 2007/05/11 09:13:07 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   UC San Diego,        raven@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_PARAM_KEYS
#define ROO_PARAM_KEYS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooParamKeysPdf : public RooAbsPdf {
public:
  enum Mirror { NoMirror, MirrorLeft, MirrorRight, MirrorBoth,
		MirrorAsymLeft, MirrorAsymLeftRight,
		MirrorAsymRight, MirrorLeftAsymRight,
		MirrorAsymBoth };
  RooParamKeysPdf() ;
  RooParamKeysPdf(
    const char *name, const char *title,
    RooAbsReal& x, RooAbsReal& deltax, 
    RooDataSet& data, Mirror mirror= NoMirror, Double_t rho=1, Int_t nPoints=1000
  );
  RooParamKeysPdf(
    const char *name, const char *title,
    RooAbsReal& x, RooAbsReal& deltax, double centralValue, RooAbsReal& multiplicativeShift,
    RooDataSet& data, Mirror mirror= NoMirror, Double_t rho=1, Int_t nPoints=1000
  );
  RooParamKeysPdf(const RooParamKeysPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const {return new RooParamKeysPdf(*this,newname); }
  virtual ~RooParamKeysPdf();
  
  void LoadDataSet( RooDataSet& data);

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  
  RooRealProxy _x ;
  RooRealProxy _deltax ;
  double _centralValue;
  RooRealProxy _multiplicativeShift;
  Double_t evaluate() const;

private:
  
  Double_t evaluateFull(Double_t x) const;

  Int_t _nEvents;
  Double_t *_dataPts;  //[_nEvents]
  Double_t *_dataWgts; //[_nEvents]
  Double_t *_weights;  //[_nEvents]
  Double_t _sumWgt ;
  mutable Double_t _normVal ;
  
  Int_t _nPoints;

  //enum { _nPoints = 1000 };
  Double_t *_lookupTable; //[_nPoints] 
  
  Double_t g(Double_t x,Double_t sigma) const;

  Bool_t _mirrorLeft, _mirrorRight;
  Bool_t _asymLeft, _asymRight;

  // cached info on variable
  Char_t _varName[128];
  Double_t _lo, _hi, _binWidth;
  Double_t _rho;
  
  ClassDef(RooParamKeysPdf,4) // One-dimensional non-parametric kernel estimation p.d.f.
};

#endif
// --- Importing file RooAbsListContainer
/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOABSLISTCONTAINER
#define ROOABSLISTCONTAINER

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooArgList.h"
#include "RooListProxy.h"
#include "RooDataHist.h"
 
class RooAbsListContainer : public RooAbsReal {
public:
  RooAbsListContainer() {} ; 
  RooAbsListContainer(const char *name, const char *title, RooDataHist& dh);
  //RooAbsListContainer(const char *name, const char *title, const RooArgList& parList);

  RooAbsListContainer(const RooAbsListContainer& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooAbsListContainer(*this,newname); }
  inline virtual ~RooAbsListContainer() { }

  Double_t parVal(const Int_t& idx) const;
  const RooArgList& paramList() const { return _parList ; }
  RooArgList& paramList() { return _parList ; }

protected:

  RooListProxy _parList ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooAbsListContainer,1) // Your description goes here...
};
 
#endif
// --- Importing file RooExpandedDataHist
// vim: ts=4:sw=4
/**********************************************************************************
 * Project: HistFitter - A ROOT-based package for statistical data analysis       *
 * Package: HistFitter                                                            *
 * Class  : RooExpandedDataHist                                                  *
 *                                                                                *
 * Description:                                                                   *
 *      Class derived from RooFitResult, to be able to add more parameters        *
 *      for error propagation (calculation & visualization)                       *
 *                                                                                *
 * Authors:                                                                       *
 *      HistFitter group, CERN, Geneva                                            *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in the file          *
 * LICENSE.                                                                       *
 **********************************************************************************/

#ifndef ROOEXPANDEDDATAHIST_H
#define ROOEXPANDEDDATAHIST_H

#include "TString.h"
#include "RooDataHist.h"
#include "RooAbsArg.h"
#include "RooArgList.h"
#include "RooListProxy.h"

#include <iostream>
#include <vector>

class RooExpandedDataHist: public RooDataHist {

 public:
  RooExpandedDataHist() {} ;
  RooExpandedDataHist(RooDataHist& dh, const char *name=0);
  RooExpandedDataHist(const RooExpandedDataHist& other, const char* name=0) ;

  virtual TObject* Clone(const char* newname=0) const { return new RooExpandedDataHist(*this,newname?newname:GetName()) ; }

  virtual ~RooExpandedDataHist();

  RooArgList& paramList() { return _container.paramList(); }
  const RooArgList& paramList() const { return _container.paramList(); }

  void updateParamList(const RooArgList& newParams) {

/*
    if (_container==0) { // || _firstCall) {
      //if (_container!=NULL) { delete _container; _container=0; }
      _container = new RooAbsListContainer(TString(GetName())+"containerName", TString(GetName())+"containerTitle", *dynamic_cast<RooDataHist*>(this) ) ;

      _firstCall = false;
    } 
*/

    //paramList().assignValueOnly(newParams) ;
    paramList().assignFast(newParams) ;

    _cache_sum_valid = kFALSE ;
  }

  void applyScaleFactor(Bool_t value) { _applyScaleFactor=value; }
  Double_t scaleFactor() { Int_t ibin = calcTreeIndex(); return _container.parVal(ibin); }

 protected:

  virtual Double_t get_wgt(const Int_t& idx)   const; 
  virtual Double_t get_curWeight()   const; 

  RooAbsListContainer _container; 

  Bool_t _firstCall;

  Bool_t   _applyScaleFactor;

 private:
  ClassDef(RooExpandedDataHist,2) // Container class for expanded fit result
};

#endif
// --- Importing file RooExpandedHistPdf
/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOEXPANDEDHISTPDF
#define ROOEXPANDEDHISTPDF

#include "RooHistPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooArgList.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooExpandedDataHist;
 
class RooExpandedHistPdf : public RooHistPdf {

 friend class RooMCHistConstraint;

 public:
  RooExpandedHistPdf() {} ; 
  RooExpandedHistPdf(const char* name, const char* title, const RooArgSet& vars, const RooExpandedDataHist& dhist, Int_t intOrder = 0);
  RooExpandedHistPdf(const char* name, const char* title, const RooArgSet& vars, const RooDataHist& dhist, RooArgList& parList, Int_t intOrder = 0);
  RooExpandedHistPdf(const char *name, RooHistPdf& histPdf, RooArgList& parList);
  RooExpandedHistPdf(const RooExpandedHistPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpandedHistPdf(*this,newname); }
  inline virtual ~RooExpandedHistPdf() { }

  double scaleFactor() { return ((RooExpandedDataHist*)_dataHist)->scaleFactor(); }

  void applyScaleFactor(Bool_t val) { ((RooExpandedDataHist*)_dataHist)->applyScaleFactor(val); }
  
  virtual TObject* Clone(const char* newname=0) const { return new RooExpandedHistPdf(*this,newname?newname:GetName()) ; }
 protected:
  RooListProxy _parList ;

  Double_t evaluate() const ; 
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  
 private:
  ClassDef(RooExpandedHistPdf,1) // Your description goes here...
};
 
#endif
// --- Importing file RooStarMomentMorph
/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOSTARMOMENTMORPH
#define ROOSTARMOMENTMORPH

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooSetProxy.h"
#include "RooListProxy.h"
#include "RooArgList.h"

#include "TMatrixD.h"

#include <vector>
#include <string>
class RooChangeTracker ;

class RooExpandedHistPdf;

class RooStarMomentMorph : public RooAbsPdf {
public:

  enum Setting { Linear, NonLinear, NonLinearPosFractions, NonLinearLinFractions, SineLinear } ;

  RooStarMomentMorph() ;
  
  RooStarMomentMorph(const char *name, const char *title, const RooArgList& parList, 
		     const RooArgList& obsList, const RooArgList& pdfList, 
		     const std::vector<int>& nnuispoints, // # of pdfs for each nuisance parameter
		     const std::vector<double>& nrefpoints, 
		     const Setting& setting = NonLinearPosFractions );

  RooStarMomentMorph(const RooStarMomentMorph& other, const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooStarMomentMorph(*this,newname); }

  virtual ~RooStarMomentMorph();

  void     setMode(const Setting& setting) { _setting = setting; }
  void     setUseBinByBinScaleFactors(Bool_t value);

  int nnuisSize() { return _nnuis.size(); }

  virtual Bool_t selfNormalized() const { 
    // P.d.f is self normalized
    return !_useBinByBinScaleFactors;
  }

  virtual Double_t getVal(const RooArgSet* set=0) const ;
  RooAbsPdf* sumPdfPos(const RooArgSet* nset) ;
  RooAbsPdf* sumPdfNeg(const RooArgSet* nset) ;

  Bool_t setBinIntegrator(RooArgSet& allVars);
  void useHorizontalMorphing(Bool_t val) { _useHorizMorph=val; }

  void fixCache() { _cacheMgr.setClearOnRedirect(kFALSE) ; }
  virtual CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised ; } ;
  

  Double_t fraction(Int_t i) const ; // returns the fraction of the i-th pdf coefficient

protected:

  class CacheElem : public RooAbsCacheElement {
  public:
  CacheElem(RooAbsPdf& sumPdfPos, RooAbsPdf& sumPdfNeg, RooChangeTracker& tracker, const RooArgList& flist) : _sumPdfPos(&sumPdfPos), _sumPdfNeg(&sumPdfNeg), _tracker(&tracker), _fractionsCalculated(false) { 
      _frac.add(flist) ; 
    } ;
    void operModeHook(RooAbsArg::OperMode) {};
    virtual ~CacheElem() ; 
    virtual RooArgList containedArgs(Action) ;
 

    RooAbsPdf* _sumPdfPos ;
    RooAbsPdf* _sumPdfNeg ;
    RooChangeTracker* _tracker ; 
    RooArgList _frac ;

    RooRealVar* frac(Int_t i ) ;
    const RooRealVar* frac(Int_t i ) const ; 
    void calculateFractions(const RooStarMomentMorph& self, Bool_t verbose=kTRUE) const;
    mutable bool _fractionsCalculated;
  } ;
  mutable RooObjCacheManager _cacheMgr ; // The cache manager
  mutable RooArgSet* _curNormSet ; //! Current normalization set

  friend class CacheElem ; // Cache needs to be able to clear _norm pointer

  Double_t evaluate() const ;

  void     initialize();
  CacheElem* getCache(const RooArgSet* nset) const ;

  Int_t ij(const Int_t& i, const Int_t& j) const;
  Int_t ijhi(const Int_t& i, const double& nval) const;
  Int_t ijlo(const Int_t& i, const double& nval) const;

  inline Int_t sij(const Int_t& i, const Int_t& j) const { return (i*_obsList.getSize()+j); }

  
  RooListProxy _parList ;
  RooSetProxy  _obsList ;
  RooListProxy _pdfList ;

  mutable std::vector<int> _nnuis; 
  mutable std::vector<double> _nref;

  TIterator* _parItr ;  //! do not persist
  TIterator* _obsItr ;  //! do not persist
  TIterator* _pdfItr ;  //!
  mutable TMatrixD* _M; //!

  Setting _setting;
  unsigned int _nnuisvar;


  Bool_t _useHorizMorph;

  Bool_t _useBinByBinScaleFactors; // this is only if the RSMM is made with RooExpandedHistPdf
  mutable RooAbsPdf* _nominalPdf; //!

  ClassDef(RooStarMomentMorph,4) // Your description goes here...
};
 
#endif


// --- Importing file RooTwoSidedCBShape
#ifndef ROOT_Hfitter_RooTwoSidedCBShape
#define ROOT_Hfitter_RooTwoSidedCBShape

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"

class RooRealVar;


class RooTwoSidedCBShape : public RooAbsPdf {
public:
  RooTwoSidedCBShape() {}
  RooTwoSidedCBShape(const char *name, const char *title, RooAbsReal& _m,
               RooAbsReal& _m0, RooAbsReal& _sigma,
               RooAbsReal& _alphaLo, RooAbsReal& _nLo,
               RooAbsReal& _alphaHi, RooAbsReal& _nHi);

  RooTwoSidedCBShape(const RooTwoSidedCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooTwoSidedCBShape(*this,newname); }

  inline virtual ~RooTwoSidedCBShape() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;
  
protected:

  Double_t gaussianIntegral(Double_t tmin, Double_t tmax) const;
  Double_t powerLawIntegral(Double_t tmin, Double_t tmax, Double_t alpha, Double_t n) const;
  Double_t evaluate() const;

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaLo;
  RooRealProxy nLo;
  RooRealProxy alphaHi;
  RooRealProxy nHi;


private:

  ClassDef(RooTwoSidedCBShape,1)
};


#endif
