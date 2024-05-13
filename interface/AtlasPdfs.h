// @(#)root/roostats:$Id: RooBSplineBases.h 873 2014-02-24 22:16:29Z adye $
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
// @(#)root/roostats:$Id: RooBSpline.h 873 2014-02-24 22:16:29Z adye $
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

//#include "RooStats/HistFactory/RooBSplinePenalty.h"
//#include "RooBSplineBases.h"

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
    //RooBSpline(const char* name, const char* title, int order, std::vector<double>& tValues,
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
/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooParamKeysPdf.h 888 2014-08-01 19:54:39Z adye $
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

#ifdef __CINT__
// Specify schema conversion rule here, rather than in LinkDef1.h, so it is included if we compile with ACLiC.
#pragma read sourceClass="RooParamKeysPdf" \
  targetClass="RooParamKeysPdf" \
  version="[-2]" \
  source="Double_t _lookupTable[1001]" \
  target="_nPoints, _lookupTable" \
  code="{ _nPoints=1000; _lookupTable=new Double_t[_nPoints]; for (Int_t i=0; i<_nPoints; i++) _lookupTable[i]= onfile._lookupTable[i]; }"
#endif

#endif
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

  int nnuisSize() { return _nnuis.size(); }

  virtual Bool_t selfNormalized() const { 
    // P.d.f is self normalized
    return kTRUE ; 
  }

  virtual Double_t getVal(const RooArgSet* set=0) const ;
  RooAbsPdf* sumPdf(const RooArgSet* nset) ;

  Bool_t setBinIntegrator(RooArgSet& allVars);
  void useHorizontalMorphing(Bool_t val) { _useHorizMorph=val; }

#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,15)
  void fixCache() { _cacheMgr.setClearOnRedirect(kFALSE) ; }
  virtual CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised ; } ;
#endif
  

protected:

  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem(RooAbsPdf& sumPdf, RooChangeTracker& tracker, const RooArgList& flist) : _sumPdf(&sumPdf), _tracker(&tracker), _fractionsCalculated(false) { 
      _frac.add(flist) ; 
    } ;
    void operModeHook(RooAbsArg::OperMode) {};
    virtual ~CacheElem() ; 
    virtual RooArgList containedArgs(Action) ;
    RooAbsPdf* _sumPdf ;
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

  ClassDef(RooStarMomentMorph,2) // Your description goes here...
};
 
#endif

#ifndef QuickStats_Shapes_RooTwoSidedCBShape
#define QuickStats_Shapes_RooTwoSidedCBShape

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

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

// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
// Modified by Hongtao Yang for xmlAnaWSBuilder: we need to have both upper and lower uncertainties as functions, not just numbers
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ResponseFunction
#define ROOSTATS_ResponseFunction

#include <TIterator.h>
#include <RooListProxy.h>
#include <RooAbsReal.h>
#include <vector>
#include <ostream>

class RooRealVar;
class RooArgList;

class ResponseFunction : public RooAbsReal
{
public:
  ResponseFunction();
  ResponseFunction(const char *name, const char *title,
                   const RooArgList &_paramList,
                   std::vector<double> nominal, const RooArgList &lowList, const RooArgList &highList, std::vector<int> code);

  ResponseFunction(const ResponseFunction &, const char *);

  void setInterpCode(RooAbsReal &param, int code);
  void setAllInterpCodes(int code);
  void setGlobalBoundary(double boundary) { _interpBoundary = boundary; }

  virtual TObject *clone(const char *newname) const { return new ResponseFunction(*this, newname); }
  virtual ~ResponseFunction();

  virtual void printMultiline(std::ostream &os, Int_t contents, Bool_t verbose = kFALSE, TString indent = "") const;
  virtual void printResponseFunctions(std::ostream &os) const;

private:
  double PolyInterpValue(int i, double x) const;
  inline double logn(double nominal, double delta, double theta) const { return pow(1 + delta / nominal, theta); };
  inline void cacheCoef(unsigned j) const;

protected:
  RooListProxy _paramList;
  std::vector<double> _nominal;
  RooListProxy _lowList;
  RooListProxy _highList;
  std::vector<int> _interpCode;
  Double_t _interpBoundary;

  mutable Bool_t _logInit;               //! flag used for chaching polynomial coefficients
  mutable std::vector<double> _polCoeff; //! cached polynomial coefficients

  Double_t evaluate() const;

public:
  double interpolationBoundary() const { return _interpBoundary; }  
  const std::vector<int>& interpolationCodes() const { return _interpCode; }
  const std::vector<double>& nominal() const { return _nominal; }
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,26,0)  
  const RooAbsCollection::Storage_t& parameters() const { return _paramList.get(); }
  const RooAbsCollection::Storage_t& low() const { return _lowList.get(); }
  const RooAbsCollection::Storage_t& high() const { return _highList.get(); }
#endif
  

  ClassDef(ResponseFunction, 2) // flexible interpolation
};

#endif

