// @(#)root/roostats:$Id: RooBSplineBases.cxx 873 2014-02-24 22:16:29Z adye $
// Author: Aaron Armbruster
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//_________________________________________________
/*
BEGIN_HTML
<p>
</p>
END_HTML
*/
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "TMath.h"

#include "HiggsAnalysis/CombinedLimit/interface/AtlasPdfs.h"

ClassImp(RooStats::HistFactory::RooBSplineBases)

using namespace RooStats;
using namespace HistFactory;

//_____________________________________________________________________________
RooBSplineBases::RooBSplineBases()
{
  // Default constructor
//   _t_ary=NULL;
//   _bin=NULL;
}


//_____________________________________________________________________________
RooBSplineBases::RooBSplineBases(const char* name, const char* title, int order, std::vector<double>& tValues,
				 RooAbsReal& t, int nrClose) :
  RooAbsReal(name, title),
  _tValues(tValues),
  _m(tValues.size()+2*order),
  //_t_ary(NULL),
  _t("t_param", "t_param", this, t),
  _n(order),
  _nrClose(nrClose)
  //_nPlusOne(order+1)//,
  //_bin(NULL)
{
  //std::cout << "in Ctor" << std::endl;
  //std::cout << "t = " << t << std::endl;
  buildTAry();

//   _bin=new double*[_n+1];
//   for (int n=0;n<_n+1;n++)
//   {
//     _bin[n] = new double[_m];
//     for (int i=0;i<_m;i++)
//     {
//       _bin[n][i] = 0;
//     }
//   }

//   _t_ary = new double[_m];
//   for (int i=_n;i<_m-_n;i++) // add the main knots
//   {
//     _t_ary[i] = tValues[i-_n];
//     //std::cout << "Adding main point:   " << i << " = " << _t_ary[i] << std::endl;
//   }

//   double firstDelta=_t_ary[_n+1]-_t_ary[_n]; // extrapolate to the lower non-closed knots
//   for (int i=nrClose;i<_n;i++) 
//   {
//     _t_ary[i] = _t_ary[_n]+firstDelta*(i-_n);
//     //std::cout << "Adding lower open:   " << i << " = " << _t_ary[i] << std::endl;
//   }

//   for (int i=0;i<nrClose;i++) // add the lower closed knots
//   {
//     _t_ary[i] = _t_ary[nrClose];
//     //std::cout << "Adding lower closed: " << i << " = " << _t_ary[i] << std::endl;
//   }


//   double lastDelta=_t_ary[_m-_n-1]-_t_ary[_m-_n-2]; //extrapolate the upper non-closed knots
//   for (int i=_m-_n;i<_m-nrClose;i++) 
//   {
//     _t_ary[i] = _t_ary[_m-_n-1]+lastDelta*(i-(_m-_n-1));
//     //std::cout << "Adding upper open:   " << i << " = " << _t_ary[i] << std::endl;
//   }

//   for (int i=_m-nrClose;i<_m;i++) // add the upper closed knots
//   {
//     _t_ary[i] = _t_ary[_m-nrClose-1];
//     //std::cout << "Adding upper closed: " << i << " = " << _t_ary[i] << std::endl;
//   }
//   //std::cout << std::endl;

//   for (int i=0;i<_m;i++)
//   {
//     if (fabs(_t_ary[i]) < pow(10., -9)) _t_ary[i] = 0;
//   }

  //std::cout << "out Ctor" << std::endl;



  _bin.resize(_n+1);
  for (int i=0;i<_n+1;i++)
  {
    _bin[i].resize(_m);
  }

}

//_____________________________________________________________________________
void RooBSplineBases::buildTAry() const
{
  //delete[] _t_ary;
  _t_ary.resize(_m);
  //std::cout << "In buildTAry. _m=" << _m << ", _n=" << _n << ", _t_ary.size()=" << _t_ary.size() << std::endl;
  //_t_ary = new double[_m];
  for (int i=_n;i<_m-_n;i++) // add the main knots
  {
    _t_ary[i] = _tValues[i-_n];
    //std::cout << "Adding main point:   " << i << " = " << _t_ary[i] << std::endl;
  }

  double firstDelta=_t_ary[_n+1]-_t_ary[_n]; // extrapolate to the lower non-closed knots
//   std::cout << "Starting loop" << std::endl;
//   std::cout << "_nrClose=" << _nrClose << std::endl;
  for (int i=_nrClose;i<_n;i++) 
  {
    _t_ary[i] = _t_ary[_n]+firstDelta*(i-_n);
    //std::cout << "Adding lower open:   " << i << " = " << _t_ary[i] << std::endl;
  }

  for (int i=0;i<_nrClose;i++) // add the lower closed knots
  {
    _t_ary[i] = _t_ary[_nrClose];
    //std::cout << "Adding lower closed: " << i << " = " << _t_ary[i] << std::endl;
  }


  double lastDelta=_t_ary[_m-_n-1]-_t_ary[_m-_n-2]; //extrapolate the upper non-closed knots
  for (int i=_m-_n;i<_m-_nrClose;i++) 
  {
    _t_ary[i] = _t_ary[_m-_n-1]+lastDelta*(i-(_m-_n-1));
    //std::cout << "Adding upper open:   " << i << " = " << _t_ary[i] << std::endl;
  }

  for (int i=_m-_nrClose;i<_m;i++) // add the upper closed knots
  {
    _t_ary[i] = _t_ary[_m-_nrClose-1];
    //std::cout << "Adding upper closed: " << i << " = " << _t_ary[i] << std::endl;
  }
  //std::cout << std::endl;

  for (int i=0;i<_m;i++)
  {
    if (fabs(_t_ary[i]) < pow(10., -9)) _t_ary[i] = 0;
  }
}

//_____________________________________________________________________________
RooBSplineBases::RooBSplineBases(const char* name, const char* title) :
  RooAbsReal(name, title)
{
  // Constructor of flat polynomial function
  //_bin=NULL;
  _bin.resize(_n+1);
  for (int i=0;i<_n+1;i++)
  {
    _bin[i].resize(_m);
  }
}

//_____________________________________________________________________________
RooBSplineBases::RooBSplineBases(const RooBSplineBases& other, const char* name) :
  RooAbsReal(other, name), 
  _tValues(other._tValues),
  _m(other._m),
  _t("t_param", this, other._t),
  _n(other._n),
  _nrClose(other._nrClose),
  //_nPlusOne(other._nPlusOne),
  _t_ary(other._t_ary),
  _bin(other._bin)
{
  // Copy constructor

  buildTAry();
  _bin.resize(_n+1);
  for (int i=0;i<_n+1;i++)
  {
    _bin[i].resize(_m);
  }
//   _t_ary = new double[_m];
//   for (int i=0;i<_m;i++)
//   {
//     _t_ary[i] = other._t_ary[i];
//   }
//   if (other._bin)
//   {
//     _bin=new double*[_n+1];
//     for (int n=0;n<_n+1;n++)
//     {
//       _bin[n] = new double[_m];
//       for (int i=0;i<_m;i++)
//       {
// 	if (other._bin[n])
// 	{
// 	  _bin[n][i] = other._bin[n][i];
// 	}
// 	else
// 	{
// 	  _bin[n][i] = 0;
// 	}
//       }
//     }
//   }
}


//_____________________________________________________________________________
RooBSplineBases::~RooBSplineBases() 
{
  // Destructor
  //delete[] _t_ary;
//   _t_ary=NULL;
//   if (_bin)
//   {
//     for (int i=0;i<_n+1;i++) 
//     {
//       delete _bin[i];
//     _bin[i]=NULL;
//     }
//     delete _bin;
//     _bin=NULL;
//   }
}




//_____________________________________________________________________________
Double_t RooBSplineBases::evaluate() const 
{
//   std::cout << "In eval, _n=" << _n << ", _m=" << _m << ", _nPlusOne=" << _nPlusOne << std::endl;
//   std::cout << "_bin=" << _bin << std::endl;
//   std::cout << "_t_ary=" << _t_ary << std::endl;
  if (!_t_ary.size()) buildTAry();

  // Calculate and return value of spline
  //std::cout << "In RooBSplineBases::evaluate()" << std::endl;
  double t = _t;
  if (t < _t_ary[_n] || t > _t_ary[_m-_n-1])
  {
    if (t > _t_ary[_m-_n-1]) t = _t_ary[_m-_n-1];
    if (t < _t_ary[_n]) t = _t_ary[_n];
  }




//build the basis splines

//   if (_bin)
//   {
//     //std::cout << "Deleting bin: " << _bin << std::endl;
//     for (int i=0;i<_n+1;i++)
//     {
//       //std::cout << "Deleting bin[" << i << "]: " << _bin[i] << std::endl;
//       delete[] _bin[i];
//       _bin[i]=NULL;
// //       for (int j=0;j<_m;j++)
// //       {
// // 	_bin[i][j]=0;
// //       }
//     }
//     delete[] _bin;
//     _bin=NULL;
//   }

//   bool remake=(_bin==NULL);
//   if (remake) _bin = new double*[_n+1];

//   if (!_bin)
//   {
//     _bin=new double*[_n+1];
//     for (int n=0;n<_n+1;n++)
//     {
//       _bin[n] = new double[_m];
//     }
//   }
  //_bin = new double*[_n+1];
  for (int n=0;n<_n+1;n++)
  {
    //std::cout << "_bin[n] = " << _bin[n] << std::endl;
    //if (remake) _bin[n] = new double[_m];
    //_bin[n] = new double[_m];
    for (int i=0;i<_m;i++)
    {
      //std::cout << "Resetting to zero" << std::endl;
      _bin[n][i] = 0;
    }
    for (int i=0;i<_m-n-1;i++)
    {
      if (n == 0)
      {
	if (t >= _t_ary[i] && t < _t_ary[i+1] && i >= _n && i <= _m-_n-1)
	{
	  //std::cout << "Setting " << i << " to 1" << std::endl;
	  _bin[n][i] = 1;
	}
      }
      else
      {
	//std::cout << "Getting term1" << std::endl;
	double term1 = 0;
	if (_t_ary[i+n] - _t_ary[i] > 0.000000000001) term1 = _bin[n-1][i] / (_t_ary[i+n] - _t_ary[i]);

	//std::cout << "Getting term2" << std::endl;
	double term2 = 0;
	if (_t_ary[i+n+1] - _t_ary[i+1] > 0.0000000000001) term2 = _bin[n-1][i+1] / (_t_ary[i+n+1] - _t_ary[i+1]);

	//std::cout << "Setting bin" << std::endl;
	_bin[n][i] = (t - _t_ary[i]) * term1 + (_t_ary[i+n+1] - t) * term2;
      }
      if (_bin[n][i] < 0.000000000001) _bin[n][i] = 0;
    }
  }
  //std::cout << "Out RooBSplineBases::evaluate()" << std::endl;
  return t;
}

Double_t RooBSplineBases::getBasisVal(int n, int i, bool rebuild) const
{
  if (rebuild) getVal();
//   if (rebuild || !_bin) getVal();
//   if (!_bin) 
//   {
//     getVal();
//   }
  if (i >= _m-_n-1) return 0.;
  //std::cout << "Getting basis for n=" << n << ", i=" << i << ", rebuild ? " << rebuild << ", order = " << _n << ", name=" << GetName() << std::endl;  
  return _bin[n][i];
}



// @(#)root/roostats:$Id: RooBSpline.cxx 873 2014-02-24 22:16:29Z adye $
// Author: Aaron Armbruster
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//_________________________________________________
/*
BEGIN_HTML
<p>
</p>
END_HTML
*/
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include <memory>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "TMath.h"

//#include "RooBSpline.h"

ClassImp(RooStats::HistFactory::RooBSpline)

using namespace RooStats;
using namespace HistFactory;

//_____________________________________________________________________________
RooBSpline::RooBSpline()
{
  // Default constructor
//   _t_ary=NULL;
}


//_____________________________________________________________________________
RooBSpline::RooBSpline(const char* name, const char* title,
		       const RooArgList& controlPoints, RooBSplineBases& bases, const RooArgSet& vars) :
  RooAbsReal(name, title),
  _controlPoints("controlPoints","List of control points",this),
  _m(bases.getTValues().size()+2*bases.getOrder()),
//   _t_ary(NULL),
//   _t("t_param", "t_param", this, *t),
  _n(bases.getOrder()),
  _weights("weights","List of weights",this),
  _bases("bases","Basis polynomials",this,bases),
  _vars("observables","List of observables",this),
  _cacheMgr(this,10)
{
  //std::cout << "in Ctor" << std::endl;
  //std::cout << "t = " << t << std::endl;

  if (_m-2*_n != controlPoints.getSize())
  {
    std::cout << "ERROR::Nr t values (" << _m-2*_n << ") != nr control points (" << controlPoints.getSize() << ")" << std::endl;
  }

  //bool even = fabs(_n/2-_n/2.) < 0.0000000001;
  bool first=1;
  TIterator* pointIter = controlPoints.createIterator() ;
  RooAbsArg* point ;
  RooAbsArg* lastPoint=NULL ;
  while((point = (RooAbsArg*)pointIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(point)) {
      coutE(InputArguments) << "RooBSpline::ctor(" << GetName() << ") ERROR: control point " << point->GetName() 
			    << " is not of type RooAbsReal" << std::endl ;
      assert(0) ;
    }
    //RooAbsReal* pointReal = (RooAbsReal*)point;
    //std::cout << "Adding control point " << point->GetName() << ", has val " << pointReal->getVal() << std::endl;
    _controlPoints.add(*point) ;
    if (first) 
    {
      for (int i=0;i<(_n)/2;i++)
      {
	_controlPoints.add(*point) ;
      }
    }
    first=false;
    lastPoint=point;
  }
  for (int i=0;i<(_n)/2;i++) _controlPoints.add(*lastPoint);
  delete pointIter ;


  TIterator* varItr = vars.createIterator();
  RooAbsArg* arg;
  while ((arg=(RooAbsArg*)varItr->Next())) {
    //std::cout << "======== Adding "<<arg->GetName()<<" to list of _vars of "<<name<<"." << std::endl;
    _vars.add(*arg);
  }
//   std::cout << "all _vars: " << std::endl;
//   _vars.Print("V");
  delete varItr;
  //std::cout << "out Ctor" << std::endl;
}

//_____________________________________________________________________________
RooBSpline::RooBSpline(const char* name, const char* title) :
  RooAbsReal(name, title),
  _controlPoints("controlPoints","List of coefficients",this),
  _cacheMgr(this,10)
{
  // Constructor of flat polynomial function

}

//_____________________________________________________________________________
RooBSpline::RooBSpline(const RooBSpline& other, const char* name) :
  RooAbsReal(other, name), 
  _controlPoints("controlPoints",this,other._controlPoints),
  _m(other._m),
//   _t_ary(NULL),
//   _t("t_param", this, other._t),
  _n(other._n),
  _weights("weights",this,other._weights),
  _bases("bases",this,other._bases),
  _vars("observables",this,other._vars),
  _cacheMgr(this,10)
{
  // Copy constructor
  
//   _t_ary = new double[_m];
//   for (int i=0;i<_m;i++)
//   {
//     _t_ary[i] = other._t_ary[i];
//   }
}


//_____________________________________________________________________________
RooBSpline::~RooBSpline() 
{
  // Destructor
//   delete _t_ary;
}




//_____________________________________________________________________________
Double_t RooBSpline::evaluate() const 
{
  //std::cout << "In RooBSpline::evaluate(): " << GetName() << std::endl;
  // Calculate and return value of spline

  //std::cout << "computing S" << std::endl;
  RooBSplineBases* bases = (RooBSplineBases*)&_bases.arg();
  bases->getVal(); // build the basis polynomials
  bool useWeight = _weights.getSize();
  double S = 0;
  //bool even = fabs(_n/2-_n/2.) < 0.0000000001;
  for (int i=0;i<_m-_n-1;i++)
  {
    //if (even && i <_m-_n-2) p=_n-1;
    double basis = bases->getBasisVal(_n,i,false);
    if (basis > 0)
    {
      int p=i;
      //if (even && i > 0) p=i-1;
      RooAbsReal* point = (RooAbsReal*)_controlPoints.at(p);
      //std::cout << "name=" << GetName() << ", point addy=" << point << std::endl;
      double weight = 1.0;
      if (useWeight)
      {
	RooAbsReal* weightVar = (RooAbsReal*)_weights.at(p);
	weight = weightVar->getVal();
      }
      S += basis * point->getVal() * weight;
    }
  }

  //std::cout << "Out RooBSpline::evaluate()" << std::endl;
  return S;
}



void RooBSpline::setWeights(const RooArgList& weights)
{
  _weights.removeAll();
  bool first=1;
  TIterator* pointIter = weights.createIterator() ;
  RooAbsArg* point ;
  RooAbsArg* lastPoint=NULL ;
  while((point = (RooAbsArg*)pointIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(point)) {
      coutE(InputArguments) << "RooBSpline::ctor(" << GetName() << ") ERROR: control point " << point->GetName() 
			    << " is not of type RooAbsReal" << std::endl ;
      assert(0) ;
    }
    _weights.add(*point) ;
    if (first) 
    {
      for (int i=0;i<_n/2;i++)
      {
	_weights.add(*point) ;
      }
    }
    first=false;
    lastPoint=point;
  }
  for (int i=0;i<(_n+1)/2;i++) _weights.add(*lastPoint);
  delete pointIter;
}




//_____________________________________________________________________________
Bool_t RooBSpline::setBinIntegrator(RooArgSet& allVars) 
{
  //std::cout << "In RooBSpline::setBinIntegrator" << std::endl;
  if(allVars.getSize()==1){
    RooAbsReal* temp = const_cast<RooBSpline*>(this);
    temp->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator")  ;
    int nbins = ((RooRealVar*) allVars.first())->numBins();
    temp->specialIntegratorConfig(kTRUE)->getConfigSection("RooBinIntegrator").setRealValue("numBins",nbins);
    return true;
  }else{
    std::cout << "Currently BinIntegrator only knows how to deal with 1-d "<<std::endl;
    return false;
  }
  return false;
}


//_____________________________________________________________________________
Int_t RooBSpline::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, 
					  const RooArgSet* normSet, const char* rangeName) const 
{
//   std::cout << "In RooBSpline["<<GetName()<<"]::getAnalyticalIntegralWN" << std::endl;
//   std::cout << "allVars:" << std::endl;
//   allVars.Print("V");
//   std::cout << "analVars:" << std::endl;
//   analVars.Print("V");
//   std::cout << "_vars:" << std::endl;
//   _vars.Print("V");
//   std::cout << "--- end ---" << std::endl;
  
  if (_forceNumInt) return 0 ;

  if (_vars.getSize()==0) return 1;
  
  if (matchArgs(allVars, analVars, *_vars.first())) {

      // From RooAddition:
      // check if we already have integrals for this combination of factors

      // we always do things ourselves -- actually, always delegate further down the line ;-)
      analVars.add(allVars);
      if( normSet ) analVars.add(*normSet);
      
      // check if we already have integrals for this combination of factors
      Int_t sterileIndex(-1);
      CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&analVars,&analVars,&sterileIndex,RooNameReg::ptr(rangeName));
      if (cache==0) {
         // we don't, so we make it right here....
         cache = new CacheElem;         
         for (int i=0;i<_m-_n-1;i++) {
            RooAbsReal* point = (RooAbsReal*)_controlPoints.at(i);
            cache->_I.addOwned( *point->createIntegral(analVars,rangeName) );
         }
      }

      Int_t code = _cacheMgr.setObj(&analVars,&analVars,(RooAbsCacheElement*)cache,RooNameReg::ptr(rangeName));
      return 2+code;
  }

  return 0;
}


//_____________________________________________________________________________
Double_t RooBSpline::analyticalIntegralWN(Int_t code, const RooArgSet* /*normSet*/,const char* rangeName) const 
{
  //std::cout << "In RooBSpline::analyticalIntegralWN" << std::endl;
  double integral = 0;
  if (code == 1)
  {
    return getVal();
  }
  else if (code >= 2)
  {
//     RooRealVar* obs = (RooRealVar*)_vars.first();
//     int nrBins = obs->getBins();
//     int initValue=obs->getBin();
// 
//     for (int ibin=0;ibin<nrBins;ibin++)
//     {
//       obs->setBin(ibin);
//       integral+=getVal();
//     }
// 
//     obs->setBin(initValue);



// 
//       // From RooAddition:
//       // check if we already have integrals for this combination of factors
// 
//       // we always do things ourselves -- actually, always delegate further down the line ;-)
//       RooArgSet analVars( _vars );
//       analVars.add(*normSet);
// 
//       Int_t sterileIndex(-1);
//       CacheElem* cache = (CacheElem*) _cacheMgr.getObj(&analVars,&analVars,&sterileIndex,RooNameReg::ptr(rangeName));
//       if (cache==0) {
//          // we don't, so we make it right here....
//          cache = new CacheElem;         
//          for (int i=0;i<_m-_n-1;i++) {
//             RooAbsReal* point = (RooAbsReal*)_controlPoints.at(i);
//             cache->_I.addOwned( *point->createIntegral(_vars,*normSet) );
//          }
//       }
// 


     // Calculate integral internally from appropriate integral cache
   
     // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
     CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-2);
     if (cache==0) {
       // cache got sterilized, trigger repopulation of this slot, then try again...
       //std::cout << "Cache got sterilized" << std::endl;
       std::auto_ptr<RooArgSet> vars( getParameters(RooArgSet()) );
       std::auto_ptr<RooArgSet> iset(  _cacheMgr.nameSet2ByIndex(code-2)->select(*vars) );
       RooArgSet dummy;
       Int_t code2 = getAnalyticalIntegral(*iset,dummy,rangeName);
       assert(code==code2); // must have revived the right (sterilized) slot...
       return analyticalIntegral(code2,rangeName);
     }
     assert(cache!=0);
   



     RooBSplineBases* bases = (RooBSplineBases*)&_bases.arg();
     bases->getVal(); // build the basis polynomials
     bool useWeight = _weights.getSize();
     double S = 0;
     //bool even = fabs(_n/2-_n/2.) < 0.0000000001;
     for (int i=0;i<_m-_n-1;i++)
     {
       //if (even && i <_m-_n-2) p=_n-1;
       double basis = bases->getBasisVal(_n,i,false);
       if (basis > 0)
       {
         int p=i;
         //if (even && i > 0) p=i-1;
         //RooAbsReal* point = (RooAbsReal*)_controlPoints.at(p);
         //std::cout << "name=" << GetName() << std::endl;
         double weight = 1.0;
         if (useWeight)
         {
           RooAbsReal* weightVar = (RooAbsReal*)_weights.at(p);
           weight = weightVar->getVal();
         }


         RooAbsReal* intReal = (RooAbsReal*)cache->_I.at(p);   //point->createIntegral(_vars,*normSet);
         S += basis * intReal->getVal() * weight;
         //std::cout << "adding "<<intReal->getVal()<<" to integral" << std::endl;
         //delete intReal;
       }
     }
     
     integral = S;
  }
  //std::cout << "Spline " << GetName() << " has integral " << integral << " obtained with code "<< code << std::endl;
  return integral;
}

// //_____________________________________________________________________________
// RooBSplinePenalty* RooBSpline::getRealPenalty(int k, RooRealVar* obs,RooRealVar* beta,const char* name) const
// {
//   if (name == "")
//   {
//     std::stringstream nameStr;
//     nameStr << GetName() << "_penalty";
//     name = nameStr.str().c_str();
//   }

//   RooArgList controlPoints;
//   for (int i=_n/2;i<_controlPoints.getSize()-(_n+1)/2;i++)
//   {
//     //std::cout << "adding point with val " << ((RooAbsReal*)_controlPoints.at(i))->getVal() << std::endl;
//     controlPoints.add(*_controlPoints.at(i));
//   }

//   std::vector<double> tValues;
//   for (int i=_n;i<_m-_n;i++)
//   {
//     tValues.push_back(_t_ary[i]);
//   }

//   RooBSplinePenalty* pen = new RooBSplinePenalty(name, name, _n, tValues, controlPoints, k, obs, beta);

//   int nrWeights = _weights.getSize();
//   if (nrWeights > 0)
//   {
//     RooArgSet weights;
//     int counter = 0;
//     for (int i=_n/2;i<nrWeights-(_n+1)/2;i++)
//     {
//       weights.add(*_weights.at(i));
//       counter++;
//     }
//     //std::cout << "added " << counter << " weights" << std::endl;
//     pen->setWeights(weights);
//   }

//   return pen;
// }






//_____________________________________________________________________________
RooArgList RooBSpline::CacheElem::containedArgs(Action)
{
  // Return list of all RooAbsArgs in cache element
  RooArgList ret(_I) ;
  return ret ;
}

RooBSpline::CacheElem::~CacheElem()
{
  // Destructor
}



/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooParamKeysPdf.cxx 888 2014-08-01 19:54:39Z adye $
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
#include "RooFit.h"

#include <math.h>
#include "Riostream.h"
#include "TMath.h"

//#include "RooParamKeysPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooDataSet.h"

using namespace std;

ClassImp(RooParamKeysPdf)


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooParamKeysPdf implements a one-dimensional kernel estimation p.d.f which model the distribution
// of an arbitrary input dataset as a superposition of Gaussian kernels, one for each data point,
// each contributing 1/N to the total integral of the p.d.f.
// <p>
// If the 'adaptive mode' is enabled, the width of the Gaussian is adaptively calculated from the
// local density of events, i.e. narrow for regions with high event density to preserve details and
// wide for regions with log event density to promote smoothness. The details of the general algorithm
// are described in the following paper: 
// <p>
// Cranmer KS, Kernel Estimation in High-Energy Physics.  
//             Computer Physics Communications 136:198-207,2001 - e-Print Archive: hep ex/0011057
// <p>
// END_HTML
//


//_____________________________________________________________________________
  RooParamKeysPdf::RooParamKeysPdf() : _centralValue(0.0), 
                             _nEvents(0), _dataPts(0), _dataWgts(0), _weights(0), _sumWgt(0), _normVal(0), _nPoints(1000),
			     _lookupTable(0), _mirrorLeft(kFALSE), _mirrorRight(kFALSE), 
			     _asymLeft(kFALSE), _asymRight(kFALSE)
{ 
  // coverity[UNINIT_CTOR]
}


//_____________________________________________________________________________
RooParamKeysPdf::RooParamKeysPdf(const char *name, const char *title,
                       RooAbsReal& x, RooAbsReal& deltax, RooDataSet& data,
		       Mirror mirror, Double_t rho, Int_t nPoints) :
  RooAbsPdf(name,title),
  _x("x","Dependent",this,x),
  _deltax("deltax","Dependent",this,deltax),
  
  _centralValue(0.0),
  
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),
  
  _nPoints(nPoints),
  _lookupTable(new Double_t[_nPoints]),

  _mirrorLeft(mirror==MirrorLeft || mirror==MirrorBoth || mirror==MirrorLeftAsymRight),
  _mirrorRight(mirror==MirrorRight || mirror==MirrorBoth || mirror==MirrorAsymLeftRight),
  _asymLeft(mirror==MirrorAsymLeft || mirror==MirrorAsymLeftRight || mirror==MirrorAsymBoth),
  _asymRight(mirror==MirrorAsymRight || mirror==MirrorLeftAsymRight || mirror==MirrorAsymBoth),
  _rho(rho)
{
  // cache stuff about x
  snprintf(_varName, 128,"%s", x.GetName());
  RooRealVar real= (RooRealVar&)(_x.arg());
  _lo = real.getMin();
  _hi = real.getMax();
  _binWidth = (_hi-_lo)/(_nPoints-1);

  // form the lookup table
  LoadDataSet(data);
}


//_____________________________________________________________________________
RooParamKeysPdf::RooParamKeysPdf(const char *name, const char *title,
                       RooAbsReal& x, RooAbsReal& deltax, double centralValue, RooAbsReal& multiplicativeShift, RooDataSet& data,
		       Mirror mirror, Double_t rho, Int_t nPoints) :
  RooAbsPdf(name,title),
  _x("x","Dependent",this,x),
  _deltax("deltax","Dependent",this,deltax),  
  
  _centralValue(centralValue),
  _multiplicativeShift("multiplicativeShift","Dependent",this,multiplicativeShift),
  
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),

  _nPoints(nPoints),
  _lookupTable(new Double_t[_nPoints]),

  _mirrorLeft(mirror==MirrorLeft || mirror==MirrorBoth || mirror==MirrorLeftAsymRight),
  _mirrorRight(mirror==MirrorRight || mirror==MirrorBoth || mirror==MirrorAsymLeftRight),
  _asymLeft(mirror==MirrorAsymLeft || mirror==MirrorAsymLeftRight || mirror==MirrorAsymBoth),
  _asymRight(mirror==MirrorAsymRight || mirror==MirrorLeftAsymRight || mirror==MirrorAsymBoth),
  _rho(rho)
{
  // cache stuff about x
  snprintf(_varName, 128,"%s", x.GetName());
  RooRealVar real= (RooRealVar&)(_x.arg());
  _lo = real.getMin();
  _hi = real.getMax();
  _binWidth = (_hi-_lo)/(_nPoints-1);

  // form the lookup table
  LoadDataSet(data);
}




//_____________________________________________________________________________
RooParamKeysPdf::RooParamKeysPdf(const RooParamKeysPdf& other, const char* name):
  RooAbsPdf(other,name), _x("x",this,other._x), _deltax("deltax",this,other._deltax), 
  _centralValue(other._centralValue),
  _multiplicativeShift("multiplicativeShift",this,other._multiplicativeShift),
  _nEvents(other._nEvents),
  _dataPts(0), _dataWgts(0), _weights(0), _sumWgt(0), _normVal(other._normVal),
  _nPoints(other._nPoints), _lookupTable(new Double_t[_nPoints]),
  _mirrorLeft( other._mirrorLeft ), _mirrorRight( other._mirrorRight ),
  _asymLeft(other._asymLeft), _asymRight(other._asymRight),
  _rho( other._rho ) {

  // cache stuff about x
  snprintf(_varName, 128, "%s", other._varName );
  _lo = other._lo;
  _hi = other._hi;
  _binWidth = other._binWidth;

  // copy over data and weights... not necessary, commented out for speed ...
//    _dataPts = new Double_t[_nEvents];
//    _dataWgts = new Double_t[_nEvents];
//    _weights = new Double_t[_nEvents];  
//    for (Int_t i= 0; i<_nEvents; i++) {
//      _dataPts[i]= other._dataPts[i];
//      _dataWgts[i]= other._dataWgts[i];
//      _weights[i]= other._weights[i];
//    }

  // copy over the lookup table
  for (Int_t i= 0; i<_nPoints; i++)
    _lookupTable[i]= other._lookupTable[i];
  
}


//_____________________________________________________________________________
RooParamKeysPdf::~RooParamKeysPdf() {
  delete[] _dataPts;
  delete[] _dataWgts;
  delete[] _weights;
  delete[] _lookupTable;
}


void

//_____________________________________________________________________________
RooParamKeysPdf::LoadDataSet( RooDataSet& data) {
  delete[] _dataPts;
  delete[] _dataWgts;
  delete[] _weights;

  // make new arrays for data and weights to fill
  _nEvents= (Int_t)data.numEntries();
  if (_mirrorLeft) _nEvents += data.numEntries();
  if (_mirrorRight) _nEvents += data.numEntries();

  _dataPts  = new Double_t[_nEvents];
  _dataWgts = new Double_t[_nEvents];
  _weights  = new Double_t[_nEvents];
  _sumWgt = 0 ;
  _normVal = 0 ;

  Double_t x0(0);
  Double_t x1(0);
  Double_t x2(0);

  Int_t i, idata=0;
  for (i=0; i<data.numEntries(); i++) {
    const RooArgSet *values= data.get(i);
    RooRealVar real= (RooRealVar&)(values->operator[](_varName));

    _dataPts[idata]= real.getVal();
    _dataWgts[idata] = data.weight() ;
    x0 += _dataWgts[idata] ; x1+=_dataWgts[idata]*_dataPts[idata]; x2+=_dataWgts[idata]*_dataPts[idata]*_dataPts[idata];
    idata++;
    _sumWgt+= data.weight() ;

    if (_mirrorLeft) {
      _dataPts[idata]= 2*_lo - real.getVal();
      _dataWgts[idata]= data.weight() ;
      _sumWgt+= data.weight() ;
      idata++;
    }

    if (_mirrorRight) {
      _dataPts[idata]  = 2*_hi - real.getVal();
      _dataWgts[idata] = data.weight() ;
      _sumWgt+= data.weight() ;
      idata++;
    }
  }

  Double_t meanv=x1/x0;
  Double_t sigmav=sqrt(x2/x0-meanv*meanv);
  Double_t h=TMath::Power(Double_t(4)/Double_t(3),0.2)*TMath::Power(_nEvents,-0.2)*_rho;
  Double_t hmin=h*sigmav*sqrt(2.)/10;
  Double_t norm=h*sqrt(sigmav)/(2.0*sqrt(3.0));

  _weights=new Double_t[_nEvents];
  for(Int_t j=0;j<_nEvents;++j) {
    _weights[j]=norm/sqrt(g(_dataPts[j],h*sigmav));
    if (_weights[j]<hmin) _weights[j]=hmin;
  }
  
  for (i=0;i<_nPoints;++i) {
    _lookupTable[i]=evaluateFull( _lo+Double_t(i)*_binWidth );
    _normVal += _lookupTable[i];
  }
  _normVal *= _binWidth;

  
}



//_____________________________________________________________________________
Double_t RooParamKeysPdf::evaluate() const {
  double deltax = _deltax-_centralValue;
  if( _multiplicativeShift.absArg() ) deltax += _centralValue*(_multiplicativeShift-1.0);

  Int_t i = (Int_t)floor((Double_t(_x-deltax)-_lo)/_binWidth);
  bool forcePositive = false; // this is just to suppress warning for values outside of range.
  if (i<0) {
//     std::cerr << "got point below lower bound:"
// 	 << Double_t(_x) << " < " << _lo
// 	 << " -- performing linear extrapolation..." << std::endl;
    i=0;
    forcePositive = true;
  }
  if (i>_nPoints-2) {
//     std::cerr << "got point above upper bound:"
// 	 << Double_t(_x) << " > " << _hi
// 	 << " -- performing linear extrapolation..." << std::endl;
    i=_nPoints-2;
    forcePositive = true;
  }
  Double_t dx = (Double_t(_x-deltax)-(_lo+i*_binWidth))/_binWidth;
  
  // for now do simple linear interpolation.
  // one day replace by splines...
  double val = (_lookupTable[i]+dx*(_lookupTable[i+1]-_lookupTable[i]));
  if( forcePositive && val<0.0 ) val = 0.0; 
  return val;
}


//_____________________________________________________________________________
Double_t RooParamKeysPdf::evaluateFull( Double_t x ) const {
  Double_t y=0;

  for (Int_t i=0;i<_nEvents;++i) {
    Double_t chi=(x-_dataPts[i])/_weights[i];
    y+=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];

    // if mirroring the distribution across either edge of
    // the range ("Boundary Kernels"), pick up the additional
    // contributions
//      if (_mirrorLeft) {
//        chi=(x-(2*_lo-_dataPts[i]))/_weights[i];
//        y+=exp(-0.5*chi*chi)/_weights[i];
//      }
    if (_asymLeft) {
      chi=(x-(2*_lo-_dataPts[i]))/_weights[i];
      y-=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];
    }
//      if (_mirrorRight) {
//        chi=(x-(2*_hi-_dataPts[i]))/_weights[i];
//        y+=exp(-0.5*chi*chi)/_weights[i];
//      }
    if (_asymRight) {
      chi=(x-(2*_hi-_dataPts[i]))/_weights[i];
      y-=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];
    }
  }
  
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));  
  return y/(sqrt2pi*_sumWgt);
}


//_____________________________________________________________________________
Double_t RooParamKeysPdf::g(Double_t x,Double_t sigmav) const {
  
  Double_t c=Double_t(1)/(2*sigmav*sigmav);

  Double_t y=0;
  for (Int_t i=0;i<_nEvents;++i) {
    Double_t r=x-_dataPts[i];
    y+=exp(-c*r*r);
  }
  
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));  
  return y/(sigmav*sqrt2pi*_nEvents);
}




//_____________________________________________________________________________
Int_t RooParamKeysPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,_x)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t RooParamKeysPdf::analyticalIntegral(Int_t code, const char* /*rangeName*/) const 
{
   assert(code==1) ;
   //if( rangeName == NULL ) return _sumWgt;

//    double sumW = 0.0;
//    for( Int_t i=0; i < _nEvents; i++ ) {
//       if( rangeName == NULL  ||
//          (_dataPts[i] > _x.min(rangeName)  && _dataPts[i] < _x.max(rangeName))
//       ) {
//          sumW += _dataWgts[i]/_weights[i];
//       }
//    }
//    return sumW;
   if (_normVal == 0.0) {
      for( Int_t i=0; i < _nPoints; i++ ) {
         _normVal += _lookupTable[i];
      }
      _normVal *= _binWidth;
//    std::cout << ClassName() << "::" << GetName() << " analyticalIntegral = " << _normVal << std::endl;
   }
  return _normVal;
}




/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

//#include "RooStarMomentMorph.h" 
#include "RooAbsCategory.h" 
#include "RooRealIntegral.h"
#include "RooRealConstant.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCustomizer.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,19)
#include "RooAbsMoment.h"
#else
#include "RooMoment.h"
#endif
#include "RooLinearVar.h"
#include "RooChangeTracker.h"
#include "RooNumIntConfig.h"
#include "RooHistPdf.h"

#include "TMath.h"
#include "TVector.h"

ClassImp(RooStarMomentMorph) 


//_____________________________________________________________________________
RooStarMomentMorph::RooStarMomentMorph() : 
  _curNormSet(0), 
  _M(0), 
  _setting(RooStarMomentMorph::Linear), 
  _nnuisvar(0),
  _useHorizMorph(true)
{

  // coverity[UNINIT_CTOR]
  _parItr    = _parList.createIterator() ;
  _obsItr    = _obsList.createIterator() ;
  _pdfItr    = _pdfList.createIterator() ; 
}



//_____________________________________________________________________________
RooStarMomentMorph::RooStarMomentMorph(const char *name, const char *title, 
				       const RooArgList& parList,
				       const RooArgList& obsList,
				       const RooArgList& pdfList, // {dn0,up0,dn1,up1,...,nominal}
				       const std::vector<int>& nnuispoints, // number of variation (not counting nominal) for each NP (2,2)
				       const std::vector<double>& nrefpoints, // (-1,1,-1,1)
				       const Setting& setting) :
  RooAbsPdf(name,title), 
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,19)
  _cacheMgr(this,10,kTRUE,kTRUE),
#else
  _cacheMgr(this,10,kFALSE),
#endif
  _parList("parList","List of fit parameters",this),
  _obsList("obsList","List of variables",this),
  _pdfList("pdfList","List of pdfs",this),
  _nnuis(nnuispoints),
  _nref(nrefpoints),
  _M(0),
  _setting(setting),
  _nnuisvar(0),
  _useHorizMorph(true)
{ 
  // CTOR

  // fit parameters
  TIterator* parItr = parList.createIterator() ;
  RooAbsArg* par ;
  for (Int_t i=0; (par = (RooAbsArg*)parItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(par)) {
      coutE(InputArguments) << "RooStarMomentMorph::ctor(" << GetName() << ") ERROR: parameter " << par->GetName() << " is not of type RooAbsReal" << endl ;
      throw std::string("RooStarMomentMorh::ctor() ERROR parameter is not of type RooAbsReal") ;
    }
    _parList.add(*par) ;
  }
  delete parItr ;

  // observables
  TIterator* obsItr = obsList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)obsItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooStarMomentMorph::ctor(" << GetName() << ") ERROR: variable " << var->GetName() << " is not of type RooAbsReal" << endl ;
      throw std::string("RooStarMomentMorh::ctor() ERROR variable is not of type RooAbsReal") ;
    }
    _obsList.add(*var) ;
  }
  delete obsItr ;

  // reference p.d.f.s
  TIterator* pdfItr = pdfList.createIterator() ;
  RooAbsPdf* pdf ;
  for (Int_t i=0; (pdf = dynamic_cast<RooAbsPdf*>(pdfItr->Next())); ++i) {
    if (!pdf) {
      coutE(InputArguments) << "RooStarMomentMorph::ctor(" << GetName() << ") ERROR: pdf " << pdf->GetName() << " is not of type RooAbsPdf" << endl ;
      throw std::string("RooStarMomentMorph::ctor() ERROR pdf is not of type RooAbsPdf") ;
    }
    _pdfList.add(*pdf) ;
  }
  delete pdfItr ;

  _parItr    = _parList.createIterator() ;
  _obsItr    = _obsList.createIterator() ;
  _pdfItr    = _pdfList.createIterator() ; 


  // initialization
  initialize();
} 


//_____________________________________________________________________________
RooStarMomentMorph::RooStarMomentMorph(const RooStarMomentMorph& other, const char* name) :  
  RooAbsPdf(other,name), 
  _cacheMgr(other._cacheMgr,this),
  _curNormSet(0),
  _parList("parList",this,other._parList),
  _obsList("obsList",this,other._obsList),
  _pdfList("pdfList",this,other._pdfList),
  _nnuis(other._nnuis),
  _nref(other._nref),
  _M(0),
  _setting(other._setting),
  _nnuisvar(other._nnuisvar),
  _useHorizMorph(other._useHorizMorph)
{ 
 

  _parItr    = _parList.createIterator() ;
  _obsItr    = _obsList.createIterator() ;
  _pdfItr    = _pdfList.createIterator() ; 

  // nref is resized in initialize, so reduce the size here
  if (_nref.size()>0) {
    _nref.resize( _nref.size()-1 );
  }

  // initialization
  initialize();
} 

//_____________________________________________________________________________
RooStarMomentMorph::~RooStarMomentMorph() 
{
  if (_parItr) delete _parItr;
  if (_obsItr) delete _obsItr;
  if (_pdfItr) delete _pdfItr;
  if (_M)      delete _M;
}



//_____________________________________________________________________________
void RooStarMomentMorph::initialize() 
{

  unsigned int nPar = _parList.getSize();
  unsigned int nPdf = _pdfList.getSize();

  // other quantities needed
  if (nPar!=_nnuis.size()) {
    coutE(InputArguments) << "RooStarMomentMorph::initialize(" << GetName() << ") ERROR: nPar != nNuis" << endl ;
    assert(0) ;
  }

  // total number of nuisance parameter variations
  _nnuisvar = 0; 
  for (unsigned int k=0; k<_nnuis.size(); ++k) { _nnuisvar += _nnuis[k]; }

  // check of minimal number inputs
  for (unsigned int k=0; k<_nnuis.size(); ++k) { 
    if (_nnuis[k]<1) {
      coutE(InputArguments) << "RooStarMomentMorph::initialize(" << GetName() << ") ERROR: nNuis < 1 for variation: " << k << endl ;
      assert(0);
    }
  }



  // check number of input pdfs
  if (nPdf!=_nnuisvar+1) {
    coutE(InputArguments) << "RooStarMomentMorph::initialize(" << GetName() << ") ERROR: nPdf != nNuis variations" << endl ;
    assert(0) ;
  }

  // bit of a hack: last element corresponds to nnuis value of nominal pdf, ala pdflist
  _nref.resize( _nref.size()+1 );
  _nref[_nnuisvar] = 0.;


  //// MB : old left-over code ...
  //   TVector* dm = new TVector(nPdf);
  //   _M = new TMatrixD(nPdf,nPdf);
  
  //   // transformation matrix for non-linear extrapolation, needed in evaluate()
  //   TMatrixD M(nPdf,nPdf);
  //   for (unsigned int i=0; i<_nref->size(); ++i) {
  //     (*dm)[i] = (*_nref)[i]-(*_nref)[0];
  //     M(i,0) = 1.;
  //     if (i>0) M(0,i) = 0.;
  //   }
  //   for (unsigned int i=1; i<_nref->size(); ++i) {
  //     for (unsigned int j=1; j<_nref->size(); ++j) {
  //       M(i,j) = TMath::Power((*dm)[i],(double)j);
  //     }
  //   }
  //   (*_M) = M.Invert();
  
  //   delete dm ;
}

//_____________________________________________________________________________
RooStarMomentMorph::CacheElem* RooStarMomentMorph::getCache(const RooArgSet* /*nset*/) const
{
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(0,(RooArgSet*)0) ;
  if (cache) {
    return cache ;
  }


  Int_t nObs = _obsList.getSize();
  Int_t nPdf = _pdfList.getSize();

  RooAbsReal* null = 0 ;
  std::vector<RooAbsReal*> meanrv(nPdf*nObs,null);
  std::vector<RooAbsReal*> sigmarv(nPdf*nObs,null); 
  std::vector<RooAbsReal*> myrms(nObs,null);      
  std::vector<RooAbsReal*> mypos(nObs,null);      
  std::vector<RooAbsReal*> slope(nPdf*nObs,null); 
  std::vector<RooAbsReal*> offsetVar(nPdf*nObs,null); 
  std::vector<RooAbsReal*> transVar(nPdf*nObs,null); 
  std::vector<RooAbsReal*> transPdf(nPdf,null);      

  RooArgSet ownedComps ;

  RooArgList fracl ;

  // fraction parameters
  RooArgList coefList("coefList");    // fractions multiplied with input pdfs
  RooArgList coefList2("coefList2");  // fractions multiplied with mean position of observable contribution
  RooArgList coefList3("coefList3");  // fractions multiplied with rms position of observable contribution

  for (Int_t i=0; i<3*nPdf; ++i) {

    std::string fracName = Form("frac_%d",i);
    double initval=(i==nPdf-1||i==2*nPdf-1||i==3*nPdf-1) ? 1. : 0.;
    RooRealVar* frac=new RooRealVar(fracName.c_str(),fracName.c_str(),initval); // to be set later 

    fracl.add(*frac); 
    if      (i<nPdf)   coefList .add(*frac) ;
    else if (i<2*nPdf) coefList2.add(*frac) ;
    else               coefList3.add(*frac) ;
    ownedComps.add(*frac) ;
  }

  if (_useHorizMorph) {
    
    // mean and sigma
    RooArgList obsList(_obsList);
    for (Int_t i=0; i<nPdf; ++i) {      
      for (Int_t j=0; j<nObs; ++j) {
	
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,19)
	RooAbsMoment* mom = nObs==1 ?
#else
	RooMoment* mom = nObs==1 ?
#endif
          ((RooAbsPdf*)_pdfList.at(i))->sigma((RooRealVar&)*obsList.at(j)) :
          ((RooAbsPdf*)_pdfList.at(i))->sigma((RooRealVar&)*obsList.at(j), obsList);

#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,19)
        mom->setLocalNoDirtyInhibit(kTRUE);
        mom->mean()->setLocalNoDirtyInhibit(kTRUE);
#endif
	
        sigmarv[sij(i,j)] = mom;
        meanrv [sij(i,j)] = mom->mean();
        ownedComps.add(*sigmarv[sij(i,j)]);
	
      }
    }    
  
    // slope and offset (to be set later, depend on nuisance parameters)
    for (Int_t j=0; j<nObs; ++j) {
      
      RooArgList meanList("meanList");
      RooArgList rmsList("rmsList");
      
      for (Int_t i=0; i<nPdf; ++i) {
	meanList.add(*meanrv [sij(i,j)]);
	rmsList. add(*sigmarv[sij(i,j)]);
      }
      
      std::string myrmsName = Form("%s_rms_%d",GetName(),j);
      std::string myposName = Form("%s_pos_%d",GetName(),j);
      
      mypos[j] = new RooAddition(myposName.c_str(),myposName.c_str(),meanList,coefList2);
      myrms[j] = new RooAddition(myrmsName.c_str(),myrmsName.c_str(),rmsList,coefList3);
      
      ownedComps.add(RooArgSet(*myrms[j],*mypos[j])) ;
    }
  }

  //for (Int_t i=0;i<_parList.getSize();i++) {    
  //RooRealVar* par=(RooRealVar*)_parList.at(i);
  //par->Print();
  //}

  // construction of unit pdfs
  _pdfItr->Reset();
  RooAbsPdf* pdf;
  RooArgList transPdfList;

  for (Int_t i=0; i<nPdf; ++i) {
    _obsItr->Reset() ;
    RooRealVar* var ;

    pdf = (RooAbsPdf*)_pdfItr->Next();

    std::string pdfName = Form("pdf_%d",i);
    RooCustomizer cust(*pdf,pdfName.c_str());

    for (Int_t j=0; j<nObs; ++j) {

      // slope and offset formulas
      std::string slopeName  = Form("%s_slope_%d_%d", GetName(),i,j);
      std::string offsetName = Form("%s_offset_%d_%d",GetName(),i,j);

      slope[sij(i,j)]  = _useHorizMorph ? 
	(RooAbsReal*)new RooFormulaVar(slopeName.c_str(),"@0/@1",
				       RooArgList(*sigmarv[sij(i,j)],*myrms[j])) : 
	(RooAbsReal*)new RooConstVar(slopeName.c_str(),slopeName.c_str(),1.0); // 
      
      
      offsetVar[sij(i,j)] = _useHorizMorph ? 
	(RooAbsReal*)new RooFormulaVar(offsetName.c_str(),"@0-(@1*@2)",
				       RooArgList(*meanrv[sij(i,j)],*mypos[j],*slope[sij(i,j)])) : 	
	(RooAbsReal*)new RooConstVar(offsetName.c_str(),offsetName.c_str(),0.0); // 
      
      ownedComps.add(RooArgSet(*slope[sij(i,j)],*offsetVar[sij(i,j)])) ;
      
      // linear transformations, so pdf can be renormalized easily
      var = (RooRealVar*)(_obsItr->Next());
      std::string transVarName = Form("%s_transVar_%d_%d",GetName(),i,j);
      transVar[sij(i,j)] = new RooLinearVar(transVarName.c_str(),transVarName.c_str(),*var,*slope[sij(i,j)],*offsetVar[sij(i,j)]);

      // *** WVE this is important *** this declares that frac effectively depends on the morphing parameters
      // This will prevent the likelihood optimizers from erroneously declaring terms constant
      transVar[sij(i,j)]->addServerList((RooAbsCollection&)_parList);

      ownedComps.add(*transVar[sij(i,j)]) ;
      cust.replaceArg(*var,*transVar[sij(i,j)]);
    }
    transPdf[i] = (RooAbsPdf*) cust.build() ;
    transPdfList.add(*transPdf[i]);

    ownedComps.add(*transPdf[i]) ;
  }
  // sum pdf
  
  std::string sumpdfName = Form("%s_sumpdf",GetName());

  RooAbsPdf* theSumPdf = new RooAddPdf(sumpdfName.c_str(),sumpdfName.c_str(),transPdfList,coefList);
  theSumPdf->setAttribute("NeverConstant") ;
  theSumPdf->addOwnedComponents(ownedComps) ;
  
  // change tracker for fraction parameters
  std::string trackerName = Form("%s_frac_tracker",GetName()) ;
  RooChangeTracker* tracker = new RooChangeTracker(trackerName.c_str(),trackerName.c_str(),_parList,kTRUE) ;

  // Store it in the cache
  cache = new CacheElem(*theSumPdf,*tracker,fracl) ;
  _cacheMgr.setObj(0,0,cache,0) ;

  //const_cast<RooStarMomentMorph*>(this)->addServer(*theSumPdf,kFALSE,kFALSE) ;
  coutI(InputArguments) << "RooStarMorphPdf::getCache("<< GetName() << ") created cache element"<< endl ;

  return cache ;
}



//_____________________________________________________________________________
RooArgList RooStarMomentMorph::CacheElem::containedArgs(Action) 
{
  return RooArgList(*_sumPdf,*_tracker) ; 
}



//_____________________________________________________________________________
RooStarMomentMorph::CacheElem::~CacheElem() 
{ 
  delete _sumPdf ; 
  delete _tracker ; 
  //if (_owner) _owner->removeServer(*_sumPdf) ;
} 



//_____________________________________________________________________________
Double_t RooStarMomentMorph::getVal(const RooArgSet* set) const 
{
  // Special version of getVal() overrides RooAbsReal::getVal() to save value of current normalization set
  _curNormSet = set ? (RooArgSet*)set : (RooArgSet*)&_obsList ;
  return RooAbsPdf::getVal(set) ;
}



//_____________________________________________________________________________
RooAbsPdf* RooStarMomentMorph::sumPdf(const RooArgSet* nset) 
{
  CacheElem* cache = getCache(nset ? nset : _curNormSet) ;
  
  if (!cache->_fractionsCalculated || cache->_tracker->hasChanged(kTRUE)) {
    cache->calculateFractions(*this,kFALSE); // verbose turned off
  } 
  
  return cache->_sumPdf ;
}


//_____________________________________________________________________________
Double_t RooStarMomentMorph::evaluate() const 
{ 
  CacheElem* cache = getCache(_curNormSet) ;

  if (!cache->_fractionsCalculated || cache->_tracker->hasChanged(kTRUE)) {
    cache->calculateFractions(*this,kFALSE); // verbose turned off
  } 

  Double_t ret = cache->_sumPdf->getVal(_pdfList.nset());
  if (ret<0.) ret=0.;
  /*std::cout<<"(RSMM) "<<this->GetName();
  _parItr->Reset();  
  for (Int_t j=0; j<_parList.getSize(); j++) {    
    RooRealVar* m = (RooRealVar*)(_parItr->Next());
    std::cout<<", par="<<m->getVal();
  }

  std::cout<<", m4l="<<_obsList.getRealValue("m4l");
  std::cout<<":  "<<ret<<std::endl;
  //exit(3);
  */

  return ret ;
}


//_____________________________________________________________________________
RooRealVar* RooStarMomentMorph::CacheElem::frac(Int_t i ) 
{ 
  return (RooRealVar*)(_frac.at(i))  ; 
}



//_____________________________________________________________________________
const RooRealVar* RooStarMomentMorph::CacheElem::frac(Int_t i ) const 
{ 
  return (RooRealVar*)(_frac.at(i))  ; 
}


//_____________________________________________________________________________
void RooStarMomentMorph::CacheElem::calculateFractions(const RooStarMomentMorph& self, Bool_t /* verbose */) const
{
  _fractionsCalculated=true;

  switch (self._setting) {
  case Linear:
  case SineLinear:
    {
      //int nObs=self._obsList.getSize();
      
      // loop over parList
      self._parItr->Reset();
      int nnuis=0;
      
      // zero all fractions
      int nPdf=self._pdfList.getSize();	
      for (Int_t i=0; i<3*nPdf; ++i) {
	double initval=(i==nPdf-1||i==2*nPdf-1||i==3*nPdf-1) ? 1. : 0.;
	((RooRealVar*)frac(i))->setVal(initval);      
      }
      for (Int_t j=0; j<self._parList.getSize(); j++) {
	
	RooRealVar* m = (RooRealVar*)(self._parItr->Next());
	double m0=m->getVal();
	
	if (m0==0.) continue;

	int imin = self.ijlo(j,m0);
	int imax = self.ijhi(j,m0);
	
	double mlo=self._nref[imin];
	double mhi=self._nref[imax];
	
	// get reference for this obs
	nnuis+=self._nnuis[j];
	       
	double mfrac = (imax>imin) ? (mhi-m0)/(mhi-mlo) : (mlo-m0)/(mhi-mlo);
	if (mfrac> 1.) mfrac= 1.;
	if (mfrac<-1.) mfrac=-1.;
	
	if (self._setting==SineLinear) {
	  mfrac = TMath::Sin( TMath::PiOver2()*mfrac ); // this gives a continuous differentiable transition between grid points. 
	}

	((RooRealVar*)frac(imin))->setVal(((RooRealVar*)frac(imin))->getVal()+mfrac); 
	((RooRealVar*)frac(nPdf+imin))->setVal(((RooRealVar*)frac(nPdf+imin))->getVal()+mfrac); 
	((RooRealVar*)frac(2*nPdf+imin))->setVal(((RooRealVar*)frac(2*nPdf+imin))->getVal()+mfrac); 

	((RooRealVar*)frac(imax))->setVal(((RooRealVar*)frac(imax))->getVal()-mfrac); 
	((RooRealVar*)frac(nPdf+imax))->setVal(((RooRealVar*)frac(nPdf+imax))->getVal()-mfrac); 
	((RooRealVar*)frac(2*nPdf+imax))->setVal(((RooRealVar*)frac(2*nPdf+imax))->getVal()-mfrac); 	
      }

      break;
    }
  default:
    {
      std::cerr << "RooStarMomentMorph::calculateFractions() ERROR: only Linear Setting implemented!" << std::endl;
      throw std::string("RooStarMomentMorph::calculateFractions() ERROR only Linear Setting implemented") ;   
      break;
    }
  }
  
  //// MB : As an example: old left-over code. 
  

  //   Int_t nPdf = self._pdfList.getSize();
  
  //   // HACK: for now, for code to compile
  //   Double_t selfm = 0.0;
  
  //   Double_t dm = selfm - (*self._nref)[0];
  
  //   // fully non-linear
  //   double sumposfrac=0.;
  //   for (Int_t i=0; i<nPdf; ++i) {
  //     double ffrac=0.;
  //     for (Int_t j=0; j<nPdf; ++j) { ffrac += (*self._M)(j,i) * (j==0?1.:TMath::Power(dm,(double)j)); }
  //     if (ffrac>=0) sumposfrac+=ffrac;
  //     // fractions for pdf
  //     ((RooRealVar*)frac(i))->setVal(ffrac);
  //     // fractions for rms and mean
  //     ((RooRealVar*)frac(nPdf+i))->setVal(ffrac);
  //     if (verbose) { std::cout << ffrac << std::endl; }
  //   }
  
  //   // various mode settings
  //   int imin = self.idxmin(selfm);
  //   int imax = self.idxmax(selfm);
  //   double mfrac = (selfm-(*self._nref)[imin])/((*self._nref)[imax]-(*self._nref)[imin]);
  //   switch (self._setting) {
  //     case NonLinear:
  //       // default already set above
  //     break;
  //     case Linear: 
  //       for (Int_t i=0; i<2*nPdf; ++i)
  //         ((RooRealVar*)frac(i))->setVal(0.);      
  //       if (imax>imin) { // m in between mmin and mmax
  //         ((RooRealVar*)frac(imin))->setVal(1.-mfrac); 
  //         ((RooRealVar*)frac(nPdf+imin))->setVal(1.-mfrac);
  //         ((RooRealVar*)frac(imax))->setVal(mfrac);
  //         ((RooRealVar*)frac(nPdf+imax))->setVal(mfrac);
  //       } else if (imax==imin) { // m outside mmin and mmax
  //         ((RooRealVar*)frac(imin))->setVal(1.);
  //         ((RooRealVar*)frac(nPdf+imin))->setVal(1.);
  //       }
  //     break;
  //     case NonLinearLinFractions:
  //       for (Int_t i=0; i<nPdf; ++i)
  //         ((RooRealVar*)frac(i))->setVal(0.);
  //       if (imax>imin) { // m in between mmin and mmax
  //         ((RooRealVar*)frac(imin))->setVal(1.-mfrac);
  //         ((RooRealVar*)frac(imax))->setVal(mfrac);
  //       } else if (imax==imin) { // m outside mmin and mmax
  //         ((RooRealVar*)frac(imin))->setVal(1.);
  //       }
  //     break;
  //     case NonLinearPosFractions:
  //       for (Int_t i=0; i<nPdf; ++i) {
  //         if (((RooRealVar*)frac(i))->getVal()<0) ((RooRealVar*)frac(i))->setVal(0.);
  //         ((RooRealVar*)frac(i))->setVal(((RooRealVar*)frac(i))->getVal()/sumposfrac);
  //       }
  //     break;
  //   }
  

}



//_____________________________________________________________________________
Int_t RooStarMomentMorph::ij(const Int_t& i, const Int_t& j) const
{
  // i = nuisance parameter
  // j = variation number

  if (i<0 || i>=static_cast<int>(_nnuis.size())) return -1;
  if (j<0 || j>=_nnuis[i]) return -1;

  Int_t idx = 0;  
  for (int k=0; k<i; ++k) { idx += _nnuis[k]; }
  idx += j;

  return idx;
}


//_____________________________________________________________________________
Int_t RooStarMomentMorph::ijhi(const Int_t& i, const double& nval) const
{
  // find higher index on right side of nval
  // if lo==hi, pass the next-to highest index in line

  // i = nuisance parameter
  if (i<0 || i>=static_cast<int>(_nnuis.size())) return -1;

  // initialize
  int ijxlo = ij(i,0);
  double nlo=_nref[ijxlo];
  if ( _nref[_nnuisvar]<nlo ) { ijxlo=_nnuisvar; }

  int ijxhi = ij(i,_nnuis[i]-1);
  double nhi=_nref[ijxhi];
  if (_nref[_nnuisvar]>nhi ) { 
    ijxhi=_nnuisvar; 
    nhi = _nref[_nnuisvar];
  }

  int ijxhiprev = ij(i,0);
  for (Int_t j=0; j<_nnuis[i]; ++j) {
    int k = ij(i,j);
    if ( _nref[k]>=nval && _nref[k]<nhi ) { 
      nhi=_nref[k]; ijxhi=k;
    } else {
      ijxhiprev=k;
      break;
    }
  }
  if ( _nref[_nnuisvar]<nhi && _nref[_nnuisvar]>=nval) { 
    ijxhiprev = ijxhi;    
    ijxhi=_nnuisvar; 
  }
  if ( ijxhi!=static_cast<int>(_nnuisvar) && _nref[_nnuisvar]<_nref[ijxhiprev] && _nref[_nnuisvar]>=_nref[ijxhi] ) {
    ijxhiprev = _nnuisvar;
  }

  // if lo==hi, pass the next hi in line
  return ( ijxhi!=ijxlo ? ijxhi : ijxhiprev );
}


//_____________________________________________________________________________
Int_t RooStarMomentMorph::ijlo(const Int_t& i, const double& nval) const
{

  // find lower index on left side of nval
  // if lo==hi, pass the next-to lower index in line

  // i = nuisance parameter
  if (i<0 || i>=static_cast<int>(_nnuis.size())) return -1;

  // initialize
  int ijxhi = ij(i,_nnuis[i]-1);
  
  double nhi=_nref[ijxhi];

  if ( _nref[_nnuisvar]>nhi ) { ijxhi=_nnuisvar; }

  int ijxlo = ij(i,0);

  double nlo= _nref[ijxlo];

  if ( _nref[_nnuisvar]<nlo ) { 
    ijxlo=_nnuisvar; 
    nlo = _nref[_nnuisvar];
  }

  int ijxloprev = _nnuis[i]-1;

  for (Int_t j=_nnuis[i]-1; j>=0; --j) {
    int k = ij(i,j);
    if ( _nref[k]>nlo && _nref[k]<=nval ) { 
      nlo=_nref[k]; ijxlo=k; 
    } else {
      ijxloprev=k;
      break;
    }
  }
  if ( _nref[_nnuisvar]>nlo && _nref[_nnuisvar]<=nval ) { 
    ijxloprev = ijxlo;
    ijxlo=_nnuisvar; 
  }
  if ( ijxlo!=static_cast<int>(_nnuisvar) && _nref[_nnuisvar]>_nref[ijxloprev] && _nref[_nnuisvar]<=_nref[ijxlo] ) {
    ijxloprev = _nnuisvar;
  }

  return ( ijxlo!=ijxhi ? ijxlo : ijxloprev  );
}

//_____________________________________________________________________________
Bool_t RooStarMomentMorph::setBinIntegrator(RooArgSet& allVars) 
{

  if(allVars.getSize()==1){
    RooAbsReal* temp = const_cast<RooStarMomentMorph*>(this);
    temp->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator")  ;
    int nbins = ((RooRealVar*) allVars.first())->numBins();
    temp->specialIntegratorConfig(kTRUE)->getConfigSection("RooBinIntegrator").setRealValue("numBins",nbins);
    return true;
  }else{
    std::cout << "Currently BinIntegrator only knows how to deal with 1-d "<<std::endl;
    return false;
  }
  return false;
}
