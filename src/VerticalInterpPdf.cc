#include "../interface/VerticalInterpPdf.h"
#include "../interface/RooCheapProduct.h"
#include "../interface/CombineMathFuncs.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,34,06)
#include "RooFit/Detail/RooNormalizedPdf.h"
#endif

#include "Riostream.h"

#include "RooRealVar.h"
#include "RooRealIntegral.h"
#include "RooMsgService.h"
#include "RooProdPdf.h"

#include <functional>

ClassImp(VerticalInterpPdf)


//_____________________________________________________________________________
VerticalInterpPdf::VerticalInterpPdf(const char *name, const char *title, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t quadraticRegion, Int_t quadraticAlgo) :
  RooAbsPdf(name,title),
  _normIntMgr(this,10),
  _funcList("!funcList","List of functions",this),
  _coefList("!coefList","List of coefficients",this),
  _quadraticRegion(quadraticRegion),
  _quadraticAlgo(quadraticAlgo)
{ 

  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "VerticalInterpPdf::VerticalInterpPdf(" << GetName() 
			  << ") number of pdfs and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl ;
    assert(0);
  }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
  // With ROOT 6.32, RooFit can do the type checking for us
  _funcList.addTyped<RooAbsReal>(inFuncList);
  _coefList.addTyped<RooAbsReal>(inCoefList);
#else
  for (RooAbsArg *func : inFuncList) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: VerticalInterpPdf::VerticalInterpPdf(" << GetName() << ") function  " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _funcList.add(*func) ;
  }

  for (RooAbsArg *coef : inCoefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      coutE(InputArguments) << "ERROR: VerticalInterpPdf::VerticalInterpPdf(" << GetName() << ") coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*coef) ;    
  }
#endif

  if (_quadraticAlgo == -1) { 
    // multiplicative morphing: no way to do analytical integrals.
    _forceNumInt = kTRUE; 
  } else if (_quadraticAlgo >= 100) {
      _quadraticAlgo -= 100;
      _forceNumInt = kTRUE; 
  }
}




//_____________________________________________________________________________
VerticalInterpPdf::VerticalInterpPdf(const VerticalInterpPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _normIntMgr(other._normIntMgr,this),
  _funcList("!funcList",this,other._funcList),
  _coefList("!coefList",this,other._coefList),
  _quadraticRegion(other._quadraticRegion),
  _quadraticAlgo(other._quadraticAlgo),
  _pdfFloorVal(other._pdfFloorVal),
  _integralFloorVal(other._integralFloorVal)
{
  // Copy constructor
}


//_____________________________________________________________________________
Double_t VerticalInterpPdf::evaluate() const 
{
  // Do running sum of coef/func pairs, calculate lastCoef.
  Double_t value = (_quadraticAlgo >= 0) ?
    RooFit::Detail::MathFuncs::opInterpolate<std::plus<Double_t>>(_coefList, _funcList, _pdfFloorVal, _quadraticRegion, _quadraticAlgo) :
    RooFit::Detail::MathFuncs::opInterpolate<std::multiplies<Double_t>>(_coefList, _funcList, _pdfFloorVal, _quadraticRegion, _quadraticAlgo);
  return value;
}




//_____________________________________________________________________________
Bool_t VerticalInterpPdf::checkObservables(const RooArgSet* nset) const 
{
  // Check if FUNC is valid for given normalization set.
  // Coeffient and FUNC must be non-overlapping, but func-coefficient 
  // pairs may overlap each other
  //
  // In the present implementation, coefficients may not be observables or derive
  // from observables
  if (_quadraticAlgo == -1) return false; // multiplicative morphing. we don't care.
  

  for (RooAbsArg *coef : _coefList) {
    if (coef->dependsOn(*nset)) {
      coutE(InputArguments) << "RooRealPdf::checkObservables(" << GetName() << "): ERROR coefficient " << coef->GetName() 
			    << " depends on one or more of the following observables" ; nset->Print("1") ;
      return true;
    }
  }

  RooAbsReal * coef = nullptr;
  for (int ifunc = 0; ifunc < _funcList.getSize(); ++ifunc) {
    RooAbsReal* func = &(RooAbsReal&)_funcList[ifunc];
    if (ifunc % 2 == 0) coef = &(RooAbsReal&)_coefList[ifunc];
    if (coef && func->observableOverlaps(nset,*coef)) {
      coutE(InputArguments) << "VerticalInterpPdf::checkObservables(" << GetName() << "): ERROR: coefficient " << coef->GetName() 
			    << " and FUNC " << func->GetName() << " have one or more observables in common" << std::endl;
      return true;
    }
    ifunc++;
  }
  
  return false;
}




//_____________________________________________________________________________
Int_t VerticalInterpPdf::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, 
					     const RooArgSet* normSet2, const char* /*rangeName*/) const 
{
  // Advertise that all integrals can be handled internally.

  // Handle trivial no-integration scenario
  if (allVars.getSize()==0) return 0 ;
  if (_forceNumInt) return 0 ;

  // Select subset of allVars that are actual dependents
  analVars.add(allVars) ;
  std::unique_ptr<RooArgSet> normSet{normSet2 ? getObservables(normSet2) : nullptr};


  // Check if this configuration was created before
  Int_t sterileIdx(-1) ;
  CacheElem* cache = (CacheElem*) _normIntMgr.getObj(normSet.get(),&analVars,&sterileIdx,(const char *)0);
  if (cache) {
    return _normIntMgr.lastIndex()+1 ;
  }
  
  // Create new cache element
  cache = new CacheElem ;

  // Make list of function projection and normalization integrals 
  for (RooAbsArg *funcAbsArg : _funcList) {
    auto func = static_cast<RooAbsReal*>(funcAbsArg);
    RooAbsReal* funcInt = nullptr;
    if (isConditionalProdPdf(func)) {
      RooProdPdf *prod = static_cast<RooProdPdf*>(func);
      funcInt = makeConditionalProdPdfIntegral(prod, analVars);
    } else {
      funcInt = func->createIntegral(analVars) ;
    }
    cache->_funcIntList.addOwned(*funcInt) ;
    if (normSet && normSet->getSize()>0) {
      RooAbsReal* funcNorm = func->createIntegral(*normSet) ;
      cache->_funcNormList.addOwned(*funcNorm) ;
    }
  }

  // Store cache element
  Int_t code = _normIntMgr.setObj(normSet.get(),&analVars,(RooAbsCacheElement*)cache,0) ;

  return code+1 ; 
}

bool VerticalInterpPdf::isConditionalProdPdf(RooAbsReal *pdf) const {
  // If pdf is not RooProdPdf, we can return false immediately
  if (!dynamic_cast<RooProdPdf*>(pdf)) {
    return false;
  }
  RooProdPdf *prod = static_cast<RooProdPdf*>(pdf);

  // Now loop through prodPdf components, and find the saved "nset" for each one
  // If this is a non-empty set, we have a conditional pdf
  for (int i = 0; i < prod->pdfList().getSize(); ++i) {
    RooAbsPdf *compPdf = static_cast<RooAbsPdf*>(prod->pdfList().at(i));
    RooArgSet* nset = prod->findPdfNSet(*compPdf);
    if (nset && TString(nset->GetName()) == "nset" && nset->getSize() > 0) {
      // We can immediately return true
      return true;
    }
  }
  return false;
}

RooAbsReal* VerticalInterpPdf::makeConditionalProdPdfIntegral(RooAbsPdf* pdf, RooArgSet const& analVars) const {
  // If pdf is not RooProdPdf, we can return false immediately
  if (!dynamic_cast<RooProdPdf*>(pdf)) {
    return nullptr;
  }
  RooProdPdf *prod = static_cast<RooProdPdf*>(pdf);

  // Make a list of component integrals
  RooArgList prodIntComps;
  for (int i = 0; i < prod->pdfList().getSize(); ++i) {
    RooAbsPdf *compPdf = static_cast<RooAbsPdf*>(prod->pdfList().at(i));
    RooArgSet* nset = prod->findPdfNSet(*compPdf);
    if (nset && TString(nset->GetName()) == "nset" && nset->getSize() > 0) {
      // We integrate the subset of analVars variables that are in nset
      std::unique_ptr<RooArgSet> selObs(static_cast<RooArgSet*>(nset->selectCommon(analVars)));
      prodIntComps.add(*compPdf->createIntegral(*selObs));
      // std::cout << "For ProdPdf=" << prod->GetName() << ", added integral of " << compPdf->GetName() << " for cond. observables: \n";
      // selObs->Print();
    } else {
      std::unique_ptr<RooArgSet> iVars(compPdf->getVariables());
      std::unique_ptr<RooArgSet> selObs(static_cast<RooArgSet*>(iVars->selectCommon(analVars)));
      prodIntComps.add(*compPdf->createIntegral(*selObs));
      // std::cout << "For ProdPdf=" << prod->GetName() << ", added integral of " << compPdf->GetName() << " for observables: \n";
      // selObs->Print();
    }
  }
  RooCheapProduct *intProd = new RooCheapProduct(TString("intProd_")+prod->GetName(), "", prodIntComps);
  intProd->addOwnedComponents(prodIntComps);
  return intProd;
}

//_____________________________________________________________________________
const RooArgList& VerticalInterpPdf::funcIntListFromCache() const
{
  // Return the list of integrals of the component functions
  auto* cache = dynamic_cast<CacheElem*>(_normIntMgr.getObjByIndex(0));
  return cache->_funcIntList;
}


//_____________________________________________________________________________
Double_t VerticalInterpPdf::analyticalIntegralWN(Int_t code, const RooArgSet* normSet2, const char* /*rangeName*/) const 
{
  // Implement analytical integrations by deferring integration of component
  // functions to integrators of components

  // Handle trivial passthrough scenario
  if (code==0) return getVal(normSet2) ;

  // WVE needs adaptation for rangeName feature
  const RooArgList& fIntL = funcIntListFromCache();
  Double_t value = RooFit::Detail::MathFuncs::opInterpolate<std::plus<Double_t>>(_coefList, fIntL, _pdfFloorVal, _quadraticRegion, _quadraticAlgo, normSet2);

  Double_t normVal(1) ;
  RooArgList& fNormL = (dynamic_cast<CacheElem*>(_normIntMgr.getObjByIndex(0)))->_funcNormList;
  if (normSet2) {
    normVal = RooFit::Detail::MathFuncs::opInterpolate<std::plus<Double_t>>(_coefList, fNormL, _pdfFloorVal, _quadraticRegion, _quadraticAlgo, normSet2);
  }

  Double_t result = 0;
  if(normVal>0.) result = value / normVal;
  return (result > 0. ? result : _integralFloorVal);
}

Double_t VerticalInterpPdf::interpolate(Double_t coeff, Double_t central, RooAbsReal *fUp, RooAbsReal *fDn) const  
{
  return RooFit::Detail::MathFuncs::interpolate(coeff, central, fUp->getVal(), fDn->getVal(), _quadraticRegion, _quadraticAlgo);
}


//_____________________________________________________________________________
void VerticalInterpPdf::setFloorVals(Double_t const& pdf_val, Double_t const& integral_val){ _pdfFloorVal = pdf_val; _integralFloorVal = integral_val; }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,34,06)

// In the new RooFit CPU backend, the computation graph is "compiled" for a fixed normalization set.
// The function RooAbsArg::compileForNormSet() can be overridden to hook into this process.
//
// The VerticalInterpPdf is unfortunately a special case, just like the
// RooRealSumPdf that served as its inspiration. It deliberately breaks the
// auto-normalization normalization of its input functions in _funcList by
// explicitly calling getVal() on them without any normalization set.
//
// Therefore, we have to override the compilation for a given normSet such that
// it doesn't propagate the normalization set to the servers.
//
// One can surely refactor the VerticalInterpPdf such that the input functions
// are evaluated with a defined normalization set, but this would take more
// effort to validate. The new RooFit CPU backend is not even used by Combine
// yet, so overriding compileForNormSet() doesn't require any validation.
std::unique_ptr<RooAbsArg>
VerticalInterpPdf::compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext &ctx) const
{
   if (normSet.empty() || selfNormalized()) {
      return RooAbsPdf::compileForNormSet({}, ctx);
   }
   std::unique_ptr<RooAbsPdf> pdfClone(static_cast<RooAbsPdf *>(this->Clone()));

   ctx.compileServers(*pdfClone, {});

   RooArgSet depList;
   pdfClone->getObservables(&normSet, depList);

   auto newArg = std::make_unique<RooFit::Detail::RooNormalizedPdf>(*pdfClone, depList);

   // The direct servers are this pdf and the normalization integral, which
   // don't need to be compiled further.
   for (RooAbsArg *server : newArg->servers()) {
      ctx.markAsCompiled(*server);
   }
   ctx.markAsCompiled(*newArg);
   newArg->addOwnedComponents(std::move(pdfClone));
   return newArg;
}
#endif
