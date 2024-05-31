#include "../interface/VerticalInterpPdf.h"
#include "../interface/RooCheapProduct.h"

#include "RooFit.h"
#include "Riostream.h"

#include "TList.h"
#include "RooRealProxy.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooAddGenContext.h"
#include "RooRealConstant.h"
#include "RooRealIntegral.h"
#include "RooMsgService.h"
#include "RooProdPdf.h"



ClassImp(VerticalInterpPdf)


//_____________________________________________________________________________
VerticalInterpPdf::VerticalInterpPdf() 
{
  // Default constructor
  // coverity[UNINIT_CTOR]
  _quadraticRegion = 0;
  _pdfFloorVal = 1e-15;
  _integralFloorVal = 1e-10;
}


//_____________________________________________________________________________
VerticalInterpPdf::VerticalInterpPdf(const char *name, const char *title, const RooArgList& inFuncList, const RooArgList& inCoefList, Double_t quadraticRegion, Int_t quadraticAlgo) :
  RooAbsPdf(name,title),
  _normIntMgr(this,10),
  _funcList("!funcList","List of functions",this),
  _coefList("!coefList","List of coefficients",this),
  _quadraticRegion(quadraticRegion),
  _quadraticAlgo(quadraticAlgo),
  _pdfFloorVal(1e-15),
  _integralFloorVal(1e-10)
{ 

  if (inFuncList.getSize()!=2*inCoefList.getSize()+1) {
    coutE(InputArguments) << "VerticalInterpPdf::VerticalInterpPdf(" << GetName() 
			  << ") number of pdfs and coefficients inconsistent, must have Nfunc=1+2*Ncoef" << std::endl ;
    assert(0);
  }

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
VerticalInterpPdf::~VerticalInterpPdf() = default;

//_____________________________________________________________________________
Double_t VerticalInterpPdf::evaluate() const 
{
  // Do running sum of coef/func pairs, calculate lastCoef.
  RooAbsReal* func = &(RooAbsReal&)_funcList[0];

  // Calculate the current value
  Double_t central = func->getVal();
  Double_t value = central;

  if (_quadraticAlgo >= 0) {
      // additive interpolation
      for (int iCoef = 0; iCoef < _coefList.getSize(); ++iCoef) {
          Double_t coefVal = static_cast<RooAbsReal&>(_coefList[iCoef]).getVal() ;
          RooAbsReal* funcUp = &(RooAbsReal&)_funcList[2 * iCoef + 1];
          RooAbsReal* funcDn = &(RooAbsReal&)_funcList[2 * iCoef + 2];
          value += interpolate(coefVal, central, funcUp, funcDn);
      }
  } else {
      // multiplicative interpolation
      for (int iCoef = 0; iCoef < _coefList.getSize(); ++iCoef) {
          Double_t coefVal = static_cast<RooAbsReal&>(_coefList[iCoef]).getVal() ;
          RooAbsReal* funcUp = &(RooAbsReal&)_funcList[2 * iCoef + 1];
          RooAbsReal* funcDn = &(RooAbsReal&)_funcList[2 * iCoef + 2];
          value *= interpolate(coefVal, central, funcUp, funcDn);
      }
  }
   
  return ( value > 0. ? value : _pdfFloorVal);
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
  RooArgSet* normSet = normSet2 ? getObservables(normSet2) : 0 ;


  // Check if this configuration was created before
  Int_t sterileIdx(-1) ;
  CacheElem* cache = (CacheElem*) _normIntMgr.getObj(normSet,&analVars,&sterileIdx,(const char *)0);
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
  Int_t code = _normIntMgr.setObj(normSet,&analVars,(RooAbsCacheElement*)cache,0) ;

  if (normSet) {
    delete normSet ;
  }

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
Double_t VerticalInterpPdf::analyticalIntegralWN(Int_t code, const RooArgSet* normSet2, const char* /*rangeName*/) const 
{
  // Implement analytical integrations by deferring integration of component
  // functions to integrators of components

  // Handle trivial passthrough scenario
  if (code==0) return getVal(normSet2) ;

  Double_t value = 0;

  // WVE needs adaptation for rangeName feature
  CacheElem* cache = (CacheElem*) _normIntMgr.getObjByIndex(code-1) ;
  RooArgList& fIntL = cache->_funcIntList;

  Double_t central = static_cast<RooAbsReal&>(fIntL[0]).getVal();
  value += central;

  for (int iCoef = 0; iCoef < _coefList.getSize(); ++iCoef) {
    Double_t coefVal = static_cast<RooAbsReal&>(_coefList[iCoef]).getVal(normSet2) ;
    RooAbsReal * funcIntUp = &(RooAbsReal&)fIntL[2 * iCoef + 1];
    RooAbsReal * funcIntDn = &(RooAbsReal&)fIntL[2 * iCoef + 2];
    value += interpolate(coefVal, central, funcIntUp, funcIntDn);
  }
  
  Double_t normVal(1) ;
  if (normSet2) {
    RooArgList& fnl = cache->_funcNormList;
    central = static_cast<RooAbsReal&>(fnl[0]).getVal(normSet2) ;
    normVal = central;

    for (int iCoef = 0; iCoef < _coefList.getSize(); ++iCoef) {
      RooAbsReal *funcNormUp = &(RooAbsReal&)fnl[2 * iCoef + 1];
      RooAbsReal *funcNormDn = &(RooAbsReal&)fnl[2 * iCoef + 2];
      Double_t coefVal = static_cast<RooAbsReal&>(_coefList[iCoef]).getVal(normSet2) ;
      normVal += interpolate(coefVal, central, funcNormUp, funcNormDn);
    }
  }

  Double_t result = 0;
  if(normVal>0.) result = value / normVal;
  return (result > 0. ? result : _integralFloorVal);
}

Double_t VerticalInterpPdf::interpolate(Double_t coeff, Double_t central, RooAbsReal *fUp, RooAbsReal *fDn) const  
{
    if (_quadraticAlgo == -1) {
        Double_t kappa = (coeff > 0 ? fUp->getVal()/central : central/fDn->getVal());
        return pow(kappa,fabs(coeff));
    }

    if (fabs(coeff) >= _quadraticRegion) {
        return coeff * (coeff > 0 ? fUp->getVal() - central : central - fDn->getVal());
    } else {
        // quadratic interpolation coefficients between the three
        if (_quadraticAlgo == 0) {
            // quadratic interpolation null at zero and continuous at boundaries, but not differentiable at boundaries
            // conditions:
            //   c_up (+_quadraticRegion) = +_quadraticRegion
            //   c_cen(+_quadraticRegion) = -_quadraticRegion
            //   c_dn (+_quadraticRegion) = 0
            //   c_up (-_quadraticRegion) = 0 
            //   c_cen(-_quadraticRegion) = -_quadraticRegion
            //   c_dn (-_quadraticRegion) = +_quadraticRegion
            //   c_up(0) = c_dn(0) = c_cen(0) = 0
            Double_t c_up  = + coeff * (_quadraticRegion + coeff) / (2 * _quadraticRegion);
            Double_t c_dn  = - coeff * (_quadraticRegion - coeff) / (2 * _quadraticRegion);
            Double_t c_cen = - coeff * coeff / _quadraticRegion;
            return c_up * fUp->getVal() + c_dn * fDn->getVal() + c_cen * central;
        } else if (_quadraticAlgo == 1) { 
            // quadratic interpolation that is everywhere differentiable, but it's not null at zero
            // conditions on the function
            //   c_up (+_quadraticRegion) = +_quadraticRegion
            //   c_cen(+_quadraticRegion) = -_quadraticRegion
            //   c_dn (+_quadraticRegion) = 0
            //   c_up (-_quadraticRegion) = 0 
            //   c_cen(-_quadraticRegion) = -_quadraticRegion
            //   c_dn (-_quadraticRegion) = +_quadraticRegion
            // conditions on the derivatives
            //   c_up '(+_quadraticRegion) = +1
            //   c_cen'(+_quadraticRegion) = -1
            //   c_dn '(+_quadraticRegion) = 0
            //   c_up '(-_quadraticRegion) = 0
            //   c_cen'(-_quadraticRegion) = +1
            //   c_dn '(-_quadraticRegion) = -1
            Double_t c_up  = (_quadraticRegion + coeff) * (_quadraticRegion + coeff) / (4 * _quadraticRegion);
            Double_t c_dn  = (_quadraticRegion - coeff) * (_quadraticRegion - coeff) / (4 * _quadraticRegion);
            Double_t c_cen = - c_up - c_dn;
            return c_up * fUp->getVal() + c_dn * fDn->getVal() + c_cen * central;
        } else/* if (_quadraticAlgo == 1)*/ {
            // P(6) interpolation that is everywhere differentiable and null at zero
            /* === how the algorithm works, in theory ===
            * let  dhi = h_hi - h_nominal
            *      dlo = h_lo - h_nominal
            * and x be the morphing parameter
            * we define alpha = x * 0.5 * ((dhi-dlo) + (dhi+dlo)*smoothStepFunc(x));
            * which satisfies:
            *     alpha(0) = 0
            *     alpha(+1) = dhi
            *     alpha(-1) = dlo
            *     alpha(x >= +1) = |x|*dhi
            *     alpha(x <= -1) = |x|*dlo
            *     alpha is continuous and has continuous first and second derivative, as smoothStepFunc has them
            * === and in practice ===
            * we already have computed the histogram for diff=(dhi-dlo) and sum=(dhi+dlo)
            * so we just do template += (0.5 * x) * (diff + smoothStepFunc(x) * sum)
            * ========================================== */
            Double_t cnorm = coeff/_quadraticRegion;
            Double_t cnorm2 = pow(cnorm, 2);
            Double_t hi = fUp->getVal() - central;
            Double_t lo = fDn->getVal() - central;
            Double_t sum = hi+lo;
            Double_t diff = hi-lo;
            Double_t a = coeff/2.; // cnorm*_quadraticRegion
            Double_t b = 0.125 * cnorm * (cnorm2 * (3.*cnorm2 - 10.) + 15.);
            Double_t result = a*(diff + b*sum);
            return result;
        }
    }
}


//_____________________________________________________________________________
void VerticalInterpPdf::setFloorVals(Double_t const& pdf_val, Double_t const& integral_val){ _pdfFloorVal = pdf_val; _integralFloorVal = integral_val; }
