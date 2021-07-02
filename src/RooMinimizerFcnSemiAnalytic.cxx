/*****************************************************************************
 * Authors:                                                                  *
 *  Andrea Carlo Marini,   MIT (MA),        andrea.carlo.marini@cern.ch        *
 *  Code based on roofit code.                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

//////////////////////////////////////////////////////////////////////////////
//
// RooMinimizerFcnSemiAnalytic is an interface class to the ROOT::Math function 
// for minization.
//                                                                                   

#include <iostream>

#include "RooFit.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMinimizerFcnSemiAnalytic.h"

#include "Riostream.h"

#include "TIterator.h"
#include "TClass.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"

#include "RooMinimizer.h"
#include <stdexcept>

using namespace std;

RooMinimizerFcnSemiAnalytic::RooMinimizerFcnSemiAnalytic(RooAbsReal *funct, RooMinimizerSemiAnalytic* context, const std::map<std::string,RooAbsReal*>& knownDerivatives, bool verbose) :
  _funct(funct), 
  _knownDerivatives (knownDerivatives),
  _context(context),
  // Reset the *largest* negative log-likelihood value we have seen so far
  _maxFCN(-1e30), _numBadNLL(0),  
  _printEvalErrors(10), _doEvalErrorWall(kTRUE),
  _nDim(0), _logfile(0),
  _verbose(verbose)
{ 

  _evalCounter = 0 ;
  
  // Examine parameter list
  RooArgSet* paramSet = _funct->getParameters(RooArgSet());
  RooArgList paramList(*paramSet);
  delete paramSet;

  _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant",kFALSE);
  if (_floatParamList->getSize()>1) {
    _floatParamList->sort();
  }
  _floatParamList->setName("floatParamList");

  _constParamList = (RooArgList*) paramList.selectByAttrib("Constant",kTRUE);
  if (_constParamList->getSize()>1) {
    _constParamList->sort();
  }
  _constParamList->setName("constParamList");

  // Remove all non-RooRealVar parameters from list (MINUIT cannot handle them)
  TIterator* pIter = _floatParamList->createIterator();
  RooAbsArg* arg;
  while ((arg=(RooAbsArg*)pIter->Next())) {
    if (!arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
      std::cout
      /*oocoutW(_context,Minimization)*/ << "RooMinimizerFcnSemiAnalytic::RooMinimizerFcnSemiAnalytic: removing parameter " 
				     << arg->GetName()
                                     << " from list because it is not of type RooRealVar" << endl;
      _floatParamList->remove(*arg);
    }
  }
  delete pIter;

  _nDim = _floatParamList->getSize();

  _verbose=true; // print the variables for which we have derivatives?
  updateFloatVec() ;
  _verbose=verbose;
  
  // Save snapshot of initial lists
  _initFloatParamList = (RooArgList*) _floatParamList->snapshot(kFALSE) ;
  _initConstParamList = (RooArgList*) _constParamList->snapshot(kFALSE) ;

}



RooMinimizerFcnSemiAnalytic::RooMinimizerFcnSemiAnalytic(const RooMinimizerFcnSemiAnalytic& other) : ROOT::Math::IMultiGradFunction(other), 
  _evalCounter(other._evalCounter),
  _funct(other._funct),
  _knownDerivatives(other._knownDerivatives),
  _context(other._context),
  _maxFCN(other._maxFCN),
  _numBadNLL(other._numBadNLL),
  _printEvalErrors(other._printEvalErrors),
  _doEvalErrorWall(other._doEvalErrorWall),
  _nDim(other._nDim),
  _logfile(other._logfile),
  _verbose(other._verbose),
  _floatParamVec(other._floatParamVec),
  _derivParamVec(other._derivParamVec)
{  
  _floatParamList = new RooArgList(*other._floatParamList) ;
  _constParamList = new RooArgList(*other._constParamList) ;
  _initFloatParamList = (RooArgList*) other._initFloatParamList->snapshot(kFALSE) ;
  _initConstParamList = (RooArgList*) other._initConstParamList->snapshot(kFALSE) ;  
}


RooMinimizerFcnSemiAnalytic::~RooMinimizerFcnSemiAnalytic()
{
  delete _floatParamList;
  delete _initFloatParamList;
  delete _constParamList;
  delete _initConstParamList;
}


ROOT::Math::IMultiGradFunction* RooMinimizerFcnSemiAnalytic::Clone() const 
{  
  return new RooMinimizerFcnSemiAnalytic(*this) ;
}


Bool_t RooMinimizerFcnSemiAnalytic::Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, 
				 Bool_t optConst, Bool_t verbose)
{

  // Internal function to synchronize TMinimizer with current
  // information in RooAbsReal function parameters
  
  Bool_t constValChange(kFALSE) ;
  Bool_t constStatChange(kFALSE) ;
  
  Int_t index(0) ;
  
  // Handle eventual migrations from constParamList -> floatParamList
  for(index= 0; index < _constParamList->getSize() ; index++) {

    RooRealVar *par= dynamic_cast<RooRealVar*>(_constParamList->at(index)) ;
    if (!par) continue ;

    RooRealVar *oldpar= dynamic_cast<RooRealVar*>(_initConstParamList->at(index)) ;
    if (!oldpar) continue ;

    // Test if constness changed
    if (!par->isConstant()) {      
    
      // Remove from constList, add to floatList
      _constParamList->remove(*par) ;
      _floatParamList->add(*par) ;
      _initFloatParamList->addClone(*oldpar) ;      
      _initConstParamList->remove(*oldpar) ;
      constStatChange=kTRUE ;
      _nDim++ ;

      if (verbose) {
	/*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: parameter " 
				     << par->GetName() << " is now floating." << endl ;
      }
    } 

    // Test if value changed
    if (par->getVal()!= oldpar->getVal()) {
      constValChange=kTRUE ;      
      if (verbose) {
	/*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: value of constant parameter " 
				       << par->GetName() 
				       << " changed from " << oldpar->getVal() << " to " 
				       << par->getVal() << endl ;
      }
    }

  }

  // Update reference list
  *_initConstParamList = *_constParamList ;
  
  // Synchronize MINUIT with function state
  // Handle floatParamList
  for(index= 0; index < _floatParamList->getSize(); index++) {
    RooRealVar *par= dynamic_cast<RooRealVar*>(_floatParamList->at(index)) ;
    
    if (!par) continue ;

    Double_t pstep(0) ;
    Double_t pmin(0) ;
    Double_t pmax(0) ;

    if(!par->isConstant()) {

      // Verify that floating parameter is indeed of type RooRealVar 
      if (!par->IsA()->InheritsFrom(RooRealVar::Class())) {
          std::cout
	      /*oocoutW(_context,Minimization)*/ << "RooMinimizerFcnSemiAnalytic::fit: Error, non-constant parameter " 
				       << par->GetName() 
				       << " is not of type RooRealVar, skipping" << endl ;
	_floatParamList->remove(*par);
	index--;
	_nDim--;
	continue ;
      }

      // Set the limits, if not infinite
      if (par->hasMin() )
	pmin = par->getMin();
      if (par->hasMax() )
	pmax = par->getMax();
      
      // Calculate step size
      pstep = par->getError();
      if(pstep <= 0) {
	// Floating parameter without error estitimate
	if (par->hasMin() && par->hasMax()) {
	  pstep= 0.1*(pmax-pmin);

	  // Trim default choice of error if within 2 sigma of limit
	  if (pmax - par->getVal() < 2*pstep) {
	    pstep = (pmax - par->getVal())/2 ;
	  } else if (par->getVal() - pmin < 2*pstep) {
	    pstep = (par->getVal() - pmin )/2 ;	    
	  }	  

	  // If trimming results in zero error, restore default
	  if (pstep==0) {
	    pstep= 0.1*(pmax-pmin);
	  }

	} else {
	  pstep=1 ;
	}						  
	if(verbose) {
        std::cout
            /*oocoutW(_context,Minimization)*/ << "RooMinimizerFcnSemiAnalytic::synchronize: WARNING: no initial error estimate available for "
					 << par->GetName() << ": using " << pstep << endl;
	}
      }       
    } else {
      pmin = par->getVal() ;
      pmax = par->getVal() ;      
    }

    // new parameter
    if (index>=Int_t(parameters.size())) {

      if (par->hasMin() && par->hasMax()) {
	parameters.push_back(ROOT::Fit::ParameterSettings(par->GetName(),
							  par->getVal(),
							  pstep,
							  pmin,pmax));
      }
      else {
	parameters.push_back(ROOT::Fit::ParameterSettings(par->GetName(),
							  par->getVal(),
							  pstep));
        if (par->hasMin() ) 
           parameters.back().SetLowerLimit(pmin);
        else if (par->hasMax() ) 
           parameters.back().SetUpperLimit(pmax);
      }
      
      continue;

    }

    Bool_t oldFixed = parameters[index].IsFixed();
    Double_t oldVar = parameters[index].Value();
    Double_t oldVerr = parameters[index].StepSize();
    Double_t oldVlo = parameters[index].LowerLimit();
    Double_t oldVhi = parameters[index].UpperLimit();

    if (par->isConstant() && !oldFixed) {

      // Parameter changes floating -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	parameters[index].SetValue(par->getVal());
	if (verbose) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: value of parameter " 
					 << par->GetName() << " changed from " << oldVar 
					 << " to " << par->getVal() << endl ;
	}
      }
      parameters[index].Fix();
      constStatChange=kTRUE ;
      if (verbose) {
	/*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: parameter " 
				       << par->GetName() << " is now fixed." << endl ;
      }

    } else if (par->isConstant() && oldFixed) {
      
      // Parameter changes constant -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	parameters[index].SetValue(par->getVal());
	constValChange=kTRUE ;

	if (verbose) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: value of fixed parameter " 
					 << par->GetName() << " changed from " << oldVar 
					 << " to " << par->getVal() << endl ;
	}
      }

    } else {
      // Parameter changes constant -> floating
      if (!par->isConstant() && oldFixed) {
	parameters[index].Release();
	constStatChange=kTRUE ;
	
	if (verbose) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: parameter " 
					 << par->GetName() << " is now floating." << endl ;
	}
      } 

      // Parameter changes constant -> floating : update all if necessary      
      if (oldVar!=par->getVal() || oldVlo!=pmin || oldVhi != pmax || oldVerr!=pstep) {
	parameters[index].SetValue(par->getVal()); 
	parameters[index].SetStepSize(pstep);
        if (par->hasMin() && par->hasMax() ) 
           parameters[index].SetLimits(pmin,pmax);  
        else if (par->hasMin() )
           parameters[index].SetLowerLimit(pmin);  
        else if (par->hasMax() )
           parameters[index].SetUpperLimit(pmax);  
      }

      // Inform user about changes in verbose mode
      if (verbose) {
	// if ierr<0, par was moved from the const list and a message was already printed

	if (oldVar!=par->getVal()) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: value of parameter " 
					 << par->GetName() << " changed from " << oldVar << " to " 
					 << par->getVal() << endl ;
	}
	if (oldVlo!=pmin || oldVhi!=pmax) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: limits of parameter " 
					 << par->GetName() << " changed from [" << oldVlo << "," << oldVhi 
					 << "] to [" << pmin << "," << pmax << "]" << endl ;
	}

	// If oldVerr=0, then parameter was previously fixed
	if (oldVerr!=pstep && oldVerr!=0) {
	  /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: error/step size of parameter " 
					 << par->GetName() << " changed from " << oldVerr << " to " << pstep << endl ;
	}
      }      
    }
  }

  if (optConst) {
    if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: set of constant parameters changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange) ;
    } else if (constValChange) {
      /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::synchronize: constant parameter values changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange) ;
    }
    
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;  

  }

  updateFloatVec() ;

  return 0 ;  

}

Double_t RooMinimizerFcnSemiAnalytic::GetPdfParamVal(Int_t index)
{
  // Access PDF parameter value by ordinal index (needed by MINUIT)

  return ((RooRealVar*)_floatParamList->at(index))->getVal() ;
}

Double_t RooMinimizerFcnSemiAnalytic::GetPdfParamErr(Int_t index)
{
  // Access PDF parameter error by ordinal index (needed by MINUIT)
  return ((RooRealVar*)_floatParamList->at(index))->getError() ;
}


void RooMinimizerFcnSemiAnalytic::SetPdfParamErr(Int_t index, Double_t value)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->setError(value) ;
}



void RooMinimizerFcnSemiAnalytic::ClearPdfParamAsymErr(Int_t index)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->removeAsymError() ;
}


void RooMinimizerFcnSemiAnalytic::SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->setAsymError(loVal,hiVal) ;
}


void RooMinimizerFcnSemiAnalytic::BackProp(const ROOT::Fit::FitResult &results)
{
  // Transfer MINUIT fit results back into RooFit objects

  for (Int_t index= 0; index < _nDim; index++) {
    Double_t value = results.Value(index);
    SetPdfParamVal(index, value);

    // Set the parabolic error    
    Double_t err = results.Error(index);
    SetPdfParamErr(index, err);
    
    Double_t eminus = results.LowerError(index);
    Double_t eplus = results.UpperError(index);

    if(eplus > 0 || eminus < 0) {
      // Store the asymmetric error, if it is available
      SetPdfParamErr(index, eminus,eplus);
    } else {
      // Clear the asymmetric error
      ClearPdfParamAsymErr(index) ;
    }
  }

}

Bool_t RooMinimizerFcnSemiAnalytic::SetLogFile(const char* inLogfile) 
{
  // Change the file name for logging of a RooMinimizer of all MINUIT steppings
  // through the parameter space. If inLogfile is null, the current log file
  // is closed and logging is stopped.

  if (_logfile) {
    /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::setLogFile: closing previous log file" << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile = 0 ;
  }
  _logfile = new ofstream(inLogfile) ;
  if (!_logfile->good()) {
    /*oocoutI(_context,Minimization)*/std::cout << "RooMinimizerFcnSemiAnalytic::setLogFile: cannot open file " << inLogfile << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile= 0;
  }  
  
  return kFALSE ;

}


void RooMinimizerFcnSemiAnalytic::ApplyCovarianceMatrix(TMatrixDSym& V) 
{
  // Apply results of given external covariance matrix. i.e. propagate its errors
  // to all RRV parameter representations and give this matrix instead of the
  // HESSE matrix at the next save() call

  for (Int_t i=0 ; i<_nDim ; i++) {
    // Skip fixed parameters
    if (_floatParamList->at(i)->isConstant()) {
      continue ;
    }
    SetPdfParamErr(i, sqrt(V(i,i))) ;		  
  }

}


Bool_t RooMinimizerFcnSemiAnalytic::SetPdfParamVal(const Int_t &index, const Double_t &value) const
{
  //RooRealVar* par = (RooRealVar*)_floatParamList->at(index);
  RooRealVar* par = (RooRealVar*)_floatParamVec[index] ;

  if (par->getVal()!=value) {
    if (_verbose) cout << par->GetName() << "=" << value << ", " ;
    
    par->setVal(value);
    return kTRUE;
  }

  return kFALSE;
}



////////////////////////////////////////////////////////////////////////////////

void RooMinimizerFcnSemiAnalytic::updateFloatVec() 
{ 
    // derivative vector and float vectors needs to be aligned
  _floatParamVec.clear() ;
  _derivParamVec.clear();

  RooFIter iter = _floatParamList->fwdIterator() ;
  RooAbsArg* arg ;
  _floatParamVec = std::vector<RooAbsArg*>(_floatParamList->getSize()) ;
  _derivParamVec = std::vector<RooAbsReal*>(_floatParamList->getSize()) ;
  Int_t i(0) ;
  while((arg=iter.next())) {
    auto derivative=_knownDerivatives.find(arg->GetName());
    if (derivative == _knownDerivatives.end()){
        if(_verbose) std::cout<<"[DEBUG]:["<<__PRETTY_FUNCTION__<<"]:"<< "derivative for param "<<arg->GetName()<<" not found. I will use numerical derivatives."<<std::endl;
        _derivParamVec[i]  = nullptr;
    }
    else{
        if(_verbose) std::cout<<"[DEBUG]:["<<__PRETTY_FUNCTION__<<"]:"<< "derivative for param "<<arg->GetName()<<" found with name "<<derivative->second->GetName()<<"."<<std::endl;
        _derivParamVec[i] = derivative->second;
    }
    _floatParamVec[i++] = arg ;
  }
}



double RooMinimizerFcnSemiAnalytic::DoEval(const double *x) const 
{
  //std::cout<<"[DEBUG]:["<<__PRETTY_FUNCTION__<<"]:"<< "DoEval"<<std::endl; // run with -v3 to get it

  // Set the parameter values for this iteration
  for (int index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters  
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = _funct->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (RooAbsReal::numEvalErrors()>0 || fvalue > 1e30) {

    if (_printEvalErrors>=0) {

      if (_doEvalErrorWall) {
          std::cout
        /*oocoutW(_context,Minimization)*/ << "RooMinimizerFcnSemiAnalytic: Minimized function has error status." << endl 
				       << "Returning maximum FCN so far (" << _maxFCN 
				       << ") to force MIGRAD to back out of this region. Error log follows" << endl ;
      } else {
          std::cout
        /*oocoutW(_context,Minimization)*/ << "RooMinimizerFcnSemiAnalytic: Minimized function has error status but is ignored" << endl ;
      } 

      Bool_t first(kTRUE) ;
      /*ooccoutW(_context,Minimization)*/std::cout << "Parameter values: " ;
      //for (const auto par : *_floatParamList) { //new loop}
      for (int ip =0 ;ip < _floatParamList->getSize();++ip){ 
        const auto par=_floatParamList->at(ip); 
        auto var = static_cast<const RooRealVar*>(par);
        if (first) { first = kFALSE ; } else /*ooccoutW(_context,Minimization)*/std::cout << ", " ;
        /*ooccoutW(_context,Minimization)*/std::cout << var->GetName() << "=" << var->getVal() ;
      }
      /*ooccoutW(_context,Minimization)*/std::cout << endl ;
      
      RooAbsReal::printEvalErrors(/*ooccoutW(_context,Minimization)*/std::cout,_printEvalErrors) ;
      /*ooccoutW(_context,Minimization)*/std::cout << endl ;
    } 

    if (_doEvalErrorWall) {
      fvalue = _maxFCN+1;
    }

    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;
  } else {
    _maxFCN = std::max(fvalue, _maxFCN);
  }
      
  // Optional logging
  if (_logfile) 
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN" << (_funct->isOffsetting()?"-offset":"") << " = " << setprecision(10) 
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;

  return fvalue;
}

double RooMinimizerFcnSemiAnalytic::DoDerivative(const double * x, unsigned int icoord) const{
    //std::cout<<"[DEBUG]:["<<__PRETTY_FUNCTION__<<"]:"<< "DoDerivative:"<<icoord<<std::endl; // run with -v3 to get it

  // Set the parameter values for this iteration
  for (int index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // compute Derivative
  if (_derivParamVec.size()<icoord) 
    {
    throw std::runtime_error(Form("Derivative vector is too small. (requested coordinate) %u < (derivative size) %lu",icoord,_derivParamVec.size()));
    }
  if (_derivParamVec[icoord] == nullptr) return DoNumericalDerivative(x,icoord);

  //verify correctness of analytical derivatives // DEBUG
  if (true){
      std::cout<<"[RooMinimizerFcnSemiAnalytic][DEBUG]"<< "Doing derivative for "<<icoord<<std::endl;
      std::cout<< "----------- NUMERICAL: "<< DoNumericalDerivative(x,icoord)<<std::endl;
      std::cout<< "----------- ANALYTICAL: "<< _derivParamVec[icoord]->getVal()<<std::endl;
  }
    
  return _derivParamVec[icoord]->getVal();

}

double RooMinimizerFcnSemiAnalytic::DoNumericalDerivative(const double * x,int icoord) const
{
    // ROOT 
    if (_useNumDerivatives==0) {
        //return ROOT::Math::MultiNumGradFunction::DoDerivative(x,icoord);
        //copied from ROOT::Math::MultiNumGradFunction::DoDerivative

        // calculate derivative using mathcore derivator class
        // step size can be changes using SetDerivPrecision()
        // this also use gsl.

        //static double kPrecision = std::sqrt ( std::numeric_limits<double>::epsilon() );
        //double x0 = std::abs(x[icoord]);
        //double step = (x0 > 0) ? kPrecision * x0 : kPrecision;
        // this seems to work better than above
        //double step = (x0>0) ? std::max( fgEps* x0, 8.0*kPrecision*(x0 + kPrecision) ) : kPrecision;
        //return ROOT::Math::Derivator::Eval(*fFunc, x, icoord, step);
        throw std::logic_error("Unimplemented");
    }
    if (_useNumDerivatives==1){
        /*
                 *
        ROOT use this has step above, maybe copy it dynamically? TODO
        double MultiNumGradFunction::fgEps = 0.001;
         
         double MultiNumGradFunction::DoDerivative (const double * x, unsigned int icoord  ) const 
               // calculate derivative using mathcore derivator class
            // step size can be changes using SetDerivPrecision()
         
            static double kPrecision = std::sqrt ( std::numeric_limits<double>::epsilon() );
            double x0 = std::abs(x[icoord]);
            //double step = (x0 > 0) ? kPrecision * x0 : kPrecision;
            // this seems to work better than above
            double step = (x0>0) ? std::max( fgEps* x0, 8.0*kPrecision*(x0 + kPrecision) ) : kPrecision;
         */

        // from gsl_deriv_central -> central_deriv
        /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
           x+h/2, x+h). Note that the central point is not used.  

           Compute the error using the difference between the 5-point and
           the 3-point rule (x-h,x,x+h). Again the central point is not
           used. */

        // smallest number representable in doubles
        static double kPrecision = std::sqrt ( std::numeric_limits<double>::epsilon() );
        static double fgEps = 0.001;

        double ax0 = std::abs(x[icoord]);
        //double step = (x0 > 0) ? kPrecision * x0 : kPrecision;
        // this seems to work better than above
        double h = (ax0>0) ? std::max( fgEps* ax0, 8.0*kPrecision*(ax0 + kPrecision) ) : kPrecision;
        //static const double h=std::numeric_limits<double>::epsilon()*8.;
        const double hhalf=h/2.;

        //double f0 = _funct->getVal();
        //double x0 = x[icoord];

        SetPdfParamVal(icoord,x[icoord]-h);
        double fm1 = _funct->getVal();
        SetPdfParamVal(icoord,x[icoord]+h);
        double fp1 = _funct->getVal();

        SetPdfParamVal(icoord,x[icoord]-hhalf);
        double fmh = _funct->getVal();
        SetPdfParamVal(icoord,x[icoord]+hhalf);
        double fph = _funct->getVal();

        double r3 = 0.5 * (fp1 - fm1);
        double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

        // don't do counts I'm not using.
        //double e3 = (fabs (fp1) + fabs (fm1)) * GSL_DBL_EPSILON;
        //double e5 = 2.0 * (fabs (fph) + fabs (fmh)) * GSL_DBL_EPSILON + e3;

        /* The next term is due to finite precision in x+h = O (eps * x) */

        //double dy = std::max (fabs (r3 / h), fabs (r5 / h)) *(fabs (x) / h) * GSL_DBL_EPSILON;

        /* The truncation error in the r5 approximation itself is O(h^4).
           However, for safety, we estimate the error from r5-r3, which is
           O(h^2).  By scaling h we will minimise this estimated error, not
           the actual truncation error in r5. */

        //*result = r5 / h;
        //*abserr_trunc = fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
        //*abserr_round = fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
        return r5 /h ;
    }
    throw std::runtime_error(Form("Numerical derivatives method un-implemented: %d",_useNumDerivatives));
}


#endif

