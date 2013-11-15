#include "../interface/RooMinimizerOpt.h"

#include <stdexcept>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooMsgService.h>

#include <Math/MinimizerOptions.h>

#include <iomanip>

using namespace std;

RooMinimizerOpt::RooMinimizerOpt(RooAbsReal& function) :
    RooMinimizer(function)
{
    delete _fcn;
    _fcn = new RooMinimizerFcnOpt(_func,this,_verbose); 
    setEps(ROOT::Math::MinimizerOptions::DefaultTolerance());
}

Double_t
RooMinimizerOpt::edm()
{
    if (_theFitter == 0) throw std::logic_error("Must have done a fit before calling edm()");
    return _theFitter->Result().Edm();    
}

Int_t RooMinimizerOpt::minimize(const char* type, const char* alg)
{
  if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
    static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  } else {
    _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  }

  _theFitter->Config().SetMinimizer(type,alg);

  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;

  bool ret = _theFitter->FitFCN(*_fcn);
  _status = ((ret) ? _theFitter->Result().Status() : -1);

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  _fcn->BackProp(_theFitter->Result());

  saveStatus("MINIMIZE",_status) ;

  return _status ;
}
//_____________________________________________________________________________
Int_t RooMinimizerOpt::improve()
{
  // Execute IMPROVE. Changes in parameter values
  // and calculated errors are automatically
  // propagated back the RooRealVars representing
  // the floating parameters in the MINUIT operation

  if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
    static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  } else {
  _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  }
  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;

  _theFitter->Config().SetMinimizer(_minimizerType.c_str(),"migradimproved");
  bool ret = _theFitter->FitFCN(*_fcn);
  _status = ((ret) ? _theFitter->Result().Status() : -1);

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  _fcn->BackProp(_theFitter->Result());

  saveStatus("IMPROVE",_status) ;

  return _status ;
}

Int_t RooMinimizerOpt::migrad()
{
  // Execute MIGRAD. Changes in parameter values
  // and calculated errors are automatically
  // propagated back the RooRealVars representing
  // the floating parameters in the MINUIT operation

  if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
    static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  } else {
  _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
            _optConst,_verbose) ;
  }

  profileStart() ;
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
  RooAbsReal::clearEvalErrorLog() ;

  _theFitter->Config().SetMinimizer(_minimizerType.c_str(),"migrad");
  bool ret = _theFitter->FitFCN(*_fcn);
  _status = ((ret) ? _theFitter->Result().Status() : -1);

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
  profileStop() ;
  _fcn->BackProp(_theFitter->Result());

  saveStatus("MIGRAD",_status) ;

  return _status ;
}



//_____________________________________________________________________________
Int_t RooMinimizerOpt::hesse()
{
  // Execute HESSE. Changes in parameter values
  // and calculated errors are automatically
  // propagated back the RooRealVars representing
  // the floating parameters in the MINUIT operation

  if (_theFitter->GetMinimizer()==0) {
    coutW(Minimization) << "RooMinimizerOpt::hesse: Error, run Migrad before Hesse!"
            << endl ;
    _status = -1;
  }
  else {

    if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
          static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    } else {
          _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    }

    profileStart() ;
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
    RooAbsReal::clearEvalErrorLog() ;

    _theFitter->Config().SetMinimizer(_minimizerType.c_str());
    bool ret = _theFitter->CalculateHessErrors();
    _status = ((ret) ? _theFitter->Result().Status() : -1);

    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
    profileStop() ;
    _fcn->BackProp(_theFitter->Result());

    saveStatus("HESSE",_status) ;
  
  }

  return _status ;

}

//_____________________________________________________________________________
Int_t RooMinimizerOpt::minos()
{
  // Execute MINOS. Changes in parameter values
  // and calculated errors are automatically
  // propagated back the RooRealVars representing
  // the floating parameters in the MINUIT operation

  if (_theFitter->GetMinimizer()==0) {
    coutW(Minimization) << "RooMinimizerOpt::minos: Error, run Migrad before Minos!"
            << endl ;
    _status = -1;
  }
  else {

    if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
          static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    } else {
          _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    }
    profileStart() ;
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
    RooAbsReal::clearEvalErrorLog() ;

    _theFitter->Config().SetMinimizer(_minimizerType.c_str());
    bool ret = _theFitter->CalculateMinosErrors();
    _status = ((ret) ? _theFitter->Result().Status() : -1);

    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
    profileStop() ;
    _fcn->BackProp(_theFitter->Result());

    saveStatus("MINOS",_status) ;

  }

  return _status ;

}


//_____________________________________________________________________________
Int_t RooMinimizerOpt::minos(const RooArgSet& minosParamList)
{
  // Execute MINOS for given list of parameters. Changes in parameter values
  // and calculated errors are automatically
  // propagated back the RooRealVars representing
  // the floating parameters in the MINUIT operation

  if (_theFitter->GetMinimizer()==0) {
    coutW(Minimization) << "RooMinimizerOpt::minos: Error, run Migrad before Minos!"
            << endl ;
    _status = -1;
  }
  else if (minosParamList.getSize()>0) {

    if (typeid(*_fcn) == typeid(RooMinimizerFcnOpt)) {
          static_cast<RooMinimizerFcnOpt*>(_fcn)->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    } else {
          _fcn->Synchronize(_theFitter->Config().ParamsSettings(),
                  _optConst,_verbose) ;
    }
    profileStart() ;
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;
    RooAbsReal::clearEvalErrorLog() ;

    // get list of parameters for Minos
    TIterator* aIter = minosParamList.createIterator() ;
    RooAbsArg* arg ;
    std::vector<unsigned int> paramInd;
    while((arg=(RooAbsArg*)aIter->Next())) {
      RooAbsArg* par = _fcn->GetFloatParamList()->find(arg->GetName());
      if (par && !par->isConstant()) {
    Int_t index = _fcn->GetFloatParamList()->index(par);
    paramInd.push_back(index);
      }
    }
    delete aIter ;

    if (paramInd.size()) {
      // set the parameter indeces
      _theFitter->Config().SetMinosErrors(paramInd);

      _theFitter->Config().SetMinimizer(_minimizerType.c_str());
      bool ret = _theFitter->CalculateMinosErrors();
      _status = ((ret) ? _theFitter->Result().Status() : -1);

    }

    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;
    profileStop() ;
    _fcn->BackProp(_theFitter->Result());

    saveStatus("MINOS",_status) ;

  }

  return _status ;
}






RooMinimizerFcnOpt::RooMinimizerFcnOpt(RooAbsReal *funct, RooMinimizer *context,  bool verbose) :
    RooMinimizerFcn(funct, context, verbose)
{
}

RooMinimizerFcnOpt::RooMinimizerFcnOpt(const RooMinimizerFcnOpt &other) :
    RooMinimizerFcn(other._funct, other._context, other._verbose),
    _vars(other._vars), _vals(other._vals), _hasOptimzedBounds(other._hasOptimzedBounds), _optimzedBounds(other._optimzedBounds)
{

}

ROOT::Math::IBaseFunctionMultiDim* 
RooMinimizerFcnOpt::Clone() const
{
      return new RooMinimizerFcnOpt(*this);
}

double
RooMinimizerFcnOpt::DoEval(const double * x) const 
{
  // Set the parameter values for this iteration
  for (int index = 0; index < _nDim; index++) {
      if (_vals[index]!=x[index]) {
          RooRealVar* par = _vars[index];
          if (_verbose) cout << par->GetName() << "=" << x[index] << ", " ;
          if (_hasOptimzedBounds[index]) {
              // boundary transformation done internally
              par->setVal( _optimzedBounds[index].transform(x[index]) );
              _vals[index] = x[index];
          } else {
              par->setVal(x[index]);
              _vals[index] = par->getVal(); // might not work otherwise if x is out of the boundary
          }
      }
  }
  if (_logfile) {
      for (int index = 0; index < _nDim; index++) {
          (*_logfile) << x[index] << " " ;
      }
  }

  // Calculate the function for these parameters
  double fvalue = _funct->getVal();
  if (RooAbsPdf::evalError() || RooAbsReal::numEvalErrors()>0) {

    if (_printEvalErrors>=0) {

      if (_doEvalErrorWall) {
        oocoutW(_context,Minimization) << "RooMinimizerFcn: Minimized function has error status." << endl 
				       << "Returning maximum FCN so far (" << _maxFCN 
				       << ") to force MIGRAD to back out of this region. Error log follows" << endl ;
      } else {
        oocoutW(_context,Minimization) << "RooMinimizerFcn: Minimized function has error status but is ignored" << endl ;
      } 

      TIterator* iter = _floatParamList->createIterator() ;
      RooRealVar* var ;
      Bool_t first(kTRUE) ;
      ooccoutW(_context,Minimization) << "Parameter values: " ;
      while((var=(RooRealVar*)iter->Next())) {
        if (first) { first = kFALSE ; } else ooccoutW(_context,Minimization) << ", " ;
        ooccoutW(_context,Minimization) << var->GetName() << "=" << var->getVal() ;
      }
      delete iter ;
      ooccoutW(_context,Minimization) << endl ;
      
      RooAbsReal::printEvalErrors(ooccoutW(_context,Minimization),_printEvalErrors) ;
      ooccoutW(_context,Minimization) << endl ;
    } 

    if (_doEvalErrorWall) {
      fvalue = _maxFCN ;
    }

    RooAbsPdf::clearEvalError() ;
    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;
  } else if (fvalue>_maxFCN) {
    _maxFCN = fvalue ;
  }
      
  // Optional logging
  if (_logfile) 
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN = " << setprecision(10) 
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }
  return fvalue;
}

Bool_t RooMinimizerFcnOpt::Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, 
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
    oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: parameter " 
                     << par->GetName() << " is now floating." << endl ;
      }
    } 

    // Test if value changed
    if (par->getVal()!= oldpar->getVal()) {
      constValChange=kTRUE ;      
      if (verbose) {
    oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: value of constant parameter " 
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
    oocoutW(_context,Minimization) << "RooMinimizerFcn::fit: Error, non-constant parameter " 
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
      oocoutW(_context,Minimization) << "RooMinimizerFcn::synchronize: WARNING: no initial error estimate available for "
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
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: value of parameter " 
                     << par->GetName() << " changed from " << oldVar 
                     << " to " << par->getVal() << endl ;
    }
      }
      parameters[index].Fix();
      constStatChange=kTRUE ;
      if (verbose) {
    oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: parameter " 
                       << par->GetName() << " is now fixed." << endl ;
      }

    } else if (par->isConstant() && oldFixed) {
      
      // Parameter changes constant -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
    parameters[index].SetValue(par->getVal());
    constValChange=kTRUE ;

    if (verbose) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: value of fixed parameter " 
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
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: parameter " 
                     << par->GetName() << " is now floating." << endl ;
    }
      } 

      // Parameter changes constant -> floating : update all if necessary      
      if (oldVar!=par->getVal() || oldVlo!=pmin || oldVhi != pmax || oldVerr!=pstep) {
    parameters[index].SetValue(par->getVal()); 
    parameters[index].SetStepSize(pstep);
        if (par->hasMin() && par->hasMax() ) {
           parameters[index].SetLimits(pmin,pmax);  
        } else if (par->hasMin() ) {
           parameters[index].SetLowerLimit(pmin);  
        } else if (par->hasMax() ) {
           parameters[index].SetUpperLimit(pmax);  
        }
      }

      // Inform user about changes in verbose mode
      if (verbose) {
    // if ierr<0, par was moved from the const list and a message was already printed

    if (oldVar!=par->getVal()) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: value of parameter " 
                     << par->GetName() << " changed from " << oldVar << " to " 
                     << par->getVal() << endl ;
    }
    if (oldVlo!=pmin || oldVhi!=pmax) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: limits of parameter " 
                     << par->GetName() << " changed from [" << oldVlo << "," << oldVhi 
                     << "] to [" << pmin << "," << pmax << "]" << endl ;
    }

    // If oldVerr=0, then parameter was previously fixed
    if (oldVerr!=pstep && oldVerr!=0) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: error/step size of parameter " 
                     << par->GetName() << " changed from " << oldVerr << " to " << pstep << endl ;
    }
      }      
    }
  }

  if (optConst) {
    if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: set of constant parameters changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange) ;
    } else if (constValChange) {
      oocoutI(_context,Minimization) << "RooMinimizerFcn::synchronize: constant parameter values changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange) ;
    }
    
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;  

  }

  updateFloatVec() ;

  initStdVects();

  return 0 ;  

}

void RooMinimizerFcnOpt::initStdVects() const {
  _vars.resize(_floatParamList->getSize());
  _vals.resize(_vars.size());
  _hasOptimzedBounds.resize(_vars.size());
  _optimzedBounds.resize(_vars.size());
  std::vector<RooRealVar *>::iterator itv = _vars.begin();
  std::vector<double>::iterator it = _vals.begin();
  RooLinkedListIter iter = _floatParamList->iterator();
  std::fill(_hasOptimzedBounds.begin(), _hasOptimzedBounds.end(), false);
  unsigned int i = 0;
  for (TObject *a = iter.Next(); a != 0; a = iter.Next(), ++itv, ++it, ++i) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      if (rrv == 0) throw std::logic_error(Form("Float param not a RooRealVar but a %s", a->ClassName()));
      *itv = rrv; 
      *it  = rrv->getVal();
      if (rrv->getAttribute("optimizeBounds")) {
        _hasOptimzedBounds[i] = true;
        OptBound &b = _optimzedBounds[i];
        b.hardMax = rrv->getMax("optimizeBoundRange");
        b.hardMin = rrv->getMin("optimizeBoundRange");
        double x0 = 0.5*(b.hardMax+b.hardMin);
        b.softMax = 0.75*b.hardMax + 0.25*x0;
        b.softMin = 0.75*b.hardMin + 0.25*x0;
      }
  }
}
