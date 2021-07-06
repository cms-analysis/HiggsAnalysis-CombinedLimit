/*****************************************************************************
 * Authors:                                                                  *
 *  Andrea Carlo Marini,   MIT (MA),        andrea.carlo.marini@cern.ch        *
 *  Code based on the equivalent roofit code/file.                                                                           *
 * Wed Oct 16 14:04:34 CEST 2019
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

#ifndef ROO_MINIMIZER_FCN_SEMIANALYTIC
#define ROO_MINIMIZER_FCN_SEMIANALYTIC

#include "Math/IFunction.h"
#include "Math/MultiNumGradFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "TMatrixDSym.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <iostream>
#include <fstream>
#include <map>

class RooMinimizerSemiAnalytic;

class RooMinimizerFcnSemiAnalytic : 
    //public ROOT::Math::MultiNumGradFunction { --> need to be constructed on a function. With the nominal Fcn?
    public ROOT::Math::IMultiGradFunction {
    //public ROOT::Math::IBaseFunctionMultiDim {

 public:

  RooMinimizerFcnSemiAnalytic(RooAbsReal *funct, RooMinimizerSemiAnalytic *context, std::map<std::string,RooAbsReal*>* knownDerivatives, bool verbose = false);
  RooMinimizerFcnSemiAnalytic(const RooMinimizerFcnSemiAnalytic& other);
  virtual ~RooMinimizerFcnSemiAnalytic();

  ROOT::Math::IMultiGradFunction* Clone() const override;
  unsigned int NDim() const override{ return _nDim; }

  RooArgList* GetFloatParamList() { return _floatParamList; }
  RooArgList* GetConstParamList() { return _constParamList; }
  RooArgList* GetInitFloatParamList() { return _initFloatParamList; }
  RooArgList* GetInitConstParamList() { return _initConstParamList; }

  void SetEvalErrorWall(Bool_t flag) { _doEvalErrorWall = flag ; }
  void SetPrintEvalErrors(Int_t numEvalErrors) { _printEvalErrors = numEvalErrors ; }
  Bool_t SetLogFile(const char* inLogfile);
  std::ofstream* GetLogFile() { return _logfile; }
  void SetVerbose(Bool_t flag=kTRUE) { _verbose = flag ; }

  Double_t& GetMaxFCN() { return _maxFCN; }
  Int_t GetNumInvalidNLL() { return _numBadNLL; }

  Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, 
		     Bool_t optConst, Bool_t verbose);
  void BackProp(const ROOT::Fit::FitResult &results);  
  void ApplyCovarianceMatrix(TMatrixDSym& V); 

  Int_t evalCounter() const { return _evalCounter ; }
  void zeroEvalCount() { _evalCounter = 0 ; }


 private:
  
  Double_t GetPdfParamVal(Int_t index);
  Double_t GetPdfParamErr(Int_t index);
  void SetPdfParamErr(Int_t index, Double_t value);
  void ClearPdfParamAsymErr(Int_t index);
  void SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal);

  inline Bool_t SetPdfParamVal(const Int_t &index, const Double_t &value) const;


  double DoEval(const double * x) const override;  
  double DoDerivative(const double * x,unsigned int icoord) const override;  
  // semi-analytic. The following function is called by DoDerivative if analytical derivative is not present.
  virtual double DoNumericalDerivative(const double * x,int icoord) const;  
  void updateFloatVec() ;

private:

  mutable Int_t _evalCounter ;
  
  RooAbsReal *_funct;
  // for name -> Derivative
  std::map<std::string,RooAbsReal*> *_knownDerivatives; // does not own them. Otherwise needs to design ad hoc copy constructors.

  RooMinimizerSemiAnalytic *_context;

  mutable double _maxFCN;
  mutable int _numBadNLL;
  mutable int _printEvalErrors;
  Bool_t _doEvalErrorWall;

  int _nDim;
  std::ofstream *_logfile;
  bool _verbose;

  RooArgList* _floatParamList;
  std::vector<RooAbsArg*> _floatParamVec ;
  std::vector<RooAbsReal*> _derivParamVec ; // vector of derivative functions - algined with the one above
  //
  RooArgList* _constParamList;
  RooArgList* _initFloatParamList;
  RooArgList* _initConstParamList;

  int _useNumDerivatives{1};// 0 = ROOT (not-impl), 1 re-implementation of gsl 5 point, 2 re-implementation of gsl 3 point

};

#endif
#endif
