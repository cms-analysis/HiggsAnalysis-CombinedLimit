/*****************************************************************************
 * Authors:                                                                  *
 *  Andrea Carlo Marini,   MIT (MA),        andrea.carlo.marini@cern.ch        *
 *  Code based on the equivalent roofit code/file.                                                                           *
 * Wed Oct 16 14:04:34 CEST 2019
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


/* ~Copy Of RooMINIMIZER, 
 * the difference is that is templated in RooMinimizerFcn. i
 * An instance of RooMinimizerCopy is present as typedef for compatibility.
 */

#ifndef __ROOFIT_NOROOMINIMIZER

#ifndef ROO_MINIMIZER_SEMIANALYTIC
#define ROO_MINIMIZER_SEMIANALYTIC

#include "TObject.h"
#include "TStopwatch.h"
#include <fstream>
#include "TMatrixDSymfwd.h"


#include "Fit/Fitter.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMinimizerFcnSemiAnalytic.h"
#include "RooMinimizerFcn.h"

class RooAbsReal ;
class RooFitResult ;
class RooArgList ;
class RooRealVar ;
class RooArgSet ;
class TH2F ;
class RooPlot ;

class RooMinimizerSemiAnalytic : public TObject {
public:

  RooMinimizerSemiAnalytic(RooAbsReal& function,std::map<std::string,RooAbsReal*>* knownDerivatives) ;
  virtual ~RooMinimizerSemiAnalytic() ;

  enum Strategy { Speed=0, Balance=1, Robustness=2 } ;
  enum PrintLevel { None=-1, Reduced=0, Normal=1, ExtraForProblem=2, Maximum=3 } ;
  void setStrategy(Int_t strat) ;
  void setErrorLevel(Double_t level) ;
  void setEps(Double_t eps) ;
  void optimizeConst(Int_t flag) ;
  void setEvalErrorWall(Bool_t flag) { fitterFcn()->SetEvalErrorWall(flag); }
  void setOffsetting(Bool_t flag) ;
  void setMaxIterations(Int_t n) ;
  void setMaxFunctionCalls(Int_t n) ; 

  RooFitResult* fit(const char* options) ;

  Int_t migrad() ;
  Int_t hesse() ;
  Int_t minos() ;
  Int_t minos(const RooArgSet& minosParamList) ;
  Int_t seek() ;
  Int_t simplex() ;
  Int_t improve() ;

  Int_t minimize(const char* type, const char* alg=0) ;

  RooFitResult* save(const char* name=0, const char* title=0) ;
  RooPlot* contour(RooRealVar& var1, RooRealVar& var2, 
		   Double_t n1=1, Double_t n2=2, Double_t n3=0,
		   Double_t n4=0, Double_t n5=0, Double_t n6=0, unsigned int npoints = 50) ;

  Int_t setPrintLevel(Int_t newLevel) ; 
  void setPrintEvalErrors(Int_t numEvalErrors) { fitterFcn()->SetPrintEvalErrors(numEvalErrors); }
  void setVerbose(Bool_t flag=kTRUE) { _verbose = flag ; fitterFcn()->SetVerbose(flag); }
  void setProfile(Bool_t flag=kTRUE) { _profile = flag ; }
  Bool_t setLogFile(const char* logf=0) { return fitterFcn()->SetLogFile(logf); }

  void setMinimizerType(const char* type) ;

  static void cleanup() ;
  static RooFitResult* lastMinuitFit(const RooArgList& varList=RooArgList()) ;

  void saveStatus(const char* label, Int_t status) { _statusHistory.push_back(std::pair<std::string,int>(label,status)) ; }

  Int_t evalCounter() const { return fitterFcn()->evalCounter() ; }
  void zeroEvalCount() { fitterFcn()->zeroEvalCount() ; }

  ROOT::Fit::Fitter* fitter() ;
  const ROOT::Fit::Fitter* fitter() const ;
  
protected:

  friend class RooAbsPdf ;
  void applyCovarianceMatrix(TMatrixDSym& V) ;

  void profileStart() ;
  void profileStop() ;

  inline Int_t getNPar() const { return fitterFcn()->NDim() ; }
  inline std::ofstream* logfile() { return fitterFcn()->GetLogFile(); }
  inline Double_t& maxFCN() { return fitterFcn()->GetMaxFCN() ; }
  
  const RooMinimizerFcnSemiAnalytic* fitterFcn() const {  return ( fitter()->GetFCN() ? (dynamic_cast<RooMinimizerFcnSemiAnalytic*> ( fitter()->GetFCN()) ): dynamic_cast<RooMinimizerFcnSemiAnalytic*>(_fcn) ) ; }
  RooMinimizerFcnSemiAnalytic* fitterFcn() { return ( fitter()->GetFCN() ? (dynamic_cast<RooMinimizerFcnSemiAnalytic*> (fitter()->GetFCN())) : dynamic_cast<RooMinimizerFcnSemiAnalytic*>(_fcn) ) ; }

private:

  Int_t       _printLevel ;
  Int_t       _status ;
  Bool_t      _optConst ;
  Bool_t      _profile ;
  RooAbsReal* _func ;

  Bool_t      _verbose ;
  TStopwatch  _timer ;
  TStopwatch  _cumulTimer ;
  Bool_t      _profileStart ;

  TMatrixDSym* _extV ;

  RooMinimizerFcnSemiAnalytic *_fcn; // Here it is the difference
  std::string _minimizerType;

  static ROOT::Fit::Fitter *_theFitter ;

  std::vector<std::pair<std::string,int> > _statusHistory ;

  RooMinimizerSemiAnalytic(const RooMinimizerSemiAnalytic&) ;

public:	
  ClassDef(RooMinimizerSemiAnalytic,0) // RooMinimizerSemiAnalytic this should be a title?
} ;


#endif

#endif
