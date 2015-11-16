#ifndef EXTENDEDMINIMIZER
#define EXTENDEDMINIMIZER

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "Math/MinimizerOptions.h"

#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooAbsArg.h"
#include "RooLinkedListIter.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

//#ifdef __MAKECINT__
//#pragma link C++ class ExtendedMinimizer;
//#endif

class ExtendedMinimizer : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  ExtendedMinimizer( std::string MinimizerName, RooAbsPdf* pdf, RooAbsData* data );
  virtual ~ExtendedMinimizer();

  void SetPdf( RooAbsPdf* Pdf ) { fPdf = Pdf; }
  RooAbsPdf* GetPdf() { return fPdf; }

  void SetData( RooAbsData* Data ) { fData = Data; }
  RooAbsData* GetData() { return fData; }

  void SetFitResult( RooFitResult* Result ) { fFitResult = Result; }
  RooFitResult* GetFitResult() { return fFitResult; }

  void SetHesseMatrix( TMatrixDSym Matrix ) { fHesseMatrix = Matrix; }
  TMatrixDSym GetHesseMatrix() { return fHesseMatrix; }

  double GetMinNll() { return fMinNll; }

  // Steering
  int minimize( const RooCmdArg& arg1 = RooCmdArg::none(), const RooCmdArg& arg2 = RooCmdArg::none(), const RooCmdArg& arg3 = RooCmdArg::none(), const RooCmdArg& arg4 = RooCmdArg::none(), const RooCmdArg& arg5 = RooCmdArg::none(), const RooCmdArg& arg6 = RooCmdArg::none(), const RooCmdArg& arg7 = RooCmdArg::none(), const RooCmdArg& arg8 = RooCmdArg::none(), const RooCmdArg& arg9 = RooCmdArg::none(), const RooCmdArg& arg10 = RooCmdArg::none(), const RooCmdArg& arg11 = RooCmdArg::none(), const RooCmdArg& arg12 = RooCmdArg::none() );
  int minimize( const RooLinkedList& cmdList );

// ____________________________________________________________________________|__________
public:

  static RooCmdArg Eigen(Bool_t flag = kTRUE) { return RooCmdArg("Eigen",flag,0,0,0,0,0,0,0); }
  static RooCmdArg NumRetryFit(Int_t retry) { return RooCmdArg("NumRetryFit",retry,0,0,0,0,0,0,0); }
  static RooCmdArg Eps(double eps) { return RooCmdArg("Eps",0,0,eps,0,0,0,0,0); }
  static RooCmdArg NSigma(double nsigma) { return RooCmdArg("NSigma",0,0,nsigma,0,0,0,0,0); }
  static RooCmdArg Scan(Bool_t flag = kTRUE) { return RooCmdArg("Scan",flag,0,0,0,0,0,0,0); }
  static RooCmdArg Scan(const RooArgSet& scanArgs) { return RooCmdArg("Scan",kTRUE,0,0,0,0,0,&scanArgs,0); }
  static RooCmdArg Cond(const RooArgSet& condArgs) { return RooCmdArg("Cond",kTRUE,0,0,0,0,0,&condArgs,0); }
  static RooCmdArg ReuseMinimizer(Bool_t flag = kFALSE) { return RooCmdArg("ReuseMinimizer",flag,0,0,0,0,0,0,0); }
  static RooCmdArg ReuseNLL(Bool_t flag = kTRUE) { return RooCmdArg("ReuseNLL",flag,0,0,0,0,0,0,0); }

// ____________________________________________________________________________|__________
protected:

  void createMinimizer();
  int parseConfig( const RooLinkedList& cmdList );
  double UseLimits( const RooRealVar* par, double val );
  void eigenAnalysis();
  virtual double findSigma( double nll_min, double val_guess, double val_mle, const RooRealVar* par, const RooLinkedList& cmdList, double nsigma = +1.0, double precision = -1.0, double fitTol = 1.0 );
  virtual void findSigma();
  virtual int robustMinimize();

// ____________________________________________________________________________|__________
private:

  ClassDef(ExtendedMinimizer,1);

  TFile* fFile;
  RooAbsPdf* fPdf;
  RooAbsData* fData;
  RooAbsReal* fNll;
  RooMinimizer* fMinimizer;
  RooFitResult* fFitResult;
  TMatrixDSym fHesseMatrix;
  TVectorD* fEigenValues;
  TMatrixD* fEigenVectors;

  RooLinkedList fFitCmdList;
  RooLinkedList fScanCmdList;
  RooLinkedList fNllCmdList;

  Int_t fOptConst;
  Int_t fVerbose;
  Int_t fSave;
  Int_t fTimer;
  Int_t fPrintLevel;
  Int_t fDefaultStrategy;
  Int_t fHesse;
  Int_t fMinos;
  Int_t fScan;
  Int_t fNumee;
  Int_t fDoEEWall;
  Int_t fRetry;
  Int_t fEigen;
  Int_t fReuseMinimizer;
  Int_t fReuseNLL;
  Double_t fEps;
  Double_t fNsigma;
  Double_t fPrecision;
  const RooArgSet* fMinosSet;
  const RooArgSet* fCondSet;
  const RooArgSet* fScanSet;
  string  fMinimizerType;
  string fMinimizerAlgo;

  double fMinNll;

// ____________________________________________________________________________|__________
protected:

  // ClassDef(ExtendedMinimizer, 1)

};

#endif
