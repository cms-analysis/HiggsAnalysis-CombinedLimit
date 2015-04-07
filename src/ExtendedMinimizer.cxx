#include "../interface/ExtendedMinimizer.h"

#include "TMath.h"
#include "TMatrixDSymEigen.h"

#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooCmdConfig.h"
#include "RooNLLVar.h"

#include "RooStats/RooStatsUtils.h"

ClassImp(ExtendedMinimizer)

// ____________________________________________________________________________|__________
// Constructor
ExtendedMinimizer::ExtendedMinimizer( std::string MinimizerName, RooAbsPdf* pdf, RooAbsData* data )
  :
  TNamed( MinimizerName.c_str(), MinimizerName.c_str() ),
  fPdf( pdf ),
  fData( data ),

  fOptConst( 2 ),
  fVerbose( 0 ),
  fSave ( 0 ),
  fTimer( 1 ),
  fPrintLevel( 1 ),
  fDefaultStrategy( 0 ),
  fHesse( 0 ),
  fMinos( 0 ),
  fScan( 0 ),
  fNumee( 5 ),
  fDoEEWall( 1 ),
  fRetry( 0 ),
  fEigen( 0 ),
  fReuseMinimizer( 0 ),
  fReuseNLL( 0 ),
  fEps( 1.0 ),
  // fNsigma( 1.959591794 ), // 1dof 95%
  // fNsigma( 2.44744765 ), // 2dof 95%
  fNsigma( 1 ), // 1sigma 1dof
  fPrecision( 0.005 ),

  fMinimizerType( "Minuit2" ),
  fMinimizerAlgo( "Migrad" )
{

  fMinosSet = NULL;
  fCondSet  = NULL;
  fScanSet  = NULL;

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(), fMinimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

  coutP(InputArguments) << "ExtendedMinimizer::ExtendedMinimizer(" << fName <<") created" << endl;
}

// ____________________________________________________________________________|__________
// Destructor
ExtendedMinimizer::~ExtendedMinimizer()
{
  // TODO
}

// ____________________________________________________________________________|__________
// Minimize function with iterative retry strategy adopted, simplified and
// extended from RooAbsPdf::fitTo()
int ExtendedMinimizer::minimize( const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
                                 const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8,
                                 const RooCmdArg& arg9, const RooCmdArg& arg10, const RooCmdArg& arg11, const RooCmdArg& arg12 )
{
  RooLinkedList l;
  l.Add((TObject*)&arg1);   l.Add((TObject*)&arg2);
  l.Add((TObject*)&arg3);   l.Add((TObject*)&arg4);
  l.Add((TObject*)&arg5);   l.Add((TObject*)&arg6);
  l.Add((TObject*)&arg7);   l.Add((TObject*)&arg8);
  l.Add((TObject*)&arg9);   l.Add((TObject*)&arg10);
  l.Add((TObject*)&arg11);  l.Add((TObject*)&arg12);
  return minimize(l);
}

// ____________________________________________________________________________|__________
void ExtendedMinimizer::createMinimizer()
{
  if (!fReuseNLL) {
    // if (fNll) delete fNll;

    fNll = fPdf->createNLL(*fData, fNllCmdList);
  }

  if (!fReuseMinimizer) {
    // if (fMinimizer) delete fMinimizer;

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(), fMinimizerAlgo.c_str());
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

    fMinimizer = new RooMinimizer(*fNll);
    fMinimizer->setPrintLevel(fPrintLevel);
    fMinimizer->optimizeConst(fOptConst);
    fMinimizer->setMinimizerType(fMinimizerType.c_str());
    fMinimizer->setEvalErrorWall(fDoEEWall);
    fMinimizer->setPrintEvalErrors(fNumee);
    fMinimizer->setVerbose(fVerbose);
    fMinimizer->setProfile(fTimer);
    fMinimizer->setStrategy(fDefaultStrategy);
    fMinimizer->setEps(fEps);
  }
}

// ____________________________________________________________________________|__________
int ExtendedMinimizer::parseConfig( const RooLinkedList& cmdList )
{
  RooCmdConfig pc(Form("ExtendedMinimizer::parseConfig(%s)", GetName()));

  fFitCmdList = cmdList;
  fScanCmdList = cmdList;
  pc.filterCmdList(fScanCmdList, "Scan,Minos,Save");
  fNllCmdList = pc.filterCmdList(fFitCmdList, "NumCPU,Constrained,Constrain,CloneData,GlobalObservables,GlobalObservablesTag,OffsetLikelihood");

  pc.defineInt("optConst",     "Optimize",         0, fOptConst);
  pc.defineInt("verbose",      "Verbose",          0, fVerbose);
  pc.defineInt("doSave",       "Save",             0, fSave);
  pc.defineInt("doTimer",      "Timer",            0, fTimer);
  pc.defineInt("plevel",       "PrintLevel",       0, fPrintLevel);
  pc.defineInt("strat",        "Strategy",         0, fDefaultStrategy);
  pc.defineInt("hesse",        "Hesse",            0, fHesse);
  pc.defineInt("minos",        "Minos",            0, fMinos);
  pc.defineInt("scan",         "Scan",             0, fScan);
  pc.defineInt("numee",        "PrintEvalErrors",  0, fNumee);
  pc.defineInt("doEEWall",     "EvalErrorWall",    0, fDoEEWall);
  pc.defineInt("retry",        "NumRetryFit",      0, fRetry);
  pc.defineInt("eigen",        "Eigen",            0, fEigen);
  pc.defineInt("reminim",      "ReuseMinimizer",   0, fReuseMinimizer);
  pc.defineInt("renll",        "ReuseNLL",         0, fReuseNLL);
  pc.defineDouble("eps",       "Eps",              0, fEps);
  pc.defineDouble("nsigma",    "NSigma",           0, fNsigma);
  pc.defineDouble("precision", "Precision",        0, fPrecision);
  pc.defineString("mintype",   "Minimizer",        0, fMinimizerType.c_str());
  pc.defineString("minalg",    "Minimizer",        1, fMinimizerAlgo.c_str());
  pc.defineObject("minosSet",  "Minos",            0, 0);
  pc.defineObject("condSet",   "Cond",             0, 0);
  pc.defineObject("scanSet",   "Scan",             0, 0);
  pc.defineSet("cPars",        "Constrain",        0, 0);
  pc.defineMutex("Scan",       "Minos");

  pc.process(fFitCmdList);
  if (!pc.ok(kTRUE)) {
    return -1;
  }

  fOptConst        = pc.getInt("optConst");
  fVerbose         = pc.getInt("verbose");
  fSave            = pc.getInt("doSave");
  fTimer           = pc.getInt("doTimer");
  fPrintLevel      = pc.getInt("plevel");
  fDefaultStrategy = pc.getInt("strat");
  fHesse           = pc.getInt("hesse");
  fMinos           = pc.getInt("minos");
  fScan            = pc.getInt("scan");
  fNumee           = pc.getInt("numee");
  fDoEEWall        = pc.getInt("doEEWall");
  fRetry           = pc.getInt("retry");
  fEigen           = pc.getInt("eigen");
  fReuseMinimizer  = pc.getInt("reminim");
  fReuseNLL        = pc.getInt("renll");
  fEps             = pc.getDouble("eps");
  fNsigma          = pc.getDouble("nsigma");
  fPrecision       = pc.getDouble("precision");
  fMinosSet        = static_cast<RooArgSet*>(pc.getObject("minosSet"));
  fCondSet         = static_cast<RooArgSet*>(pc.getObject("condSet"));
  fScanSet         = static_cast<RooArgSet*>(pc.getObject("scanSet"));
  fMinimizerType   = string(pc.getString("mintype", "Minuit2"));
  fMinimizerAlgo   = string(pc.getString("minalg", "Migrad"));

  return 1;
}

// ____________________________________________________________________________|__________
void ExtendedMinimizer::eigenAnalysis()
{
  int n = fHesseMatrix.GetNrows();

  // Construct reduced Hessian matrix
  TMatrixDSym Gred(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Gred(i,j) = fHesseMatrix(i,j) / sqrt( fHesseMatrix(i,i) * fHesseMatrix(j,j) );
    }
  }

  // Perform eigenvalue analysis using ROOT standard tools
  TMatrixDSymEigen Geigen(Gred);
  fEigenValues = new TVectorD(Geigen.GetEigenValues());
  fEigenVectors = new TMatrixD(Geigen.GetEigenVectors());

  // Simple printing of eigenvalues and eigenvectors
  fEigenValues->Print();
  fEigenVectors->Print();
}

// ____________________________________________________________________________|__________
// Robust minimization, using a iterative retry strategy
int ExtendedMinimizer::robustMinimize()
{
  int strat = fDefaultStrategy;
  int retry = fRetry;
  int status = fMinimizer->minimize(fMinimizerType.c_str(), fMinimizerAlgo.c_str());

  while (status != 0 && status != 1 && strat < 2 && retry > 0) {
    strat++;
    retry--;
    coutW(ObjectHandling) << "ExtendedMinimizer::robustMinimize(" << fName << ") fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    fMinimizer->setStrategy(strat);
    status =  fMinimizer->minimize(fMinimizerType.c_str(), fMinimizerAlgo.c_str());
  }

  if (status != 0 && status != 1) {
    coutE(ObjectHandling) << "ExtendedMinimizer::robustMinimize(" << fName << ") fit failed with status " << status << endl;
  }

  fMinimizer->setStrategy(fDefaultStrategy);

  return status;
}

// ____________________________________________________________________________|__________
// Minimize function  adopted, simplified and extended from RooAbsPdf::fitTo()
int ExtendedMinimizer::minimize( const RooLinkedList& cmdList )
{
  parseConfig(cmdList);

  if (fCondSet) {
    RooArgSet* attachedSet = fNll->getVariables();
    for (RooLinkedListIter it = fCondSet->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      if (RooRealVar* var = dynamic_cast<RooRealVar*>( attachedSet->find(v->GetName()))) {
        var->setVal(v->getVal());
        var->setConstant(1);
      }
    }
  }

  createMinimizer();
  int status = robustMinimize();

  // Evaluate errors with Hesse
  if (fHesse) {
    fMinimizer->hesse();
  }

  // Obtain Hessian matrix either from patched Minuit or after inversion
  // TMatrixDSym G = Minuit2::MnHesse::lastHessian();
  TMatrixDSym G = ((TMatrixDSym)fMinimizer->lastMinuitFit()->covarianceMatrix()).Invert();
  int n = G.GetNrows();
  fHesseMatrix.ResizeTo(n,n);
  fHesseMatrix = G;

  // Eigenvalue and eigenvector analysis
  if (fEigen) {
    eigenAnalysis();
  }

  // Evaluate errors with Minos
  if (fMinos) {
    if (fMinosSet) {
      fMinimizer->minos(*fMinosSet);
    } else {
      fMinimizer->minos();
    }
  }

  fMinNll = fNll->getVal();

  if (fScan) {
    findSigma();
  }

  // Return fit result
  if (fSave) {
    string name = Form("fitresult_%s_%s",GetName(),fData->GetName());
    string title = Form("Result of fit of p.d.f. %s to dataset %s",GetName(),fData->GetName());
    coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName << ") saving results as " << name << endl;
    fFitResult = fMinimizer->save(name.c_str(),title.c_str());
  }

  if (fCondSet) {
    RooArgSet* attachedSet = fNll->getVariables();
    for (RooLinkedListIter it = fCondSet->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      if (RooRealVar* var = dynamic_cast<RooRealVar*>( attachedSet->find(v->GetName()))) {
        var->setVal(v->getVal());
        var->setConstant(v->isConstant());
      }
    }
  }

  return status;
}

// ____________________________________________________________________________|__________
Double_t ExtendedMinimizer::UseLimits(const RooRealVar* par, Double_t val)
{
  if (val < par->getMin()) {
    coutP(ObjectHandling) << "ExtendedMinimizer::UseLimits(" << fName << ") " << Form("%s = %g limited by minimum at %g", par->GetName(), val, par->getMin()) << endl;
    return par->getMin();
  } else if (val > par->getMax()) {
    coutP(ObjectHandling) << "ExtendedMinimizer::UseLimits(" << fName << ") " << Form("%s = %g limited by maximum at %g", par->GetName(), val, par->getMax()) << endl;
    return par->getMax();
  }
  return val;
}

// ____________________________________________________________________________|__________
void ExtendedMinimizer::findSigma()
{
  if (!fScanSet && fScanSet->getSize() == 0) {
    RooArgSet* attachedSet = fNll->getVariables();
    RemoveConstantParameters(attachedSet);
    fScanSet = attachedSet;
  }

  for (RooLinkedListIter it = fScanSet->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
    RooArgSet* vars = fPdf->getVariables();
    RooStats::RemoveConstantParameters(vars);
    vars->add(*v, kTRUE);
    RooArgSet* snap = dynamic_cast<RooArgSet*>(vars->snapshot());

    fNsigma = fabs(fNsigma);
    Double_t val = v->getVal();
    Double_t err = fNsigma * v->getError();

    double sMinNll = fMinNll;

    Double_t shi = findSigma(sMinNll, val+err, val, v, fScanCmdList, fNsigma, fPrecision, fEps);
    *vars = *snap;
    Double_t slo = findSigma(sMinNll, val-err, val, v, fScanCmdList, -fNsigma, fPrecision, fEps);
    *vars = *snap;

    fMinNll = sMinNll;

    v->setAsymError (std::isnan(slo)?1.0:slo, std::isnan(shi)?-1.0:shi);

    coutI(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName << ") ";
    v->Print();

    delete vars;
    delete snap;
  }
}

// ____________________________________________________________________________|__________
// Find the value of sigma evaluated at a specified nsigma, assuming NLL -2logL is roughly parabolic in par.
// The precision is specified as a fraction of the error, or based on the Minuit default tolerance.
// Based on https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/HiggsPhys/HSG3/WWDileptonAnalysisCode/HWWStatisticsCode/trunk/macros/findSigma.C
// by Aaron Armbruster <aaron.james.armbruster@cern.ch> and adopted by Tim Adye <T.J.Adye@rl.ac.uk>.
Double_t ExtendedMinimizer::findSigma( double nll_min, double val_guess, double val_mle, const RooRealVar* par, const RooLinkedList& cmdList, double nsigma, double precision, double fitTol )
{
  const int maxiter = 25;
  val_guess = UseLimits (par, val_guess);
  int direction = nsigma>=0.0 ? +1 : -1;
  int nDamping = 1;
  double damping_factor = 1.0;
  std::map<double,double> guess_to_corr;
  double tmu = TMath::QuietNaN();

  RooLinkedList sFitCmdList(fFitCmdList);
  RooLinkedList sScanCmdList(fScanCmdList);
  RooLinkedList sNllCmdList(fNllCmdList);
  RooLinkedList scanCmdList(cmdList);

  if (precision <= 0.0) {
    // RooFit default tolerance is 1.0.
    double eps = 0.001 * (fitTol > 0.0 ? fitTol : 1.0);
    precision = 5.0*eps / (nsigma*nsigma);
  }

  int iter = 0;
  for (; iter < maxiter; iter++) {
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") " << Form("Parameter %s %+gsigma iteration %d: start %g (MLE%+g)", par->GetName(), nsigma, iter+1, val_guess, val_guess-val_mle) << endl;
    double val_pre = val_guess;

    RooArgSet* poiSet = dynamic_cast<RooArgSet*>(RooArgSet(*par).snapshot());
    dynamic_cast<RooRealVar*>(poiSet->first())->setVal(val_guess);

    const RooCmdArg& arg1 = Cond(*poiSet); scanCmdList.Add((TObject*)&arg1);
    const RooCmdArg& arg2 = ReuseNLL(1); scanCmdList.Add((TObject*)&arg2);
    const RooCmdArg& arg3 = Scan(0); scanCmdList.Add((TObject*)&arg3);

    minimize(scanCmdList);
    double nll = fMinNll;
    delete poiSet;

    tmu = 2.0 * (nll-nll_min);
    double sigma_guess = fabs(val_guess-val_mle);
    if (tmu > 0.01) sigma_guess /= sqrt(tmu);
    else sigma_guess *= 10.0;  // protect against tmu<=0, and also don't move too far

    double corr = damping_factor*(val_pre - val_mle - nsigma*sigma_guess);

    for (std::map<double,double>::iterator iguess= guess_to_corr.begin(); iguess != guess_to_corr.end(); ++iguess) {
      if (fabs(iguess->first - val_pre) < direction*val_pre*0.02) {
        damping_factor *= 0.8;
        coutW(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") Changing damping factor to " << damping_factor << endl;
        if (nDamping++ > 10) {
          nDamping = 1;
          damping_factor = 1.0;
        }
        corr *= damping_factor;
        break;
      }
    }

    // subtract off the difference in the new and damped correction
    val_guess -= corr;
    guess_to_corr[val_pre] = corr;
    val_guess= UseLimits (par, val_guess);
    double relprecision = precision*fabs(val_guess-val_mle);
    double delta = val_guess-val_pre;

    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") " << Form("%s %.3f (MLE%+.3f) -> %.3f (MLE%+.3f), change %+.3f, precision %.3f, -2lnL %.4f, sigma(guess) %.3f", par->GetName(), val_pre, val_pre-val_mle, val_guess, val_guess-val_mle, delta, relprecision, tmu, sigma_guess) << endl;

    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") NLL:                 " << nll << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") delta(NLL):          " << nll-nll_min << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") nsigma*sigma(pre):   " << fabs(val_pre-val_mle) << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") sigma(guess):        " << sigma_guess << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") par(guess):          " << val_guess+corr << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") best-fit val:        " << val_mle << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") tmu:                 " << tmu << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") Precision:           " << direction*val_guess*precision << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") Correction:          " << -corr << endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") nsigma*sigma(guess): " << fabs(val_guess-val_mle) << endl;
    cout << endl;

    if (fabs(delta) <= relprecision) break;
  }
  if (iter >= maxiter) {
    cerr << "findSigma failed after " << iter << " iterations" << endl;
    return TMath::QuietNaN();
  }

  double err= val_guess-val_mle;
  coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") " << Form("%s %+gsigma = %.3f at -2lnL = %.4f after %d iterations", par->GetName(), nsigma, err, tmu, iter+1) << endl;

  RooLinkedList fFitCmdList(sFitCmdList);
  RooLinkedList fScanCmdList(sScanCmdList);
  RooLinkedList fNllCmdList(sNllCmdList);

  return err;
}
