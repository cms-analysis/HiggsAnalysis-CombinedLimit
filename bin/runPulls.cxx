// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@nikhef.nl
// Date        : 2014-02-10
// Description : Compute pulls and impact on the POI

#include <string>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "../interface/AtlasPdfs.h"
#include "RooProdPdf.h"   
#include "../interface/RooMultiPdf.h"
#include "RooGaussian.h"
#include "RooRealSumPdf.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"

#include "boost/program_options.hpp"
#include "boost/program_options/cmdline.hpp"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/variables_map.hpp"

#include "../interface/ExtendedMinimizer.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components );
bool AlmostEqualUlpsAndAbs( float A, float B, float maxDiff, int maxUlpsDiff );
vector<string> parseString(string str, string sep);

// ____________________________________________________________________________|__________
int main(int argc, char** argv)
{
  TStopwatch timer_pulls;
  timer_pulls.Start();

  string inFileName      = "path/to/workspace.root";
  string poiName         = "mu";
  string profileName     = "";
  string wsName          = "combined";
  string modelConfigName = "ModelConfig";
  string dataName        = "combData";
  string folder          = "test";
  string variable        = "";

  string minimizerType  = "Minuit2";
  string minimizerAlgo  = "Migrad";
  int defaultStrategy   = 0;
  int binnedLikelihood  = 1;
  int fixCache          = 1;
  int fixMulti          = 1;
  int offsetting        = 1;
  int constOpt          = 2;
  double eps            = 1.0;
  int numCPU            = 1;
  double precision      = 0.01;

  using namespace boost;
  namespace po = boost::program_options;
  po::options_description desc( "Program options" );
  desc.add_options()
    ( "help,h"        , "Print this help message")
    ( "input,i"       , po::value<string>( &inFileName )->default_value( inFileName )              , "File to run over." )
    ( "poi,p"         , po::value<string>( &poiName )->default_value( poiName )                    , "POIs to use." )
    ( "profile,o"     , po::value<string>( &profileName )->default_value( profileName )            , "POIs to profile." )
    ( "workspace,w"   , po::value<string>( &wsName )->default_value( wsName )                      , "WS to grab." )
    ( "modelconfig,m" , po::value<string>( &modelConfigName )->default_value( modelConfigName )    , "MC to load." )
    ( "data,d"        , po::value<string>( &dataName )->default_value( dataName )                  , "Data to use." )
    ( "folder,f"      , po::value<string>( &folder )->default_value( folder )                      , "Folder name for output." )
    ( "parameter,n"   , po::value<string>( &variable )->default_value( variable )                  , "Parameter to rank." )
    ( "precision,r"   , po::value<double>( &precision )->default_value( precision )                , "Precision to scan." )
    ( "minimizerType" , po::value<string>( &minimizerType )->default_value( minimizerType )        , "Minimizer type." )
    ( "minimizerAlgo" , po::value<string>( &minimizerAlgo )->default_value( minimizerAlgo )        , "Minimizer algorithm." )
    ( "strategy"      , po::value<int>( &defaultStrategy )->default_value( defaultStrategy )       , "Default strategy." )
    ( "numCPU"        , po::value<int>( &numCPU )->default_value( numCPU )                         , "Number of CPUs." )
    ( "binned"        , po::value<int>( &binnedLikelihood )->default_value( binnedLikelihood )     , "Binned likelihood." )
    ( "offset"        , po::value<int>( &offsetting )->default_value( offsetting )                 , "Offset LLH." )
    ( "optimize"      , po::value<int>( &constOpt )->default_value( constOpt )                     , "Optimize constant terms." )
    ( "starfix"       , po::value<int>( &fixCache )->default_value( fixCache )                     , "Fix StarMomentMorph cache." )
    ( "multifix"      , po::value<int>( &fixMulti )->default_value( fixMulti )                     , "Fix MultiPdf level 2." )
    ( "eps"           , po::value<double>( &eps )->default_value( eps )                            , "Convergence criterium." )
    ;

  po::variables_map vm0;

  try {
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm0 );
    po::notify( vm0 );
  }

  catch ( std::exception& ex ) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Use ./a.out --help to print the allowed program options" << endl;
    return -1;
  }

  catch ( ... ) {
    cerr << "Unidentified error parsing program options." << endl;
    return -1;
  }

  // if help, print help
  if ( vm0.count( "help" ) ) {
    cout << "Usage: ./a.out [PROGRAMOPTIONS]\n";
    cout << desc;
    return 0;
  }

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  cout << "Running over workspace: " << inFileName << endl;
  system(("mkdir -vp root-files/" + string(folder) + "_pulls").c_str());
  vector<string> parsed = parseString(poiName, ",");
  vector<string> parsedProfile;
  if (profileName != "none") parsedProfile = parseString(profileName, ",");

  TFile* file = new TFile(inFileName.c_str());

  RooWorkspace* ws = (RooWorkspace*)file->Get(wsName.c_str());
  if (!ws) {
    cout << "Workspace: " << wsName << " doesn't exist!" << endl;
    exit(1);
  }

  if (binnedLikelihood) {
    RooFIter iter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        cout << "Activating binned likelihood for " << arg->GetName() << endl;
      }
    }
  }

  if (fixCache) {
    RooFIter iter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooStarMomentMorph::Class()) {
        ((RooStarMomentMorph*)arg)->fixCache();
        cout << "Fixing cache for " << arg->GetName() << endl;
      }
    }
  }

  if (fixMulti) {
    RooFIter iter = ws->components().fwdIterator();
    RooAbsArg* arg;
    while ((arg = iter.next())) {
      if (arg->IsA() == RooMultiPdf::Class()) {
        arg->setAttribute("NOCacheAndTrack");
        cout << "De-activation level 2 constant term optimization for " << arg->GetName() << endl;
      }
    }
  }

  ModelConfig* mc = (ModelConfig*)(ws->obj(modelConfigName.c_str()));
  if (!mc) {
    cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
    exit(1);
  }

  RooAbsPdf* pdf = (RooAbsPdf*)mc->GetPdf();
  if (!pdf) {
    cout << "ERROR::PDF not found!" << endl;
    exit(1);
  }

  RooAbsData* data = (RooAbsData*)(ws->data(dataName.c_str()));
  if (!data) {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    exit(1);
  }

  RooArgSet* nuis = (RooArgSet*)mc->GetNuisanceParameters();
  if (!nuis) {
    cout << "ERROR::Nuisance parameter set doesn't exist!" << endl;
    exit(1);
  }

  // nuis->add(*ws->var("muF_gg"));
  // nuis->add(*ws->var("muV_gg"));
  // nuis->add(*ws->var("mu_ZZ"));

  RooArgSet* globs = (RooArgSet*)mc->GetGlobalObservables();
  if (!globs) {
    cout << "ERROR::Global observables don't exist!" << endl;
    exit(1);
  }


  ws->loadSnapshot("ucmles");

  vector<RooRealVar*> pois;
  for (size_t i = 0; i < parsed.size(); i++) {
    RooRealVar* poi = (RooRealVar*)ws->var(parsed[i].c_str());
    cout << "Getting POI " << poi->GetName() << endl;
    if (!poi) {
      cout << "POI: " << poi->GetName() << " doesn't exist!" << endl;
      exit(1);
    }

    poi->Print();
    poi->setRange(-1000, 200000);
    poi->setError(0.5);
    poi->setConstant(0);
    pois.push_back(poi);
  }

  for (size_t i = 0; i < parsedProfile.size(); i++) {
    RooRealVar* poi = (RooRealVar*)ws->var(parsedProfile[i].c_str());
    cout << "Getting POI " << poi->GetName() << endl;
    if (!poi) {
      cout << "POI: " << poi->GetName() << " doesn't exist!" << endl;
      exit(1);
    }

    poi->Print();
    poi->setRange(-1000, 200000);
    poi->setError(0.5);
    poi->setConstant(0);
  }

  // set all POIs to consider floating
  for (size_t i = 0; i < pois.size(); i++) {
    pois[i]->setConstant(0);
  }

  // TMP!!!
  // ws->var("mu")->setConstant(1);


  if (variable != "") {
    RooRealVar* nuip = (RooRealVar*)ws->var(variable.c_str());
    nuip->setConstant(0);

    cout << "Computing error for var " << nuip->GetName() << " at " << nuip->getVal() << endl;

    ExtendedMinimizer minimizer("minimizer", pdf, data);
    minimizer.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                       Constrain(*nuis), GlobalObservables(*globs),
                       NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                       ExtendedMinimizer::Scan(RooArgSet(*nuip)), Precision(precision), Save());

    vector<double> pois_hat;
    for (size_t i = 0; i < pois.size(); i++) {
      pois_hat.push_back(pois[i]->getVal());
      cout << pois[i]->GetName() << " " << pois[i]->getVal() << endl;
    }

    // save snapshot
    ws->saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));
    cout << "Made unconditional snapshot" << endl;

    double nuip_hat = nuip->getVal();
    double nuip_errup = nuip->getErrorHi(), nuip_errdown = nuip->getErrorLo();
    double prefitvariation = 1.0;

    // loop over all constraint of pdf to determine
    // constraint type of nuisance parameter

    // Load information needed to determine attributes from ModelConfig
    RooAbsPdf* tmpPdf = (RooAbsPdf*)mc->GetPdf();
    RooArgSet* tmpAllNuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
    RooArgSet* tmpAllObservables = (RooArgSet*)mc->GetObservables();
    RooArgSet* tmpAllGlobalObservables = (RooArgSet*)mc->GetGlobalObservables();

    // Copies, to keep original sets in place after getAllconstraints call
    RooArgSet tmpAllNuisanceParameters2 = *tmpAllNuisanceParameters;
    RooArgSet tmpAllObservables2 = *tmpAllObservables;
    RooArgSet* AllConstraints = tmpPdf->getAllConstraints(tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);

    // Take care of the case where we have a product of constraint terms
    TIterator* ConstraintItrAll = AllConstraints->createIterator();
    RooAbsArg* nextConstraint;
    RooArgSet* tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
    while ((nextConstraint = (RooAbsArg*)ConstraintItrAll->Next())) {
      if (nextConstraint->IsA() == RooProdPdf::Class()) {
        RooArgSet thisComponents;
        FindUniqueProdComponents((RooProdPdf*)nextConstraint, thisComponents);
        tmpAllConstraints->add(thisComponents);
      } else {
        tmpAllConstraints->add(*nextConstraint);
      }
    }

    TIterator* ConstraintItr = tmpAllConstraints->createIterator();
    bool foundConstraint = kFALSE;
    while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next()) && !foundConstraint) {
      if (nextConstraint->dependsOn(*nuip)) {
        foundConstraint = kTRUE;

        // loop over global observables to match nuisance parameter and
        // global observable in case of a constrained nuisnace parameter
        TIterator* GlobsItr = tmpAllGlobalObservables->createIterator();
        RooRealVar* nextGlobalObservable;
        bool foundGlobalObservable = kFALSE;
        while ((nextGlobalObservable = (RooRealVar*)GlobsItr->Next()) && !foundGlobalObservable) {
          if (nextConstraint->dependsOn(*nextGlobalObservable)) {
            foundGlobalObservable = kTRUE;

            // find constraint width in case of a Gaussian
            if (nextConstraint->IsA() == RooGaussian::Class()) {
              double oldSigmaVal = 1.0;
              TIterator* ServerItr = nextConstraint->serverIterator();
              RooRealVar* nextServer;
              bool foundSigma = kFALSE;
              while ((nextServer = (RooRealVar*)ServerItr->Next()) && !foundSigma) {
                if (nextServer != nextGlobalObservable && nextServer != nuip) {
                  oldSigmaVal = nextServer->getVal();
                  foundSigma = kTRUE;
                }
              }

              if (AlmostEqualUlpsAndAbs(oldSigmaVal, 1.0, 0.001, 4)) {
                oldSigmaVal = 1.0;
              }

              if (!foundSigma) {
                cout << "Sigma for pdf " << nextConstraint->GetName() << " not found. Using 1.0." << endl;
              } else {
                cout << "Using " << oldSigmaVal << " for sigma of pdf " << nextConstraint->GetName() << endl;
              }

              prefitvariation = oldSigmaVal;
            }
          }
        }
        delete GlobsItr;
      }
    }
    delete ConstraintItr;

    if (!foundConstraint) {
      cout << "Not a constrained parameter. Stop here." << endl;
      exit(-1);
    }

    // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat+fabs(nuip_errup));
    nuip->setConstant(1);

    // SetAllConstant(*nuis);
    ExtendedMinimizer minimizer1("minimizer", pdf, data);
    minimizer1.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                        Constrain(*nuis), GlobalObservables(*globs),
                        NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                        Precision(precision), Save());

    vector<double> pois_up;
    for (size_t i = 0; i < pois.size(); i++) {
      pois_up.push_back(pois[i]->getVal());
    }

    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat-fabs(nuip_errdown));
    nuip->setConstant(1);

    // SetAllConstant(*nuis);
    ExtendedMinimizer minimizer2("minimizer", mc->GetPdf(), data);
    minimizer2.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                        Constrain(*nuis), GlobalObservables(*globs),
                        NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                        Precision(precision), Save());


    vector<double> pois_down;
    for (size_t i = 0; i < pois.size(); i++) {
      pois_down.push_back(pois[i]->getVal());
    }

    // fix theta at the MLE value +/- prefit uncertainty and minimize again to estimate the change in the POI
    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat+prefitvariation);
    nuip->setConstant(1);

    // SetAllConstant(*nuis);
    ExtendedMinimizer minimizer3("minimizer", mc->GetPdf(), data);
    minimizer3.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                        Constrain(*nuis), GlobalObservables(*globs),
                        NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                        Precision(precision), Save());

    vector<double> pois_nom_up;
    for (size_t i = 0; i < pois.size(); i++) {
      pois_nom_up.push_back(pois[i]->getVal());
    }

    ws->loadSnapshot("tmp_snapshot");
    nuip->setVal(nuip_hat-prefitvariation);
    nuip->setConstant(1);

    // SetAllConstant(*nuis);
    ExtendedMinimizer minimizer4("minimizer", mc->GetPdf(), data);
    minimizer4.minimize(Minimizer(minimizerType.c_str(), minimizerAlgo.c_str()), Strategy(defaultStrategy), ExtendedMinimizer::Eps(eps),
                        Constrain(*nuis), GlobalObservables(*globs),
                        NumCPU(numCPU, 3), Offset(offsetting), Optimize(constOpt),
                        Precision(precision), Save());

    vector<double> pois_nom_down;
    for (size_t i = 0; i < pois.size(); i++) {
      pois_nom_down.push_back(pois[i]->getVal());
    }

    for (size_t i = 0; i < pois.size(); i++) {
      cout << pois[i]->GetName() << " = " << pois_hat[i] << ", " << variable << " = " << nuip_hat << " +" << nuip_errup << " -" << nuip_errdown << endl;
      cout << "Variation of " << pois[i]->GetName() << " = " << pois_up[i] << " (" << pois_nom_up[i] << ") / " << pois_down[i] << " (" << pois_nom_down[i] << ")" << endl;
    }

    // store result in root file
    stringstream fileName;
    fileName << "root-files/" << folder << "_pulls/" << variable << ".root";
    TFile fout(fileName.str().c_str(), "recreate");

    TH1D* h_out = new TH1D(variable.c_str(), variable.c_str(), 3 + 5 * pois.size(), 0, 3 + 5 * pois.size());

    h_out->SetBinContent(1, nuip_hat);
    h_out->SetBinContent(2, fabs(nuip_errup));
    h_out->SetBinContent(3, fabs(nuip_errdown));

    h_out->GetXaxis()->SetBinLabel(1, "nuip_hat");
    h_out->GetXaxis()->SetBinLabel(2, "nuip_up");
    h_out->GetXaxis()->SetBinLabel(3, "nuip_down");

    int bin = 4;
    for (size_t i = 0; i < pois.size(); i++) {
      h_out->SetBinContent(bin, pois_hat[i]);
      h_out->SetBinContent(bin+1, pois_up[i]);
      h_out->SetBinContent(bin+2, pois_down[i]);
      h_out->SetBinContent(bin+3, pois_nom_up[i]);
      h_out->SetBinContent(bin+4, pois_nom_down[i]);

      h_out->GetXaxis()->SetBinLabel(bin, pois[i]->GetName());
      h_out->GetXaxis()->SetBinLabel(bin+1, "poi_up");
      h_out->GetXaxis()->SetBinLabel(bin+2, "poi_down");
      h_out->GetXaxis()->SetBinLabel(bin+3, "poi_nom_up");
      h_out->GetXaxis()->SetBinLabel(bin+4, "poi_nom_down");

      bin += 5;
    }

    fout.Write();
    fout.Close();
  }

  timer_pulls.Stop();
  timer_pulls.Print();

  return 0;
}

// ____________________________________________________________________________|__________
vector<string> parseString(string str, string sep)
{
  vector<string> parsed;
  int pos = 0;
  bool first = true;
  if (str.size() == 0) return parsed;
  if (str.find(sep) == string::npos) {
    parsed.push_back(str);
    return parsed;
  }

  while (true) {
    int newPos = str.find(sep, pos);
    if (str.find(sep, pos) == string::npos) {
      if (!first) parsed.push_back(str.substr(pos, newPos-pos));
      break;
    }

    string sub = str.substr(pos, newPos-pos);
    parsed.push_back(sub);
    pos = newPos+1;
    first = false;
  }

  return parsed;
}

// ____________________________________________________________________________|__________
void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components )
{
  static int counter = 0;
  counter++;

  if (counter > 50) {
    cout << "FindUniqueProdComponents detected infinite loop. Please check." << endl;
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    cout << "FindUniqueProdComponents " << pdfList.at(0)->GetName() << " is fundamental." << endl;
    Components.add(pdfList);
  } else {
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      RooProdPdf* Pdf = (RooProdPdf*)nextArg;
      if (string(Pdf->ClassName()) != "RooProdPdf") {
        cout << "FindUniqueProdComponents " << Pdf->GetName() << " is no RooProdPdf. Adding it." << endl;
        Components.add(*Pdf);
        continue;
      }
      FindUniqueProdComponents(Pdf, Components);
    }
    delete pdfItr;
  }
  counter = 0;
}

// ____________________________________________________________________________|__________
// Compare two floating point numbers, see http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
union MyFloat_t
{
  MyFloat_t(float num = 0.0f) : f(num) {}

  // Portable extraction of components.
  bool Negative() const { return (i >> 31) != 0; }
  int32_t RawMantissa() const { return i & ((1 << 23) - 1); }
  int32_t RawExponent() const { return (i >> 23) & 0xFF; }

  int32_t i;
  float f;
};

bool AlmostEqualUlpsAndAbs( float A, float B, float maxDiff, int maxUlpsDiff )
{
  // Check if the numbers are really close -- needed  when comparing numbers near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff)
    return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
    return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff)
    return true;

  return false;
}
