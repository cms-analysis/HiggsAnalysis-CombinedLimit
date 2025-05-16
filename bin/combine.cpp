#include "../interface/Combine.h"
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <RooRandom.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>
#include "../interface/Significance.h"
#include "../interface/HybridNew.h"
#include "../interface/BayesianFlatPrior.h"
#include "../interface/BayesianToyMC.h"
#include "../interface/MarkovChainMC.h"
#include "../interface/FeldmanCousins.h"
#include "../interface/FitDiagnostics.h"
#include "../interface/AsymptoticLimits.h"
#include "../interface/GoodnessOfFit.h"
#include "../interface/ChannelCompatibilityCheck.h"
#include "../interface/MultiDimFit.h"
#include "../interface/CascadeMinimizer.h"
#include "../interface/ProfilingTools.h"
#include "../interface/GenerateOnly.h"
#include "../interface/CombineLogger.h"
#include <map>

using namespace std;

// Update whenever we have a new Tag
std::string combineTagString = "v10.2.0";
// 

int main(int argc, char **argv) {
  using namespace boost;
  namespace po = boost::program_options;

  string name;
  string datacard, dataset, dataMapName;
  float iMass;
  string whichMethod, whichHintMethod;
  int runToys;
  int    seed;
  string toysFile;

  vector<string> librariesToLoad;
  vector<string> runtimeDefines;
  vector<string> modelPoints;
  vector<string> modelParamNameVector_;
  vector<string> modelParamValVector_;

  Combine combiner;

  // Set name of Log file (ideally would make this similar format to out file)
  CombineLogger::instance().setName("combine_logger.out");

  map<string, LimitAlgo *> methods;
  algo = new Significance(); methods.insert(make_pair(algo->name(), algo));
  algo = new BayesianFlatPrior(); methods.insert(make_pair(algo->name(), algo));
  algo = new BayesianToyMC(); methods.insert(make_pair(algo->name(), algo));
  algo = new MarkovChainMC();  methods.insert(make_pair(algo->name(), algo));
  algo = new HybridNew();  methods.insert(make_pair(algo->name(), algo));
  algo = new FeldmanCousins();  methods.insert(make_pair(algo->name(), algo));
  algo = new FitDiagnostics();  methods.insert(make_pair(algo->name(), algo));
  algo = new AsymptoticLimits();  methods.insert(make_pair(algo->name(), algo));
  algo = new GoodnessOfFit();  methods.insert(make_pair(algo->name(), algo));
  algo = new ChannelCompatibilityCheck();  methods.insert(make_pair(algo->name(), algo));
  algo = new MultiDimFit();  methods.insert(make_pair(algo->name(), algo));
  algo = new GenerateOnly();  methods.insert(make_pair(algo->name(), algo));
  
  CascadeMinimizer::initOptions();

  string methodsDesc("Method to extract upper limit. Supported methods are: ");
  for(map<string, LimitAlgo *>::const_iterator i = methods.begin(); i != methods.end(); ++i) {
    if(i != methods.begin()) methodsDesc += ", ";
    methodsDesc += i->first;
  }
  
  po::options_description desc("Main options");
  desc.add_options()
    ("datacard,d", po::value<string>(&datacard), "Datacard file (can also be specified directly without the -d or --datacard)")
    ("method,M",      po::value<string>(&whichMethod)->default_value("AsymptoticLimits"), methodsDesc.c_str())
    ("verbose,v",  po::value<int>(&verbose)->default_value(0), "Verbosity level (-1 = very quiet; 0 = quiet, 1 = verbose, 2+ = debug)")
    ("help,h", "Produce help message")
    ;
  combiner.statOptions().add_options()
    ("toys,t", po::value<int>(&runToys)->default_value(0), "Number of Toy MC extractions")
    ("seed,s", po::value<int>(&seed)->default_value(123456), "Toy MC random seed")
    ("hintMethod,H",  po::value<string>(&whichHintMethod)->default_value(""), "First run this method to provide a hint on the result")
    ;
  combiner.ioOptions().add_options()
    ("name,n",     po::value<string>(&name)->default_value("Test"), "Name of the job, affects the name of the output tree")
    ("mass,m",     po::value<float>(&iMass)->default_value(120.), "Mass to store in the output tree")
    ("dataset,D",  po::value<string>(&dataset)->default_value("data_obs"), "Name of the dataset for observed limit - use this to replace dataset in workspace for example with a toy dataset. Format as file:workspace:object or file:object")
    ("dataMapName",  po::value<string>(&dataMapName)->default_value("data_obs"), "Name of the dataset for observed limit pattern in the datacard")
    ("toysFile",   po::value<string>(&toysFile)->default_value(""), "Read toy mc or other intermediate results from this file")
    ;
  combiner.miscOptions().add_options()
    ("igpMem", "Setup support for memory profiling using IgProf")
    ("perfCounters", "Dump performance counters at end of job")
    ("LoadLibrary,L", po::value<vector<string> >(&librariesToLoad), "Load library through gSystem->Load(...). Can specify multiple libraries using this option multiple times")
    ("keyword-value",  po::value<vector<string> >(&modelPoints), "Set keyword values with 'WORD=VALUE', will replace $WORD with VALUE in datacards. Filename will also be extended with 'WORDVALUE'. Can specify multiple times")
    ("X-rtd",  po::value<vector<string> >(&runtimeDefines), "Define some constants to be used at runtime (for debugging purposes). The syntax is --X-rtd identifier[=value], where value is an integer and defaults to 1. Can specify multiple times")
    ("X-fpeMask", po::value<int>(), "Set FPE mask: 1=NaN, 2=Div0, 4=Overfl, 8=Underf, 16=Inexact; 7=default")
    ;
  desc.add(combiner.statOptions());
  desc.add(combiner.ioOptions());
  desc.add(CascadeMinimizer::options());
  desc.add(combiner.miscOptions());
  po::positional_options_description p;
  p.add("datacard", -1);
  po::variables_map vm, vm0;
 
  
  // parse the first time, using only common options and allow unregistered options 
  try{
    po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm0);
    po::notify(vm0);
  } catch(std::exception &ex) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Invalid options: " << ex.what() << endl;
    cout << "Use combine --help to get a list of all the allowed options"  << endl;
    return 999;
  } catch(...) {
    cerr << "Unidentified error parsing options." << endl;
    return 1000;
  }

  // if help, print help
  if(vm0.count("help")) {
    cout << "Usage: combine datacard [options]\n";
    cout << desc;
    map<string, LimitAlgo *>::const_iterator i;
    for(i = methods.begin(); i != methods.end(); ++i) {
        if ( (!vm0["method"].defaulted()) && whichMethod.compare(i->second->name())!=0 ) continue;
        cout << i->second->options() << "\n";
    }
    return 0;
  }

  if (verbose>1){
   std::ifstream splashFile; splashFile.open("data/splash.txt");
   while(splashFile.good()) { 
     std::cout << (char)splashFile.get() ; 
   }
   std::cout << std::endl;
   splashFile.close(); 
  } else {
   std::cout << " <<< Combine >>> " << std::endl;
   // UPDATE THIS TO THE LATEST TAG WHENEVER RELEASED 
   std::cout << Form(" <<< %s >>>",combineTagString.c_str() ) << std::endl;
  }

  // now search for algo, and add option
  map<string, LimitAlgo *>::const_iterator it_algo = methods.find(whichMethod);
  if (it_algo == methods.end()) {
    cerr << "Unsupported method: " << whichMethod << endl;
    cout << "Use combine --help to get a list of all the allowed methods and options"  << endl;
    return 1003;
  } 
  desc.add(it_algo->second->options());  

  // parse the first time, now include options of the algo but not unregistered ones
  try{
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);
  } catch(std::exception &ex) {
    cerr << "Invalid options: " << ex.what() << endl;
    cout << "Invalid options: " << ex.what() << endl;
    cout << "Use combine --help to get a list of all the allowed options"  << endl;
    return 999;
  } catch(...) {
    cerr << "Unidentified error parsing options." << endl;
    return 1000;
  }

  if(datacard == "") {
    cerr << "Missing datacard file" << endl;
    cout << "Usage: combine [options]\n";
    cout << "Use combine --help to get a list of all the allowed methods and options"  << endl;
    return 1002;
  }

  if (seed == -1) {
    if (verbose > 0) std::cout << ">>> Using OpenSSL to get a truly random seed " << std::endl;
    FILE *rpipe = popen("openssl rand 8", "r");
    if (rpipe == 0) { std::cout << "Error when running 'openssl rand 8'" << std::endl; return 2101; }
    if (fread(&seed, sizeof(int), 1, rpipe) != 1) {
        std::cout << "Error when reading from 'openssl rand 8'" << std::endl; return 2102; 
    }
    std::cout << ">>> Used OpenSSL to get a truly random seed " << seed << std::endl;
  } else {
    std::cout << ">>> Random number generator seed is " << seed << std::endl;
  }
  RooRandom::randomGenerator()->SetSeed(seed); 

  TString massName = TString::Format("mH%g.", iMass);
  TString toyName  = "";  if (runToys > 0 || seed != 123456 || vm.count("saveToys")) toyName  = TString::Format("%d.", seed);
  if (vm.count("expectedFromGrid") && !vm["expectedFromGrid"].defaulted()) toyName += TString::Format("quant%.3f.", vm["expectedFromGrid"].as<float>());
  if (vm.count("expected")         && !vm["expected"].defaulted())         toyName += TString::Format("quant%.3f.", vm["expected"].as<float>());

  if (vm.count("keyword-value") ) {
    for (vector<string>::const_iterator rmp = modelPoints.begin(), endrmp = modelPoints.end(); rmp != endrmp; ++rmp) {
      std::string::size_type idx = rmp->find('=');
      if (idx == std::string::npos) {
     	cerr << "No value found for keyword :\n\t" << *rmp << " use --keyword-value WORD=VALUE " << std::endl;
      } else {
        std::string name   = rmp->substr(0, idx);
        std::string svalue = rmp->substr(idx+1);
        if (verbose > 0) std::cout << "Setting keyword " << name << " to " << svalue << std::endl;
        modelParamNameVector_.push_back(name);
        modelParamValVector_.push_back(svalue);
        massName += TString(name.c_str())+TString(svalue.c_str())+".";
     }
    }
  }

  //store derived values in variable map as defaults for new keys
  po::options_description derived_desc("Derived values");
  derived_desc.add_options()
    ("massName", po::value<string>()->default_value(string(massName.View())), "massName for output file names")
    ("toyName", po::value<string>()->default_value(string(toyName.View())), "toyName for output file names")
  ;
  stringstream dummy;
  po::store(po::parse_config_file(dummy, derived_desc), vm);

  try {
    combiner.applyOptions(vm);
    CascadeMinimizer::applyOptions(vm);
  } catch (std::exception &ex) {
    cerr << "Error when configuring combine:\n\t" << ex.what() << std::endl;
    return 2001;
  }

  algo = it_algo->second;
  try {
      algo->applyOptions(vm);
  } catch (std::exception &ex) {
      cerr << "Error when configuring the algorithm " << whichMethod << ":\n\t" << ex.what() << std::endl;
      return 2002;
  }
  cout << ">>> Method used is " << whichMethod << endl;

  if (!whichHintMethod.empty()) {
      map<string, LimitAlgo *>::const_iterator it_hint = methods.find(whichHintMethod);
      if (it_hint == methods.end()) {
          cerr << "Unsupported hint method: " << whichHintMethod << endl;
          cout << "Usage: combine [options]\n";
          cout << "Use combine --help to get a list of all the allowed methods and options"  << endl;
          return 1003;
      } 
      hintAlgo = it_hint->second;
      hintAlgo->applyDefaultOptions();
      cout << ">>> Method used to hint where the upper limit is " << whichHintMethod << endl;
  }
  
  TString fileName = "higgsCombine" + name + "."+whichMethod+"."+massName+toyName+"root";

  TFile *test = new TFile(fileName, "RECREATE"); outputFile = test;
  TTree *t = new TTree("limit", "limit");
  int syst, iToy, iSeed, iChannel; 
  double mass, limit, limitErr; 
  t->Branch("limit",&limit,"limit/D");
  t->Branch("limitErr",&limitErr,"limitErr/D");
  t->Branch("mh",   &mass, "mh/D");
  t->Branch("syst", &syst, "syst/I");
  t->Branch("iToy", &iToy, "iToy/I");
  t->Branch("iSeed", &iSeed, "iSeed/I");
  t->Branch("iChannel", &iChannel, "iChannel/I");
  t->Branch("t_cpu",   &t_cpu_,  "t_cpu/F");
  t->Branch("t_real",  &t_real_, "t_real/F");
  t->Branch("quantileExpected",  &g_quantileExpected_, "quantileExpected/F");
  for (unsigned int mpi=0;mpi<modelParamNameVector_.size();++mpi){
	std::string name = modelParamNameVector_[mpi];
  	t->Branch(Form("%s",name.c_str()),  &modelParamValVector_[mpi]);
  }
  
  writeToysHere = test->mkdir("toys","toys"); 
  if (toysFile != "") readToysFromHere = TFile::Open(toysFile.c_str());
  
  syst = withSystematics;
  mass = iMass;
  iSeed = seed;
  iChannel = 0;

  // if you have libraries, it's time to load them now
  for (vector<string>::const_iterator lib = librariesToLoad.begin(), endlib = librariesToLoad.end(); lib != endlib; ++lib) {
    gSystem->Load(lib->c_str());
  }

  if (vm.count("igpMem")) setupIgProfDumpHook();

  if (vm.count("X-fpeMask")) gSystem->SetFPEMask(vm["X-fpeMask"].as<int>());

  // CMSDAS Defaults (you can turn off with --X-rtd <name>=0
  runtimedef::set("OPTIMIZE_BOUNDS", 1);
  runtimedef::set("ADDNLL_RECURSIVE", 1);
  runtimedef::set("ADDNLL_GAUSSNLL", 1);
  runtimedef::set("ADDNLL_HISTNLL", 1);
  runtimedef::set("ADDNLL_CBNLL", 1);
  runtimedef::set("TMCSO_AdaptivePseudoAsimov", 1);
  // Optimization for bare RooFit likelihoods (--optimizeSimPdf=0)
  runtimedef::set("MINIMIZER_optimizeConst", 2); 
  runtimedef::set("MINIMIZER_rooFitOffset", 1); 
  // Optimization for ATLAS HistFactory likelihoods
  runtimedef::set("ADDNLL_ROOREALSUM_FACTOR",1);
  runtimedef::set("ADDNLL_ROOREALSUM_NONORM",1);
  runtimedef::set("ADDNLL_ROOREALSUM_BASICINT",1);
  runtimedef::set("ADDNLL_ROOREALSUM_KEEPZEROS",1);
  runtimedef::set("ADDNLL_PRODNLL",1);
  runtimedef::set("ADDNLL_HFNLL",1);
  runtimedef::set("ADDNLL_HISTFUNCNLL",1);
  runtimedef::set("ADDNLL_ROOREALSUM_CHEAPPROD",1);
 


  for (vector<string>::const_iterator rtdp = runtimeDefines.begin(), endrtdp = runtimeDefines.end(); rtdp != endrtdp; ++rtdp) {
    std::string::size_type idx = rtdp->find('=');
    if (idx == std::string::npos) {
        runtimedef::set(*rtdp, 1); 
        if (verbose > 0) std::cout << "Turning on runtime-define " << *rtdp << std::endl;
    } else {
        std::string name  = rtdp->substr(0, idx);
        std::string svalue = rtdp->substr(idx+1);
        int ivalue = atoi( svalue.c_str() );
        if (verbose > 0) std::cout << "Setting runtime-define " << name << " to " << ivalue << std::endl;
        runtimedef::set(name, ivalue);
    }
  }

  try {
     combiner.run(datacard, dataset, limit, limitErr, iToy, t, runToys);
     if (verbose>0) CombineLogger::instance().printLog(); 
  } catch (std::exception &ex) {
     cerr << "Error when running the combination:\n\t" << ex.what() << std::endl;
     test->Close();
     return 3001;
  }
  
  test->WriteTObject(t);
  test->Close();

  for(map<string, LimitAlgo *>::const_iterator i = methods.begin(); i != methods.end(); ++i)
    delete i->second;

  if (vm.count("perfCounters")) PerfCounter::printAll();
}


