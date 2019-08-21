#ifndef HiggsAnalysis_CombinedLimit_Combine_h
#define HiggsAnalysis_CombinedLimit_Combine_h
#include <TString.h>
#include <TFile.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

class TDirectory;
class TTree;
class LimitAlgo;
class RooWorkspace;
class RooAbsData;
namespace RooStats { class ModelConfig; }

extern Float_t t_cpu_, t_real_, g_quantileExpected_; 
extern bool g_fillTree_; 
//RooWorkspace *writeToysHere = 0;
extern TDirectory *outputFile;
extern TDirectory *writeToysHere;
extern TDirectory *readToysFromHere;
extern LimitAlgo * algo, * hintAlgo ;
extern int verbose;
extern bool withSystematics;
extern bool doSignificance_, lowerLimit_;
extern float cl;
extern bool bypassFrequentistFit_;
extern  std::string setPhysicsModelParameterExpression_;
extern  std::string setPhysicsModelParameterRangeExpression_;
extern  std::string defineBackgroundOnlyModelParameterExpression_;

namespace { 
    struct ToCleanUp {
        TFile *tfile; std::string file, path;
        ToCleanUp() : tfile(0), file(""), path("") {}
        ~ToCleanUp() {
              if (tfile) { tfile->Close(); delete tfile; }
              if (!file.empty()) {  
                 unlink(file.c_str());  // FIXME, we should check that the file deleted safely but currently when running HybridNew, we get a status of -1 even though the file is in fact removed?!
		 //if (unlink(file.c_str()) == -1) std::cerr << "Failed to delete temporary file " << file << ": " << strerror(errno) << std::endl;
              }
              if (!path.empty()) {  boost::filesystem::remove_all(path); }
        }
    };
}

class Combine {
public:
  Combine() ;
  
  boost::program_options::options_description & statOptions() { return statOptions_; }    
  boost::program_options::options_description & ioOptions() { return ioOptions_; }    
  boost::program_options::options_description & miscOptions() { return miscOptions_; }    
  void applyOptions(const boost::program_options::variables_map &vm) ;
  
  void run(TString hlfFile, const std::string &dataset, double &limit, double &limitErr, int &iToy, TTree *tree, int nToys);
 
  /// Stop combine from fillint the tree (some algos need control)
  static void toggleGlobalFillTree(bool flag=false);

  /// Save a point into the output tree. Usually if expected = false, quantile should be set to -1 (except e.g. for saveGrid option of HybridNew)
  static void commitPoint(bool expected, float quantile);

  /// Add a branch to the output tree (for advanced use or debugging only)
  static void addBranch(const char *name, void *address, const char *leaflist) ;
private:
  bool mklimit(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr) ;
 
  void addDiscreteNuisances(RooWorkspace *);
  void addNuisances(const RooArgSet *);
  void addFloatingParameters(const RooArgSet &);
  void addPOI(const RooArgSet *);

  boost::program_options::options_description statOptions_, ioOptions_, miscOptions_;

  // statistics-related variables
  bool unbinned_, generateBinnedWorkaround_, newGen_, guessGenMode_; 
  std::string genAsBinned_, genAsUnbinned_;
  float rMin_, rMax_;
  std::string prior_;
  bool hintUsesStatOnly_;
  bool toysNoSystematics_;
  bool toysFrequentist_;
  float expectSignal_;
  bool expectSignalSet_;  // keep track of whether or not expectSignal was defaulted
  float expectSignalMass_;
  std::string redefineSignalPOIs_;
  std::string freezeNuisances_;
  std::string floatNuisances_;
  std::string freezeNuisanceGroups_;
  std::string freezeWithAttributes_;

  // input-output related variables
  bool saveWorkspace_;
  std::string workspaceName_;
  std::string snapshotName_;
  std::string modelConfigName_, modelConfigNameB_;
  bool overrideSnapshotMass_;
  bool validateModel_;
  bool saveToys_;
  double mass_;

  // implementation-related variables
  bool compiledExpr_;
  bool makeTempDir_;
  bool rebuildSimPdf_;
  bool optSimPdf_;
  bool noMCbonly_;
  bool noDefaultPrior_;
  bool makeToyGenSnapshot_;
  bool floatAllNuisances_;
  bool freezeAllGlobalObs_;
  std::vector<std::string> librariesToLoad_;
  std::vector<std::string> modelPoints_;
  
  static TTree *tree_;

  static std::vector<std::pair<RooAbsReal*,float> > trackedParametersMap_;
  static std::string  trackParametersNameString_;
  static std::string  textToWorkspaceString_;

};

#endif
