#ifndef HiggsAnalysis_CombinedLimit_MultiDimFit_h
#define HiggsAnalysis_CombinedLimit_MultiDimFit_h
/** \class MultiDimFit
 *
 * Do a ML fit with multiple POI 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "HiggsAnalysis/CombinedLimit/interface/FitterAlgoBase.h"
#include <RooRealVar.h>
#include <vector>

class MultiDimFit : public FitterAlgoBase {
public:
  MultiDimFit() ;
  virtual const std::string & name() const {
    static const std::string name("MultiDimFit");
    return name;
  }
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;

protected:
  virtual bool runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);

  enum Algo { None, Singles, Cross, Grid, RandomPoints, Contour2D, Stitch2D, FixedPoint, Impact };
  static Algo algo_;

  enum GridType { G1x1, G3x3 };
  static GridType gridType_;

  static std::vector<std::string>  poi_;
  static std::vector<RooRealVar*>  poiVars_;
  static std::vector<float>        poiVals_;
  static RooArgList                poiList_; 
  static unsigned int              nOtherFloatingPoi_; // keep a count of other POIs that we're ignoring, for proper chisquare normalization
  static float                     deltaNLL_;

  static std::string name_;
  std::auto_ptr<TFile> fitOut;

  // options    
  static unsigned int points_, firstPoint_, lastPoint_;
  static bool floatOtherPOIs_;
  static bool squareDistPoiStep_;
  static bool skipInitialFit_;
  static bool fastScan_;
  static bool hasMaxDeltaNLLForProf_;
  static bool loadedSnapshot_,  savingSnapshot_;
  static float maxDeltaNLLForProf_;
  static float autoRange_;
  static bool  startFromPreFit_;
  static bool  alignEdges_;
  static bool  saveFitResult_;
  static std::string fixedPointPOIs_;
  static float centeredRange_;

  static bool robustHesse_;
  static std::string robustHesseLoad_;
  static std::string robustHesseSave_;

  static std::string saveSpecifiedFuncs_;
  static std::string saveSpecifiedNuis_;
  static std::string saveSpecifiedIndex_;
  static std::vector<std::string>  specifiedFuncNames_;
  static std::vector<RooAbsReal*> specifiedFunc_;
  static std::vector<float>        specifiedFuncVals_;
  static RooArgList                specifiedFuncList_;
  static std::vector<std::string>  specifiedCatNames_;
  static std::vector<RooCategory*> specifiedCat_;
  static std::vector<int>        specifiedCatVals_;
  static RooArgList                specifiedCatList_;
  static std::vector<std::string>  specifiedNuis_;
  static std::vector<RooRealVar *> specifiedVars_;
  static std::vector<float>        specifiedVals_;
  static RooArgList                specifiedList_;
  static bool saveInactivePOI_;
  // initialize variables
  void initOnce(RooWorkspace *w, RooStats::ModelConfig *mc_s) ;

  // variables
  void doSingles(RooFitResult &res) ;
  void doGrid(RooWorkspace *w, RooAbsReal &nll) ;
  void doRandomPoints(RooWorkspace *w, RooAbsReal &nll) ;
  void doFixedPoint(RooWorkspace *w, RooAbsReal &nll) ;
  void doContour2D(RooWorkspace *w, RooAbsReal &nll) ;
  void doStitch2D(RooWorkspace *w, RooAbsReal &nll) ;
  void doImpact(RooFitResult &res, RooAbsReal &nll) ;

  // utilities
  /// for each RooRealVar, set a range 'box' from the PL profiling all other parameters
  void doBox(RooAbsReal &nll, double cl, const char *name="box", bool commitPoints=true) ;
  /// save a file with the RooFitResult inside
  void saveResult(RooFitResult &res);
};


#endif
