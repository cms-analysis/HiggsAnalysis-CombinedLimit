#ifndef HiggsAnalysis_CombinedLimit_FitDiagnostics_h
#define HiggsAnalysis_CombinedLimit_FitDiagnostics_h
/** \class FitDiagnostics
 *
 * Do a ML fit of the data with background and signal+background hypothesis and print out diagnostics plots 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "HiggsAnalysis/CombinedLimit/interface/FitterAlgoBase.h"
#include <TTree.h>
#include <RooArgList.h>
#include <RooFitResult.h>
#include <boost/utility.hpp>
#include <map>

class FitDiagnostics : public FitterAlgoBase {
public:
  FitDiagnostics() ;
  virtual const std::string & name() const {
    static const std::string name("FitDiagnostics");
    return name;
  }
  ~FitDiagnostics();
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;
  virtual void setToyNumber(const int) ;
  virtual void setNToys(const int);

protected:
  virtual bool runSpecific(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);

  static std::string name_;

  static std::string minos_;

  static bool justFit_,  skipBOnlyFit_, noErrors_;
  static std::string out_; 
  static bool        makePlots_;
  static float       rebinFactor_;
  static int         numToysForShapes_;
  static std::string signalPdfNames_, backgroundPdfNames_;
  static std::string filterString_;
  static bool        saveNormalizations_;
  static bool        savePredictionsPerToy_;
  static bool        oldNormNames_;
  static bool        saveShapes_;
  static bool        saveOverallShapes_;
  static bool        saveWithUncertainties_;
  static bool	     saveWorkspace_;
  static bool        reuseParams_;
  static bool        customStartingPoint_;
  static bool       robustHesse_;
  static bool        saveWithUncertsRequested_;
  static bool        ignoreCovWarning_;
  int currentToy_, nToys;
  int overallBins_,overallNorms_,overallNuis_,overallCons_;
  int fitStatus_, numbadnll_;
  double mu_, muErr_, muLoErr_, muHiErr_, nll_nll0_, nll_bonly_, nll_sb_;
  std::auto_ptr<TFile> fitOut;
  double* globalObservables_;
  double* nuisanceParameters_;
  double* processNormalizations_;
  double* processNormalizationsShapes_;

  TTree *t_fit_b_, *t_fit_sb_, *t_prefit_;
   
  void getNormalizationsSimple(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out);
  void createFitResultTrees(const RooStats::ModelConfig &,bool);
  void resetFitResultTrees(bool);
  void setFitResultTrees(const RooArgSet *, double *);
  void setNormsFitResultTrees(const RooArgSet *, double *);

  struct ShapeAndNorm {
    bool        signal;
    std::string process;
    std::string channel;
    RooArgList  obs;
    const RooAbsReal *norm;
    const RooAbsReal  *pdf;
    bool isfunc;
  };
  void setShapesFitResultTrees(std::map<std::string,ShapeAndNorm> &snm, double *);
  void getShapesAndNorms(RooAbsPdf *pdf, const RooArgSet &obs, std::map<std::string, ShapeAndNorm> &shapesAndNorms, const std::string &channel);

  class NuisanceSampler { 
    public:
        virtual ~NuisanceSampler() {}
        virtual void  generate(int ntoys) = 0;
        virtual const RooAbsCollection & get(int itoy) = 0;
        virtual const RooAbsCollection & centralValues() = 0;
  };
  void getNormalizations(RooAbsPdf *pdf, const RooArgSet &obs, RooArgSet &out, NuisanceSampler &sampler, TDirectory *fOut, const std::string &postfix,RooAbsData &data);

  class CovarianceReSampler : public NuisanceSampler {
    public:
        CovarianceReSampler(RooFitResult *res) : res_(res) {}
        virtual void  generate(int ntoys) {}
        virtual const RooAbsCollection & get(int) { return res_->randomizePars(); }
        virtual const RooAbsCollection & centralValues() { return res_->floatParsFinal(); }
    protected:
        RooFitResult *res_;
  };
  class ToySampler : public NuisanceSampler, boost::noncopyable {
    public:
        ToySampler(RooAbsPdf *pdf, const RooArgSet *nuisances) ;    
        virtual ~ToySampler() ;
        virtual void  generate(int ntoys);
        virtual const RooAbsCollection & get(int itoy);
        virtual const RooAbsCollection & centralValues();
    private:
        RooAbsPdf  *pdf_;
        RooAbsData *data_;
        RooArgSet  snapshot_; 
  };
};


#endif
