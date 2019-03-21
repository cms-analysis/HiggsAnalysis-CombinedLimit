#ifndef CombineTools_CombineHarvester_h
#define CombineTools_CombineHarvester_h
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <list>
#include <cmath>
#include <set>
#include <functional>
#include "boost/range/algorithm_ext/erase.hpp"
#include "TFile.h"
#include "TH1.h"
#include "RooWorkspace.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/HistMapping.h"


namespace ch {

// Define some useful CombineHarvester-specific typedefs
typedef std::vector<std::pair<int, std::string>> Categories;

class CombineHarvester {
 public:
  /**
   * \name Constructors, destructors and copying
   *
   * \details Distinct shallow and deep copying methods are provided. A deep
   * copy creates a completely independent CombineHarvester instace: all
   * Observation, Process, Systematic and Parameter objects are cloned, as
   * well as any histogram objects they may contain and any attached
   * RooWorkspaces. In a shallow copy these objects are shared with the
   * original CombineHarvester instance.
   */
  /**@{*/
  CombineHarvester();
  ~CombineHarvester();

  /**
   * Copy constructor (makes a shallow copy)
   *
   * When copying a CombineHarvester instance it is important to remember that
   * the stored Observation, Process, Systematic and Parameter objects
   * themselves are not duplicated, rather their pointers are simply copied.
   * This is called making a *shallow-copy*, and means modifying the contents
   * of one of these objects in the copied CombineHarvester also modifies it
   * for the original. However, filtering methods only affect the instance
   * they are called from. For example, if all the signal Process entries are
   * filtered from a copied CombineHarvester instance they will remain in the
   * original instance.
   */
  CombineHarvester(CombineHarvester const& other);
  CombineHarvester(CombineHarvester&& other);
  CombineHarvester& operator=(CombineHarvester other);

  /**
   * Set a named flag value
   *
   * Creates the flag if it doesn't already
   */
  void SetFlag(std::string const& flag, bool const& value);
  bool GetFlag(std::string const& flag) const;

  /**
   * Creates and returns a shallow copy of the CombineHarvester instance
   */
  CombineHarvester cp();

  /**
   * Creates and retunrs a deep copy of the CombineHarvester instance
   *
   * Unlike the shallow copy, a deep copy will duplicate every internal
   * object, including any attached RooWorkspaces. This makes it completely
   * independent of the original instance.
   */
  CombineHarvester deep();
  /**@}*/

  /**
   * \name Logging and Printing
   */
  /**@{*/
  CombineHarvester& PrintAll();
  CombineHarvester& PrintObs();
  CombineHarvester& PrintProcs();
  CombineHarvester& PrintSysts();
  CombineHarvester& PrintParams();
  inline void SetVerbosity(unsigned verbosity) { verbosity_ = verbosity; }
  inline unsigned Verbosity() { return verbosity_; }
  /**@}*/

  /**
   * \name Datacards
   *
   * \brief Methods for the reading and writing of datacards.
   *
   * \details As well as reading or writing the plain-text datacard files
   * these methods also handle the automatic loading and saving of any TH1
   * objects or RooWorkspaces.
   */
  /**@{*/
  int ParseDatacard(std::string const& filename,
      std::string const& analysis,
      std::string const& era,
      std::string const& channel,
      int bin_id,
      std::string const& mass);
  int ParseDatacard(std::string const& filename,
      std::string parse_rule = "");

  void WriteDatacard(std::string const& name, std::string const& root_file);
  void WriteDatacard(std::string const& name, TFile & root_file);
  void WriteDatacard(std::string const& name);
  /**@}*/

  /**
   * \name Filters
   * \anchor CH-Filters
   * \brief A collection of methods to filter the Observation, Process and
   * Systematic objects
   *
   * \details Each method here will remove objects from the Observation,
   * Process and Systematic collections where the value of the relevant
   * property is not equal to one of those specified in the list vec. Each
   * methods returns a reference to the class instance, meaning it is simple
   * to chain multiple filters. Note that filters only act on a set of objects
   * where that property is defined. For example, the method systs_name will
   * only filter the Systematic collection, and process will only filter the
   * Process and Systematic collections.
   *
   * \param vec A vector containing the possible values of the given property.
   */
  /**@{*/
  CombineHarvester& bin(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& bin_id(std::vector<int> const& vec, bool cond = true);
  CombineHarvester& process(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& analysis(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& era(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& channel(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& mass(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& attr(std::vector<std::string> const& vec,std::string attr_label, bool cond = true);
  CombineHarvester& syst_name(std::vector<std::string> const& vec, bool cond = true);
  CombineHarvester& syst_type(std::vector<std::string> const& vec, bool cond = true);

  CombineHarvester& process_rgx(std::vector<std::string> const& vec, bool cond = true);

  CombineHarvester& signals();
  CombineHarvester& backgrounds();
  CombineHarvester& histograms();
  CombineHarvester& pdfs();
  CombineHarvester& data();

  template<typename Function>
  CombineHarvester& FilterAll(Function func);
  template<typename Function>
  CombineHarvester& FilterObs(Function func);
  template<typename Function>
  CombineHarvester& FilterProcs(Function func);
  template<typename Function>
  CombineHarvester& FilterSysts(Function func);
  /**@}*/


  /**
   * \name Set producers
   * \anchor CH-Set-Producers
   * \brief Methods that extract sets of properties from the Observation,
   * Process and Systematic entries
   */
  /**@{*/
  // Set generation
  std::set<std::string> bin_set();
  std::set<int> bin_id_set();
  std::set<std::string> process_set();
  std::set<std::string> analysis_set();
  std::set<std::string> era_set();
  std::set<std::string> channel_set();
  std::set<std::string> mass_set();
  std::set<std::string> syst_name_set();
  std::set<std::string> syst_type_set();

  /**
   * Fill an std::set with the return values from an arbitrary function
   *
   * This method will loop through all ch::Observation, ch::Process and
   * ch::Systematic entries and call the user-supplied function `func`. The
   * return value is then inserted into the set. 
   *
   * @tparam T A function (or other callable) that must have a single
   * `ch::Object const*` argument.
   * @tparam R The return type of the function, which is deduced by using
   * `std::result_of`, then `std::decay`. The latter is needed to handle
   * functions that return by reference, i.e. turning a type R& into a type R.
   */
  template <typename T,
            typename R = typename std::decay<
                typename std::result_of<T(Object const*)>::type>::type>
  std::set<R> SetFromAll(T func);

  /**
   * Fill an std::set using only the Observation entries
   *
   * \sa SetFromAll
   */
  template <typename T,
            typename R = typename std::decay<
                typename std::result_of<T(Observation const*)>::type>::type>
  std::set<R> SetFromObs(T func);

  /**
   * Fill an std::set using only the Process entries
   *
   * \sa SetFromAll
   */
  template <typename T,
            typename R = typename std::decay<
                typename std::result_of<T(Process const*)>::type>::type>
  std::set<R> SetFromProcs(T func);

  /**
   * Fill an std::set using only the Systematic entries
   *
   * \sa SetFromAll
   */
  template <typename T,
            typename R = typename std::decay<
                typename std::result_of<T(Systematic const*)>::type>::type>
  std::set<R> SetFromSysts(T func);
  /**@}*/

  /**
   * \name Modification
   *
   * \brief Methods to modify existing objects.
   *
   * \details See the documentation of each method for details
   */
  /**@{*/
  ch::Parameter const* GetParameter(std::string const& name) const;
  ch::Parameter* GetParameter(std::string const& name);

  void UpdateParameters(std::vector<ch::Parameter> const& params);

  /**
   * Update the parameters to the post-fit values in a RooFitResult
   *
   * \deprecated For consistency and compatibility with the python interface,
   * please use the UpdateParameters(RooFitResult const&) version of this
   * method instead - this method will be removed in an upcoming release
   */
  void UpdateParameters(RooFitResult const* fit);
  void UpdateParameters(RooFitResult const& fit);

  std::vector<ch::Parameter> GetParameters() const;
  void RenameParameter(std::string const& oldname, std::string const& newname);

  template<typename Function>
  void ForEachObj(Function func);

  template<typename Function>
  void ForEachProc(Function func);

  template<typename Function>
  void ForEachObs(Function func);

  template<typename Function>
  void ForEachSyst(Function func);

  void VariableRebin(std::vector<double> bins);
  void ZeroBins(double min, double max);
  void SetPdfBins(unsigned nbins);

  /**
   * Add parameters to a given group
   *
   * A parameter will be added to the given group if its name passes a regex
   * match with any of the input patterns. A parameter can belong to any
   * number of groups.
   *
   * @param name The group name to set
   * @param patterns A vector of regex pattern strings
   */
  void SetGroup(std::string const& name, std::vector<std::string> const& patterns);

  /**
   * Remove parameters to a given group
   *
   * A parameter will be removed from the given group if its name passes a
   * regex match with any of the input patterns.
   *
   * @param name The group name to remove
   * @param patterns A vector of regex pattern strings
   */
  void RemoveGroup(std::string const& name, std::vector<std::string> const& patterns);

  /**
   * Rename a nuisance parameter group
   *
   * @param oldname The current name
   * @param newname The new name
   */
  void RenameGroup(std::string const& oldname, std::string const& newname);

  /**
   * Add a line of text at the end of all datacards
   *
   * @param line Line of text to add
   */
  void AddDatacardLineAtEnd(std::string const& line);

  /**
   * Clear all added datacard lines
   */
  void ClearDatacardLinesAtEnd();
  /**@}*/

  /**
   * \name Rate, shape and uncertainty evaluation
   *
   * \brief Methods to calculate total yields and shapes
   *
   * \details All of these methods are greedy, meaning they will sum over all
   * available objects and evaluate the effect of all uncertainties. They
   * should be used at the end of a chain of filter methods to give the
   * desired yield, shape or uncertainty.
   */
  /**@{*/
  double GetRate();
  double GetObservedRate();
  double GetUncertainty();

  /**
   * Sum all Process yields and evaluate uncertainty by sampling from the fit
   * convariance matrix
   *
   * \deprecated For consistency and compatibility with the python interface,
   * please use the GetUncertainty(RooFitResult const&, unsigned) version of
   * this method instead - this method will be removed in an upcoming release
   */
  double GetUncertainty(RooFitResult const* fit, unsigned n_samples);
  double GetUncertainty(RooFitResult const& fit, unsigned n_samples);
  TH1F GetShape();
  TH1F GetShapeWithUncertainty();

  /**
   * Sum all Process shapes and evaluate bin-wise uncertainty by sampling from
   * the fit convariance matrix
   *
   * \deprecated For consistency and compatibility with the python interface,
   * please use the GetShapeWithUncertainty(RooFitResult const&, unsigned)
   * version of this method instead - this method will be removed in an
   * upcoming release
   */
  TH1F GetShapeWithUncertainty(RooFitResult const* fit, unsigned n_samples);
  TH1F GetShapeWithUncertainty(RooFitResult const& fit, unsigned n_samples);
  TH1F GetObservedShape();

  TH2F GetRateCovariance(RooFitResult const& fit, unsigned n_samples);
  TH2F GetRateCorrelation(RooFitResult const& fit, unsigned n_samples);
  /**@}*/

  /**
   * \name Creating new entries
   *
   * \brief Methods to create new Observation, Process, Systematic and
   * Parameter objects
   *
   * \details The general order in which these are run is given below.
   *
   * 1. The AddObservations and AddProcesses methods build the desired
   *    Observation and ch::Process entries. The rate entries will default to
   *    zero at this stage, and are determined automatically in step 3.
   * 2. The AddSyst method will iterate through each Process and add a
   *    matching Systematic entry, with the name, type and value specfied in
   *    the method arguments. This method can be used at the end of a chain of
   *    filters to only create Systematic objects for a subset of processes.
   * 3. The method CombineHarvester::ExtractShapes opens a specified ROOT file
   *    and uses the provided pattern rules to load the TH1 objects for all
   *    entries. The TH1 integrals are used to set the event rates at this
   *    stage.
   * 4. Optionally, methods such as AddBinByBin may be used to prepare
   *      additional entries.
   */
  /**@{*/
  void AddObservations(std::vector<std::string> mass,
                       std::vector<std::string> analysis,
                       std::vector<std::string> era,
                       std::vector<std::string> channel,
                       ch::Categories bin);

  void AddProcesses(std::vector<std::string> mass,
                    std::vector<std::string> analysis,
                    std::vector<std::string> era,
                    std::vector<std::string> channel,
                    std::vector<std::string> procs,
                    ch::Categories bin, bool signal);

  void AddSystFromProc(Process const& proc, std::string const& name,
                       std::string const& type, bool asymm, double val_u,
                       double val_d, std::string const& formula,
                       std::string const& args);

  template <class Map>
  void AddSyst(CombineHarvester & target, std::string const& name,
               std::string const& type, Map const& valmap);

  void ExtractShapes(std::string const& file, std::string const& rule,
                     std::string const& syst_rule);
  void ExtractPdfs(CombineHarvester& target, std::string const& ws_name,
                   std::string const& rule, std::string norm_rule = "");
  void ExtractData(std::string const& ws_name, std::string const& rule);

  void AddWorkspace(RooWorkspace const& ws, bool can_rename = false);

  void InsertObservation(ch::Observation const& obs);
  void InsertProcess(ch::Process const& proc);
  void InsertSystematic(ch::Systematic const& sys);

  /**
   * Rename a systematic from 'old_name' to 'new_name' and add a parameter
   * 'new_name' to CH instance 'target' if that parameter doesn't exist yet. 
   * Usage similar to AddSyst()
   */
  void RenameSystematic(CombineHarvester& target, std::string const& old_name, std::string const& new_name);

  void CreateParameterIfEmpty(std::string const& name);

  /**
   * Create bin-by-bin uncertainties
   *
   * \deprecated Please use the standalone ch::BinByBinFactory tool, defined
   * in CombineTools/interface/BinByBin.h
   */
  void AddBinByBin(double threshold, bool fixed_norm, CombineHarvester* other);


  /**
   * Create bin-by-bin uncertainties
   *
   * \deprecated Please use the standalone ch::BinByBinFactory tool, defined
   * in CombineTools/interface/BinByBin.h
   */
  void AddBinByBin(double threshold, bool fixed_norm, CombineHarvester & other);


  /**
   * Merge bin errors within a bin
   *
   * \deprecated Please use the standalone ch::BinByBinFactory tool, defined
   * in CombineTools/interface/BinByBin.h
   */
  void MergeBinErrors(double bbb_threshold, double merge_threshold);
  /**@}*/

  void SetAutoMCStats(CombineHarvester &target, double thresh, bool sig=false, int mode=1);
  void RenameAutoMCStatsBin(std::string const& oldname, std::string const& newname);
  std::set<std::string> GetAutoMCStatsBins() const;

  void AddExtArgValue(std::string const& name, double const& value);
 private:
  friend void swap(CombineHarvester& first, CombineHarvester& second);

  // ---------------------------------------------------------------
  // Main data members
  // ---------------------------------------------------------------
  std::vector<std::shared_ptr<Observation>> obs_;
  std::vector<std::shared_ptr<Process>> procs_;
  std::vector<std::shared_ptr<Systematic>> systs_;
  std::map<std::string, std::shared_ptr<Parameter>> params_;
  std::map<std::string, std::shared_ptr<RooWorkspace>> wspaces_;

  std::unordered_map<std::string, bool> flags_;

  struct AutoMCStatsSettings {
    double event_threshold;
    bool include_signal;
    int hist_mode;

    AutoMCStatsSettings(double thresh, bool sig=false, int mode=1) {
      event_threshold = thresh;
      include_signal = sig;
      hist_mode = mode;
    }

    AutoMCStatsSettings() : AutoMCStatsSettings(0.) {}
  };

  std::map<std::string, AutoMCStatsSettings> auto_stats_settings_;
  std::vector<std::string> post_lines_;

  // ---------------------------------------------------------------
  // typedefs
  // ---------------------------------------------------------------
  typedef std::pair<std::string, std::string> StrPair;
  typedef std::vector<StrPair> StrPairVec;
  typedef std::vector<std::string> StrVec;


  // ---------------------------------------------------------------
  // Logging
  // ---------------------------------------------------------------
  unsigned verbosity_;
  std::ostream * log_;
  std::ostream& log() const { return *log_; }

  // ---------------------------------------------------------------
  // Private methods for the shape extraction routines
  // --> implementation in src/CombineHarvester.cc
  // ---------------------------------------------------------------
  void LoadShapes(Observation* entry,
                     std::vector<HistMapping> const& mappings);
  void LoadShapes(Process* entry,
                     std::vector<HistMapping> const& mappings);
  void LoadShapes(Systematic* entry,
                     std::vector<HistMapping> const& mappings);

  HistMapping const& ResolveMapping(std::string const& process,
                                    std::string const& bin,
                                    std::vector<HistMapping> const& mappings);

  StrPairVec GenerateShapeMapAttempts(std::string process,
      std::string category);

  std::shared_ptr<RooWorkspace> SetupWorkspace(RooWorkspace const& ws,
                                    bool can_rename = false);

  void ImportParameters(RooArgSet *vars);

  RooAbsData const* FindMatchingData(Process const* proc);

  ch::Parameter * SetupRateParamVar(std::string const& name, double val, bool is_ext_arg = false);
  void SetupRateParamFunc(std::string const& name, std::string const& formula,
                          std::string const& pars);
  void SetupRateParamWspObj(std::string const& name, std::string const& obj, bool is_ext_arg = false);
  // ---------------------------------------------------------------
  // Private methods for the shape writing routines
  // ---------------------------------------------------------------
  void WriteHistToFile(
      TH1 * hist,
      TFile * file,
      std::vector<HistMapping> const& mappings,
      std::string const& bin,
      std::string const& process,
      std::string const& mass,
      std::string const& nuisance,
      unsigned type);

void FillHistMappings(std::vector<HistMapping> & mappings);

  // ---------------------------------------------------------------
  // Private methods for shape/yield evaluation
  // --> implementation in src/CombineHarvester_Evaluate.cc
  // ---------------------------------------------------------------
  typedef std::vector<std::vector<Systematic const*>> ProcSystMap;
  ProcSystMap GenerateProcSystMap();

  double GetRateInternal(ProcSystMap const& lookup,
    std::string const& single_sys = "");

  TH1F GetShapeInternal(ProcSystMap const& lookup,
    std::string const& single_sys = "");

  inline double smoothStepFunc(double x) const {
    if (std::fabs(x) >= 1.0/*_smoothRegion*/) return x > 0 ? +1 : -1;
    double xnorm = x/1.0;/*_smoothRegion*/
    double xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }

  double logKappaForX(double x, double k_low, double k_high) const;

  void ShapeDiff(double x, TH1F* target, TH1 const* nom, TH1 const* low,
                 TH1 const* high, bool linear);
  void ShapeDiff(double x, TH1F* target, RooDataHist const* nom,
                 RooDataHist const* low, RooDataHist const* high);
};


// ---------------------------------------------------------------
// Template method implementation
// ---------------------------------------------------------------
template <typename T, typename R>
std::set<R> CombineHarvester::SetFromAll(T func) {
  std::set<R> ret;
  for (auto const& item : obs_) ret.insert(func(item.get()));
  for (auto const& item : procs_) ret.insert(func(item.get()));
  for (auto const& item : systs_) ret.insert(func(item.get()));
  return ret;
};

template <typename T, typename R>
std::set<R> CombineHarvester::SetFromObs(T func) {
  std::set<R> ret;
  for (auto const& item : obs_) ret.insert(func(item.get()));
  return ret;
};

template <typename T, typename R>
std::set<R> CombineHarvester::SetFromProcs(T func) {
  std::set<R> ret;
  for (auto const& item : procs_) ret.insert(func(item.get()));
  return ret;
};

template <typename T, typename R>
std::set<R> CombineHarvester::SetFromSysts(T func) {
  std::set<R> ret;
  for (auto const& item : systs_) ret.insert(func(item.get()));
  return ret;
};

template<typename Function>
void CombineHarvester::ForEachObj(Function func) {
  ForEachObs(func);
  ForEachProc(func);
  ForEachSyst(func);
}

template<typename Function>
void CombineHarvester::ForEachProc(Function func) {
  for (auto & item: procs_) func(item.get());
}

template<typename Function>
void CombineHarvester::ForEachObs(Function func) {
  for (auto & item: obs_) func(item.get());
}

template<typename Function>
void CombineHarvester::ForEachSyst(Function func) {
  for (auto & item: systs_) func(item.get());
}

template<typename Function>
CombineHarvester& CombineHarvester::FilterAll(Function func) {
  FilterObs(func);
  FilterProcs(func);
  FilterSysts(func);
  return *this;
}

template<typename Function>
CombineHarvester& CombineHarvester::FilterObs(Function func) {
  boost::remove_erase_if(
      obs_, [&](std::shared_ptr<Observation> ptr) { return func(ptr.get());
  });
  return *this;
}

template<typename Function>
CombineHarvester& CombineHarvester::FilterProcs(Function func) {
  boost::remove_erase_if(
      procs_, [&](std::shared_ptr<Process> ptr) { return func(ptr.get());
  });
  return *this;
}
template<typename Function>
CombineHarvester& CombineHarvester::FilterSysts(Function func) {
  boost::remove_erase_if(
      systs_, [&](std::shared_ptr<Systematic> ptr) { return func(ptr.get());
  });
  return *this;
}

template <class Map>
void CombineHarvester::AddSyst(CombineHarvester& target,
                               std::string const& name, std::string const& type,
                               Map const& valmap) {
  // Keep track of which Process entries get a Systematic assigned and which
  // don't. If verbosity is on we'll print lists of these processes at the end.
  std::vector<ch::Process *> not_added_procs;
  std::vector<ch::Process *> added_procs;
  // Also track which tuples in the map did not get used. Do this by getting the
  // full map here and then removing elements as they are used to create a
  // Systematic.
  auto tuples = valmap.GetTupleSet();
  if (verbosity_ >= 1) {
    log() << (name + ":" + type) << "\n";
  }
  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (!valmap.Contains(procs_[i].get())) {
      not_added_procs.push_back(procs_[i].get());
      continue;
    }
    tuples.erase(valmap.GetTuple(procs_[i].get()));
    added_procs.push_back(procs_[i].get());
    double val_u = valmap.ValU(procs_[i].get());
    double val_d = valmap.ValD(procs_[i].get());
    std::string formula = valmap.Formula(procs_[i].get());
    std::string args = valmap.Args(procs_[i].get());
    target.AddSystFromProc(*(procs_[i]), name, type, valmap.IsAsymm(),
                           val_u, val_d, formula, args);
  }
  if (tuples.size() && verbosity_ >= 1) {
    log() << ">> Map keys that were not used to create a Systematic:\n";
    for (auto s : tuples) {
      log() << ch::Tuple2String(s) << "\n";
    }
  }
  if (verbosity_ >= 2) {
    Process::PrintHeader(log());
    log() << ">> Process entries that did not get a Systematic:\n";
    for (auto p : not_added_procs) {
      log() << *p << "\n";
    }
    log() << ">> Process entries that did get a Systematic:\n";
    for (auto p : added_procs) {
      log() << *p << "\n";
    }
  }
}
}

#endif
