#include "CombineTools/interface/CombineHarvester.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <list>
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Systematic.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/Logging.h"
#include "CombineTools/interface/TFileIO.h"

namespace ch {

CombineHarvester::CombineHarvester() : verbosity_(0), log_(&(std::cout)) {
  flags_["zero-negative-bins-on-import"] = true;
}

CombineHarvester::~CombineHarvester() { }

void swap(CombineHarvester& first, CombineHarvester& second) {
  using std::swap;
  swap(first.obs_, second.obs_);
  swap(first.procs_, second.procs_);
  swap(first.systs_, second.systs_);
  swap(first.params_, second.params_);
  swap(first.wspaces_, second.wspaces_);
  swap(first.verbosity_, second.verbosity_);
  swap(first.flags_, second.flags_);
  swap(first.log_, second.log_);
}

CombineHarvester::CombineHarvester(CombineHarvester const& other)
    : obs_(other.obs_.size(), nullptr),
      procs_(other.procs_.size(), nullptr),
      systs_(other.systs_.size(), nullptr),
      flags_(other.flags_),
      verbosity_(other.verbosity_),
      log_(other.log_) {
  if (verbosity_ >= 3) {
    log() << "[CombineHarvester] Copy-constructor called\n";
  }

  // Build a map of workspace object pointers
  std::map<RooAbsData const*, RooAbsData *> data_map;
  std::map<RooAbsPdf const*, RooAbsPdf *> pdf_map;
  std::map<RooRealVar const*, RooRealVar *> var_map;
  std::map<RooAbsReal const*, RooAbsReal *> real_map;

  for (auto const& it : other.wspaces_) {
    if (it.second) {
      wspaces_.insert({it.first, std::make_shared<RooWorkspace>(*(it.second))});
      std::list<RooAbsData *> o_data = it.second->allData();
      RooArgSet const& o_pdf = it.second->allPdfs();
      RooArgSet const& o_var = it.second->allVars();
      RooArgSet const& o_real = it.second->allFunctions();

      std::list<RooAbsData *> n_data = wspaces_.at(it.first)->allData();
      RooArgSet const& n_pdf = wspaces_.at(it.first)->allPdfs();
      RooArgSet const& n_var = wspaces_.at(it.first)->allVars();
      RooArgSet const& n_real = wspaces_.at(it.first)->allFunctions();
      auto o_iter = o_data.begin();
      auto n_iter = n_data.begin();
      for (; o_iter != o_data.end(); ++o_iter, ++n_iter) {
        data_map[*o_iter] = *n_iter;
      }
      auto o_pdf_it = o_pdf.createIterator();
      auto n_pdf_it = n_pdf.createIterator();
      do {
        RooAbsPdf *o_pdf_ptr = dynamic_cast<RooAbsPdf*>(**o_pdf_it);
        RooAbsPdf *n_pdf_ptr = dynamic_cast<RooAbsPdf*>(**n_pdf_it);
        if (o_pdf_ptr && n_pdf_ptr) pdf_map[o_pdf_ptr] = n_pdf_ptr;
        n_pdf_it->Next();
      } while (o_pdf_it->Next());

      auto o_var_it = o_var.createIterator();
      auto n_var_it = n_var.createIterator();
      do {
        RooRealVar *o_var_ptr = dynamic_cast<RooRealVar*>(**o_var_it);
        RooRealVar *n_var_ptr = dynamic_cast<RooRealVar*>(**n_var_it);
        if (o_var_ptr && n_var_ptr) var_map[o_var_ptr] = n_var_ptr;
        n_var_it->Next();
      } while (o_var_it->Next());

      auto o_real_it = o_real.createIterator();
      auto n_real_it = n_real.createIterator();
      do {
        RooAbsReal *o_real_ptr = dynamic_cast<RooAbsReal*>(**o_real_it);
        RooAbsReal *n_real_ptr = dynamic_cast<RooAbsReal*>(**n_real_it);
        if (o_real_ptr && n_real_ptr) real_map[o_real_ptr] = n_real_ptr;
        n_real_it->Next();
      } while (o_real_it->Next());
    } else {
      wspaces_.insert({it.first, nullptr});
    }
  }


  for (std::size_t i = 0; i < obs_.size(); ++i) {
    if (other.obs_[i]) {
      obs_[i] = std::make_shared<Observation>(*(other.obs_[i]));
      if (obs_[i]->data()) obs_[i]->set_data(data_map.at(obs_[i]->data()));
    }
  }

  for (std::size_t i = 0; i < procs_.size(); ++i) {
    if (other.procs_[i]) {
      procs_[i] = std::make_shared<Process>(*(other.procs_[i]));
      if (procs_[i]->pdf()) procs_[i]->set_pdf(pdf_map.at(procs_[i]->pdf()));
    }
  }

  // Need to update RooAbsPdf pointers here

  for (std::size_t i = 0; i < systs_.size(); ++i) {
    if (other.systs_[i]) {
      systs_[i] = std::make_shared<Systematic>(*(other.systs_[i]));
    }
  }
  for (auto const& it : other.params_) {
    if (it.second) {
      params_.insert({it.first, std::make_shared<Parameter>(*(it.second))});
    } else {
      params_.insert({it.first, nullptr});
    }
  }

  // Need to update RooRealVar pointers here

}

CombineHarvester::CombineHarvester(CombineHarvester&& other) {
  swap(*this, other);
  if (verbosity_ >= 3) {
    log() << "[CombineHarvester] Move-constructor called\n";
  }
}

CombineHarvester& CombineHarvester::operator=(CombineHarvester other) {
  swap(*this, other);
  return (*this);
}

CombineHarvester CombineHarvester::cp() {
  if (verbosity_ >= 2) {
    log() << "[CombineHarvester] Shallow copy method cp() called\n";
  }
  CombineHarvester cpy;
  cpy.obs_      = obs_;
  cpy.procs_    = procs_;
  cpy.systs_    = systs_;
  cpy.params_   = params_;
  cpy.wspaces_  = wspaces_;
  cpy.log_      = log_;
  cpy.verbosity_ = verbosity_;
  cpy.flags_     = flags_;
  return cpy;
}

CombineHarvester & CombineHarvester::PrintAll() {
  std::cout << Observation::PrintHeader;
  for (unsigned i = 0; i < obs_.size(); ++i)
      std::cout << *(obs_[i]) << std::endl;
  std::cout << Process::PrintHeader;
  for (unsigned i = 0; i < procs_.size(); ++i)
      std::cout << *(procs_[i]) << std::endl;
  std::cout << Systematic::PrintHeader;
  for (unsigned i = 0; i < systs_.size(); ++i)
      std::cout << *(systs_[i]) << std::endl;
  std::cout << Parameter::PrintHeader;
  for (auto const& it : params_) std::cout << *(it.second) << std::endl;
  return *this;
}

/**
 * \brief Resolve a HistMapping object for the given Observation and load the
 *TH1 or RooAbsData from the TFile and naming pattern this HistMapping provides
 *
 * \pre
 *  1. Input Observation must not already store a shape (TH1 or RooAbsData)
 *  2. It must be possible to resolve this input to a single HistMapping
 *     object. See the documentation of CombineHarvester::ResolveMapping for
 *     further details.
 *  3. The resolved HistMapping object must provide a valid TFile. This file
 *     must contain the TH1 or RooAbsData object as directed by the pattern in
 *     the HistMapping after the substitution of the Observation properties.
 *
 * **If any of these conditions is unmet an exception will be thrown
 *and the calling CombineHarvester instance will remain unmodified.**
 *
 * \post
 *  1. The Observation object will contain either a TH1 or a pointer to a
 *     RooAbsData object.
 *  2. If a TH1 is loaded, the Observation rate will be set to the Integral of
 *     the histogram, discarding any existing value.
 */
void CombineHarvester::LoadShapes(Observation* entry,
                                     std::vector<HistMapping> const& mappings) {
  // Pre-condition #1
  if (entry->shape() || entry->data()) {
    throw std::runtime_error(FNERROR("Observation already contains a shape"));
  }

  if (verbosity_ >= 2) {
    LOGLINE(log(), "Extracting shapes for Observation:");
    log() << Observation::PrintHeader << *entry << "\n";
    LOGLINE(log(), "Mappings:");
    for (HistMapping const& m : mappings) log() << m << "\n";
  }

  // Pre-condition #2
  // ResolveMapping will throw if this fails
  HistMapping mapping =
      ResolveMapping(entry->process(), entry->bin(), mappings);
  boost::replace_all(mapping.pattern, "$CHANNEL", entry->bin());
  boost::replace_all(mapping.pattern, "$BIN", entry->bin());
  boost::replace_all(mapping.pattern, "$PROCESS", entry->process());
  boost::replace_all(mapping.pattern, "$MASS", entry->mass());

  if (verbosity_ >= 2) {
    LOGLINE(log(), "Resolved Mapping:");
    log() << mapping << "\n";
  }

  if (mapping.IsHist()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type in TH1");
    // Pre-condition #3
    // GetClonedTH1 will throw if this fails
    std::unique_ptr<TH1> h = GetClonedTH1(mapping.file.get(), mapping.pattern);
    // Post-conditions #1 and #2
    entry->set_shape(std::move(h), true);
  } else if (mapping.IsData()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is RooAbsData");
    // Pre-condition #3
    // SetupWorkspace will throw if workspace not found
    StrPair res = SetupWorkspace(mapping);
    RooAbsData* dat = wspaces_[res.first]->data(res.second.c_str());
    if (!dat) {
      throw std::runtime_error(FNERROR("RooAbsData not found in workspace"));
    }
    // Post-condition #1
    entry->set_data(dat);
    entry->set_rate(dat->sumEntries());
  }
}

/**
 * \brief Resolve a HistMapping object for the given Process and load the TH1 or
 *RooAbsPdf from the TFile and naming pattern this HistMapping provides
 *
 * \pre
 *  1. Input Process must not already store a shape (TH1 or RooAbsPdf) or
 *     normalisation (RooAbsReal) object
 *  2. It must be possible to resolve this input to a single HistMapping
 *     object. See the documentation of CombineHarvester::ResolveMapping for
 *     further details.
 *  3. The resolved HistMapping object must provide a valid TFile. This file
 *     must contain the TH1 or RooAbsPdf object as directed by the pattern in
 *     the HistMapping after the substitution of the Process properties.
 *
 * **If any of these conditions is unmet an exception will be thrown
 *and the calling CombineHarvester instance will remain unmodified.**
 *
 * \post
 *  1. The Process object will contain either a TH1 or a pointer to a RooAbsPdf
 *     object.
 *  2. If a TH1 is loaded, the Process rate will be set to the Integral of the
 *     histogram, discarding any existing value.
 *  3. It will also contain a pointer to a RooAbsReal object for the
 *     normalisation, but only when the mapping resolved to a RooAbsPdf and a
 *     RooAbsReal object can be found by appending the  suffix "_norm" to the
 *     same pattern.
 *  4. If a RooAbsPdf is loaded, any dependent parameters will be added to the
 *     CombineHarvester instance.
 */
void CombineHarvester::LoadShapes(Process* entry,
                                     std::vector<HistMapping> const& mappings) {
  // Pre-condition #1
  if (entry->shape() || entry->pdf()) {
    throw std::runtime_error(FNERROR("Process already contains a shape"));
  }

  if (verbosity_ >= 2) {
    LOGLINE(log(), "Extracting shapes for Process:");
    log() << Process::PrintHeader << *entry << "\n";
    LOGLINE(log(), "Mappings:");
    for (HistMapping const& m : mappings) log() << m << "\n";
  }

  // Pre-condition #2
  // ResolveMapping will throw if this fails
  HistMapping mapping =
      ResolveMapping(entry->process(), entry->bin(), mappings);
  boost::replace_all(mapping.pattern, "$CHANNEL", entry->bin());
  boost::replace_all(mapping.pattern, "$BIN", entry->bin());
  boost::replace_all(mapping.pattern, "$PROCESS", entry->process());
  boost::replace_all(mapping.pattern, "$MASS", entry->mass());

  if (verbosity_ >= 2) {
    LOGLINE(log(), "Resolved Mapping:");
    log() << mapping << "\n";
  }

  if (mapping.IsHist()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type in TH1");
    // Pre-condition #3
    // GetClonedTH1 will throw if this fails
    std::unique_ptr<TH1> h = GetClonedTH1(mapping.file.get(), mapping.pattern);

    if (flags_.at("zero-negative-bins-on-import")) {
      if (HasNegativeBins(h.get())) {
        LOGLINE(log(), "Process shape has negative bin => fixing to zero");
        log() << Process::PrintHeader << *entry << "\n";
        ZeroNegativeBins(h.get());
      }
    }
    // Post-conditions #1 and #2
    entry->set_shape(std::move(h), true);
  } else if (mapping.IsPdf()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is RooAbsPdf");
    // Pre-condition #3
    // SetupWorkspace will throw if workspace not found
    StrPair res = SetupWorkspace(mapping);
    RooAbsPdf* pdf = wspaces_[res.first]->pdf(res.second.c_str());
    // Pre-condition #3
    if (!pdf) {
      throw std::runtime_error(FNERROR("RooAbsPdf not found in workspace"));
    }
    if (verbosity_ >= 2) {
      pdf->printStream(log(), pdf->defaultPrintContents(0),
                     pdf->defaultPrintStyle(0), "[LoadShapes] ");
    }
    // Post-condition #1
    entry->set_pdf(pdf);

    HistMapping norm_mapping = mapping;
    norm_mapping.pattern += "_norm";
    StrPair norm_res = SetupWorkspace(norm_mapping);
    RooAbsReal* norm =
        wspaces_[norm_res.first]->function(norm_res.second.c_str());
    if (norm) {
      // Post-condition #3
      entry->set_norm(norm);
      if (verbosity_ >= 2) {
        LOGLINE(log(), "Normalisation RooAbsReal found");
        norm->printStream(log(), norm->defaultPrintContents(0),
                          norm->defaultPrintStyle(0), "[LoadShapes] ");
      }
      // If we can upcast norm to a RooRealVar then we can interpret
      // it as a free parameter that should be added to the list
      RooRealVar* norm_var = dynamic_cast<RooRealVar*>(norm);
      if (norm_var) {
        RooArgSet tmp_set(*norm_var);
        ImportParameters(&tmp_set);
      }
    }

    // Post-condition #4
    // Import any paramters of the RooAbsPdf and the RooRealVar
    RooAbsData const* data_obj = FindMatchingData(entry);
    if (data_obj) {
      if (verbosity_ >= 2) LOGLINE(log(), "Matching RooAbsData has been found");
      ImportParameters(pdf->getParameters(data_obj));
      if (norm) ImportParameters(norm->getParameters(data_obj));
    } else {
      if (verbosity_ >= 2)
        LOGLINE(log(), "No RooAbsData found, assume observable CMS_th1x");
      RooRealVar mx("CMS_th1x" , "CMS_th1x", 0, 1);
      RooArgSet tmp_set(mx);
      ImportParameters(pdf->getParameters(&tmp_set));
      if (norm) ImportParameters(norm->getParameters(&tmp_set));
    }
  }
}

void CombineHarvester::LoadShapes(Systematic* entry,
                                     std::vector<HistMapping> const& mappings) {
  if (entry->shape_u() || entry->shape_d()) {
    throw std::runtime_error(FNERROR("Systematic already contains a shape"));
  }

  if (verbosity_ >= 2) {
    LOGLINE(log(), "Extracting shapes for Systematic:");
    log() << Systematic::PrintHeader << *entry << "\n";
    LOGLINE(log(), "Mappings:");
    for (HistMapping const& m : mappings) log() << m << "\n";
  }

  // Pre-condition #2
  // ResolveMapping will throw if this fails
  HistMapping const& mapping =
      ResolveMapping(entry->process(), entry->bin(), mappings);

  if (!mapping.IsHist()) {
    throw std::runtime_error(
        FNERROR("Resolved mapping is not of histogram type"));
  }

  std::string p = mapping.pattern;
  boost::replace_all(p, "$CHANNEL", entry->bin());
  boost::replace_all(p, "$BIN", entry->bin());
  boost::replace_all(p, "$PROCESS", entry->process());
  boost::replace_all(p, "$MASS", entry->mass());
  std::string p_s = mapping.syst_pattern;
  boost::replace_all(p_s, "$CHANNEL", entry->bin());
  boost::replace_all(p_s, "$BIN", entry->bin());
  boost::replace_all(p_s, "$PROCESS", entry->process());
  boost::replace_all(p_s, "$MASS", entry->mass());
  std::string p_s_hi = p_s;
  std::string p_s_lo = p_s;
  boost::replace_all(p_s_hi, "$SYSTEMATIC", entry->name() + "Up");
  boost::replace_all(p_s_lo, "$SYSTEMATIC", entry->name() + "Down");

  std::unique_ptr<TH1> h = GetClonedTH1(mapping.file.get(), p);
  std::unique_ptr<TH1> h_u = GetClonedTH1(mapping.file.get(), p_s_hi);
  std::unique_ptr<TH1> h_d = GetClonedTH1(mapping.file.get(), p_s_lo);

  if (flags_.at("zero-negative-bins-on-import")) {
    if (HasNegativeBins(h.get())) {
      LOGLINE(log(), "Systematic shape has negative bins => fixing to zero");
      log() << Systematic::PrintHeader << *entry << "\n";
      ZeroNegativeBins(h.get());
    }

    if (HasNegativeBins(h_u.get())) {
      LOGLINE(log(), "Systematic shape_u has negative bins => fixing to zero");
      log() << Systematic::PrintHeader << *entry << "\n";
      ZeroNegativeBins(h_u.get());
    }

    if (HasNegativeBins(h_d.get())) {
      LOGLINE(log(), "Systematic shape_d has negative bins => fixing to zero");
      log() << Systematic::PrintHeader << *entry << "\n";
      ZeroNegativeBins(h_d.get());
    }
  }

  entry->set_shapes(std::move(h_u), std::move(h_d), h.get());
}

/**
 * Determines the best-matched HistMapping for a given process
 *
 * Given a process `p` and and bin `b`, this function will first call
 * GenerateShapeMapAttempts to obtain a list of process,bin pairs to search
 * for in the given vector of HistMapping object. The fist attempt will be for
 * a mapping that matches `p,b` exactly, followed by the substitution of each
 * with the wildcard `*`, making a total of four attempts: `p,b`, `*,b`,
 * `p,*`, `*,*`
 *
 *  \pre One of the generated mapping attempts must match exactly with one of
 *  the supplied HistMapping objects, otherwise an exception will be thrown.
 *
 * @param process Process name to match
 * @param bin Bin name to match
 * @param mappings Vector of HistMapping objects
 *
 * @return A reference to the selected HistMapping object
 */
HistMapping const& CombineHarvester::ResolveMapping(
    std::string const& process, std::string const& bin,
    std::vector<HistMapping> const& mappings) {
  StrPairVec attempts = GenerateShapeMapAttempts(process, bin);
  for (unsigned a = 0; a < attempts.size(); ++a) {
    for (unsigned m = 0; m < mappings.size(); ++m) {
      if ((attempts[a].first == mappings[m].process) &&
        (attempts[a].second == mappings[m].category)) {
        return mappings[m];
      }
    }
  }
  // If we get this far then we didn't find a valid mapping
  throw std::runtime_error("ResolveMapping: valid mapping not found");
}

CombineHarvester::StrPairVec CombineHarvester::GenerateShapeMapAttempts(
    std::string process, std::string category) {
  CombineHarvester::StrPairVec result;
  result.push_back(std::make_pair(process   , category));
  result.push_back(std::make_pair("*"       , category));
  result.push_back(std::make_pair(process   , "*"));
  result.push_back(std::make_pair("*"       , "*"));
  return result;
}

std::pair<std::string, std::string> CombineHarvester::SetupWorkspace(
    HistMapping const& mapping) {
  std::string p = mapping.pattern;
  std::pair<std::string, std::string> res;
  std::size_t colon = p.find_last_of(':');
  if (colon != p.npos) {
    res.first = p.substr(0, colon);
    res.second = p.substr(colon+1);
  } else {
    std::cout << "Something went wrong here\n";
    return res;
  }
  if (!wspaces_.count(res.first)) {
    if (mapping.file) {
      mapping.file->cd();
      RooWorkspace* w =
          dynamic_cast<RooWorkspace*>(gDirectory->Get(res.first.c_str()));
      AddWorkspace(w);
      delete w;
    }
  }
  return res;
}

void CombineHarvester::ImportParameters(RooArgSet *vars) {
  auto x = vars->createIterator();
  do {
    RooRealVar *y = dynamic_cast<RooRealVar*>(**x);
    if (y) {
      if (!params_.count(y->GetName())) {
        if (verbosity_ >= 1) {
          log() << "[ImportParameters] Creating parameter from RooRealVar:\n";
          y->printStream(log(), y->defaultPrintContents(0),
                           y->defaultPrintStyle(0), "[ImportParameters] ");
        }
        Parameter par;
        par.set_name(y->GetName());
        par.set_val(y->getVal());
        par.set_err_d(0.0);
        par.set_err_u(0.0);
        params_[y->GetName()] = std::make_shared<Parameter>(par);
      } else {
        if (verbosity_ >= 1)
          log() << "[ImportParameters] Parameter \"" << y->GetName()
                << "\" already exists\n";
      }
      Parameter *param = params_[y->GetName()].get();
      std::vector<RooRealVar *> & p_vars = param->vars();
      if (std::find(p_vars.begin(), p_vars.end(), y) == p_vars.end()) {
        p_vars.push_back(y);
       if (verbosity_ >= 1)
         log() << "[ImportParameters] Parameter now stores " << p_vars.size()
               << " link(s) to RooRealVar objects\n";
      } else {
        if (verbosity_ >= 1)
          log() << "[ImportParameters] Parameter already stores link to "
                   "RooRealVar object\n";
      }
    }
  } while (x->Next());
}

RooAbsData const* CombineHarvester::FindMatchingData(Process const* proc) {
  RooAbsData const* data_obj = nullptr;
  for (unsigned i = 0; i < obs_.size(); ++i) {
    if (proc->bin() == obs_[i]->bin() &&
        proc->bin_id() == obs_[i]->bin_id()) {
      data_obj = obs_[i]->data();
    }
  }
  return data_obj;
}
}
