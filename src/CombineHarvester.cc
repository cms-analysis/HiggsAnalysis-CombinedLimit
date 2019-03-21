#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <list>
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"

namespace ch {

CombineHarvester::CombineHarvester() : verbosity_(0), log_(&(std::cout)) {
  // if (verbosity_ >= 3) {
    // log() << "[CombineHarvester] Constructor called: " << this << "\n";
  // }
  flags_["check-negative-bins-on-import"] = true;
  flags_["zero-negative-bins-on-import"] = false;
  flags_["allow-missing-shapes"] = true;
  flags_["workspaces-use-clone"] = false;
  flags_["workspace-uuid-recycle"] = true;
  flags_["import-parameter-err"] = true;
  flags_["filters-use-regex"] = false;
  // std::cout << "[CombineHarvester] Constructor called for " << this << "\n";
}

CombineHarvester::~CombineHarvester() {
  // std::cout << "[CombineHarvester] Destructor called for " << this << "\n";
}

void swap(CombineHarvester& first, CombineHarvester& second) {
  using std::swap;
  // std::cout << "[CombineHarvester] Swap " << &first << " <-> "
  //           << &second << "\n";
  swap(first.obs_, second.obs_);
  swap(first.procs_, second.procs_);
  swap(first.systs_, second.systs_);
  swap(first.params_, second.params_);
  swap(first.wspaces_, second.wspaces_);
  swap(first.verbosity_, second.verbosity_);
  swap(first.flags_, second.flags_);
  swap(first.post_lines_, second.post_lines_);
  swap(first.log_, second.log_);
  swap(first.auto_stats_settings_, second.auto_stats_settings_);
}

CombineHarvester::CombineHarvester(CombineHarvester const& other)
    : obs_(other.obs_),
      procs_(other.procs_),
      systs_(other.systs_),
      params_(other.params_),
      wspaces_(other.wspaces_),
      flags_(other.flags_),
      auto_stats_settings_(other.auto_stats_settings_),
      post_lines_(other.post_lines_),
      verbosity_(other.verbosity_),
      log_(other.log_) {
  // std::cout << "[CombineHarvester] Copy-constructor called " << &other
  //     << " -> " << this << "\n";
}


void CombineHarvester::SetFlag(std::string const& flag, bool const& value) {
  auto it = flags_.find(flag);
  if (it != flags_.end()) {
    FNLOG(std::cout) << "Changing value of flag \"" << it->first << "\" from " << it->second << " to " << value << "\n";
    it->second = value;
  } else {
    FNLOG(std::cout) << "Created new flag \"" << flag << "\" with value " << value << "\n";
    flags_[flag] = value;
  }
}

bool CombineHarvester::GetFlag(std::string const& flag) const {
  auto it = flags_.find(flag);
  if (it != flags_.end()) {
    return it->second;
  } else {
    throw std::runtime_error(FNERROR("Flag " + flag + " is not defined"));
  }
}

CombineHarvester CombineHarvester::deep() {
  CombineHarvester cpy;
  // std::cout << "[CombineHarvester] Deep copy called " << this
  //     << " -> " << &cpy << "\n";
  cpy.obs_.resize(obs_.size());
  cpy.procs_.resize(procs_.size());
  cpy.systs_.resize(systs_.size());
  cpy.flags_ = flags_;
  cpy.verbosity_ = verbosity_;
  cpy.post_lines_ = post_lines_;
  cpy.log_ = log_;

  // Build a map of workspace object pointers
  std::map<RooAbsData const*, RooAbsData *> dat_map;
  std::map<RooAbsReal  const*, RooAbsReal  *> pdf_map;
  std::map<RooRealVar const*, RooRealVar *> var_map;
  std::map<RooAbsReal const*, RooAbsReal *> fun_map;

  // Loop through each workspace and make a full copy
  for (auto const& it : wspaces_) {
    std::string ws_name = it.first;
    RooWorkspace const& o_wsp = *(it.second.get());
    cpy.wspaces_.insert({ws_name, std::make_shared<RooWorkspace>(o_wsp)});
    RooWorkspace & n_wsp = *(cpy.wspaces_.at(ws_name).get());

    std::list<RooAbsData *> o_dat = o_wsp.allData();
    RooArgSet const&        o_pdf = o_wsp.allPdfs();
    RooArgSet const&        o_var = o_wsp.allVars();
    RooArgSet const&        o_fun = o_wsp.allFunctions();

    std::list<RooAbsData *> n_dat = n_wsp.allData();
    RooArgSet const&        n_pdf = n_wsp.allPdfs();
    RooArgSet const&        n_var = n_wsp.allVars();
    RooArgSet const&        n_fun = n_wsp.allFunctions();

    auto o_dat_it = o_dat.begin();
    auto n_dat_it = n_dat.begin();
    for (; o_dat_it != o_dat.end(); ++o_dat_it, ++n_dat_it) {
      dat_map[*o_dat_it] = *n_dat_it;
    }

    auto o_pdf_it = o_pdf.createIterator();
    auto n_pdf_it = n_pdf.createIterator();
    do {
      RooAbsPdf *o_pdf_ptr = static_cast<RooAbsPdf*>(**o_pdf_it);
      RooAbsPdf *n_pdf_ptr = static_cast<RooAbsPdf*>(**n_pdf_it);
      if (o_pdf_ptr && n_pdf_ptr) pdf_map[o_pdf_ptr] = n_pdf_ptr;
      n_pdf_it->Next();
    } while (o_pdf_it->Next());

    auto o_var_it = o_var.createIterator();
    auto n_var_it = n_var.createIterator();
    do {
      RooRealVar *o_var_ptr = static_cast<RooRealVar*>(**o_var_it);
      RooRealVar *n_var_ptr = static_cast<RooRealVar*>(**n_var_it);
      if (o_var_ptr && n_var_ptr) var_map[o_var_ptr] = n_var_ptr;
      n_var_it->Next();
    } while (o_var_it->Next());

    auto o_fun_it = o_fun.createIterator();
    auto n_fun_it = n_fun.createIterator();
    do {
      RooAbsReal *o_fun_ptr = static_cast<RooAbsReal*>(**o_fun_it);
      RooAbsReal *n_fun_ptr = static_cast<RooAbsReal*>(**n_fun_it);
      if (o_fun_ptr && n_fun_ptr) fun_map[o_fun_ptr] = n_fun_ptr;
      n_fun_it->Next();
    } while (o_fun_it->Next());
  }


  for (std::size_t i = 0; i < cpy.obs_.size(); ++i) {
    if (obs_[i]) {
      cpy.obs_[i] = std::make_shared<Observation>(*(obs_[i]));
      if (obs_[i]->data())
        cpy.obs_[i]->set_data(dat_map.at(obs_[i]->data()));
    }
  }

  for (std::size_t i = 0; i < cpy.procs_.size(); ++i) {
    if (procs_[i]) {
      cpy.procs_[i] = std::make_shared<Process>(*(procs_[i]));
      if (procs_[i]->pdf())
        cpy.procs_[i]->set_pdf(pdf_map.at(procs_[i]->pdf()));
      if (procs_[i]->observable())
        cpy.procs_[i]->set_observable(var_map.at(procs_[i]->observable()));
      if (procs_[i]->data())
        cpy.procs_[i]->set_data(dat_map.at(procs_[i]->data()));
      if (procs_[i]->norm())
        cpy.procs_[i]->set_norm(fun_map.at(procs_[i]->norm()));
    }
  }

  for (std::size_t i = 0; i < cpy.systs_.size(); ++i) {
    if (systs_[i]) {
      cpy.systs_[i] = std::make_shared<Systematic>(*(systs_[i]));
      if (systs_[i]->data_u() || systs_[i]->data_d()) {
        cpy.systs_[i]->set_data(
            static_cast<RooDataHist*>(dat_map.at(systs_[i]->data_u())),
            static_cast<RooDataHist*>(dat_map.at(systs_[i]->data_d())),
            nullptr);
      }
    }
  }

  for (auto const& it : params_) {
    if (it.second) {
      cpy.params_.insert({it.first, std::make_shared<Parameter>(*(it.second))});
    } else {
      params_.insert({it.first, nullptr});
    }
  }

  for (auto it : cpy.params_) {
    for (unsigned i = 0; i < it.second->vars().size(); ++i) {
      it.second->vars()[i] = var_map.at(it.second->vars()[i]);
    }
  }
  return cpy;
}

CombineHarvester::CombineHarvester(CombineHarvester&& other) {
  // std::cout << "[CombineHarvester] Move-constructor called " <<
  //     &other << " -> " << this << "\n";
  swap(*this, other);
}

CombineHarvester& CombineHarvester::operator=(CombineHarvester other) {
  // std::cout << "[CombineHarvester] Assignment operator called " << &other <<
  // " -> " << this << "\n";
  swap(*this, other);
  return (*this);
}

CombineHarvester CombineHarvester::cp() {
  // std::cout << "[CombineHarvester] Shallow copy method cp() called from " <<
  // this << "\n";
  return CombineHarvester(*this);
}

CombineHarvester & CombineHarvester::PrintAll() {
  return PrintObs().PrintProcs().PrintSysts().PrintParams();
}

CombineHarvester & CombineHarvester::PrintObs() {
  std::cout << Observation::PrintHeader;
  for (unsigned i = 0; i < obs_.size(); ++i)
      std::cout << *(obs_[i]) << "\n";
  return *this;
}

CombineHarvester & CombineHarvester::PrintProcs() {
  std::cout << Process::PrintHeader;
  for (unsigned i = 0; i < procs_.size(); ++i)
      std::cout << *(procs_[i]) << "\n";
  return *this;
}

CombineHarvester & CombineHarvester::PrintSysts() {
  std::cout << Systematic::PrintHeader;
  for (unsigned i = 0; i < systs_.size(); ++i)
      std::cout << *(systs_[i]) << "\n";
  return *this;
}

CombineHarvester & CombineHarvester::PrintParams() {
  std::cout << Parameter::PrintHeader;
  for (auto const& it : params_) std::cout << *(it.second) << "\n";
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
 *     further details. If no mappings have been provided, assume this is a
 *     counting experiment and return.
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
  if (mappings.size() == 0) {
    if (verbosity_ >= 2) LOGLINE(log(), "No mappings in this card, assuming FAKE");
    return;
  }
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

  if (mapping.is_fake) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is FAKE");
  } else if (mapping.IsHist()) {
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
    RooAbsData* dat = mapping.ws->data(mapping.WorkspaceObj().c_str());
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
  if (mappings.size() == 0) {
    if (verbosity_ >= 2) LOGLINE(log(), "No mappings in this card, assuming FAKE");
    return;
  }
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

  if (mapping.is_fake) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is FAKE");
  } else if (mapping.IsHist()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is TH1");
    // Pre-condition #3
    // GetClonedTH1 will throw if this fails
    std::unique_ptr<TH1> h = GetClonedTH1(mapping.file.get(), mapping.pattern);

    if (flags_.at("check-negative-bins-on-import")) {
      if (HasNegativeBins(h.get())) {
        LOGLINE(log(), "Warning: process shape has negative bins");
        log() << Process::PrintHeader << *entry << "\n";
        if (flags_.at("zero-negative-bins-on-import")) {
          ZeroNegativeBins(h.get());
        }
      }
    }
    // Post-conditions #1 and #2
    entry->set_shape(std::move(h), true);
  } else if (mapping.IsPdf()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is RooAbsPdf/RooAbsData");
    // Pre-condition #3
    // Try and get this as RooAbsData first. If this doesn't work try pdf
    RooAbsData* data = mapping.ws->data(mapping.WorkspaceObj().c_str());
    RooAbsReal* pdf = mapping.ws->function(mapping.WorkspaceObj().c_str());
    if (data) {
      if (verbosity_ >= 2) {
        data->printStream(log(), data->defaultPrintContents(0),
                       data->defaultPrintStyle(0), "[LoadShapes] ");
      }
      entry->set_data(data);
      entry->set_rate(data->sumEntries());
    } else if (pdf) {
      if (verbosity_ >= 2) {
        pdf->printStream(log(), pdf->defaultPrintContents(0),
                       pdf->defaultPrintStyle(0), "[LoadShapes] ");
      }
      // Post-condition #1
      entry->set_pdf(pdf);
    } else { // Pre-condition #3
      if (flags_.at("allow-missing-shapes")) {
        LOGLINE(log(), "Warning, shape missing:");
        log() << Process::PrintHeader << *entry << "\n";
      } else {
        throw std::runtime_error(FNERROR("RooAbsPdf not found in workspace"));
      }
    }

    HistMapping norm_mapping = mapping;
    if (!data && norm_mapping.syst_pattern != "") {
      // For when we're not parsing a datacard, syst_pattern is being using to
      // note a different mapping for the normalisation term
      norm_mapping.pattern = norm_mapping.syst_pattern;
      boost::replace_all(norm_mapping.pattern, "$CHANNEL", entry->bin());
      boost::replace_all(norm_mapping.pattern, "$BIN", entry->bin());
      boost::replace_all(norm_mapping.pattern, "$PROCESS", entry->process());
      boost::replace_all(norm_mapping.pattern, "$MASS", entry->mass());
    } else {
      norm_mapping.pattern += "_norm";
    }
    RooAbsReal* norm =
        norm_mapping.ws->function(norm_mapping.WorkspaceObj().c_str());
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
      if (pdf&&!data) {
        RooArgSet argset = ParametersByName(pdf, data_obj->get());
        ImportParameters(&argset);
        if (!entry->observable()) {
          std::string var_name;
          if (data_obj) var_name = data_obj->get()->first()->GetName();
          entry->set_observable(
              (RooRealVar*)entry->pdf()->findServer(var_name.c_str()));
        }
      }
      if (norm) {
        RooArgSet argset = ParametersByName(norm, data_obj->get());
        ImportParameters(&argset);
      }
    } else {
      if (verbosity_ >= 2)
        LOGLINE(log(), "No RooAbsData found, assume observable CMS_th1x");
      RooRealVar mx("CMS_th1x" , "CMS_th1x", 0, 1);
      RooArgSet tmp_set(mx);
      if (pdf) {
        RooArgSet argset = ParametersByName(pdf, &tmp_set);
        ImportParameters(&argset);
      }
      if (norm) {
        RooArgSet argset = ParametersByName(norm, &tmp_set);
        ImportParameters(&argset);
      }
    }
  }
}

void CombineHarvester::LoadShapes(Systematic* entry,
                                     std::vector<HistMapping> const& mappings) {
  if (entry->shape_u() || entry->shape_d() ||
      entry->data_u() || entry->data_d()) {
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
  HistMapping mapping =
      ResolveMapping(entry->process(), entry->bin(), mappings);
  boost::replace_all(mapping.pattern, "$CHANNEL", entry->bin());
  boost::replace_all(mapping.pattern, "$BIN", entry->bin());
  boost::replace_all(mapping.pattern, "$PROCESS", entry->process());
  boost::replace_all(mapping.pattern, "$MASS", entry->mass());
  std::string p_s =
      mapping.IsPdf() ? mapping.SystWorkspaceObj() : mapping.syst_pattern;
  boost::replace_all(p_s, "$CHANNEL", entry->bin());
  boost::replace_all(p_s, "$BIN", entry->bin());
  boost::replace_all(p_s, "$PROCESS", entry->process());
  boost::replace_all(p_s, "$MASS", entry->mass());
  std::string p_s_hi = p_s;
  std::string p_s_lo = p_s;
  boost::replace_all(p_s_hi, "$SYSTEMATIC", entry->name() + "Up");
  boost::replace_all(p_s_lo, "$SYSTEMATIC", entry->name() + "Down");
  if (mapping.IsHist()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is TH1");
    std::unique_ptr<TH1> h = GetClonedTH1(mapping.file.get(), mapping.pattern);
    std::unique_ptr<TH1> h_u = GetClonedTH1(mapping.file.get(), p_s_hi);
    std::unique_ptr<TH1> h_d = GetClonedTH1(mapping.file.get(), p_s_lo);

    if (flags_.at("check-negative-bins-on-import")) {
      if (HasNegativeBins(h.get())) {
        LOGLINE(log(), "Warning: Systematic shape has negative bins");
        log() << Systematic::PrintHeader << *entry << "\n";
        if (flags_.at("zero-negative-bins-on-import")) {
          ZeroNegativeBins(h.get());
        }
      }

      if (HasNegativeBins(h_u.get())) {
        LOGLINE(log(), "Warning: Systematic shape_u has negative bins");
        log() << Systematic::PrintHeader << *entry << "\n";
        if (flags_.at("zero-negative-bins-on-import")) {
          ZeroNegativeBins(h_u.get());
        }
      }

      if (HasNegativeBins(h_d.get())) {
        LOGLINE(log(), "Warning: Systematic shape_d has negative bins");
        log() << Systematic::PrintHeader << *entry << "\n";
        if (flags_.at("zero-negative-bins-on-import")) {
          ZeroNegativeBins(h_d.get());
        }
      }
    }
    entry->set_shapes(std::move(h_u), std::move(h_d), h.get());
  } else if (mapping.IsPdf()) {
    if (verbosity_ >= 2) LOGLINE(log(), "Mapping type is RooDataHist");
    // Try and get this as RooAbsData first. If this doesn't work try pdf
    RooDataHist* h =
        dynamic_cast<RooDataHist*>(mapping.sys_ws->data(mapping.WorkspaceObj().c_str()));
    RooDataHist* h_u =
        dynamic_cast<RooDataHist*>(mapping.sys_ws->data(p_s_hi.c_str()));
    RooDataHist* h_d =
        dynamic_cast<RooDataHist*>(mapping.sys_ws->data(p_s_lo.c_str()));
    if (!h || !h_u || !h_d) {
      if (flags_.at("allow-missing-shapes")) {
        LOGLINE(log(), "Warning, shape missing:");
        log() << Systematic::PrintHeader << *entry << "\n";
      } else {
        throw std::runtime_error(
            FNERROR("All shapes must be of type RooDataHist"));
      }
    } else {
      entry->set_data(h_u, h_d, h);
    }
  } else {
    throw std::runtime_error(
        FNERROR("Resolved mapping is not of TH1 / RooAbsData type"));
  }
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
  FNLOG(log()) << "Searching for mapping for (" << bin << "," << process << ")\n";
  FNLOG(log()) << "Avaiable mappings:\n";
  for (auto const& m : mappings) {
    FNLOG(log()) << m << "\n";
  }
  throw std::runtime_error(FNERROR("Valid mapping not found"));
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

std::shared_ptr<RooWorkspace> CombineHarvester::SetupWorkspace(
    RooWorkspace const& ws, bool can_rename) {
  bool name_in_use = false;

  // 1) Does the UUID of this ws match any other ws in the map?
  for (auto const& it : wspaces_) {
    // - Yes: just return a ptr to the matched ws
    if (it.second->uuid() == ws.uuid() && flags_["workspace-uuid-recycle"]) {
      FNLOGC(log(), verbosity_ >= 1)
          << "Workspace with name " << it.second->GetName()
          << " has the same UUID, will use this one\n";
      return it.second;
    }
    if (!name_in_use && strcmp(it.second->GetName(), ws.GetName()) == 0) {
      name_in_use = true;
    }
  }

  // 2) No: Ok will clone it in. Is the ws name already in use?
  if (!name_in_use) {
    // - No: clone with same name and return
    // IMPORTANT: Don't used RooWorkspace::Clone(), it seems to introduce
    // bugs
    if (GetFlag("workspaces-use-clone")) {
      wspaces_[std::string(ws.GetName())] = std::shared_ptr<RooWorkspace>(
          reinterpret_cast<RooWorkspace*>(ws.Clone()));      
    } else {
      wspaces_[std::string(ws.GetName())] =
          std::make_shared<RooWorkspace>(RooWorkspace(ws));  
    }
    return wspaces_.at(ws.GetName());
  }

  // 3) Am I allowed to rename (default no)?
  if (!can_rename) {
    // - No: throw runtime error
    throw std::runtime_error(FNERROR("A different workspace with name " +
                                     std::string(ws.GetName()) +
                                     " already exists"));
  }

  // 4) Yes: determine first available nameX, clone, return
  std::set<int> used_ints = {0};
  std::string src_name(ws.GetName());
  for (auto const& it : wspaces_) {
    std::string test_name(it.second->GetName());
    if (test_name.find(src_name) == 0) {
      std::string postfix = test_name.substr(src_name.size());
      try {
        int number = boost::lexical_cast<int>(postfix);
        used_ints.insert(number);
      } catch (boost::bad_lexical_cast & e) {
      }
    }
  }
  std::string new_name =
      src_name + boost::lexical_cast<std::string>(*(used_ints.rbegin()) + 1);
  FNLOGC(log(), verbosity_ >= 1) << "Workspace with name " << src_name
                                 << " already defined, renaming to " << new_name
                                 << "\n";

  std::shared_ptr<RooWorkspace> new_wsp;
  if (GetFlag("workspaces-use-clone")) {
    new_wsp = std::shared_ptr<RooWorkspace>(
        reinterpret_cast<RooWorkspace*>(ws.Clone(new_name.c_str())));    
  } else {
    new_wsp = std::make_shared<RooWorkspace>(RooWorkspace(ws));
  }
  new_wsp->SetName(new_name.c_str());
  wspaces_[new_name] = new_wsp;
  return wspaces_.at(new_name);
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
        if ((y->hasError() || y->hasAsymError()) &&
            flags_["import-parameter-err"]) {
          par.set_err_d(y->getErrorLo());
          par.set_err_u(y->getErrorHi());          
        } else {
          par.set_err_d(0.);
          par.set_err_u(0.);
        }
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

ch::Parameter* CombineHarvester::SetupRateParamVar(std::string const& name,
                                                   double val, bool is_ext_arg) {
  RooWorkspace *ws = nullptr;
  if (!wspaces_.count("_rateParams")) {
    ws = this->SetupWorkspace(RooWorkspace("_rateParams","_rateParams")).get();
  } else {
    ws = wspaces_.at("_rateParams").get();
  }
  // Parameter doesn't exist in the workspace - let's create it
  RooRealVar *var = ws->var(name.c_str());
  if (!var) {
    // The value doesn't matter, our Parameter object will set it later
    RooRealVar tmp_var(name.c_str(), name.c_str(), 0);
    ws->import(tmp_var);
    var = ws->var(name.c_str());
    FNLOGC(log(), verbosity_ > 1)
        << "Created new RooRealVar for rateParam: " << name << "\n";
    if (verbosity_ > 1) var->Print();
  } else {
    FNLOGC(log(), verbosity_ > 1)
        << "Reusing existing RooRealVar for rateParam: " << name << "\n";
  }
  if (is_ext_arg) var->setAttribute("extArg");
  Parameter * param = nullptr;
  if (!params_.count(name)) {
    params_[name] = std::make_shared<Parameter>(Parameter());
    param = params_.at(name).get();
    param->set_name(name);
    param->set_err_u(0.);
    param->set_err_d(0.);
  } else {
    param = params_.at(name).get();
  }
  // If the RooRealVar in the workpsace isn't in the list, add it
  bool var_in_par = false;
  for (auto const& ptr : param->vars()) {
    if (ptr == var) {
      var_in_par = true;
      break;
    }
  }
  if (!var_in_par) {
    param->vars().push_back(var);
  }
  // Then this propagates the value to the RooRealVar
  FNLOGC(log(), verbosity_ > 1) << "Updating parameter value from "
                                << param->val();
  param->set_val(val);
  if (verbosity_ > 1) log() << " to " << param->val() << "\n";
  return param;
}

void CombineHarvester::SetupRateParamFunc(std::string const& name,
                                          std::string const& formula,
                                          std::string const& pars) {
  RooWorkspace *ws = nullptr;
  if (!wspaces_.count("_rateParams")) {
    ws = this->SetupWorkspace(RooWorkspace("_rateParams","_rateParams")).get();
  } else {
    ws = wspaces_.at("_rateParams").get();
  }
  RooAbsReal *form = ws->function(name.c_str());
  // No parameter to make, this is a formula
  if (!form) {
    RooFormulaVar formularvar(name.c_str(), name.c_str(),
                          formula.c_str(),
                          RooArgList(ws->argSet(pars.c_str())));
    ws->import(formularvar);
    form = ws->function(name.c_str());
    FNLOGC(log(), verbosity_ > 1)
        << "Created new RooFormulaVar for rateParam: " << name << "\n";
    if (verbosity_ > 1) form->Print();
  }
}

void CombineHarvester::SetupRateParamWspObj(std::string const& name,
                                          std::string const& obj, bool is_ext_arg) {
  RooWorkspace *ws = nullptr;
  if (!wspaces_.count("_rateParams")) {
    ws = this->SetupWorkspace(RooWorkspace("_rateParams","_rateParams")).get();
  } else {
    ws = wspaces_.at("_rateParams").get();
  }
  ws->import((obj+":"+name).c_str(), RooFit::RecycleConflictNodes());
  ws->arg(name.c_str())->setStringAttribute("wspSource", obj.c_str());
  if (is_ext_arg) ws->arg(name.c_str())->setAttribute("extArg");
}

void CombineHarvester::SetAutoMCStats(CombineHarvester &target, double thresh, bool sig, int mode) {
  for (auto const& bin : this->bin_set()) {
    target.auto_stats_settings_[bin] = AutoMCStatsSettings(thresh, sig, mode);
  }
}

void CombineHarvester::RenameAutoMCStatsBin(std::string const& oldname, std::string const& newname) {
  auto it = auto_stats_settings_.find(oldname);
  if (it != auto_stats_settings_.end()) {
    auto_stats_settings_[newname] = it->second;
    auto_stats_settings_.erase(it);
  }
}

std::set<std::string> CombineHarvester::GetAutoMCStatsBins() const {
  std::set<std::string> result;
  for (auto const& it : auto_stats_settings_) {
    result.insert(it.first);
  }
  return result;
}

void CombineHarvester::AddExtArgValue(std::string const& name, double const& value) {
  ch::Parameter* param = SetupRateParamVar(name, value, true);
  param->set_err_u(0.);
  param->set_err_d(0.);
}
}
