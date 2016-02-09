#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <utility>
#include <set>
#include <fstream>
#include <sstream>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "TDirectory.h"
#include "TH1.h"
#include "RooRealVar.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/MakeUnique.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/GitVersion.h"

namespace ch {

// Extract info from filename using parse rule like:
// ".*{MASS}/{ANALYSIS}_{CHANNEL}_{BINID}_{ERA}.txt"
int CombineHarvester::ParseDatacard(std::string const& filename,
    std::string parse_rules) {
  boost::replace_all(parse_rules, "$ANALYSIS",  "(?<ANALYSIS>[\\w\\.]+)");
  boost::replace_all(parse_rules, "$ERA",       "(?<ERA>[\\w\\.]+)");
  boost::replace_all(parse_rules, "$CHANNEL",   "(?<CHANNEL>[\\w\\.]+)");
  boost::replace_all(parse_rules, "$BINID",     "(?<BINID>[\\w\\.]+)");
  boost::replace_all(parse_rules, "$MASS",      "(?<MASS>[\\w\\.]+)");
  boost::regex rgx(parse_rules);
  boost::smatch matches;
  boost::regex_search(filename, matches, rgx);
  this->ParseDatacard(filename,
    matches.str("ANALYSIS"),
    matches.str("ERA"),
    matches.str("CHANNEL"),
    matches.str("BINID").length() ?
      boost::lexical_cast<int>(matches.str("BINID")) : 0,
    matches.str("MASS"));
  return 0;
}

int CombineHarvester::ParseDatacard(std::string const& filename,
    std::string const& analysis,
    std::string const& era,
    std::string const& channel,
    int bin_id,
    std::string const& mass) {
  // Load the entire datacard into memory as a vector of strings
  std::vector<std::string> lines = ch::ParseFileLines(filename);
  // Loop through lines, trimming whitespace at the beginning or end
  // then splitting each line into a vector of words (using any amount
  // of whitespace as the separator).  We skip any line of zero length
  // or which starts with a "#" or "-" character.
  std::vector<std::vector<std::string>> words;
  for (unsigned i = 0; i < lines.size(); ++i) {
    boost::trim(lines[i]);
    if (lines[i].size() == 0) continue;
    if (lines[i].at(0) == '#' || lines[i].at(0) == '-') continue;
    words.push_back(std::vector<std::string>());
    boost::split(words.back(), lines[i], boost::is_any_of("\t "),
        boost::token_compress_on);
  }

  std::vector<HistMapping> hist_mapping;
  // std::map<std::string, RooAbsData*> data_map;
  std::map<std::string, std::shared_ptr<TFile>> file_store;
  std::map<std::string, std::shared_ptr<RooWorkspace>> ws_store;

  bool start_nuisance_scan = false;
  unsigned r = 0;

  // We will allow cards that describe a single bin to have an "observation"
  // line without a "bin" line above it. We probably won't know the bin name
  // when we parse this line, so we'll store it here and fix it later
  std::shared_ptr<ch::Observation> single_obs = nullptr;
  std::set<std::string> bin_names;

  // Store the groups that we encounter
  std::map<std::string, std::set<std::string>> groups;

  // Loop through the vector of word vectors
  for (unsigned i = 0; i < words.size(); ++i) {
    // Ignore line if it only has one word
    if (words[i].size() <= 1) continue;

    // If the line begins "shapes" then we've
    // found process --> TH1 mapping information
    if (boost::iequals(words[i][0], "shapes") && words[i].size() >= 5) {
      hist_mapping.push_back(HistMapping());
      HistMapping &mapping = hist_mapping.back();
      mapping.process = words[i][1];
      mapping.category = words[i][2];
      // The root file path given in the datacard is relative to the datacard
      // path, so we join the path to the datacard with the path to the file
      std::string dc_path;
      std::size_t slash = filename.find_last_of('/');
      if (slash != filename.npos) {
        dc_path = filename.substr(0, slash) + "/" + words[i][3];
      } else {
        dc_path = words[i][3];
      }
      if (!file_store.count(dc_path))
        file_store[dc_path] = std::make_shared<TFile>(dc_path.c_str());
      mapping.file = file_store.at(dc_path);
      mapping.pattern = words[i][4];
      if (words[i].size() > 5) mapping.syst_pattern = words[i][5];

      if (mapping.IsPdf()) {
        std::string store_key =
            mapping.file->GetName() + mapping.WorkspaceName();
        if (!ws_store.count(store_key)) {
          mapping.file->cd();
          std::shared_ptr<RooWorkspace> ptr(dynamic_cast<RooWorkspace*>(
              gDirectory->Get(mapping.WorkspaceName().c_str())));
          if (!ptr) {
            throw std::runtime_error(FNERROR("Workspace not found in file"));
          }
          ws_store[store_key] = ptr;
        }
        mapping.ws = SetupWorkspace(*(ws_store[store_key]), true);
      }
      if (mapping.IsPdf() && mapping.syst_pattern != "") {
        std::string store_key =
            mapping.file->GetName() + mapping.SystWorkspaceName();
        if (!ws_store.count(store_key)) {
          mapping.file->cd();
          std::shared_ptr<RooWorkspace> ptr(dynamic_cast<RooWorkspace*>(
              gDirectory->Get(mapping.SystWorkspaceName().c_str())));
          if (!ptr) {
            throw std::runtime_error(FNERROR("Workspace not found in file"));
          }
          ws_store[store_key] = ptr;
        }
        mapping.sys_ws = SetupWorkspace(*(ws_store[store_key]), true);
      }
    }

    // We can also have a "FAKE" shape directive
    // Must be four words long: shapes * * FAKE
    if (boost::iequals(words[i][0], "shapes") && words[i].size() == 4 &&
        boost::iequals(words[i][3], "FAKE")) {
      hist_mapping.push_back(HistMapping());
      HistMapping &mapping = hist_mapping.back();
      mapping.process = words[i][1];
      mapping.category = words[i][2];
      mapping.is_fake = true;
    }

    // Want to check this line and the previous one, so need i >= 1.
    // If the first word on this line is "observation" and "bin" on
    // the previous line then we've found the entries for data, and
    // can add Observation objects
    if (i >= 1) {
      if (  boost::iequals(words[i][0], "observation") &&
            boost::iequals(words[i-1][0], "bin") &&
            words[i].size() == words[i-1].size()) {
        for (unsigned p = 1; p < words[i].size(); ++p) {
          auto obs = std::make_shared<Observation>();
          obs->set_bin(words[i-1][p]);
          obs->set_rate(boost::lexical_cast<double>(words[i][p]));
          obs->set_analysis(analysis);
          obs->set_era(era);
          obs->set_channel(channel);
          obs->set_bin_id(bin_id);
          obs->set_mass(mass);

          LoadShapes(obs.get(), hist_mapping);

          obs_.push_back(obs);
        }
      }
    }

    if (boost::iequals(words[i][0], "observation") &&
        !boost::iequals(words[i-1][0], "bin") &&
        words[i].size() == 2 &&
        single_obs.get() == nullptr) {
      for (unsigned p = 1; p < words[i].size(); ++p) {
        single_obs = std::make_shared<Observation>();
        single_obs->set_bin("");
        single_obs->set_rate(boost::lexical_cast<double>(words[i][p]));
        single_obs->set_analysis(analysis);
        single_obs->set_era(era);
        single_obs->set_channel(channel);
        single_obs->set_bin_id(bin_id);
        single_obs->set_mass(mass);
      }
    }

    // Similarly look for the lines indicating the different signal
    // and background processes
    // Once these are found save in line index for the rate line as r
    // to we can refer back to these later, then assume that every
    // line that follows is a nuisance parameter

    if (i >= 3) {
      if (  boost::iequals(words[i][0], "rate") &&
            boost::iequals(words[i-1][0], "process") &&
            boost::iequals(words[i-2][0], "process") &&
            boost::iequals(words[i-3][0], "bin") &&
            words[i].size() == words[i-1].size() &&
            words[i].size() == words[i-2].size() &&
            words[i].size() == words[i-3].size()) {
        for (unsigned p = 1; p < words[i].size(); ++p) {
          auto proc = std::make_shared<Process>();
          proc->set_bin(words[i-3][p]);
          bin_names.insert(words[i-3][p]);
          try {
            int process_id = boost::lexical_cast<int>(words[i-2][p]);
            proc->set_signal(process_id <= 0);
            proc->set_process(words[i-1][p]);
          } catch(boost::bad_lexical_cast &) {
            int process_id = boost::lexical_cast<int>(words[i-1][p]);
            proc->set_signal(process_id <= 0);
            proc->set_process(words[i-2][p]);
          }
          proc->set_rate(boost::lexical_cast<double>(words[i][p]));
          proc->set_analysis(analysis);
          proc->set_era(era);
          proc->set_channel(channel);
          proc->set_bin_id(bin_id);
          proc->set_mass(mass);

          LoadShapes(proc.get(), hist_mapping);

          procs_.push_back(proc);
        }
        r = i;
        start_nuisance_scan = true;
      }
    }

    if (start_nuisance_scan && words[i].size() >= 4) {
      if (boost::iequals(words[i][1], "param")) {
        std::string param_name = words[i][0];
        if (!params_.count(param_name))
          params_[param_name] = std::make_shared<Parameter>(Parameter());
        Parameter * param = params_.at(param_name).get();
        param->set_name(param_name);
        param->set_val(boost::lexical_cast<double>(words[i][2]));
        std::size_t slash_pos = words[i][3].find("/");
        if (slash_pos != words[i][3].npos) {
          param->set_err_d(
            boost::lexical_cast<double>(words[i][3].substr(0, slash_pos)));
          param->set_err_u(
            boost::lexical_cast<double>(words[i][3].substr(slash_pos+1)));
        } else {
          param->set_err_u(+1.0 * boost::lexical_cast<double>(words[i][3]));
          param->set_err_d(-1.0 * boost::lexical_cast<double>(words[i][3]));
        }
        if (words[i].size() >= 5) {
          // We have a range
          std::vector<std::string> tokens;
          boost::split(tokens, words[i][4], boost::is_any_of("[],"));
          if (tokens.size() == 4) {
            param->set_range_d(boost::lexical_cast<double>(tokens[1]));
            param->set_range_u(boost::lexical_cast<double>(tokens[2]));
          }
        }
        continue;  // skip the rest of this now
      }
    }

    if (start_nuisance_scan && words[i].size() >= 5 &&
        boost::iequals(words[i][1], "rateParam")) {
      if (verbosity_ > 1) {
        FNLOG(log()) << "Processing rateParam line:\n";
        for (auto const& str : words[i]) {
          log() << str << "\t";
        }
        log() << "\n";
      }

      bool has_range = words[i].size() == 6 && words[i][5][0] == '[';
      std::string param_name = words[i][0];
      // If this is a free param may need to create a Parameter object
      // If the line has 5 words this must be a floating param, otherwise
      // if it has 6 then it's a floating param
      if (words[i].size() == 5 || has_range) {
        ch::Parameter* param = SetupRateParamVar(
            param_name, boost::lexical_cast<double>(words[i][4]));
        param->set_err_u(0.);
        param->set_err_d(0.);
        if (has_range) {
          std::vector<std::string> tokens;
          boost::split(tokens, words[i][5], boost::is_any_of("[],"));
          if (tokens.size() == 4) {
            param->set_range_d(boost::lexical_cast<double>(tokens[1]));
            param->set_range_u(boost::lexical_cast<double>(tokens[2]));
            FNLOGC(log(), verbosity_ > 1) << "Setting parameter range to " << words[i][5];
          }
        }
      } else if (words[i].size() == 6 && !has_range) {
        SetupRateParamFunc(param_name, words[i][4], words[i][5]);
      }
      for (unsigned p = 1; p < words[r].size(); ++p) {
        bool matches_bin = false;
        bool matches_proc = false;
        int process_id;
        std::string process;
        std::string bin = words[r-3][p];
        // Not great that we repeat this below
        try {
          process_id = boost::lexical_cast<int>(words[r-2][p]);
          process = words[r-1][p];
        } catch(boost::bad_lexical_cast &) {
          process_id = boost::lexical_cast<int>(words[r-1][p]);
          process = words[r-2][p];
        }
        if (words[i][2] == "*" || words[i][2] == bin) {
          matches_bin = true;
        }
        if (words[i][3] == "*" || words[i][3] == process) {
          matches_proc = true;
        }
        if (!matches_bin || !matches_proc) continue;
        auto sys = std::make_shared<Systematic>();
        sys->set_bin(bin);
        sys->set_signal(process_id <= 0);
        sys->set_process(process);
        sys->set_name(param_name);
        sys->set_type("rateParam");
        sys->set_analysis(analysis);
        sys->set_era(era);
        sys->set_channel(channel);
        sys->set_bin_id(bin_id);
        sys->set_mass(mass);
        systs_.push_back(sys);
      }
      continue;
    }

    if (start_nuisance_scan && words[i].size() >= 4 &&
        boost::iequals(words[i][1], "group")) {
      if (!groups.count(words[i][0])) {
        groups[words[i][0]] = std::set<std::string>();
      }
      for (unsigned ig = 3; ig < words[i].size(); ++ig) {
        groups[words[i][0]].insert(words[i][ig]);
      }
      continue;
    }

    if (start_nuisance_scan && words[i].size()-1 == words[r].size()) {
      for (unsigned p = 2; p < words[i].size(); ++p) {
        if (words[i][p] == "-") continue;
        auto sys = std::make_shared<Systematic>();
        sys->set_bin(words[r-3][p-1]);
        try {
          int process_id = boost::lexical_cast<int>(words[r-2][p-1]);
          sys->set_signal(process_id <= 0);
          sys->set_process(words[r-1][p-1]);
        } catch(boost::bad_lexical_cast &) {
          int process_id = boost::lexical_cast<int>(words[r-1][p-1]);
          sys->set_signal(process_id <= 0);
          sys->set_process(words[r-2][p-1]);
        }
        sys->set_name(words[i][0]);
        std::string type = words[i][1];
        if (!contains(std::vector<std::string>{"shape", "shape?", "shapeN2", "lnN", "lnU"},
                      type)) {
          throw std::runtime_error(
              FNERROR("Systematic type " + type + " not supported"));
        }
        sys->set_type(words[i][1]);
        sys->set_analysis(analysis);
        sys->set_era(era);
        sys->set_channel(channel);
        sys->set_bin_id(bin_id);
        sys->set_mass(mass);
        sys->set_scale(1.0);
        std::size_t slash_pos = words[i][p].find("/");
        if (slash_pos != words[i][p].npos) {
          // Assume asymmetric of form kDown/kUp
          sys->set_value_d(
            boost::lexical_cast<double>(words[i][p].substr(0, slash_pos)));
          sys->set_value_u(
            boost::lexical_cast<double>(words[i][p].substr(slash_pos+1)));
          sys->set_asymm(true);
        } else {
          sys->set_value_u(boost::lexical_cast<double>(words[i][p]));
          sys->set_asymm(false);
        }
        if (sys->type() == "shape" || sys->type() == "shapeN2") {
          sys->set_scale(boost::lexical_cast<double>(words[i][p]));
          LoadShapes(sys.get(), hist_mapping);
        } else if (sys->type() == "shape?") {
          // This might fail, so we have to "try"
          try {
            LoadShapes(sys.get(), hist_mapping);
          } catch (std::exception & e) {
            LOGLINE(log(), "Systematic with shape? did not resolve to a shape");
            if (verbosity_ > 0) log() << e.what();
          }
          if (!sys->shape_u() || !sys->shape_d()) {
            sys->set_type("lnN");
          } else {
            sys->set_type("shape");
            sys->set_scale(boost::lexical_cast<double>(words[i][p]));
          }
        }
        if (sys->type() == "shape" || sys->type() == "shapeN2")
          sys->set_asymm(true);

        CombineHarvester::CreateParameterIfEmpty(sys->name());
        if (sys->type() == "lnU") {
          params_.at(sys->name())->set_err_d(0.);
          params_.at(sys->name())->set_err_u(0.);
        }
        systs_.push_back(sys);
      }
    }
  }
  if (single_obs) {
    if (bin_names.size() == 1) {
      single_obs->set_bin(*(bin_names.begin()));
      LoadShapes(single_obs.get(), hist_mapping);
      obs_.push_back(single_obs);
    } else {
      throw std::runtime_error(FNERROR(
          "Input card specifies a single observation entry without a bin, but "
          "multiple bins defined elsewhere"));
    }
  }

  // Finally populate the groups
  for (auto const& grp : groups) {
    this->SetGroup(grp.first, ch::Set2Vec(grp.second));
  }
  return 0;
}

void CombineHarvester::WriteDatacard(std::string const& name,
                                     std::string const& root_file) {
  TFile file(root_file.c_str(), "RECREATE");
  CombineHarvester::WriteDatacard(name, file);
  file.Close();
}

void CombineHarvester::WriteDatacard(std::string const& name) {
  TFile dummy;
  CombineHarvester::WriteDatacard(name, dummy);
  dummy.Close();
}

void CombineHarvester::FillHistMappings(std::vector<HistMapping> & mappings) {
  // Build maps of
  //  RooAbsData --> RooWorkspace
  //  RooAbsPdf --> RooWorkspace
  std::map<RooAbsData const*, RooWorkspace*> data_ws_map;
  std::map<RooAbsPdf const*, RooWorkspace*> pdf_ws_map;
  for (auto const& iter : wspaces_) {
    auto dat = iter.second->allData();
    for (auto d : dat) {
      data_ws_map[d] = iter.second.get();
    }
    RooArgSet vars = iter.second->allPdfs();
    auto v = vars.createIterator();
    do {
      RooAbsPdf *y = dynamic_cast<RooAbsPdf*>(**v);
      if (y) pdf_ws_map[iter.second->pdf(y->GetName())] = iter.second.get();
    } while (v->Next());
  }

  // For writing TH1s we will hard code a set of patterns for each bin
  // This assumes that the backgrounds will not depend on "mass" but the
  // signal will. Will probably want to change this in the future
  std::set<std::string> hist_bins;
  auto bins = this->bin_set();
  for (auto bin : bins) {
    unsigned shape_count = std::count_if(procs_.begin(), procs_.end(),
        [&](std::shared_ptr<ch::Process> p) {
          return (p->bin() == bin && p->shape() && (!p->signal()));
        });
    shape_count += std::count_if(obs_.begin(), obs_.end(),
        [&](std::shared_ptr<ch::Observation> p) {
          return (p->bin() == bin && p->shape());
        });
    unsigned counting = std::count_if(
        procs_.begin(), procs_.end(), [&](std::shared_ptr<ch::Process> p) {
          return (p->bin() == bin && p->shape() == nullptr &&
                  p->pdf() == nullptr && p->data() == nullptr);
        });
    counting += std::count_if(
        obs_.begin(), obs_.end(), [&](std::shared_ptr<ch::Observation> p) {
          return (p->bin() == bin && p->shape() == nullptr &&
                  p->data() == nullptr);
        });

    if (shape_count > 0) {
      mappings.emplace_back("*", bin, bin + "/$PROCESS",
                            bin + "/$PROCESS_$SYSTEMATIC");
      hist_bins.insert(bin);
    } else if (counting > 0) {
      mappings.emplace_back("*", bin, "", "");
      mappings.back().is_fake = true;
    }

    CombineHarvester ch_signals =
        std::move(this->cp().bin({bin}).signals().histograms());
    auto sig_proc_set =
        ch_signals.SetFromProcs(std::mem_fn(&ch::Process::process));
    for (auto sig_proc : sig_proc_set) {
      // should only add this mapping if the signal process has a numeric mass
      // value, otherwise we will write it using the background rule above
      auto masses = Set2Vec(ch_signals.cp()
                                .process({sig_proc})
                                .SetFromProcs(std::mem_fn(&ch::Process::mass)));
      if (masses.size() != 1) {
        throw std::runtime_error(FNERROR("Process " + sig_proc + " in bin " +
                                         bin +
                                         " has multiple entries with multiple "
                                         "mass values, this is not supported"));

      }
      if (is_float(masses[0])) {
        mappings.emplace_back(sig_proc, bin, bin + "/" + sig_proc + "$MASS",
                              bin + "/" + sig_proc + "$MASS_$SYSTEMATIC");
      }
    }
  }

  // Generate mappings for RooFit objects
  for (auto bin : bins) {
    CombineHarvester ch_bin = std::move(this->cp().bin({bin}));
    for (auto obs : ch_bin.obs_) {
      if (!obs->data()) continue;
      std::string obj_name = std::string(data_ws_map[obs->data()]->GetName()) +
                             ":" + std::string(obs->data()->GetName());
      mappings.emplace_back("data_obs", obs->bin(), obj_name, "");
    }

    bool prototype_ok = false;
    HistMapping prototype;
    std::vector<HistMapping> full_list;
    auto pmap = ch_bin.GenerateProcSystMap();
    for (unsigned i = 0; i < ch_bin.procs_.size(); ++i) {
      ch::Process * proc = ch_bin.procs_[i].get();
      if (!proc->data() && !proc->pdf()) continue;
      std::string obj_name;
      std::string obj_sys_name;
      if (proc->data()) {
        obj_name = std::string(data_ws_map[proc->data()]->GetName()) + ":" +
                   std::string(proc->data()->GetName());
      }
      if (proc->pdf()) {
        obj_name = std::string(pdf_ws_map[proc->pdf()]->GetName()) +
                   ":" + std::string(proc->pdf()->GetName());
      }
      for (unsigned j = 0; j < pmap[i].size(); ++j) {
        ch::Systematic const* sys = pmap[i][j];
        if (sys->data_u()) {
          obj_sys_name = std::string(data_ws_map[sys->data_u()]->GetName()) +
                         ":" + std::string(sys->data_u()->GetName());
          boost::replace_all(obj_sys_name, sys->name() + "Up", "$SYSTEMATIC");
          boost::replace_all(obj_sys_name, sys->process(), "$PROCESS");
          break;
        }
      }
      boost::replace_all(obj_name, proc->process(), "$PROCESS");

      // If the prototype pattern is already filled, but doesn't equal this
      // new pattern - then we can't use the prototype
      if (prototype.pattern.size() && prototype.pattern != obj_name) {
        prototype_ok = false;
      }

      if (prototype.syst_pattern.size() && obj_sys_name.size() &&
          prototype.syst_pattern != obj_sys_name) {
        prototype_ok = false;
      }

      if (!prototype.pattern.size()) {
        prototype_ok = true;
        prototype.process = "*";
        prototype.category = bin;
        prototype.pattern = obj_name;
      }
      if (!prototype.syst_pattern.size()) {
        prototype.syst_pattern = obj_sys_name;
      }

      full_list.emplace_back(proc->process(), proc->bin(), obj_name,
                             obj_sys_name);
    }
    // There are two reasons we won't want to write a generic mapping
    // for the processes in this bin:
    //    1) One or more processes is not described by the prototype mapping
    //    2) We also have a generic histogram mapping for this bin
    if (!prototype_ok || hist_bins.count(bin)) {
      for (auto m : full_list) {
        mappings.push_back(m);
      }
    } else {
      mappings.push_back(prototype);
    }
  }
}


void CombineHarvester::WriteDatacard(std::string const& name,
                                     TFile& root_file) {
  using boost::format;

  // First figure out if this is a counting-experiment only
  bool is_counting = true;
  this->ForEachObs([&](ch::Observation *obs) {
    if (obs->shape() != nullptr || obs->data() != nullptr) {
      is_counting = false;
    }
  });
  if (is_counting) {
    this->ForEachProc([&](ch::Process *proc) {
      if (proc->shape() != nullptr || proc->data() != nullptr ||
          proc->pdf() != nullptr) {
        is_counting = false;
      }
    });
  }

  // Allow a non-open ROOT file if this is purely a counting experiment
  if (!root_file.IsOpen() && !is_counting) {
    throw std::runtime_error(FNERROR(
        std::string("Output ROOT file is not open: ") + root_file.GetName()));
  }
  std::ofstream txt_file;
  txt_file.open(name);
  if (!txt_file.is_open()) {
    throw std::runtime_error(FNERROR("Unable to create file: " + name));
  }

  txt_file << "# Datacard produced by CombineHarvester with git status: "
           << ch::GitVersion() << "\n";

  std::string dashes(80, '-');

  auto bin_set = this->SetFromObs(std::mem_fn(&ch::Observation::bin));
  auto proc_set = this->SetFromProcs(std::mem_fn(&ch::Process::process));
  std::set<std::string> sys_set;
  std::set<std::string> rateparam_set;
  this->ForEachSyst([&](ch::Systematic const* sys) {
    if (sys->type() == "rateParam") {
      rateparam_set.insert(sys->name());
    } else {
      sys_set.insert(sys->name());
    }
  });
  txt_file << "imax    " << bin_set.size()
            << " number of bins\n";
  txt_file << "jmax    " << proc_set.size() - 1
           << " number of processes minus 1\n";
  txt_file << "kmax    " << "*" << " number of nuisance parameters\n";
  txt_file << dashes << "\n";


  std::vector<HistMapping> mappings;
  FillHistMappings(mappings);

  auto bins = this->SetFromObs(std::mem_fn(&ch::Observation::bin));

  auto proc_sys_map = this->GenerateProcSystMap();

  std::set<std::string> all_dependents_pars;
  for (auto proc : procs_) {
    if (!proc->pdf()) continue;
    // The rest of this is building the list of dependents
    /**
     * \todo This reproduces quite a lot of CombineHarvester::ImportParameters
     */
    RooAbsData const* data_obj = FindMatchingData(proc.get());
    RooArgSet par_list;
    RooArgSet norm_par_list;

    if (data_obj) {
      par_list.add(ParametersByName(proc->pdf(), data_obj->get()));
      if (proc->norm()) {
        norm_par_list.add(ParametersByName(proc->norm(), data_obj->get()));
      }
    } else {
      RooRealVar mx("CMS_th1x" , "CMS_th1x", 0, 1);
      RooArgSet tmp_set(mx);
      par_list.add(ParametersByName(proc->pdf(), &tmp_set));
      if (proc->norm()) {
        norm_par_list.add(ParametersByName(proc->norm(), &tmp_set));
      }
    }
    RooFIter par_list_it = par_list.fwdIterator();
    RooAbsArg *par_list_var = nullptr;
    while ((par_list_var = par_list_it.next())) {
      if (dynamic_cast<RooRealVar*>(par_list_var)) {
        all_dependents_pars.insert(par_list_var->GetName());
      }
    }
    if (proc->norm()) {
      RooFIter nm_list_it = norm_par_list.fwdIterator();
      RooAbsArg *nm_list_var = nullptr;
      while ((nm_list_var = nm_list_it.next())) {
        if (dynamic_cast<RooRealVar*>(nm_list_var)) {
          all_dependents_pars.insert(nm_list_var->GetName());
        }
      }
    }
  }

  // The ROOT file mapping should be given as a relative path
  std::string file_name = root_file.GetName();
  // Get the full path to the output root file
  // NOTE: was using canonical here instead of absolute, but not
  // supported in boost 1.47
  boost::filesystem::path root_file_path =
      boost::filesystem::absolute(file_name);
  // Get the full path to the directory containing the txt file
  boost::filesystem::path txt_file_path =
      boost::filesystem::absolute(name).parent_path();
  // Compute the relative path from the txt file to the root file
  file_name = make_relative(txt_file_path, root_file_path).string();

  for (auto const& mapping : mappings) {
    if (!mapping.is_fake) {
      txt_file << format("shapes %s %s %s %s %s\n")
        % mapping.process
        % mapping.category
        % file_name
        % mapping.pattern
        % mapping.syst_pattern;
    } else {
      txt_file << format("shapes %s %s %s\n")
        % mapping.process
        % mapping.category
        % "FAKE";
    }
  }
  txt_file << dashes << "\n";

  if (!is_counting) {
    for (auto ws_it : wspaces_) {
      if (ws_it.first == "_rateParams") continue; // don't write this one
      ch::WriteToTFile(ws_it.second.get(), &root_file, ws_it.second->GetName());
    }
  }

  // Writing observations
  txt_file << "bin          ";
  for (auto const& obs : obs_) {
    txt_file << format("%-15s ") % obs->bin();
    if (obs->shape()) {
      std::unique_ptr<TH1> h((TH1*)(obs->shape()->Clone()));
      h->Scale(obs->rate());
      WriteHistToFile(h.get(), &root_file, mappings, obs->bin(), "data_obs",
                      obs->mass(), "", 0);
    }
  }
  txt_file << "\n";
  txt_file << "observation  ";
  // On the precision of the observation yields: .1f is not sufficient for
  // combine to be happy if we have some asimov dataset with non-integer values.
  // We could just always give .4f but this doesn't look nice for the majority
  // of cards that have real data. Instead we'll check...
  std::string obs_fmt_int = "%-15.1f ";
  std::string obs_fmt_flt = "%-15.4f ";
  for (auto const& obs : obs_) {
    bool is_float =
        std::fabs(obs->rate() - std::round(obs->rate())) > 1E-4;
    txt_file << format(is_float ? obs_fmt_flt : obs_fmt_int)
        % obs->rate();
  }
  txt_file << "\n";
  txt_file << dashes << "\n";

  unsigned sys_str_len = 14;
  for (auto const& sys : sys_set) {
    if (sys.length() > sys_str_len) sys_str_len = sys.length();
  }
  for (auto const& sys : all_dependents_pars) {
    if (sys.length() > sys_str_len) sys_str_len = sys.length();
  }
  std::string sys_str_short = boost::lexical_cast<std::string>(sys_str_len);
  std::string sys_str_long = boost::lexical_cast<std::string>(sys_str_len+9);

  txt_file << format("%-"+sys_str_long+"s") % "bin";
  for (auto const& proc : procs_) {
    if (proc->shape()) {
      std::unique_ptr<TH1> h = proc->ClonedScaledShape();
      WriteHistToFile(h.get(), &root_file, mappings, proc->bin(),
                      proc->process(), proc->mass(), "", 0);
    }
    txt_file << format("%-15s ") % proc->bin();
  }
  txt_file << "\n";

  txt_file << format("%-"+sys_str_long+"s") % "process";

  for (auto const& proc : procs_) {
    txt_file << format("%-15s ") % proc->process();
  }
  txt_file << "\n";

  txt_file << format("%-"+sys_str_long+"s") % "process";

  // Setup process_ids first
  std::map<std::string, int> p_ids;
  unsigned sig_id = 0;
  unsigned bkg_id = 1;
  for (unsigned p = 0; p < procs_.size(); ++p) {
    if (!procs_[p]->signal() && p_ids.count(procs_[p]->process()) == 0) {
      p_ids[procs_[p]->process()] = bkg_id;
      ++bkg_id;
    }
    unsigned q = procs_.size() - (p + 1);
    if (procs_[q]->signal() && p_ids.count(procs_[q]->process()) == 0) {
      p_ids[procs_[q]->process()] = sig_id;
      --sig_id;
    }
  }
  for (auto const& proc : procs_) {
    txt_file << format("%-15s ") % p_ids[proc->process()];
  }
  txt_file << "\n";


  txt_file << format("%-"+sys_str_long+"s") % "rate";
  for (auto const& proc : procs_) {
    txt_file << format("%-15.6g ") % proc->no_norm_rate();
  }
  txt_file << "\n";
  txt_file << dashes << "\n";

  // Need to write parameters here that feature both in the list of pdf
  // dependents and sys_set.
  for (auto par : params_) {
    Parameter const* p = par.second.get();
    if (p->err_d() != 0.0 && p->err_u() != 0.0 &&
        all_dependents_pars.count(p->name()) && sys_set.count(p->name())) {
      txt_file << format((format("%%-%is param %%g %%g") % sys_str_len).str()) %
                      p->name() % p->val() % ((p->err_u() - p->err_d()) / 2.0);
      if (p->range_d() != std::numeric_limits<double>::lowest() &&
          p->range_u() != std::numeric_limits<double>::max()) {
        txt_file << format(" [%.4g,%.4g]") % p->range_d() % p->range_u();
      }
      txt_file << "\n";
    }
  }

  for (auto const& sys : sys_set) {
    std::vector<std::string> line(procs_.size() + 2);
    line[0] = sys;
    bool seen_lnN = false;
    bool seen_lnU = false;
    bool seen_shape = false;
    bool seen_shapeN2 = false;
    for (unsigned p = 0; p < procs_.size(); ++p) {
      line[p + 2] = "-";
      for (unsigned n = 0; n < proc_sys_map[p].size(); ++n) {
        ch::Systematic const* ptr = proc_sys_map[p][n];
        // Only interested into Systematics with name "sys"
        if (ptr->name() != sys) continue;
        std::string tp = ptr->type();
        if (tp == "lnN" || tp == "lnU") {
          if (tp == "lnN") seen_lnN = true;
          if (tp == "lnU") seen_lnU = true;
          line[p + 2] =
              ptr->asymm()
                  ? (format("%g/%g") % ptr->value_d() % ptr->value_u()).str()
                  : (format("%g") % ptr->value_u()).str();
          break;
        }
        if (tp == "shape" || tp == "shapeN2") {
          if (tp == "shape") seen_shape = true;
          if (tp == "shapeN2") seen_shapeN2 = true;
          line[p + 2] = (format("%g") % ptr->scale()).str();
          if (ptr->shape_u() && ptr->shape_d()) {
            bool add_dir = TH1::AddDirectoryStatus();
            TH1::AddDirectory(false);
            std::unique_ptr<TH1> h_d = ptr->ClonedShapeD();
            h_d->Scale(procs_[p]->rate() * ptr->value_d());
            WriteHistToFile(h_d.get(), &root_file, mappings, ptr->bin(),
                            ptr->process(), ptr->mass(), ptr->name(), 1);
            std::unique_ptr<TH1> h_u = ptr->ClonedShapeU();
            h_u->Scale(procs_[p]->rate() * ptr->value_u());
            WriteHistToFile(h_u.get(), &root_file, mappings, ptr->bin(),
                            ptr->process(), ptr->mass(), ptr->name(), 2);
            TH1::AddDirectory(add_dir);
            break;
          } else if (ptr->data_u() && ptr->data_d()) {
          } else {
            if (!flags_.at("allow-missing-shapes")) {
              std::stringstream err;
              err << "Trying to write shape uncertainty with missing "
                     "shapes:\n";
              err << Systematic::PrintHeader << *ptr;
              throw std::runtime_error(FNERROR(err.str()));
            }
          }
        }
      }
    }
    if (seen_shapeN2) {
      line[1] = "shapeN2";
    } else if (seen_lnU) {
      line[1] = "lnU";
    } else if (seen_lnN && !seen_shape) {
      line[1] = "lnN";
    } else if (!seen_lnN && seen_shape) {
      line[1] = "shape";
    } else if (seen_lnN && seen_shape) {
      line[1] = "shape?";
    } else {
      throw std::runtime_error(FNERROR("Systematic type could not be deduced"));
    }
    txt_file << format("%-" + sys_str_short + "s %-7s ") % line[0] % line[1];
    for (unsigned p = 0; p < procs_.size(); ++p) {
      txt_file << format("%-15s ") % line[p + 2];
    }
    txt_file << "\n";
  }
  // write param line for any parameter which has a non-zero error
  // and which doesn't appear in list of nuisances
  CombineHarvester ch_rp = this->cp().syst_type({"rateParam"});
  // These are all the processes with a rateParam. We'll look at each one and
  // see if every bin with this process also has this rateParam.
  std::vector<std::vector<std::string> > floating_params;
  for (auto const& sys : rateparam_set) {
    FNLOGC(log(), verbosity_ > 1) << "Doing rateParam: " << sys << "\n";
    auto procs_rp = ch_rp.cp().syst_name({sys}).SetFromSysts(
        std::mem_fn(&ch::Systematic::process));
    FNLOGC(log(), verbosity_ > 1) << "Procs for this rateParam: \n";
    for (auto const& proc : procs_rp) {
      FNLOGC(log(), verbosity_ > 1) << "  - " << proc << "\n";
      // Get the set of bins
      auto bins_rp = ch_rp.cp().process({proc}).syst_name({sys}).SetFromSysts(
          std::mem_fn(&ch::Systematic::bin));
      auto bins_all = this->cp().process({proc}).bin_set();
      if (bins_rp.size() > 0 && bins_rp.size() == bins_all.size()) {
        floating_params.push_back({sys, "*", proc});
      } else {
        for (auto const& bin : bins_rp) {
          floating_params.push_back({sys, bin, proc});
        }
      }
      FNLOGC(log(), verbosity_ > 1)
          << bins_rp.size() << "/" << bins_all.size()
          << " bins with this process have a rateParam\n";
    }
  }
  for (auto const& rp : floating_params) {
    if (params_.count(rp[0])) {
      ch::Parameter const* par = params_.at(rp[0]).get();
      txt_file << format("%-" + sys_str_short + "s %-10s %-10s %-10s %g") %
                      rp[0] % "rateParam" % rp[1] % rp[2] % par->val();
      // Check if we have to set the range too
      // TODO: should package this into some function in ch::Parameter
      if (par->range_d() != std::numeric_limits<double>::lowest() &&
          par->range_u() != std::numeric_limits<double>::max()) {
        txt_file << format(" [%.4g,%.4g]") % par->range_d() % par->range_u();
      }
      txt_file << "\n";
    } else {
      // If we don't have a ch::Parameter with this name, then we'll assume
      // this is a function
      RooWorkspace *rp_ws = wspaces_.at("_rateParams").get();
      RooFormulaVar* form =
          dynamic_cast<RooFormulaVar*>(rp_ws->function(rp[0].c_str()));
      // Want to extract the formula string, apparently "printMetaArgs" is the
      // only way to do it
      std::stringstream ss;
      form->printMetaArgs(ss);
      std::string form_str = ss.str();
      // The string will be like 'formula="blah"', strip everything on either
      // side of the double quotes (including the quotes themselves)
      std::size_t first = form_str.find_first_of('"');
      std::size_t last = form_str.find_last_of('"');
      form_str = form_str.substr(first + 1, last - first - 1);
      // unique_ptr because we own the RooArgSet that is generated
      unsigned npars =
          std::unique_ptr<RooArgSet>(form->getParameters(RooArgList()))->getSize();
      std::string args = "";
      for (unsigned i = 0; i < npars; ++i) {
        args += std::string(form->getParameter(i)->GetName());
        if (i < (npars-1)) {
          args += ",";
        }
      }
      txt_file << format("%-" + sys_str_short +
                         "s %-10s %-10s %-10s %-20s %s\n") %
                      rp[0] % "rateParam" % rp[1] % rp[2] % form_str % args;
    }
  }

  std::set<std::string> ws_vars;
  for (auto iter : wspaces_) {
    RooArgSet vars = iter.second->allVars();
    auto v = vars.createIterator();
    do {
      RooRealVar *y = dynamic_cast<RooRealVar*>(**v);
      if (y) {
        ws_vars.insert(y->GetName());
      }
    } while (v->Next());
  }

  // How to check for params we need to write:
  // Get a list of all pdf dependents for the processes we've just written
  // Also need a list of all Systematic entries we've just written
  // Then we write param if:
  //  - it has non-zero errors
  //  - it appears in the first list
  //  - it doesn't appear in the second list
  for (auto par : params_) {
    Parameter const* p = par.second.get();
    if (p->err_d() != 0.0 && p->err_u() != 0.0 &&
        all_dependents_pars.count(p->name()) && !sys_set.count(p->name())) {
      txt_file << format((format("%%-%is param %%g %%g") % sys_str_len).str()) %
                      p->name() % p->val() % ((p->err_u() - p->err_d()) / 2.0);
      if (p->range_d() != std::numeric_limits<double>::lowest() &&
          p->range_u() != std::numeric_limits<double>::max()) {
        txt_file << format(" [%.4g,%.4g]") % p->range_d() % p->range_u();
      }
      txt_file << "\n";
    }
  }

  // Finally write the NP groups
  // We should only consider parameters are either:
  // - in the list of pdf dependents
  // - in the list of systematics
  std::map<std::string, std::set<std::string>> groups;
  for (auto const& par : params_) {
    Parameter * p = par.second.get();
    if (!all_dependents_pars.count(p->name()) && !sys_set.count(p->name())) {
      continue;
    }
    for (auto const& str : p->groups()) {
      if (!groups.count(str)) {
        groups[str] = std::set<std::string>();
      }
      groups[str].insert(p->name());
    }
  }
  for (auto const& gr : groups) {
    txt_file << format("%s group =") % gr.first;
    for (auto const& np : gr.second)  {
      txt_file << " " << np;
    }
    txt_file << "\n";
  }

  txt_file.close();
}

void CombineHarvester::WriteHistToFile(
    TH1 * hist,
    TFile * file,
    std::vector<HistMapping> const& mappings,
    std::string const& bin,
    std::string const& process,
    std::string const& mass,
    std::string const& nuisance,
    unsigned type) {
  StrPairVec attempts = this->GenerateShapeMapAttempts(process, bin);
  for (unsigned a = 0; a < attempts.size(); ++a) {
    for (unsigned m = 0; m < mappings.size(); ++m) {
      if ((attempts[a].first == mappings[m].process) &&
        (attempts[a].second == mappings[m].category)) {
        std::string p = (type == 0 ?
            mappings[m].pattern : mappings[m].syst_pattern);
        boost::replace_all(p, "$CHANNEL", bin);
        boost::replace_all(p, "$PROCESS", process);
        boost::replace_all(p, "$MASS", mass);
        if (type == 1) boost::replace_all(p, "$SYSTEMATIC", nuisance+"Down");
        if (type == 2) boost::replace_all(p, "$SYSTEMATIC", nuisance+"Up");
        WriteToTFile(hist, file, p);
        return;
      }
    }
  }
}
}
