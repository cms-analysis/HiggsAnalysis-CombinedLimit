#include "CombineTools/interface/CombineHarvester.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <utility>
#include <set>
#include <fstream>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/range/algorithm_ext/erase.hpp"
#include "boost/range/algorithm/find.hpp"
#include "boost/format.hpp"
#include "boost/regex.hpp"
#include "TDirectory.h"
#include "TH1.h"
#include "RooRealVar.h"
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/MakeUnique.h"
#include "CombineTools/interface/Utilities.h"
#include "CombineTools/interface/TFileIO.h"
#include "CombineTools/interface/algorithm.h"

namespace ch {

// Extract info from filename using parse rule like:
// ".*{MASS}/{ANALYSIS}_{CHANNEL}_{BINID}_{ERA}.txt"
int CombineHarvester::ParseDatacard(std::string const& filename,
    std::string parse_rules) {
  boost::replace_all(parse_rules, "$ANALYSIS",  "(?<ANALYSIS>\\w+)");
  boost::replace_all(parse_rules, "$ERA",       "(?<ERA>\\w+)");
  boost::replace_all(parse_rules, "$CHANNEL",   "(?<CHANNEL>\\w+)");
  boost::replace_all(parse_rules, "$BINID",     "(?<BINID>\\w+)");
  boost::replace_all(parse_rules, "$MASS",      "(?<MASS>\\w+)");
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
  std::map<std::string, RooAbsData*> data_map;
  bool start_nuisance_scan = false;
  unsigned r = 0;

  // Loop through the vector of word vectors
  for (unsigned i = 0; i < words.size(); ++i) {
    // Ignore line if it only has one word
    if (words[i].size() <= 1) continue;

    // If the line begins "shapes" then we've
    // found process --> TH1 mapping information
    if (words[i][0] == "shapes" && words[i].size() >= 5) {
      hist_mapping.push_back(HistMapping());
      hist_mapping.back().process = words[i][1];
      hist_mapping.back().category = words[i][2];
      // The root file path given in the datacard is relative to the datacard
      // path, so we join the path to the datacard with the path to the file
      std::string dc_path;
      std::size_t slash = filename.find_last_of('/');
      if (slash != filename.npos) {
        dc_path = filename.substr(0, slash) + "/" + words[i][3];
      } else {
        dc_path = words[i][3];
      }
      hist_mapping.back().file = std::make_shared<TFile>(dc_path.c_str());
      hist_mapping.back().pattern = words[i][4];
      if (words[i].size() > 5) hist_mapping.back().syst_pattern = words[i][5];
    }

    // Want to check this line and the previous one, so need i >= 1.
    // If the first word on this line is "observation" and "bin" on
    // the previous line then we've found the entries for data, and
    // can add Observation objects
    if (i >= 1) {
      if (  words[i][0]     == "observation" &&
            words[i-1][0]   == "bin" &&
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

    // Similarly look for the lines indicating the different signal
    // and background processes
    // Once these are found save in line index for the rate line as r
    // to we can refer back to these later, then assume that every
    // line that follows is a nuisance parameter
    if (i >= 3) {
      if (  words[i][0]   == "rate" &&
            words[i-1][0] == "process" &&
            words[i-2][0] == "process" &&
            words[i-3][0] == "bin" &&
            words[i].size() == words[i-1].size() &&
            words[i].size() == words[i-2].size() &&
            words[i].size() == words[i-3].size()) {
        for (unsigned p = 1; p < words[i].size(); ++p) {
          auto proc = std::make_shared<Process>();
          proc->set_bin(words[i-3][p]);
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
      if (words[i][1] == "param") {
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
        continue;  // skip the rest of this now
      }
    }

    if (start_nuisance_scan && words[i].size()-1 == words[r].size()) {
      for (unsigned p = 2; p < words[i].size(); ++p) {
        if (words[i][p] == "-") continue;
        auto nus = std::make_shared<Nuisance>();
        nus->set_bin(words[r-3][p-1]);
        try {
          int process_id = boost::lexical_cast<int>(words[r-2][p-1]);
          nus->set_signal(process_id <= 0);
          nus->set_process(words[r-1][p-1]);
        } catch(boost::bad_lexical_cast &) {
          int process_id = boost::lexical_cast<int>(words[r-1][p-1]);
          nus->set_signal(process_id <= 0);
          nus->set_process(words[r-2][p-1]);
        }
        nus->set_name(words[i][0]);
        std::string type = words[i][1];
        if (!contains(std::vector<std::string>{"shape", "shape?", "lnN"},
                      type)) {
          throw std::runtime_error(
              "[CombineHarvester::ParseDatacard] Systematic type " + type +
              " not supported");
        }
        nus->set_type(words[i][1]);
        nus->set_analysis(analysis);
        nus->set_era(era);
        nus->set_channel(channel);
        nus->set_bin_id(bin_id);
        nus->set_mass(mass);
        nus->set_scale(1.0);
        std::size_t slash_pos = words[i][p].find("/");
        if (slash_pos != words[i][p].npos) {
          // Assume asymmetric of form kDown/kUp
          nus->set_value_d(
            boost::lexical_cast<double>(words[i][p].substr(0, slash_pos)));
          nus->set_value_u(
            boost::lexical_cast<double>(words[i][p].substr(slash_pos+1)));
          nus->set_asymm(true);
        } else {
          nus->set_value_u(boost::lexical_cast<double>(words[i][p]));
          nus->set_asymm(false);
        }
        if (nus->type() == "shape") {
          nus->set_scale(boost::lexical_cast<double>(words[i][p]));
          LoadShapes(nus.get(), hist_mapping);
        } else if (nus->type() == "shape?") {
          // This might fail, so we have to "try"
          try {
            LoadShapes(nus.get(), hist_mapping);
          } catch (std::exception & e) {
            LOGLINE(log(), "Nuisance with shape? did not resolve to a shape");
            if (verbosity_ > 0) log() << e.what();
          }
          if (!nus->shape_u() || !nus->shape_d()) {
            nus->set_type("lnN");
          } else {
            nus->set_type("shape");
            nus->set_scale(boost::lexical_cast<double>(words[i][p]));
          }
        }

        CombineHarvester::CreateParameterIfEmpty(this, nus->name());
        nus_.push_back(nus);
      }
    }
  }
  return 0;
}

void CombineHarvester::WriteDatacard(std::string const& name,
                                     std::string const& root_file) {
  TFile file(root_file.c_str(), "RECREATE");
  CombineHarvester::WriteDatacard(name, root_file);
  file.Close();
}

void CombineHarvester::WriteDatacard(std::string const& name,
                                     TFile& root_file) {
  std::ofstream txt_file;
  txt_file.open(name);

  std::string dashes(80, '-');

  auto bin_set =
      this->GenerateSetFromObs<std::string>(std::mem_fn(&ch::Observation::bin));
  auto proc_set =
      this->GenerateSetFromProcs<std::string>(std::mem_fn(&ch::Process::process));
  auto nus_set =
      this->GenerateSetFromNus<std::string>(std::mem_fn(&ch::Nuisance::name));
  txt_file << "imax    " << bin_set.size()
            << " number of bins\n";
  txt_file << "jmax    " << proc_set.size() - 1
           << " number of processes minus 1\n";
  txt_file << "kmax    " << "*" << " number of nuisance parameters\n";
  txt_file << dashes << "\n";


  std::vector<HistMapping> mappings;
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
      if (y) {
        pdf_ws_map[iter.second->pdf(y->GetName())] = iter.second.get();
      }
    } while (v->Next());
  }

  auto bins =
      this->GenerateSetFromObs<std::string>(std::mem_fn(&ch::Observation::bin));
  for (auto bin : bins) {
    unsigned shape_count = std::count_if(procs_.begin(), procs_.end(),
        [&](std::shared_ptr<ch::Process> p) {
          return (p->bin() == bin && p->shape() && (!p->signal()));
        });
    shape_count += std::count_if(obs_.begin(), obs_.end(),
        [&](std::shared_ptr<ch::Observation> p) {
          return (p->bin() == bin && p->shape());
        });

    if (shape_count > 0) {
      mappings.push_back({
        "*", bin, nullptr, bin+"/$PROCESS", bin+"/$PROCESS_$SYSTEMATIC"
      });
    }

    CombineHarvester ch_signals =
        std::move(this->cp().bin({bin}).signals().histograms());
    auto sig_proc_set = ch_signals.GenerateSetFromProcs<std::string>(
      std::mem_fn(&ch::Process::process));
    for (auto sig_proc : sig_proc_set) {
      mappings.push_back({sig_proc, bin, nullptr,
                          bin + "/" + sig_proc + "$MASS",
                          bin + "/" + sig_proc + "$MASS_$SYSTEMATIC"});
    }
  }

  for (auto obs : obs_) {
    if (obs->data()) {
      std::string obj_name = std::string(data_ws_map[obs->data()]->GetName()) +
                             ":" + std::string(obs->data()->GetName());
      mappings.push_back({
        "data_obs", obs->bin(), nullptr, obj_name, ""
      });
    }
  }
  for (auto proc : procs_) {
    if (proc->pdf()) {
      std::string obj_name = std::string(pdf_ws_map[proc->pdf()]->GetName()) +
                             ":" + std::string(proc->pdf()->GetName());
      mappings.push_back({
        proc->process(), proc->bin(), nullptr, obj_name, ""
      });
    }
  }

  std::string file_name = root_file.GetName();
  std::size_t file_slash_pos = file_name.find_last_of('/');
  if (file_slash_pos != file_name.npos) {
    file_name = file_name.substr(file_slash_pos+1);
  }

  for (auto const& mapping : mappings) {
    txt_file << boost::format("shapes %s %s %s %s %s\n")
      % mapping.process
      % mapping.category
      % file_name
      % mapping.pattern
      % mapping.syst_pattern;
  }
  txt_file << dashes << "\n";

  for (auto ws_it : wspaces_) {
    ch::WriteToTFile(ws_it.second.get(), &root_file, ws_it.second->GetName());
  }

  // Writing observations
  txt_file << "bin          ";
  for (auto const& obs : obs_) {
    txt_file << boost::format("%-15s ") % obs->bin();
    if (obs->shape()) {
      std::unique_ptr<TH1> h((TH1*)(obs->shape()->Clone()));
      h->Scale(obs->rate());
      WriteHistToFile(h.get(), &root_file, mappings, obs->bin(), "data_obs",
                      obs->mass(), "", 0);
    }
  }
  txt_file << "\n";
  txt_file << "observation  ";
  for (auto const& obs : obs_) {
    txt_file << boost::format("%-15.1f ") % obs->rate();
  }
  txt_file << "\n";
  txt_file << dashes << "\n";

  auto proc_nus_map = this->GenerateProcNusMap();
  unsigned nus_str_len = 14;
  for (auto const& nus : nus_set) {
    if (nus.length() > nus_str_len) nus_str_len = nus.length();
  }
  std::string nus_str_short = boost::lexical_cast<std::string>(nus_str_len);
  std::string nus_str_long = boost::lexical_cast<std::string>(nus_str_len+9);

  txt_file << boost::format("%-"+nus_str_long+"s") % "bin";
  for (auto const& proc : procs_) {
    if (proc->shape()) {
      std::unique_ptr<TH1> h((TH1*)(proc->shape()->Clone()));
      h->Scale(proc->rate());
      WriteHistToFile(h.get(), &root_file, mappings, proc->bin(),
                      proc->process(), proc->mass(), "", 0);
    }
    txt_file << boost::format("%-15s ") % proc->bin();
  }
  txt_file << "\n";

  txt_file << boost::format("%-"+nus_str_long+"s") % "process";

  for (auto const& proc : procs_) {
    txt_file << boost::format("%-15s ") % proc->process();
  }
  txt_file << "\n";

  txt_file << boost::format("%-"+nus_str_long+"s") % "process";

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
    txt_file << boost::format("%-15s ") % p_ids[proc->process()];
  }
  txt_file << "\n";


  txt_file << boost::format("%-"+nus_str_long+"s") % "rate";
  for (auto const& proc : procs_) {
    txt_file << boost::format("%-15.4f ") % proc->no_norm_rate();
  }
  txt_file << "\n";
  txt_file << dashes << "\n";

  for (auto const& nus : nus_set) {
    std::vector<std::string> line(procs_.size() + 2);
    line[0] = nus;
    bool seen_lnN = false;
    bool seen_shape = false;
    for (unsigned p = 0; p < procs_.size(); ++p) {
      line[p+2] = "-";
      for (unsigned n = 0; n < proc_nus_map[p].size(); ++n) {
        ch::Nuisance const* nus_ptr = proc_nus_map[p][n];
        if (nus_ptr->name() == nus) {
          if (nus_ptr->type() == "lnN") {
            seen_lnN = true;
            line[p + 2] =
                nus_ptr->asymm()
                    ? (boost::format("%g/%g") % nus_ptr->value_d() %
                       nus_ptr->value_u()).str()
                    : (boost::format("%g") % nus_ptr->value_u()).str();
            break;
          }
          if (nus_ptr->type() == "shape") {
            if (!nus_ptr->shape_u() || !nus_ptr->shape_d()) {
              std::stringstream err;
              err << "Trying to write shape uncertainty with missing shapes:\n";
              err << Nuisance::PrintHeader << *nus_ptr;
              throw std::runtime_error(FNERROR(err.str()));
            }
            bool add_dir = TH1::AddDirectoryStatus();
            TH1::AddDirectory(false);
            std::unique_ptr<TH1> h_d(
                static_cast<TH1*>(nus_ptr->shape_d()->Clone()));
            h_d->Scale(procs_[p]->rate()*nus_ptr->value_d());
            WriteHistToFile(h_d.get(), &root_file, mappings, nus_ptr->bin(),
                            nus_ptr->process(), nus_ptr->mass(),
                            nus_ptr->name(), 1);
            std::unique_ptr<TH1> h_u(
                static_cast<TH1*>(nus_ptr->shape_u()->Clone()));
            h_u->Scale(procs_[p]->rate()*nus_ptr->value_u());
            WriteHistToFile(h_u.get(), &root_file, mappings, nus_ptr->bin(),
                            nus_ptr->process(), nus_ptr->mass(),
                            nus_ptr->name(), 2);
            seen_shape = true;
            line[p+2] = (boost::format("%g") % nus_ptr->scale()).str();
            TH1::AddDirectory(add_dir);
            break;
          }
        }
      }
    }
    if (seen_lnN && !seen_shape) line[1] = "lnN";
    if (!seen_lnN && seen_shape) line[1] = "shape";
    if (seen_lnN && seen_shape) line[1] = "shape?";
    txt_file << boost::format(
      "%-"+nus_str_short+"s %-7s ")
      % line[0] % line[1];
    for (unsigned p = 0; p < procs_.size(); ++p) {
      txt_file << boost::format("%-15s ") % line[p+2];
    }
    txt_file << "\n";
  }
  // write param line for any parameter which has a non-zero error
  // and which doesn't appear in list of nuisances


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


  for (auto par : params_) {
    // if (par.second->err_d() != 0.0 && par.second->err_u() != 0.0 &&
    //     !nus_set.count(par.second->name())) {
    if (par.second->err_d() != 0.0 && par.second->err_u() != 0.0 &&
        ws_vars.count(par.first)) {
      txt_file << boost::format(
        "%-"+boost::lexical_cast<std::string>(nus_str_len)+"s param %g %g\n")
        % par.second->name() % par.second->val() %
        ((par.second->err_u()-par.second->err_d())/2.0);
    }
  }
  txt_file.close();
}

void CombineHarvester::WriteHistToFile(
    TH1 const* hist,
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
