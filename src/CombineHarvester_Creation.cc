#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <iostream>
#include <utility>
#include <algorithm>
#include "TDirectory.h"
#include "TH1.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"

namespace ch {
void CombineHarvester::AddObservations(
    std::vector<std::string> mass,
    std::vector<std::string> analysis,
    std::vector<std::string> era,
    std::vector<std::string> channel,
    ch::Categories bin) {
  std::vector<unsigned> lengths = {
      unsigned(mass.size()),
      unsigned(analysis.size()),
      unsigned(era.size()),
      unsigned(channel.size()),
      unsigned(bin.size())};
  auto comb = ch::GenerateCombinations(lengths);
  for (auto const& c : comb) {
    auto obs = std::make_shared<Observation>();
    obs->set_mass(mass[c[0]]);
    obs->set_analysis(analysis[c[1]]);
    obs->set_era(era[c[2]]);
    obs->set_channel(channel[c[3]]);
    obs->set_bin_id(bin[c[4]].first);
    obs->set_bin(bin[c[4]].second);
    obs_.push_back(obs);
  }
}

void CombineHarvester::AddProcesses(
    std::vector<std::string> mass,
    std::vector<std::string> analysis,
    std::vector<std::string> era,
    std::vector<std::string> channel,
    std::vector<std::string> procs,
    ch::Categories bin,
    bool signal) {
  std::vector<unsigned> lengths = {
      unsigned(mass.size()),
      unsigned(analysis.size()),
      unsigned(era.size()),
      unsigned(channel.size()),
      unsigned(bin.size())};
  auto comb = ch::GenerateCombinations(lengths);
  for (auto const& c : comb) {
    for (unsigned i = 0; i < procs.size(); ++i) {
      auto proc = std::make_shared<Process>();
      proc->set_mass(mass[c[0]]);
      proc->set_analysis(analysis[c[1]]);
      proc->set_era(era[c[2]]);
      proc->set_channel(channel[c[3]]);
      proc->set_bin_id(bin[c[4]].first);
      proc->set_bin(bin[c[4]].second);
      proc->set_process(procs[i]);
      proc->set_signal(signal);
      procs_.push_back(proc);
    }
  }
}

void CombineHarvester::AddSystFromProc(Process const& proc,
                                       std::string const& name,
                                       std::string const& type, bool asymm,
                                       double val_u, double val_d,
                                       std::string const& formula,
                                       std::string const& args) {
  std::string subbed_name = name;
  boost::replace_all(subbed_name, "$BINID",
                     boost::lexical_cast<std::string>(proc.bin_id()));
  boost::replace_all(subbed_name, "$BIN", proc.bin());
  boost::replace_all(subbed_name, "$PROCESS", proc.process());
  boost::replace_all(subbed_name, "$MASS", proc.mass());
  boost::replace_all(subbed_name, "$ERA", proc.era());
  boost::replace_all(subbed_name, "$CHANNEL", proc.channel());
  boost::replace_all(subbed_name, "$ANALYSIS", proc.analysis());
  std::map<std::string,std::string> attrs = proc.all_attributes();
  for( const auto it : attrs){
      boost::replace_all(subbed_name, "$ATTR("+it.first+")",proc.attribute(it.first));
  }
  auto sys = std::make_shared<Systematic>();
  ch::SetProperties(sys.get(), &proc);
  sys->set_name(subbed_name);
  sys->set_type(type);
  if (type == "lnN" || type == "lnU") {
    sys->set_asymm(asymm);
    sys->set_value_u(val_u);
    sys->set_value_d(val_d);
    CreateParameterIfEmpty(sys->name());
  } else if (type == "shape" || type == "shapeN2" || type == "shapeU") {
    sys->set_asymm(true);
    sys->set_value_u(1.0);
    sys->set_value_d(1.0);
    sys->set_scale(val_u);
    CreateParameterIfEmpty(sys->name());
  } else if (type == "rateParam") {
    sys->set_asymm(false);
    if (formula == "" && args == "") {
      SetupRateParamVar(subbed_name, val_u);
    } else {
      std::string subbed_args = args;
      boost::replace_all(subbed_args, "$BINID",
                         boost::lexical_cast<std::string>(proc.bin_id()));
      boost::replace_all(subbed_args, "$BIN", proc.bin());
      boost::replace_all(subbed_args, "$PROCESS", proc.process());
      boost::replace_all(subbed_args, "$MASS", proc.mass());
      boost::replace_all(subbed_args, "$ERA", proc.era());
      boost::replace_all(subbed_args, "$CHANNEL", proc.channel());
      boost::replace_all(subbed_args, "$ANALYSIS", proc.analysis());
      for( const auto it : attrs){
          boost::replace_all(subbed_args, "$ATTR("+it.first+")",proc.attribute(it.first));
      }
      SetupRateParamFunc(subbed_name, formula, subbed_args);
    }
  }
  if (sys->type() == "lnU" || sys->type() == "shapeU") {
    params_.at(sys->name())->set_err_d(0.);
    params_.at(sys->name())->set_err_u(0.);
  }
  systs_.push_back(sys);
}

void CombineHarvester::AddSystVar(std::string const& name,
                                  double val, double err) {
  Parameter * param = SetupRateParamVar(name,val);
  param->set_val(val);
  param->set_err_u(+1.0 *err);
  param->set_err_d(-1.0 *err);
  auto sys = std::make_shared<Systematic>();
  sys->set_name(name);
  sys->set_type("param");
  systs_.push_back(sys);
}

void CombineHarvester::RenameSystematic(CombineHarvester &target, std::string const& old_name,
                                        std::string const& new_name) {
 for(unsigned i = 0; i<systs_.size(); ++i){
    if(systs_[i]->name()==old_name){
      systs_[i]->set_name(new_name);
      target.CreateParameterIfEmpty(systs_[i]->name());
    }
  }
}

void CombineHarvester::ExtractShapes(std::string const& file,
                                     std::string const& rule,
                                     std::string const& syst_rule) {
  std::vector<HistMapping> mapping(1);
  mapping[0].process = "*";
  mapping[0].category = "*";
  mapping[0].file = std::make_shared<TFile>(file.c_str());
  mapping[0].pattern = rule;
  mapping[0].syst_pattern = syst_rule;

  // Note that these LoadShapes calls will fail if we encounter
  // any object that already has shapes
  for (unsigned  i = 0; i < obs_.size(); ++i) {
    if (obs_[i]->shape() || obs_[i]->data()) continue;
    LoadShapes(obs_[i].get(), mapping);
  }
  for (unsigned  i = 0; i < procs_.size(); ++i) {
    if (procs_[i]->shape() || procs_[i]->pdf()) continue;
    LoadShapes(procs_[i].get(), mapping);
  }
  if (syst_rule == "") return;
  for (unsigned  i = 0; i < systs_.size(); ++i) {
    if (systs_[i]->type() != "shape" && systs_[i]->type() != "shapeN2" &&
        systs_[i]->type() != "shapeU")
      continue;
    LoadShapes(systs_[i].get(), mapping);
  }
}

void CombineHarvester::AddWorkspace(RooWorkspace const& ws,
                                    bool can_rename) {
  SetupWorkspace(ws, can_rename);
}

void CombineHarvester::ExtractPdfs(CombineHarvester& target,
                                   std::string const& ws_name,
                                   std::string const& rule,
                                   std::string norm_rule) {
  std::vector<HistMapping> mapping(1);
  mapping[0].process = "*";
  mapping[0].category = "*";
  mapping[0].pattern = ws_name+":"+rule;
  if (norm_rule != "") mapping[0].syst_pattern = ws_name + ":" + norm_rule;
  if (!wspaces_.count(ws_name)) return;
  mapping[0].ws = wspaces_.at(ws_name);
  for (unsigned  i = 0; i < procs_.size(); ++i) {
    if (!procs_[i]->pdf()) {
      target.LoadShapes(procs_[i].get(), mapping);
    }
  }
}

void CombineHarvester::ExtractData(std::string const &ws_name,
                                   std::string const &rule) {
  std::vector<HistMapping> mapping(1);
  mapping[0].process = "*";
  mapping[0].category = "*";
  mapping[0].pattern = ws_name+":"+rule;
  if (!wspaces_.count(ws_name)) return;
  mapping[0].ws = wspaces_.at(ws_name);
  for (unsigned  i = 0; i < obs_.size(); ++i) {
    if (!obs_[i]->data()) {
      LoadShapes(obs_[i].get(), mapping);
    }
  }
}

void CombineHarvester::AddBinByBin(double threshold, bool fixed_norm,
                                   CombineHarvester &other) {
  AddBinByBin(threshold, fixed_norm, &other);
}

void CombineHarvester::AddBinByBin(double threshold, bool fixed_norm,
                                   CombineHarvester *other) {
  auto bbb_factory = ch::BinByBinFactory()
      .SetAddThreshold(threshold)
      .SetFixNorm(fixed_norm)
      .SetVerbosity(verbosity_);
  bbb_factory.AddBinByBin(*this, *other);
}

void CombineHarvester::CreateParameterIfEmpty(std::string const &name) {
  if (!params_.count(name)) {
    auto param = std::make_shared<Parameter>(Parameter());
    param->set_name(name);
    params_.insert({name, param});
  }
}

void CombineHarvester::MergeBinErrors(double bbb_threshold,
                                      double merge_threshold) {
  auto bbb_factory = ch::BinByBinFactory()
      .SetAddThreshold(bbb_threshold)
      .SetMergeThreshold(merge_threshold)
      .SetVerbosity(verbosity_);
  bbb_factory.MergeBinErrors(*this);
}

void CombineHarvester::InsertObservation(ch::Observation const& obs) {
  obs_.push_back(std::make_shared<ch::Observation>(obs));
}

void CombineHarvester::InsertProcess(ch::Process const& proc) {
  procs_.push_back(std::make_shared<ch::Process>(proc));
}

void CombineHarvester::InsertSystematic(ch::Systematic const& sys) {
  systs_.push_back(std::make_shared<ch::Systematic>(sys));
}
}
