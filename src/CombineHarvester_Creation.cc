#include "CombineTools/interface/CombineHarvester.h"
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <iostream>
#include <utility>
#include <algorithm>
#include "TDirectory.h"
#include "TH1.h"
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/Utilities.h"

namespace ch {
void CombineHarvester::AddObservations(
    std::vector<std::string> mass,
    std::vector<std::string> analysis,
    std::vector<std::string> era,
    std::vector<std::string> channel,
    std::vector<std::pair<int, std::string>> bin) {
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
    std::vector<std::pair<int, std::string>> bin,
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
  for (unsigned  i = 0; i < nus_.size(); ++i) {
    if (nus_[i]->type() != "shape") continue;
    LoadShapes(nus_[i].get(), mapping);
  }
}

void CombineHarvester::AddWorkspace(RooWorkspace const *ws) {
  if (wspaces_.count(ws->GetName())) return;
  if (verbosity_ >= 1)
    log() << "[AddWorkspace] Cloning workspace \"" << ws->GetName() << "\"\n";
  wspaces_[ws->GetName()] = std::shared_ptr<RooWorkspace>(
      reinterpret_cast<RooWorkspace *>(ws->Clone()));
}

void CombineHarvester::ExtractPdfs(std::string const &ws_name,
                                   std::string const &rule,
                                   CombineHarvester *other) {
  CombineHarvester *target = other ? other : this;
  std::vector<HistMapping> mapping(1);
  mapping[0].process = "*";
  mapping[0].category = "*";
  mapping[0].pattern = ws_name+":"+rule;
  if (!wspaces_.count(ws_name)) return;
  for (unsigned  i = 0; i < procs_.size(); ++i) {
    if (!procs_[i]->pdf()) {
      target->LoadShapes(procs_[i].get(), mapping);
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
  for (unsigned  i = 0; i < obs_.size(); ++i) {
    if (!obs_[i]->data()) {
      LoadShapes(obs_[i].get(), mapping);
    }
  }
}

void CombineHarvester::AddBinByBin(double threshold, bool fixed_norm,
                                   CombineHarvester *other) {
  unsigned bbb_added = 0;
  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (!procs_[i]->shape()) continue;
    TH1 const* h = procs_[i]->shape();
    unsigned n_pop_bins = 0;
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      if (h->GetBinContent(j) > 0.0) ++n_pop_bins;
    }
    if (n_pop_bins <= 1 && fixed_norm) {
      if (verbosity_ >= 1) {
        std::cout << "Requested fixed_norm but template has <= 1 populated "
                     "bins, skipping\n";
        std::cout << Process::PrintHeader << *(procs_[i]) << "\n";
      }
      continue;
    }
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      bool do_bbb = false;
      double val = h->GetBinContent(j);
      double err = h->GetBinError(j);
      if (val == 0. && err > 0.) do_bbb = true;
      if (val > 0. && (err / val) > threshold) do_bbb = true;
      // if (h->GetBinContent(j) <= 0.0) {
      //   if (h->GetBinError(j) > 0.0) {
      //     std::cout << *(procs_[i]) << "\n";
      //     std::cout << "Bin with content <= 0 and error > 0 found, skipping\n";
      //   }
      //   continue;
      // }
      if (do_bbb) {
        ++bbb_added;
        auto nus = std::make_shared<Nuisance>();
        ch::SetProperties(nus.get(), procs_[i].get());
        nus->set_type("shape");
        // nus->set_name("CMS_" + nus->bin() + "_" + nus->process() + "_bin_" +
        //               boost::lexical_cast<std::string>(j));
        nus->set_name("CMS_" + nus->analysis() + "_" + nus->channel() + "_" +
                      nus->bin() + "_" + nus->era() + "_" + nus->process() +
                      "_bin_" + boost::lexical_cast<std::string>(j));
        nus->set_asymm(true);
        TH1 *h_d = static_cast<TH1 *>(h->Clone());
        TH1 *h_u = static_cast<TH1 *>(h->Clone());
        h_d->SetBinContent(j, val - err);
        if (h_d->GetBinContent(j) < 0.) h_d->SetBinContent(j, 0.);
        h_u->SetBinContent(j, val + err);
        if (fixed_norm) {
          nus->set_value_d(1.0);
          nus->set_value_u(1.0);
        } else {
          nus->set_value_d(h_d->Integral()/h->Integral());
          nus->set_value_u(h_u->Integral()/h->Integral());
        }
        if (h_d->Integral() > 0.0) h_d->Scale(1.0/h_d->Integral());
        if (h_u->Integral() > 0.0) h_u->Scale(1.0/h_u->Integral());
        if (h_d) nus->set_shape_d(std::unique_ptr<TH1>(h_d));
        if (h_u) nus->set_shape_u(std::unique_ptr<TH1>(h_u));
        CombineHarvester::CreateParameterIfEmpty(other ? other : this,
                                                 nus->name());
        if (other) {
          other->nus_.push_back(nus);
        } else {
          nus_.push_back(nus);
        }
      }
    }
  }
  // std::cout << "bbb added: " << bbb_added << std::endl;
}

void CombineHarvester::CreateParameterIfEmpty(CombineHarvester *cmb,
                                              std::string const &name) {
  if (!params_.count(name)) {
    auto param = std::make_shared<Parameter>(Parameter());
    param->set_name(name);
    (*cmb).params_.insert({name, param});
  }
}

void CombineHarvester::MergeBinErrors(double bbb_threshold,
                                      double merge_threshold) {
  // Reduce merge_threshold very slightly to avoid numerical issues
  // E.g. two backgrounds each with bin error 1.0. merge_threshold of
  // 0.5 should not result in merging - but can do depending on
  // machine and compiler
  merge_threshold -= 1E-9 * merge_threshold;
  auto bins = this->bin_set();
  for (auto const& bin : bins) {
    unsigned bbb_added = 0;
    unsigned bbb_removed = 0;
    CombineHarvester tmp = std::move(this->cp().bin({bin}).histograms());
    if (tmp.procs_.size() == 0) continue;

    std::vector<TH1 *> h_copies(tmp.procs_.size(), nullptr);
    for (unsigned i = 0; i < h_copies.size(); ++i) {
      h_copies[i] = static_cast<TH1 *>(tmp.procs_[i]->shape()->Clone());
      h_copies[i]->Scale(tmp.procs_[i]->rate());
    }

    for (int i = 1; i <= h_copies[0]->GetNbinsX(); ++i) {
      double tot_bbb_added = 0.0;
      std::vector<std::pair<double, TH1 *>> result;
      for (unsigned j = 0; j < h_copies.size(); ++j) {
        double val = h_copies[j]->GetBinContent(i);
        double err = h_copies[j]->GetBinError(i);
        if (val == 0.0 &&  err == 0.0) continue;
        if (val == 0 || (err/val) > bbb_threshold) {
          bbb_added += 1;
          tot_bbb_added += (err * err);
          result.push_back(std::make_pair(err*err, h_copies[j]));
        }
      }
      if (tot_bbb_added == 0.0) continue;
      std::sort(result.begin(), result.end());
      double removed = 0.0;
      for (unsigned r = 0; r < result.size(); ++r) {
        if ((result[r].first + removed) < (merge_threshold * tot_bbb_added) &&
            r < (result.size() - 1)) {
          bbb_removed += 1;
          removed += result[r].first;
          result[r].second->SetBinError(i, 0.0);
        }
      }
      double expand = std::sqrt(1. / (1. - (removed / tot_bbb_added)));
      for (unsigned r = 0; r < result.size(); ++r) {
        result[r]
            .second->SetBinError(i, result[r].second->GetBinError(i) * expand);
      }
    }
    for (unsigned i = 0; i < h_copies.size(); ++i) {
      if (h_copies[i]->Integral() > 0.) {
        h_copies[i]->Scale(1. /  h_copies[i]->Integral());
      }
      tmp.procs_[i]->set_shape(std::unique_ptr<TH1>(h_copies[i]));
    }
    // std::cout << "BIN: " << bin << std::endl;
    // std::cout << "Total bbb added:    " << bbb_added << "\n";
    // std::cout << "Total bbb removed:  " << bbb_removed << "\n";
    // std::cout << "Total bbb =======>: " << bbb_added-bbb_removed << "\n";
  }
}
}
