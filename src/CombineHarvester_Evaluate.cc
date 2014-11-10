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
#include "TDirectory.h"
#include "TH1.h"
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/MakeUnique.h"
#include "CombineTools/interface/Utilities.h"
#include "CombineTools/interface/Algorithm.h"

// #include "TMath.h"
// #include "boost/format.hpp"
// #include "Utilities/interface/FnPredicates.h"
// #include "Math/QuantFuncMathCore.h"

namespace ch {

CombineHarvester::ProcNusMap CombineHarvester::GenerateProcNusMap() {
  ProcNusMap lookup(procs_.size());
  for (unsigned i = 0; i < nus_.size(); ++i) {
    for (unsigned j = 0; j < procs_.size(); ++j) {
      if (MatchingProcess(*(nus_[i]), *(procs_[j]))) {
        lookup[j].push_back(nus_[i].get());
      }
    }
  }
  return lookup;
}

double CombineHarvester::GetUncertainty() {
  auto lookup = GenerateProcNusMap();
  double err_sq = 0.0;
  for (auto param_it : params_) {
    double backup = param_it.second->val();
    param_it.second->set_val(backup+param_it.second->err_d());
    double rate_d = this->GetRateInternal(lookup, param_it.first);
    param_it.second->set_val(backup+param_it.second->err_u());
    double rate_u = this->GetRateInternal(lookup, param_it.first);
    double err = std::fabs(rate_u-rate_d) / 2.0;
    err_sq += err * err;
    param_it.second->set_val(backup);
  }
  return std::sqrt(err_sq);
}

double CombineHarvester::GetUncertainty(RooFitResult const* fit, unsigned n_samples) {
  auto lookup = GenerateProcNusMap();
  double rate = GetRateInternal(lookup);
  double err_sq = 0.0;
  auto backup = GetParameters();
  for (unsigned i = 0; i < n_samples; ++i) {
    // std::cout << "Sampling " << i << "/" << n_samples << std::endl;
    this->UpdateParameters(ch::ExtractSampledFitParameters(*(fit)));
    double rand_rate = this->GetRateInternal(lookup);
    double err = std::fabs(rand_rate-rate);
    err_sq += (err*err);
  }
  this->UpdateParameters(backup);
  return std::sqrt(err_sq/double(n_samples));
}

TH1F CombineHarvester::GetShapeWithUncertainty() {
  auto lookup = GenerateProcNusMap();
  TH1F shape = GetShape();
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, 0.0);
  }
  for (auto param_it : params_) {
    double backup = param_it.second->val();
    param_it.second->set_val(backup+param_it.second->err_d());
    TH1F shape_d = this->GetShapeInternal(lookup, param_it.first);
    param_it.second->set_val(backup+param_it.second->err_u());
    TH1F shape_u = this->GetShapeInternal(lookup, param_it.first);
    for (int i = 1; i <= shape.GetNbinsX(); ++i) {
      double err = std::fabs(shape_u.GetBinContent(i)-shape_d.GetBinContent(i)) / 2.0;
      shape.SetBinError(i, err*err + shape.GetBinError(i));
    }
    param_it.second->set_val(backup);
  }
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, std::sqrt(shape.GetBinError(i)));
  }
  return shape;
}

TH1F CombineHarvester::GetShapeWithUncertainty(RooFitResult const* fit, unsigned n_samples) {
  auto lookup = GenerateProcNusMap();
  TH1F shape = GetShapeInternal(lookup);
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, 0.0);
  }
  auto backup = GetParameters();
  for (unsigned i = 0; i < n_samples; ++i) {
    // std::cout << "Sampling " << i << "/" << n_samples << std::endl;
    this->UpdateParameters(ch::ExtractSampledFitParameters(*(fit)));
    TH1F rand_shape = this->GetShapeInternal(lookup);
    for (int i = 1; i <= shape.GetNbinsX(); ++i) {
      double err = std::fabs(rand_shape.GetBinContent(i)-shape.GetBinContent(i));
      shape.SetBinError(i, err*err + shape.GetBinError(i));
    }
  }
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, std::sqrt(shape.GetBinError(i)/double(n_samples)));
  }
  this->UpdateParameters(backup);
  return shape;
}

double CombineHarvester::GetRate() {
  auto lookup = GenerateProcNusMap();
  return GetRateInternal(lookup);
}

TH1F CombineHarvester::GetShape() {
  auto lookup = GenerateProcNusMap();
  return GetShapeInternal(lookup);
}


double CombineHarvester::GetRateInternal(ProcNusMap const& lookup,
    std::string const& single_nus) {
  double rate = 0.0;
  // TH1F tot_shape
  for (unsigned i = 0; i < procs_.size(); ++i) {
    double p_rate = procs_[i]->rate();
    // If we are evaluating the effect of a single parameter
    // check the list of associated nuisances and skip if
    // this "single_nus" is not in the list
    // However - we can't skip if the process has a pdf, as
    // we haven't checked what the parameters are
    if (single_nus != "" && !procs_[i]->pdf()) {
      if (!ch::any_of(lookup[i], [&](Nuisance const* nus) {
        return nus->name() == single_nus;
      })) continue;
    }
    for (auto nus_it : lookup[i]) {
      double x = params_[nus_it->name()]->val();
      if (nus_it->asymm()) {
        if (x >= 0) {
          p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
        } else {
          p_rate *= std::pow(nus_it->value_d(), -1.0 * x * nus_it->scale());
        }
      } else {
        p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
      }
    }
    rate += p_rate;
  }
  return rate;
}

TH1F CombineHarvester::GetShapeInternal(ProcNusMap const& lookup,
    std::string const& single_nus) {
  TH1F shape;
  bool shape_init = false;

  for (unsigned i = 0; i < procs_.size(); ++i) {
    // Might be able to skip if only interested in one nuisance
    // However - we can't skip if the process has a pdf, as
    // we haven't checked what the parameters are
    if (single_nus != "" && !procs_[i]->pdf()) {
      if (!ch::any_of(lookup[i], [&](Nuisance const* nus) {
        return nus->name() == single_nus;
      })) continue;
    }

    double p_rate = procs_[i]->rate();
    if (procs_[i]->shape()) {
      TH1F proc_shape;
      TH1F const* test_f = dynamic_cast<TH1F const*>(procs_[i]->shape());
      TH1D const* test_d = dynamic_cast<TH1D const*>(procs_[i]->shape());
      if (test_f) {
        test_f->Copy(proc_shape);
      } else if (test_d) {
        test_d->Copy(proc_shape);
      } else {
        throw std::runtime_error(
            "[CombineHarvester] TH1 is not a TH1F or a TH1D");
      }

      for (auto nus_it : lookup[i]) {
        double x = params_[nus_it->name()]->val();
        if (nus_it->asymm()) {
          if (x >= 0) {
            p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
          } else {
            p_rate *= std::pow(nus_it->value_d(), -1.0 * x * nus_it->scale());
          }
          if (nus_it->type() == "shape") {
            ShapeDiff(x * nus_it->scale(), &proc_shape, procs_[i]->shape(),
                      nus_it->shape_d(), nus_it->shape_u());
          }
        } else {
          p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
        }
      }
      for (int b = 1; b <= proc_shape.GetNbinsX(); ++b) {
        if (proc_shape.GetBinContent(b) < 0.) proc_shape.SetBinContent(b, 0.);
      }
      proc_shape.Scale(p_rate);
      if (!shape_init) {
        proc_shape.Copy(shape);
        shape.Reset();
        shape_init = true;
      }
      shape.Add(&proc_shape);
    } else if (procs_[i]->pdf()) {
      RooAbsData const* data_obj = FindMatchingData(procs_[i].get());
      std::string var_name = "CMS_th1x";
      if (data_obj) var_name = data_obj->get()->first()->GetName();
      TH1::AddDirectory(false);
      TH1F proc_shape = *(dynamic_cast<TH1F*>(procs_[i]->pdf()->createHistogram(
                             var_name.c_str())));
      for (auto nus_it : lookup[i]) {
        double x = params_[nus_it->name()]->val();
        if (nus_it->asymm()) {
          if (x >= 0) {
            p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
          } else {
            p_rate *= std::pow(nus_it->value_d(), -1.0 * x * nus_it->scale());
          }
        } else {
          p_rate *= std::pow(nus_it->value_u(), x * nus_it->scale());
        }
      }
      proc_shape.Scale(p_rate);
      if (!shape_init) {
        proc_shape.Copy(shape);
        shape.Reset();
        shape_init = true;
      }
      shape.Add(&proc_shape);
    }
  }
  return shape;
}

double CombineHarvester::GetObservedRate() {
  double rate = 0.0;
  for (unsigned i = 0; i < obs_.size(); ++i) {
    rate += obs_[i]->rate();
  }
  return rate;
}

TH1F CombineHarvester::GetObservedShape() {
  TH1F shape;
  bool shape_init = false;

  for (unsigned i = 0; i < obs_.size(); ++i) {
    TH1F proc_shape;
    double p_rate = obs_[i]->rate();
    if (obs_[i]->shape()) {
      TH1F const* test_f = dynamic_cast<TH1F const*>(obs_[i]->shape());
      TH1D const* test_d = dynamic_cast<TH1D const*>(obs_[i]->shape());
      if (test_f) {
        test_f->Copy(proc_shape);
      } else if (test_d) {
        test_d->Copy(proc_shape);
      } else {
        throw std::runtime_error(
            "[CombineHarvester] TH1 is not a TH1F or a TH1D");
      }
    } else if (obs_[i]->data()) {
      std::string var_name = obs_[i]->data()->get()->first()->GetName();
      proc_shape = *(dynamic_cast<TH1F*>(obs_[i]->data()->createHistogram(
                             var_name.c_str())));
      proc_shape.Scale(1. / proc_shape.Integral());
    }
    proc_shape.Scale(p_rate);
    if (!shape_init) {
      proc_shape.Copy(shape);
      shape.Reset();
      shape_init = true;
    }
    shape.Add(&proc_shape);
  }
  return shape;
}

void CombineHarvester::ShapeDiff(double x,
    TH1F * target,
    TH1 const* nom,
    TH1 const* low,
    TH1 const* high) {
  double fx = smoothStepFunc(x);
  for (int i = 1; i <= target->GetNbinsX(); ++i) {
    float h = high->GetBinContent(i);
    float l = low->GetBinContent(i);
    float n = nom->GetBinContent(i);
    target->SetBinContent(
        i, target->GetBinContent(i) + 0.5 * ((h - l) + (h + l - 2. * n) * fx));
  }
}

void CombineHarvester::SetParameters(std::vector<ch::Parameter> params) {
  params_.clear();
  for (unsigned i = 0; i < params.size(); ++i) {
    params_[params[i].name()] = std::make_shared<ch::Parameter>(params[i]);
  }
}

void CombineHarvester::RenameParameter(std::string const& oldname,
                                       std::string const& newname) {
  auto it = params_.find(oldname);
  if (it != params_.end()) {
    params_[newname] = it->second;
    params_[newname]->set_name(newname);
    params_.erase(it);
  }
}

void CombineHarvester::UpdateParameters(std::vector<ch::Parameter> params) {
  for (unsigned i = 0; i < params.size(); ++i) {
    auto it = params_.find(params[i].name());
    if (it != params_.end()) {
      it->second->set_val(params[i].val());
      it->second->set_err_d(params[i].err_d());
      it->second->set_err_u(params[i].err_u());
    } else {
      if (verbosity_ >= 1) {
        log() << "[UpdateParameters] Parameter \"" << params[i].name()
              << "\" is not defined\n";
      }
    }
  }
}

std::vector<ch::Parameter> CombineHarvester::GetParameters() const {
  std::vector<ch::Parameter> params;
  for (auto const& it : params_) {
    params.push_back(*(it.second));
  }
  return params;
}

void CombineHarvester::VariableRebin(std::vector<double> bins) {
  std::vector<double> proc_scaling(procs_.size(), 0.0);
  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (procs_[i]->shape()) {
      TH1 *copy = static_cast<TH1*>(procs_[i]->shape()->Clone());
      double rate_before = copy->Integral();
      TH1 *copy2 = copy->Rebin(bins.size()-1, "", &(bins[0]));
      double rate_after = copy2->Integral();
      if (rate_after > 0.) copy2->Scale(1.0/rate_after);
      procs_[i]->set_shape(std::unique_ptr<TH1>(copy2));
      procs_[i]->set_rate(procs_[i]->rate() * (rate_after/rate_before));
      proc_scaling[i] = rate_after/rate_before;
      delete copy;
    }
  }
  for (unsigned i = 0; i < obs_.size(); ++i) {
    if (obs_[i]->shape()) {
      TH1 *copy = static_cast<TH1*>(obs_[i]->shape()->Clone());
      double rate_before = copy->Integral();
      TH1 *copy2 = copy->Rebin(bins.size()-1, "", &(bins[0]));
      double rate_after = copy2->Integral();
      if (rate_after > 0.) copy2->Scale(1.0/rate_after);
      obs_[i]->set_shape(std::unique_ptr<TH1>(copy2));
      obs_[i]->set_rate(obs_[i]->rate() * (rate_after/rate_before));
      delete copy;
    }
  }
  for (unsigned i = 0; i < nus_.size(); ++i) {
    double proc_scale = 1.0;
    for (unsigned j = 0; j < procs_.size(); ++j) {
      if (MatchingProcess(*(procs_[j]), *(nus_[i].get()))) {
        proc_scale = proc_scaling[j];
      }
    }
    if (nus_[i]->shape_u()) {
      TH1 *copy = static_cast<TH1*>(nus_[i]->shape_u()->Clone());
      double rate_before = copy->Integral();
      TH1 *copy2 = copy->Rebin(bins.size()-1, "", &(bins[0]));
      double rate_after = copy2->Integral();
      if (rate_after > 0.) copy2->Scale(1.0/rate_after);
      nus_[i]->set_shape_u(std::unique_ptr<TH1>(copy2));
      nus_[i]->set_value_u((rate_after / rate_before) * nus_[i]->value_u() /
                           proc_scale);
      delete copy;
    }
    if (nus_[i]->shape_d()) {
      TH1 *copy = static_cast<TH1*>(nus_[i]->shape_d()->Clone());
      double rate_before = copy->Integral();
      TH1 *copy2 = copy->Rebin(bins.size()-1, "", &(bins[0]));
      double rate_after = copy2->Integral();
      if (rate_after > 0.) copy2->Scale(1.0/rate_after);
      nus_[i]->set_shape_d(std::unique_ptr<TH1>(copy2));
      nus_[i]->set_value_d((rate_after / rate_before) * nus_[i]->value_d() /
                           proc_scale);
      delete copy;
    }
  }
}

void CombineHarvester::SetPdfBins(unsigned nbins) {
  for (unsigned i = 0; i < procs_.size(); ++i) {
    std::set<std::string> binning_vars;
    if (procs_[i]->pdf()) {
      RooAbsData const* data_obj = FindMatchingData(procs_[i].get());
      std::string var_name = "CMS_th1x";
      if (data_obj) var_name = data_obj->get()->first()->GetName();
      binning_vars.insert(var_name);
    }
    for (auto & it : wspaces_) {
      for (auto & var : binning_vars) {
        RooRealVar* avar =
            dynamic_cast<RooRealVar*>(it.second->var(var.c_str()));
        if (avar) avar->setBins(nbins);
      }
    }
  }

}

}
