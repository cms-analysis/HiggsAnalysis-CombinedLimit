#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
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
#include "TH2.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/MakeUnique.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"

// #include "TMath.h"
// #include "boost/format.hpp"
// #include "Utilities/interface/FnPredicates.h"
// #include "Math/QuantFuncMathCore.h"

namespace ch {

CombineHarvester::ProcSystMap CombineHarvester::GenerateProcSystMap() {
  ProcSystMap lookup(procs_.size());
  for (unsigned i = 0; i < systs_.size(); ++i) {
    for (unsigned j = 0; j < procs_.size(); ++j) {
      if (MatchingProcess(*(systs_[i]), *(procs_[j]))) {
        lookup[j].push_back(systs_[i].get());
      }
    }
  }
  return lookup;
}

double CombineHarvester::GetUncertainty() {
  auto lookup = GenerateProcSystMap();
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

double CombineHarvester::GetUncertainty(RooFitResult const* fit,
                                        unsigned n_samples) {
  return GetUncertainty(*fit, n_samples);
}

double CombineHarvester::GetUncertainty(RooFitResult const& fit,
                                        unsigned n_samples) {
  auto lookup = GenerateProcSystMap();
  double rate = GetRateInternal(lookup);
  double err_sq = 0.0;

  // Create a backup copy of the current parameter values
  auto backup = GetParameters();

  // Calling randomizePars() ensures that the RooArgList of sampled parameters
  // is already created within the RooFitResult
  RooArgList const& rands = fit.randomizePars();

  // Now create two aligned vectors of the RooRealVar parameters and the
  // corresponding ch::Parameter pointers
  int n_pars = rands.getSize();
  std::vector<RooRealVar const*> r_vec(n_pars, nullptr);
  std::vector<ch::Parameter*> p_vec(n_pars, nullptr);
  for (unsigned n = 0; n < p_vec.size(); ++n) {
    r_vec[n] = dynamic_cast<RooRealVar const*>(rands.at(n));
    p_vec[n] = GetParameter(r_vec[n]->GetName());
  }

  for (unsigned i = 0; i < n_samples; ++i) {
    // Randomise and update values
    fit.randomizePars();
    for (int n = 0; n < n_pars; ++n) {
      if (p_vec[n]) p_vec[n]->set_val(r_vec[n]->getVal());
    }

    double rand_rate = this->GetRateInternal(lookup);
    double err = std::fabs(rand_rate-rate);
    err_sq += (err*err);
  }
  this->UpdateParameters(backup);
  return std::sqrt(err_sq/double(n_samples));
}

TH1F CombineHarvester::GetShapeWithUncertainty() {
  auto lookup = GenerateProcSystMap();
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
      double err =
          std::fabs(shape_u.GetBinContent(i) - shape_d.GetBinContent(i)) / 2.0;
      shape.SetBinError(i, err*err + shape.GetBinError(i));
    }
    param_it.second->set_val(backup);
  }
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, std::sqrt(shape.GetBinError(i)));
  }
  return shape;
}

TH1F CombineHarvester::GetShapeWithUncertainty(RooFitResult const* fit,
                                               unsigned n_samples) {
  return GetShapeWithUncertainty(*fit, n_samples);
}

TH1F CombineHarvester::GetShapeWithUncertainty(RooFitResult const& fit,
                                               unsigned n_samples) {
  auto lookup = GenerateProcSystMap();
  TH1F shape = GetShapeInternal(lookup);
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, 0.0);
  }
  // Create a backup copy of the current parameter values
  auto backup = GetParameters();

  // Calling randomizePars() ensures that the RooArgList of sampled parameters
  // is already created within the RooFitResult
  RooArgList const& rands = fit.randomizePars();

  // Now create two aligned vectors of the RooRealVar parameters and the
  // corresponding ch::Parameter pointers
  int n_pars = rands.getSize();
  std::vector<RooRealVar const*> r_vec(n_pars, nullptr);
  std::vector<ch::Parameter*> p_vec(n_pars, nullptr);
  for (unsigned n = 0; n < p_vec.size(); ++n) {
    r_vec[n] = dynamic_cast<RooRealVar const*>(rands.at(n));
    p_vec[n] = GetParameter(r_vec[n]->GetName());
  }

  // Main loop through n_samples
  for (unsigned i = 0; i < n_samples; ++i) {
    // Randomise and update values
    fit.randomizePars();
    for (int n = 0; n < n_pars; ++n) {
      if (p_vec[n]) p_vec[n]->set_val(r_vec[n]->getVal());
    }

    TH1F rand_shape = this->GetShapeInternal(lookup);
    for (int i = 1; i <= shape.GetNbinsX(); ++i) {
      double err =
          std::fabs(rand_shape.GetBinContent(i) - shape.GetBinContent(i));
      shape.SetBinError(i, err*err + shape.GetBinError(i));
    }
  }
  for (int i = 1; i <= shape.GetNbinsX(); ++i) {
    shape.SetBinError(i, std::sqrt(shape.GetBinError(i)/double(n_samples)));
  }
  this->UpdateParameters(backup);
  return shape;
}

TH2F CombineHarvester::GetRateCovariance(RooFitResult const& fit,
                                         unsigned n_samples) {
  auto lookup = GenerateProcSystMap();

  unsigned n = procs_.size();
  TH1F nom("nominal", "nominal", n, 0, n);
  TH2F res("covariance", "covariance", n, 0, n, n, 0, n);

  std::vector<CombineHarvester> ch_procs;
  std::vector<std::string> labels;
  unsigned nbins = this->bin_set().size();

  for (unsigned i = 0; i < procs_.size(); ++i) {
    ch_procs.push_back(
        this->cp().bin({procs_[i]->bin()}).process({procs_[i]->process()}));
    if (nbins > 1) {
      labels.push_back(procs_[i]->bin() + "," + procs_[i]->process());
    } else {
      labels.push_back(procs_[i]->process());
    }
  }
  for (unsigned i = 0; i < procs_.size(); ++i) {
    nom.SetBinContent(i + 1, ch_procs[i].GetRate());
  }
  auto backup = GetParameters();

  // Calling randomizePars() ensures that the RooArgList of sampled parameters
  // is already created within the RooFitResult
  RooArgList const& rands = fit.randomizePars();

  // Now create two aligned vectors of the RooRealVar parameters and the
  // corresponding ch::Parameter pointers
  int n_pars = rands.getSize();
  std::vector<RooRealVar const*> r_vec(n_pars, nullptr);
  std::vector<ch::Parameter*> p_vec(n_pars, nullptr);
  for (unsigned n = 0; n < p_vec.size(); ++n) {
    r_vec[n] = dynamic_cast<RooRealVar const*>(rands.at(n));
    p_vec[n] = GetParameter(r_vec[n]->GetName());
  }

  // Main loop through n_samples
  for (unsigned rnd = 0; rnd < n_samples; ++rnd) {
    // Randomise and update values
    fit.randomizePars();
    for (int n = 0; n < n_pars; ++n) {
      if (p_vec[n]) p_vec[n]->set_val(r_vec[n]->getVal());
    }

    for (int i = 1; i <= nom.GetNbinsX(); ++i) {
      for (int j = 1; j <= nom.GetNbinsX(); ++j) {
        int x = j;
        int y = nom.GetNbinsX() - (i - 1);
        res.GetXaxis()->SetBinLabel(x, labels[j - 1].c_str());
        res.GetYaxis()->SetBinLabel(y, labels[i - 1].c_str());
        res.SetBinContent(
            x, y, res.GetBinContent(x, y) +
                      (ch_procs[i - 1].GetRate() - nom.GetBinContent(i)) *
                      (ch_procs[j - 1].GetRate() - nom.GetBinContent(j)));
      }
    }
  }

  for (int i = 1; i <= nom.GetNbinsX(); ++i) {
    for (int j = 1; j <= nom.GetNbinsX(); ++j) {
      int x = j;
      int y = nom.GetNbinsX() - (i - 1);
      res.SetBinContent(x, y, res.GetBinContent(x, y) / double(n_samples));
    }
  }
  this->UpdateParameters(backup);
  return res;
}

TH2F CombineHarvester::GetRateCorrelation(RooFitResult const& fit,
                                          unsigned n_samples) {
  TH2F cov = GetRateCovariance(fit, n_samples);
  TH2F res = cov;
  res.SetName("correlation");
  res.SetTitle("correlation");
  for (int i = 1; i <= cov.GetNbinsX(); ++i) {
    for (int j = 1; j <= cov.GetNbinsX(); ++j) {
      int x = j;
      int y = cov.GetNbinsX() - (i - 1);
      res.SetBinContent(
          x, y,
          cov.GetBinContent(x, y) /
              (std::sqrt(cov.GetBinContent(i, cov.GetNbinsX() - (i - 1))) *
               std::sqrt(cov.GetBinContent(j, cov.GetNbinsX() - (j - 1)))));
    }
  }
  return res;
}

double CombineHarvester::GetRate() {
  auto lookup = GenerateProcSystMap();
  return GetRateInternal(lookup);
}

TH1F CombineHarvester::GetShape() {
  auto lookup = GenerateProcSystMap();
  return GetShapeInternal(lookup);
}

double CombineHarvester::GetRateInternal(ProcSystMap const& lookup,
    std::string const& single_sys) {
  double rate = 0.0;
  // TH1F tot_shape
  for (unsigned i = 0; i < procs_.size(); ++i) {
    double p_rate = procs_[i]->rate();
    // If we are evaluating the effect of a single parameter
    // check the list of associated nuisances and skip if
    // this "single_sys" is not in the list
    // However - we can't skip if the process has a pdf, as
    // we haven't checked what the parameters are
    if (single_sys != "" && !procs_[i]->pdf()) {
      if (!ch::any_of(lookup[i], [&](Systematic const* sys) {
        return sys->name() == single_sys;
      })) continue;
    }
    for (auto sys_it : lookup[i]) {
      if (sys_it->type() == "rateParam") {
        continue;  // don't evaluate this for now
      }
      double x = params_[sys_it->name()]->val();
      if (sys_it->asymm()) {
        p_rate *= logKappaForX(x * sys_it->scale(), sys_it->value_d(),
                               sys_it->value_u());
      } else {
        p_rate *= std::pow(sys_it->value_u(), x * sys_it->scale());
      }
    }
    rate += p_rate;
  }
  return rate;
}

TH1F CombineHarvester::GetShapeInternal(ProcSystMap const& lookup,
    std::string const& single_sys) {
  TH1F shape;
  bool shape_init = false;

  for (unsigned i = 0; i < procs_.size(); ++i) {
    // Might be able to skip if only interested in one nuisance
    // However - we can't skip if the process has a pdf, as
    // we haven't checked what the parameters are
    if (single_sys != "" && !procs_[i]->pdf()) {
      if (!ch::any_of(lookup[i], [&](Systematic const* sys) {
        return sys->name() == single_sys;
      })) continue;
    }

    double p_rate = procs_[i]->rate();
    if (procs_[i]->shape() || procs_[i]->data()) {
      TH1F proc_shape = procs_[i]->ShapeAsTH1F();
      for (auto sys_it : lookup[i]) {
        if (sys_it->type() == "rateParam") {
          continue;  // don't evaluate this for now
        }
        auto param_it = params_.find(sys_it->name());
        if (param_it == params_.end()) {
          throw std::runtime_error(
              FNERROR("Parameter " + sys_it->name() +
                      " not found in CombineHarvester instance"));
        }
        double x = param_it->second->val();
        if (sys_it->asymm()) {
          p_rate *= logKappaForX(x * sys_it->scale(), sys_it->value_d(),
                                 sys_it->value_u());
          if (sys_it->type() == "shape" || sys_it->type() == "shapeN2" ||
              sys_it->type() == "shapeU") {
            bool linear = true;
            if (sys_it->type() == "shapeN2") linear = false;
            if (sys_it->shape_u() && sys_it->shape_d()) {
              ShapeDiff(x * sys_it->scale(), &proc_shape, procs_[i]->shape(),
                        sys_it->shape_d(), sys_it->shape_u(), linear);
            }
            if (sys_it->data_u() && sys_it->data_d()) {
              RooDataHist const* nom =
                  dynamic_cast<RooDataHist const*>(procs_[i]->data());
              if (nom) {
                ShapeDiff(x * sys_it->scale(), &proc_shape, nom,
                          sys_it->data_d(), sys_it->data_u());
              }
            }
          }
        } else {
          p_rate *= std::pow(sys_it->value_u(), x * sys_it->scale());
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
      if (!procs_[i]->observable()) {
        RooAbsData const* data_obj = FindMatchingData(procs_[i].get());
        std::string var_name = "CMS_th1x";
        if (data_obj) var_name = data_obj->get()->first()->GetName();
        procs_[i]->set_observable((RooRealVar *)procs_[i]->pdf()->findServer(var_name.c_str()));
      }
      TH1::AddDirectory(false);
      TH1F* tmp = (TH1F*)procs_[i]->observable()->createHistogram("");
      for (int b = 1; b <= tmp->GetNbinsX(); ++b) {
        procs_[i]->observable()->setVal(tmp->GetBinCenter(b));
        tmp->SetBinContent(b, tmp->GetBinWidth(b) * procs_[i]->pdf()->getVal());
      }
      TH1F proc_shape = *tmp;
      delete tmp;
      RooAbsPdf const* aspdf = dynamic_cast<RooAbsPdf const*>(procs_[i]->pdf());
      if ((aspdf && !aspdf->selfNormalized()) || (!aspdf)) {
        // LOGLINE(log(), "Have a pdf that is not selfNormalized");
        // std::cout << "Integral: " << proc_shape.Integral() << "\n";
        if (proc_shape.Integral() > 0.) {
          proc_shape.Scale(1. / proc_shape.Integral());
        }
      }
      for (auto sys_it : lookup[i]) {
        if (sys_it->type() == "rateParam") {
          continue;  // don't evaluate this for now
        }
        double x = params_[sys_it->name()]->val();
        if (sys_it->asymm()) {
          p_rate *= logKappaForX(x * sys_it->scale(), sys_it->value_d(),
                                 sys_it->value_u());
        } else {
          p_rate *= std::pow(sys_it->value_u(), x * sys_it->scale());
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
      proc_shape = obs_[i]->ShapeAsTH1F();
    } else if (obs_[i]->data()) {
      TH1F* tmp = dynamic_cast<TH1F*>(obs_[i]->data()->createHistogram(
          "", *(RooRealVar*)obs_[i]->data()->get()->first()));
      tmp->Sumw2(false);
      tmp->SetBinErrorOption(TH1::kPoisson);
      proc_shape = *tmp;
      delete tmp;
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
    TH1 const* high,
    bool linear) {
  double fx = smoothStepFunc(x);
  for (int i = 1; i <= target->GetNbinsX(); ++i) {
    float h = high->GetBinContent(i);
    float l = low->GetBinContent(i);
    float n = nom->GetBinContent(i);
    if (!linear) {
      float t = target->GetBinContent(i);
      target->SetBinContent(i, t > 0. ? std::log(t) : -999.);
      h = (h > 0. && n > 0.) ? std::log(h/n) : 0.;
      l = (l > 0. && n > 0.) ? std::log(l/n) : 0.;
      target->SetBinContent(i, target->GetBinContent(i) +
                                   0.5 * x * ((h - l) + (h + l) * fx));
      target->SetBinContent(i, std::exp(target->GetBinContent(i)));
    } else {
      target->SetBinContent(i, target->GetBinContent(i) +
                                   0.5 * x * ((h - l) + (h + l - 2. * n) * fx));
    }
  }
}

void CombineHarvester::ShapeDiff(double x,
    TH1F * target,
    RooDataHist const* nom,
    RooDataHist const* low,
    RooDataHist const* high) {
  double fx = smoothStepFunc(x);
  for (int i = 1; i <= target->GetNbinsX(); ++i) {
    high->get(i-1);
    low->get(i-1);
    nom->get(i-1);
    // The RooDataHists are not scaled to unity (unlike in the TH1 case above)
    // so we have to normalise the bin contents on the fly
    float h = high->weight() / high->sumEntries();
    float l = low->weight() / low->sumEntries();
    float n = nom->weight() / nom->sumEntries();
    target->SetBinContent(i, target->GetBinContent(i) +
                                 0.5 * x * ((h - l) + (h + l - 2. * n) * fx));
  }
}

// void CombineHarvester::SetParameters(std::vector<ch::Parameter> params) {
//   params_.clear();
//   for (unsigned i = 0; i < params.size(); ++i) {
//     params_[params[i].name()] = std::make_shared<ch::Parameter>(params[i]);
//   }
// }

void CombineHarvester::RenameParameter(std::string const& oldname,
                                       std::string const& newname) {
  auto it = params_.find(oldname);
  if (it != params_.end()) {
    params_[newname] = it->second;
    params_[newname]->set_name(newname);
    params_.erase(it);
  }
}

ch::Parameter const* CombineHarvester::GetParameter(
    std::string const& name) const {
  auto it = params_.find(name);
  if (it != params_.end()) {
    return it->second.get();
  } else {
    return nullptr;
  }
}

ch::Parameter* CombineHarvester::GetParameter(std::string const& name) {
  auto it = params_.find(name);
  if (it != params_.end()) {
    return it->second.get();
  } else {
    return nullptr;
  }
}

void CombineHarvester::UpdateParameters(
    std::vector<ch::Parameter> const& params) {
  for (unsigned i = 0; i < params.size(); ++i) {
    auto it = params_.find(params[i].name());
    if (it != params_.end()) {
      it->second->set_val(params[i].val());
      it->second->set_err_d(params[i].err_d());
      it->second->set_err_u(params[i].err_u());
    } else {
      if (verbosity_ >= 1) {
        LOGLINE(log(), "Parameter " + params[i].name() + " is not defined");
      }
    }
  }
}

void CombineHarvester::UpdateParameters(RooFitResult const& fit) {
  for (int i = 0; i < fit.floatParsFinal().getSize(); ++i) {
    RooRealVar const* var =
        dynamic_cast<RooRealVar const*>(fit.floatParsFinal().at(i));
    // check for failed cast here
    auto it = params_.find(std::string(var->GetName()));
    if (it != params_.end()) {
      it->second->set_val(var->getVal());
      it->second->set_err_d(var->getErrorLo());
      it->second->set_err_u(var->getErrorHi());
    } else {
      if (verbosity_ >= 1) {
        LOGLINE(log(),
                "Parameter " + std::string(var->GetName()) + " is not defined");
      }
    }
  }
}

void CombineHarvester::UpdateParameters(RooFitResult const* fit) {
  UpdateParameters(*fit);
}

std::vector<ch::Parameter> CombineHarvester::GetParameters() const {
  std::vector<ch::Parameter> params;
  for (auto const& it : params_) {
    params.push_back(*(it.second));
  }
  return params;
}

void CombineHarvester::VariableRebin(std::vector<double> bins) {
  // We need to keep a record of the Process rates before we rebin. The
  // reasoning comes from the following scenario: the user might choose a new
  // binning which excludes some of the existing bins - thus changing the
  // process normalisation. This is fine, but we also need to adjust the shape
  // Systematic entries - both the rebinning and the adjustment of the value_u
  // and value_d shifts.
  std::vector<double> prev_proc_rates(procs_.size());

  // Also hold on the scaled Process hists *after* the rebinning is done - these
  // are needed to update the associated Systematic entries
  std::vector<std::unique_ptr<TH1>> scaled_procs(procs_.size());

  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (procs_[i]->shape()) {
      // Get the scaled shape here
      std::unique_ptr<TH1> copy(procs_[i]->ClonedScaledShape());
      // shape norm should only be "no_norm_rate"
      prev_proc_rates[i] = procs_[i]->no_norm_rate();
      std::unique_ptr<TH1> copy2(copy->Rebin(bins.size()-1, "", &(bins[0])));
      // The process shape & rate will be reset here
      procs_[i]->set_shape(std::move(copy2), true);
      scaled_procs[i] = procs_[i]->ClonedScaledShape();
    }
  }
  for (unsigned i = 0; i < obs_.size(); ++i) {
    if (obs_[i]->shape()) {
      std::unique_ptr<TH1> copy(obs_[i]->ClonedScaledShape());
      std::unique_ptr<TH1> copy2(copy->Rebin(bins.size()-1, "", &(bins[0])));
      obs_[i]->set_shape(std::move(copy2), true);
    }
  }
  for (unsigned i = 0; i < systs_.size(); ++i) {
    TH1 const* proc_hist = nullptr;
    double prev_rate = 0.;
    for (unsigned j = 0; j < procs_.size(); ++j) {
      if (MatchingProcess(*(procs_[j]), *(systs_[i].get()))) {
        proc_hist = scaled_procs[j].get();
        prev_rate = prev_proc_rates[j];
      }
    }
    if (systs_[i]->shape_u() && systs_[i]->shape_d()) {
      // These hists will be normalised to unity
      std::unique_ptr<TH1> copy_u(systs_[i]->ClonedShapeU());
      std::unique_ptr<TH1> copy_d(systs_[i]->ClonedShapeD());

      // If we found a matching Process we will scale this back up to their
      // initial rates
      if (proc_hist) {
        copy_u->Scale(systs_[i]->value_u() * prev_rate);
        copy_d->Scale(systs_[i]->value_d() * prev_rate);
      }
      std::unique_ptr<TH1> copy2_u(
          copy_u->Rebin(bins.size() - 1, "", &(bins[0])));
      std::unique_ptr<TH1> copy2_d(
          copy_d->Rebin(bins.size() - 1, "", &(bins[0])));
      // If we have proc_hist != nullptr, set_shapes will re-calculate value_u
      // and value_d for us, before scaling the new hists back to unity
      systs_[i]->set_shapes(std::move(copy2_u), std::move(copy2_d), proc_hist);
    }
  }
}

void CombineHarvester::ZeroBins(double min, double max) {
  // We need to keep a record of the Process rates before we set bins to 0. The
  // This is necessary because we need to make sure the process normalisation
  // and shape systematic entries are correctly adjusted
  std::vector<double> prev_proc_rates(procs_.size());

  // Also hold on the scaled Process hists *after* some of the bins have
  // been set to 0 - these
  // are needed to update the associated Systematic entries
  std::vector<std::unique_ptr<TH1>> scaled_procs(procs_.size());

  for (unsigned i = 0; i < procs_.size(); ++i) {
    if (procs_[i]->shape()) {
      // Get the scaled shape here
      std::unique_ptr<TH1> copy(procs_[i]->ClonedScaledShape());
      // shape norm should only be "no_norm_rate"
      prev_proc_rates[i] = procs_[i]->no_norm_rate();
      for (int j = 1; j <= copy->GetNbinsX();j ++){
        if(copy->GetBinLowEdge(j) >= min && copy->GetBinLowEdge(j+1) <= max){
          copy->SetBinContent(j,0.);
          copy->SetBinError(j,0.);
        }
      }
      // The process shape & rate will be reset here
      procs_[i]->set_shape(std::move(copy), true);
      scaled_procs[i] = procs_[i]->ClonedScaledShape();
    }
  }
  for (unsigned i = 0; i < obs_.size(); ++i) {
    if (obs_[i]->shape()) {
      std::unique_ptr<TH1> copy(obs_[i]->ClonedScaledShape());
       for (int j = 1; j <= copy->GetNbinsX();j ++){
        if(copy->GetBinLowEdge(j) >= min && copy->GetBinLowEdge(j+1) <= max){
          copy->SetBinContent(j,0.);
          copy->SetBinError(j,0.);
        }
      }
      obs_[i]->set_shape(std::move(copy), true);
    }
  }
  for (unsigned i = 0; i < systs_.size(); ++i) {
    TH1 const* proc_hist = nullptr;
    double prev_rate = 0.;
    for (unsigned j = 0; j < procs_.size(); ++j) {
      if (MatchingProcess(*(procs_[j]), *(systs_[i].get()))) {
        proc_hist = scaled_procs[j].get();
        prev_rate = prev_proc_rates[j];
      }
    }
    if (systs_[i]->shape_u() && systs_[i]->shape_d()) {
      // These hists will be normalised to unity
      std::unique_ptr<TH1> copy_u(systs_[i]->ClonedShapeU());
      std::unique_ptr<TH1> copy_d(systs_[i]->ClonedShapeD());

      // If we found a matching Process we will scale this back up to their
      // initial rates
      if (proc_hist) {
        copy_u->Scale(systs_[i]->value_u() * prev_rate);
        copy_d->Scale(systs_[i]->value_d() * prev_rate);
      }
       for (int j = 1; j <= copy_u->GetNbinsX();j ++){
        if(copy_u->GetBinLowEdge(j) >= min && copy_u->GetBinLowEdge(j+1) <= max){
          copy_u->SetBinContent(j,0.);
          copy_u->SetBinError(j,0.);
        }
        if(copy_d->GetBinLowEdge(j) >= min && copy_d->GetBinLowEdge(j+1) <= max){
          copy_d->SetBinContent(j,0.);
          copy_d->SetBinError(j,0.);
        }
      }
      // If we have proc_hist != nullptr, set_shapes will re-calculate value_u
      // and value_d for us, before scaling the new hists back to unity
      systs_[i]->set_shapes(std::move(copy_u), std::move(copy_d), proc_hist);
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

// This implementation "borowed" from
// HiggsAnalysis/CombinedLimit/src/ProcessNormalization.cc
double CombineHarvester::logKappaForX(double x, double k_low,
                                      double k_high) const {
  if (k_high == 0. || k_low == 0.) {
    if (verbosity_ >= 1) {
      LOGLINE(log(), "Have kappa=0.0 (scaling ill-defined), returning 1.0");
    }
    return 1.;
  }
  if (fabs(x) >= 0.5)
    return (x >= 0 ? std::pow(k_high, x) : std::pow(k_low, -1.0 * x));
  // interpolate between log(kappaHigh) and -log(kappaLow)
  //    logKappa(x) = avg + halfdiff * h(2x)
  // where h(x) is the 3th order polynomial
  //    h(x) = (3 x^5 - 10 x^3 + 15 x)/8;
  // chosen so that h(x) satisfies the following:
  //      h (+/-1) = +/-1
  //      h'(+/-1) = 0
  //      h"(+/-1) = 0
  double logKhi = std::log(k_high);
  double logKlo = -std::log(k_low);
  double avg = 0.5 * (logKhi + logKlo), halfdiff = 0.5 * (logKhi - logKlo);
  double twox = x + x, twox2 = twox * twox;
  double alpha = 0.125 * twox * (twox2 * (3 * twox2 - 10.) + 15.);
  double ret = avg + alpha * halfdiff;
  return std::exp(ret * x);
}

void CombineHarvester::SetGroup(std::string const& name,
                                std::vector<std::string> const& patterns) {
  std::vector<boost::regex> rgx;
  for (auto const& pt : patterns) rgx.emplace_back(pt);
  for (auto it = params_.begin(); it != params_.end(); ++it) {
    std::string par = it->first;
    auto & groups = it->second->groups();
    if (groups.count(name)) continue;
    if (ch::contains_rgx(rgx, par)) {
      groups.insert(name);
    };
  }
}

void CombineHarvester::RemoveGroup(std::string const& name,
                                   std::vector<std::string> const& patterns) {
  std::vector<boost::regex> rgx;
  for (auto const& pt : patterns) rgx.emplace_back(pt);
  for (auto it = params_.begin(); it != params_.end(); ++it) {
    std::string par = it->first;
    auto & groups = it->second->groups();
    if (!groups.count(name)) continue;
    if (ch::contains_rgx(rgx, par)) {
      groups.erase(name);
    };
  }
}

void CombineHarvester::RenameGroup(std::string const& oldname,
                                   std::string const& newname) {
  for (auto it = params_.begin(); it != params_.end(); ++it) {
    auto & groups = it->second->groups();
    if (groups.count(oldname)) {
      groups.erase(oldname);
      groups.insert(newname);
    }
  }
}

void CombineHarvester::AddDatacardLineAtEnd(std::string const& line) {
  post_lines_.push_back(line);
}

void CombineHarvester::ClearDatacardLinesAtEnd() {
  post_lines_.clear();
}
}
