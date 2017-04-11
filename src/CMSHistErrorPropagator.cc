#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include <stdexcept>
#include <vector>
#include <ostream>
#include <memory>
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "vectorized.h"

CMSHistErrorPropagator::CMSHistErrorPropagator() : v(0), initialized_(false) {}

CMSHistErrorPropagator::CMSHistErrorPropagator(const char* name,
                                               const char* title,
                                               RooRealVar& x,
                                               RooArgList const& funcs,
                                               RooArgList const& coeffs)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      funcs_("funcs", "", this),
      coeffs_("coeffs", "", this),
      binpars_("binpars", "", this),
      sentry_(TString(name) + "_sentry", ""),
      binsentry_(TString(name) + "_binsentry", ""),
      v(0),
      initialized_(false) {
  funcs_.add(funcs);
  coeffs_.add(coeffs);
  // binpars_.add(binpars);
  // initialize();
}

CMSHistErrorPropagator::CMSHistErrorPropagator(
    CMSHistErrorPropagator const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      funcs_("funcs", this, other.funcs_),
      coeffs_("coeffs", this, other.coeffs_),
      binpars_("binpars", this, other.binpars_),
      bintypes_(other.bintypes_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.sentry_.GetName()), ""),
      binsentry_(name ? TString(name) + "_binsentry" : TString(other.binsentry_.GetName()), ""),
      v(other.v),
      initialized_(false) {
  // initialize();
}

void CMSHistErrorPropagator::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
  FNLOGC(std::cout, v) << "Initialising vectors\n";
  unsigned nf = funcs_.getSize();
  vfuncs_.resize(nf);
  vcoeffs_.resize(nf);
  for (unsigned i = 0; i < nf; ++i) {
    vfuncs_[i] = dynamic_cast<CMSHistFunc const*>(funcs_.at(i));
    vcoeffs_[i] = dynamic_cast<RooAbsReal const*>(coeffs_.at(i));
    auto sargs = vfuncs_[i]->getSentryArgs();
    // std::cout << sargs.get() << "\n";
    // sargs->Print();
    sentry_.addVars(*sargs);
  }
  unsigned nb = vfuncs_[0]->getCacheHisto().size();
  vbinpars_.resize(nb);
  if (bintypes_.size()) {
    for (unsigned j = 0, r = 0; j < nb; ++j) {
      vbinpars_[j].resize(bintypes_[j].size(), nullptr);
      for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
        if (bintypes_[j][i] >= 1 && bintypes_[j][i] < 4) {
          vbinpars_[j][i] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
          ++r;
        }
      }
    }
  }
  valsum_ = vfuncs_[0]->getCacheHisto();
  valsum_.Clear();
  cache_ = vfuncs_[0]->getCacheHisto();
  cache_.Clear();
  err2sum_.resize(nb, 0.);
  toterr_.resize(nb, 0.);
  binmods_.resize(nf, std::vector<double>(nb, 0.));
  scaledbinmods_.resize(nf, std::vector<double>(nb, 0.));
  coeffvals_.resize(nf, 0.);

  sentry_.addVars(coeffs_);
  binsentry_.addVars(binpars_);

  // sentry_.Print("v");

  sentry_.setValueDirty();
  binsentry_.setValueDirty();

  initialized_ = true;
}


void CMSHistErrorPropagator::updateCache(int eval) const {
  // IMPORTANT: Need to check if we are transitioning from eval=1 to eval=0


  // FNLOGC(std::cout, v) << "Start of function\n";
  initialize();

  // FNLOGC(std::cout, v) << "Sentry: " << sentry_.good() << "\n";
  if (!sentry_.good()) {
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      // FNLOGC(std::cout, v) << "Triggering updateCache() of function " << i << "\n";
      vfuncs_[i]->updateCache();
      coeffvals_[i] = vcoeffs_[i]->getVal();
    }

    valsum_.Clear();
    std::fill(err2sum_.begin(), err2sum_.end(), 0.);
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      vectorized::mul_add(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->getCacheHisto()[0]), &valsum_[0]);
      vectorized::mul_add_sqr(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->getBinErrors()[0]), &err2sum_[0]);
    }
    vectorized::sqrt(valsum_.size(), &err2sum_[0], &toterr_[0]);
    cache_ = valsum_;

    if (eval == 0) {
      for (unsigned j = 0; j < valsum_.size(); ++j) {
        if (bintypes_[j][0] == 1) {
          // if (v > 1) std::cout << "Bin " << j << "\n";
          // if (v > 1) printf(" | %.6f/%.6f/%.6f\n", valsum_[j], err2sum_[j], toterr_[j]);
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            if (err2sum_[j] > 0. && coeffvals_[i] > 0.) {
              double e =  vfuncs_[i]->getBinErrors()[j] * coeffvals_[i];
              binmods_[i][j] = (toterr_[j] *  e * e) / (err2sum_[j] * coeffvals_[i]);
            } else {
              binmods_[i][j] = 0.;
            }
            // if (v > 1) printf("%.6f   ", binmods_[i][j]);
          }
          // if (v > 1 ) printf("\n");
        }
      }
    }


    sentry_.reset();
    binsentry_.setValueDirty();
  }


  if (!binsentry_.good()) {
    // bintypes might have size == 0 if we never ran setupBinPars()
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      cache_[j] = valsum_[j];
      if (bintypes_[j][0] == 0) {
        continue;
      } else if (bintypes_[j][0] == 1) {
        if (data_.size() && toterr_[j] > 0. && vbinpars_[j][0]->isConstant()) {
          double err = toterr_[j] / valsum_[j];
          double b = toterr_[j] * err - 1.;
          double c = - data_[j] * err * err;
          double tmp = -0.5 * (b + copysign(1.0, b) * std::sqrt(b * b - 4. * c));
          double x1 = tmp;
          double x2 = c /tmp;

          // if (std::isnan(x1) || std::isnan(x2) || !((x1 > 0.) != (x2 > 0.))) {
          //   std::cout << "something went wrong!: " << "val: " << val << "\terr: " << err << "\tdat: " << data_[j] << "\n";
          //   std::cout << "something went wrong!: " << x1 << "\t" << x2 << "\n";
          // }

          // This part is a massive hack, but we can save time by only propagating the dirty flag to the constraint terms
          RooAbsArg::setDirtyInhibit(true);
          ((RooRealVar*)vbinpars_[j][0])->setVal((std::max(x1, x2) - 1.) / err);
          RooAbsArg::setDirtyInhibit(false);
          // Rely on the first value client being the constraint
          vbinpars_[j][0]->valueClientMIterator().next()->setValueDirty();
        }
        double x = vbinpars_[j][0]->getVal();
        cache_[j] += toterr_[j] * x;
        // Only fill the scaledbinmods if we're in eval == 0 mode (i.e. need to
        // propagate to wrappers)
        if (eval == 0) {
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            scaledbinmods_[i][j] = binmods_[i][j] * x;
          }
        }
      } else {
        for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
          if (bintypes_[j][i] == 2) {
            // Poisson: this is a multiplier on the process yield
            scaledbinmods_[i][j] = ((vbinpars_[j][i]->getVal() - 1.) *
                 vfuncs_[i]->getCacheHisto()[j] * coeffvals_[i]);
            cache_[j] += scaledbinmods_[i][j];
          } else if (bintypes_[j][i] == 3) {
            // Gaussian This is the addition of the scaled error
            scaledbinmods_[i][j] = vbinpars_[j][i]->getVal() * vfuncs_[i]->getBinErrors()[j] * coeffvals_[i];
            cache_[j] += scaledbinmods_[i][j];
          }
        }
      }
    }
    RooAbsArg::setDirtyInhibit(false);
    cache_.CropUnderflows();
    binsentry_.reset();
  }
}

RooArgList * CMSHistErrorPropagator::setupBinPars(double poissonThreshold) {
  RooArgList * res = new RooArgList();
  if (bintypes_.size()) {
    std::cout << "setupBinPars() already called for " << this->GetName() << "\n";
    return res;
  }

  // First initialize all the storage
  initialize();
  // Now fill the bin contents and errors
  updateCache(1); // the arg (1) forces updateCacheor to fill the caches for all bins

  bintypes_.resize(valsum_.size(), std::vector<unsigned>(1, 0));


  std::cout << std::string(60, '=') << "\n";
  std::cout << "Analysing bin errors for: " << this->GetName() << "\n";
  std::cout << "Poisson cut-off: " << poissonThreshold << "\n";
  std::set<unsigned> skip_idx;
  std::vector<std::string> skipped_procs;
  for (unsigned i = 0; i < vfuncs_.size(); ++i) {
    if (vfuncs_[i]->attributes().count("skipForErrorSum")) {
      skipped_procs.push_back(vfuncs_[i]->getStringAttribute("combine.process"));
      skip_idx.insert(i);
    }
  }
  if (skipped_procs.size()) {
    std::cout << "Processes excluded for sums:";
    for (auto &s: skipped_procs) std::cout << " " << s;
    std::cout << "\n";
  }
  std::cout << std::string(60, '=') << "\n";
  std::cout << TString::Format("%-10s %-15s %-15s %-30s\n", "Bin", "Contents", "Error", "Notes");

  for (unsigned j = 0; j < valsum_.size(); ++j) {
    std::cout << TString::Format("%-10i %-15f %-15f %-30s\n", j, valsum_[j], toterr_[j], "total sum");
    double sub_sum = 0.;
    double sub_err = 0.;
    // Check using a possible sub-set of bins
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      if (skip_idx.count(i)) {
        continue;
      }
      sub_sum += vfuncs_[i]->getCacheHisto()[j] * coeffvals_[i];
      sub_err += std::pow(vfuncs_[i]->getBinErrors()[j] * coeffvals_[i], 2.);;
    }
    sub_err = std::sqrt(sub_err);
    if (skipped_procs.size()) {
      std::cout << TString::Format("%-10i %-15f %-15f %-30s\n", j, sub_sum, sub_err, "excluding marked processes");
    }

    if (sub_err <= 0.) {
      std::cout << TString::Format("  %-30s\n", "=> Error is zero, ignore");
      std::cout << std::string(60, '-') << "\n";
      continue;
    }

    // Now check if we are below the poisson threshold
    double n = int(0.5 + ((sub_sum * sub_sum) / (sub_err * sub_err)));
    double alpha = valsum_[j] / n;
    std::cout << TString::Format(
        "%-10i %-15f %-15f %-30s\n", j, n, std::sqrt(n),
        TString::Format("Unweighted events, alpha=%f", alpha).Data());

    if (n <= poissonThreshold) {
      std::cout << TString::Format("  %-30s\n", "=> Number of weighted events is below poisson threshold");

      bintypes_[j].resize(vfuncs_.size(), 4);

      for (unsigned i = 0; i < vfuncs_.size(); ++i) {
        std::string proc =
            vfuncs_[i]->stringAttributes().count("combine.process")
                ? vfuncs_[i]->getStringAttribute("combine.process")
                : vfuncs_[i]->GetName();
        double v_p = vfuncs_[i]->getCacheHisto()[j];
        double e_p = vfuncs_[i]->getBinErrors()[j];
        std::cout << TString::Format("    %-20s %-15f %-15f %-30s\n", proc.c_str(), v_p, e_p, "");
        // relax the condition of v_p >= e_p slightly due to numerical rounding...
        // Possibilities:
        //    v_p = any   e_p <= 0       : skip
        //    v_p < 0     e_p > 0        : for now, skip but technically we should be able to handle this in the future
        //    v_p >= 0    e_p > v_p      : Create an additive gaussian constraint for this bin
        //    v_p > 0     0 < e_p <= v_p : do the poisson
        if (e_p <= 0.) {
          std::cout << TString::Format("      %-30s\n", "=> Error is zero, ignore");
          bintypes_[j][i] = 4;
        } else if (v_p < 0. && e_p > 0.) {
          std::cout << TString::Format("      %-30s\n", "=> Cannot handle negative content, ignore");
          bintypes_[j][i] = 4;
        } else if (v_p > 0. && e_p > 0. && v_p >= (e_p*0.999)) {
          double n_p_r = int(0.5 + ((v_p * v_p) / (e_p * e_p)));
          double alpha_p_r = v_p / n_p_r;
          std::cout << TString::Format(
              "    %-20s %-15f %-15f %-30s\n", "", n_p_r, std::sqrt(n_p_r),
              TString::Format("Unweighted events, alpha=%f", alpha_p_r).Data());
          if (n_p_r <= poissonThreshold) {
            double sigma = 7.;
            double rmin = 0.5*ROOT::Math::chisquared_quantile(ROOT::Math::normal_cdf_c(sigma), n_p_r * 2.);
            double rmax = 0.5*ROOT::Math::chisquared_quantile(1. - ROOT::Math::normal_cdf_c(sigma), n_p_r * 2. + 2.);
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 1, rmin/n_p_r, rmax/n_p_r);
            RooConstVar *cvar = new RooConstVar(TString::Format("%g", n_p_r), "", n_p_r);
            RooProduct *prod = new RooProduct(TString::Format("%s_prod", var->GetName()), "", RooArgList(*var, *cvar));
            prod->addOwnedComponents(RooArgSet(*var, *cvar));
            prod->setAttribute("createPoissonConstraint");
            res->addOwned(*prod);
            binpars_.add(*var);

            std::cout << TString::Format(
                "      => Product of %s[%.2f,%.2f,%.2f] and const [%.0f] to be poisson constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax(), cvar->getVal());
            bintypes_[j][i] = 2;
          } else {
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
            std::cout << TString::Format(
                "      => Parameter %s[%.2f,%.2f,%.2f] to be gaussian constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax());
            var->setAttribute("createGaussianConstraint");
            res->addOwned(*var);
            binpars_.add(*var);
            bintypes_[j][i] = 3;
          }
        } else if (v_p >= 0 && e_p > v_p) {
          RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
          std::cout << TString::Format(
              "      => Poisson not viable, %s[%.2f,%.2f,%.2f] to be gaussian constrained\n",
              var->GetName(), var->getVal(), var->getMin(), var->getMax());
          var->setAttribute("createGaussianConstraint");
          res->addOwned(*var);
          binpars_.add(*var);
          bintypes_[j][i] = 3;
        } else{
          std::cout << "      => ERROR: shouldn't be here\n";
        }
        std::cout << "  " << std::string(58, '-') << "\n";

      }
    } else if (toterr_[j] > 0.) {
      bintypes_[j][0] = 1;
      RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i", this->GetName(), j), "", 0, -7, 7);
      std::cout << TString::Format(
          "  => Total parameter %s[%.2f,%.2f,%.2f] to be gaussian constrained\n",
          var->GetName(), var->getVal(), var->getMin(), var->getMax());
      var->setAttribute("createGaussianConstraint");
      var->setAttribute("forBarlowBeeston");
      res->addOwned(*var);
      binpars_.add(*var);
    }
    std::cout << std::string(60, '-') << "\n";
  }

  // binpars_.add(*res);
  binsentry_.addVars(binpars_);
  binsentry_.setValueDirty();

  for (unsigned j = 0, r = 0; j < valsum_.size(); ++j) {
    vbinpars_[j].resize(bintypes_[j].size());
    for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
      if (bintypes_[j][i] >= 1) {
        vbinpars_[j][i] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
        ++r;
      }
    }
  }
  return res;
}


void CMSHistErrorPropagator::applyErrorShifts(unsigned idx,
                                              FastHisto const& nominal,
                                              FastHisto& result) {
  // We can skip the whole evaluation if there's nothing to evaluate
  // if (bintypes_.size() == 0) return;
  FNLOGC(std::cout, v) << "Start of function\n";
  updateCache(0);
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = nominal[i] + scaledbinmods_[idx][i];
  }
}

std::unique_ptr<RooArgSet> CMSHistErrorPropagator::getSentryArgs() const {
  // We can do this without initialising because we're going to hand over
  // the sentry directly
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
  std::unique_ptr<RooArgSet> args(new RooArgSet(sentry_, binsentry_));
  return args;
}

Double_t CMSHistErrorPropagator::evaluate() const {
  updateCache(1);
  return cache_.GetAt(x_);
}


void CMSHistErrorPropagator::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                    TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  updateCache();
  if (bintypes_.size()) {
    std::cout << ">> Bin types are:\n";
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      std::cout << ">> " << j << ":";
      for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
        std::cout << " " << bintypes_[j][i];
      }
      std::cout << "\n";
    }
  }
  std::cout << ">> Current cache:\n";
  valsum_.Dump();
  std::cout << ">> Current cache (bin scaled):\n";
  cache_.Dump();
  std::cout << ">> Sentry: " << sentry_.good() << "\n";
  sentry_.Print("v");

}

Int_t CMSHistErrorPropagator::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistErrorPropagator::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache(1);
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

void CMSHistErrorPropagator::setData(RooAbsData const& data) const {
  updateCache(1);
  data_.clear();
  data_.resize(cache_.fullsize(), 0.); // fullsize is important here is we've used activeBins
  RooArgSet obs(x_.arg());
  const RooRealVar& x = static_cast<const RooRealVar&>(*obs.first());
  for (int i = 0, n = data.numEntries(); i < n; ++i) {
    obs = *data.get(i);
    int idx = cache_.FindBin(x.getVal());
    data_[idx] = data.weight();
  }
}
