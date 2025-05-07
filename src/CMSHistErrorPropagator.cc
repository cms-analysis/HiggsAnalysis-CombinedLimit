#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistFuncWrapper.h"
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
#include "RooGaussian.h"
#include "RooProduct.h"
#include "vectorized.h"

#define HFVERBOSE 0

CMSHistErrorPropagator::CMSHistErrorPropagator() : initialized_(false) {}

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
      initialized_(false),
      last_eval_(-1) {
  funcs_.add(funcs);
  coeffs_.add(coeffs);
}

CMSHistErrorPropagator::CMSHistErrorPropagator(
    CMSHistErrorPropagator const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      funcs_("funcs", this, other.funcs_),
      coeffs_("coeffs", this, other.coeffs_),
      binpars_("binpars", this, other.binpars_),
      bintypes_(other.bintypes_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      binsentry_(name ? TString(name) + "_binsentry" : TString(other.GetName())+"_binsentry", ""),
      initialized_(false),
      last_eval_(-1) {
}

void CMSHistErrorPropagator::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
#if HFVERBOSE > 0
  std::cout << "Initializing vectors\n";
#endif
  unsigned nf = funcs_.getSize();
  vfuncs_.resize(nf);
  vcoeffs_.resize(nf);
  for (unsigned i = 0; i < nf; ++i) {
    vfuncs_[i] = dynamic_cast<CMSHistFunc const*>(funcs_.at(i));
    vcoeffs_[i] = dynamic_cast<RooAbsReal const*>(coeffs_.at(i));
    auto sargs = vfuncs_[i]->getSentryArgs();
    sentry_.addVars(*sargs);
  }
  unsigned nb = vfuncs_[0]->cache().size();
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
  valsum_ = vfuncs_[0]->cache();
  valsum_.Clear();
  cache_ = vfuncs_[0]->cache();
  cache_.Clear();
  err2sum_.resize(nb, 0.);
  toterr_.resize(nb, 0.);
  binmods_.resize(nf, std::vector<double>(nb, 0.));
  scaledbinmods_.resize(nf, std::vector<double>(nb, 0.));
  coeffvals_.resize(nf, 0.);

  sentry_.addVars(coeffs_);
  binsentry_.addVars(binpars_);

  sentry_.setValueDirty();
  binsentry_.setValueDirty();

  initialized_ = true;
}


void CMSHistErrorPropagator::updateCache(int eval) const {
  initialize();

#if HFVERBOSE > 0
  std::cout << "Sentry: " << sentry_.good() << "\n";
#endif
  if (!sentry_.good() || eval != last_eval_) {
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      vfuncs_[i]->updateCache();
      coeffvals_[i] = vcoeffs_[i]->getVal();
    }

    valsum_.Clear();
    std::fill(err2sum_.begin(), err2sum_.end(), 0.);
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      vectorized::mul_add(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->cache()[0]), &valsum_[0]);
      vectorized::mul_add_sqr(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->errors()[0]), &err2sum_[0]);
    }
    vectorized::sqrt(valsum_.size(), &err2sum_[0], &toterr_[0]);
    cache_ = valsum_;

    if (eval == 0 && bintypes_.size()) {
      for (unsigned j = 0; j < valsum_.size(); ++j) {
        if (bintypes_[j][0] == 1) {
#if HFVERBOSE > 1
          std::cout << "Bin " << j << "\n";
          printf(" | %.6f/%.6f/%.6f\n", valsum_[j], err2sum_[j], toterr_[j]);
#endif
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            if (err2sum_[j] > 0. && coeffvals_[i] > 0.) {
              double e =  vfuncs_[i]->errors()[j] * coeffvals_[i];
              binmods_[i][j] = (toterr_[j] *  e * e) / (err2sum_[j] * coeffvals_[i]);
            } else {
              binmods_[i][j] = 0.;
            }
#if HFVERBOSE > 1
            printf("%.6f   ", binmods_[i][j]);
#endif
          }
#if HFVERBOSE > 1
          printf("\n");
#endif
        }
      }
    }


    sentry_.reset();
    binsentry_.setValueDirty();
  }


  if (!binsentry_.good() || eval != last_eval_) {
    runBarlowBeeston();
    // bintypes might have size == 0 if we never ran setupBinPars()
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      cache_[j] = valsum_[j];
      if (bintypes_[j][0] == 0) {
        continue;
      } else if (bintypes_[j][0] == 1) {
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
                 vfuncs_[i]->cache()[j]);
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          } else if (bintypes_[j][i] == 3) {
            // Gaussian This is the addition of the scaled error
            scaledbinmods_[i][j] = vbinpars_[j][i]->getVal() * vfuncs_[i]->errors()[j];
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          }
        }
      }
    }
    cache_.CropUnderflows();
    binsentry_.reset();
  }

  last_eval_ = eval;
}

void CMSHistErrorPropagator::runBarlowBeeston() const {
  if (!bb_.init) return;
  RooAbsArg::setDirtyInhibit(true);

  const unsigned n = bb_.use.size();
  for (unsigned j = 0; j < n; ++j) {
    bb_.dat[j] = data_[bb_.use[j]];
    bb_.valsum[j] = valsum_[bb_.use[j]] * cache_.GetWidth(bb_.use[j]);
    bb_.toterr[j] = toterr_[bb_.use[j]] * cache_.GetWidth(bb_.use[j]);
  }
  // This pragma statement tells (modern) gcc that loop can be safely
  // vectorized
  #pragma GCC ivdep
  for (unsigned j = 0; j < n; ++j) {
    bb_.b[j] = bb_.toterr[j] + (bb_.valsum[j] / bb_.toterr[j]) - bb_.gobs[j];
    bb_.c[j] = bb_.valsum[j] - bb_.dat[j] - (bb_.valsum[j] / bb_.toterr[j]) * bb_.gobs[j];
    bb_.tmp[j] = -0.5 * (bb_.b[j] + copysign(1.0, bb_.b[j]) * std::sqrt(bb_.b[j] * bb_.b[j] - 4. * bb_.c[j]));
    bb_.x1[j] = bb_.tmp[j];
    bb_.x2[j] = bb_.c[j] / bb_.tmp[j];
    bb_.res[j] = std::max(bb_.x1[j], bb_.x2[j]);
  }
  for (unsigned j = 0; j < n; ++j) {
    if (toterr_[bb_.use[j]] > 0.) bb_.push_res[j]->setVal(bb_.res[j]);
  }
  RooAbsArg::setDirtyInhibit(false);
  for (RooAbsArg *arg : bb_.dirty_prop) {
    arg->setValueDirty();
  }
}

void CMSHistErrorPropagator::setAnalyticBarlowBeeston(bool flag) const {
  // Clear it if it's already initialised
  if (bb_.init && flag) return;
  if (bb_.init && !flag) {
    for (unsigned i = 0; i < bb_.push_res.size(); ++i) {
      bb_.push_res[i]->setConstant(false);
    }
    bb_.use.clear();
    bb_.dat.clear();
    bb_.valsum.clear();
    bb_.toterr.clear();
    bb_.err.clear();
    bb_.b.clear();
    bb_.c.clear();
    bb_.tmp.clear();
    bb_.x1.clear();
    bb_.x2.clear();
    bb_.res.clear();
    bb_.gobs.clear();
    bb_.dirty_prop.clear();
    bb_.push_res.clear();
    bb_.init = false;
  }
  if (flag && data_.size()) {
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      if (bintypes_[j][0] == 1 && !vbinpars_[j][0]->isConstant()) {
        bb_.use.push_back(j);
        double gobs_val = 0.;
        for (RooAbsArg * arg : vbinpars_[j][0]->valueClients()) {
          if (arg == this || arg == &binsentry_) {
            // std::cout << "Skipping " << this << " " << this->GetName() << "\n";
          } else {
            // std::cout << "Adding " << arg << " " << arg->GetName() << "\n";
            bb_.dirty_prop.insert(arg);
            auto as_gauss = dynamic_cast<RooGaussian*>(arg);
            if (as_gauss) {
              auto gobs = dynamic_cast<RooAbsReal*>(as_gauss->findServer(TString(vbinpars_[j][0]->GetName())+"_In"));
              if (gobs) gobs_val = gobs->getVal();
            }
          }
        }
        bb_.gobs.push_back(gobs_val);
        bb_.push_res.push_back((RooRealVar*)vbinpars_[j][0]);
        bb_.push_res.back()->setConstant(true);
      }
    }
    unsigned n = bb_.use.size();
    bb_.dat.resize(n);
    bb_.valsum.resize(n);
    bb_.toterr.resize(n);
    bb_.err.resize(n);
    bb_.b.resize(n);
    bb_.c.resize(n);
    bb_.tmp.resize(n);
    bb_.x1.resize(n);
    bb_.x2.resize(n);
    bb_.res.resize(n);
    bb_.init = true;
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
  updateCache(1); // the arg (1) forces updateCache to fill the caches for all bins

  bintypes_.resize(valsum_.size(), std::vector<unsigned>(1, 0));


  std::cout << std::string(60, '=') << "\n";
  std::cout << "Analyzing bin errors for: " << this->GetName() << "\n";
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
      sub_sum += vfuncs_[i]->cache()[j] * coeffvals_[i];
      sub_err += std::pow(vfuncs_[i]->errors()[j] * coeffvals_[i], 2.);;
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
    double n = std::floor(0.5 + ((sub_sum * sub_sum) / (sub_err * sub_err)));
    double alpha = valsum_[j] / n;
    std::cout << TString::Format(
        "%-10i %-15f %-15f %-30s\n", j, n, std::sqrt(n),
        TString::Format("Unweighted events, alpha=%f", alpha).Data());

    if (n <= poissonThreshold) {
      std::cout << TString::Format("  %-30s\n", "=> Number of weighted events is below Poisson threshold");

      bintypes_[j].resize(vfuncs_.size(), 4);

      for (unsigned i = 0; i < vfuncs_.size(); ++i) {
        std::string proc =
            vfuncs_[i]->stringAttributes().count("combine.process")
                ? vfuncs_[i]->getStringAttribute("combine.process")
                : vfuncs_[i]->GetName();
        double v_p = vfuncs_[i]->cache()[j];
        double e_p = vfuncs_[i]->errors()[j];
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
          double n_p_r = std::floor(0.5 + ((v_p * v_p) / (e_p * e_p)));
          double alpha_p_r = v_p / n_p_r;
          std::cout << TString::Format(
              "    %-20s %-15f %-15f %-30s\n", "", n_p_r, std::sqrt(n_p_r),
              TString::Format("Unweighted events, alpha=%f", alpha_p_r).Data());
          if (n_p_r <= poissonThreshold) {
            double sigma = 7.;
            double rmin = 0.5*ROOT::Math::chisquared_quantile(ROOT::Math::normal_cdf_c(sigma), n_p_r * 2.);
            double rmax = 0.5*ROOT::Math::chisquared_quantile(1. - ROOT::Math::normal_cdf_c(sigma), n_p_r * 2. + 2.);
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", n_p_r, rmin, rmax);
            RooConstVar *cvar = new RooConstVar(TString::Format("%g", 1. / n_p_r), "", 1. / n_p_r);
            RooProduct *prod = new RooProduct(TString::Format("%s_prod", var->GetName()), "", RooArgList(*var, *cvar));
	    RooArgSet ownedComps;
	    ownedComps.add(*prod);
	    ownedComps.add(*cvar);
            var->addOwnedComponents(ownedComps);
            var->setAttribute("createPoissonConstraint");
            res->addOwned(*var);
            binpars_.add(*prod);

            std::cout << TString::Format(
                "      => Product of %s[%.2f,%.2f,%.2f] and const [%.4f] to be Poisson constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax(), cvar->getVal());
            bintypes_[j][i] = 2;
          } else {
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
            std::cout << TString::Format(
                "      => Parameter %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax());
            var->setAttribute("createGaussianConstraint");
            res->addOwned(*var);
            binpars_.add(*var);
            bintypes_[j][i] = 3;
          }
        } else if (v_p >= 0 && e_p > v_p) {
          RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
          std::cout << TString::Format(
              "      => Poisson not viable, %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
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
          "  => Total parameter %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
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
      if (bintypes_[j][i] >= 1 && bintypes_[j][i] < 4) {
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
  // std::cout << "Start of function\n";
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
  return cache().GetAt(x_);
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
  if (allVars.find(x_.arg().GetName()) == nullptr) allVars.add(x_.arg());
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistErrorPropagator::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache(1);
      return cache().IntegralWidth();
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

RooArgList CMSHistErrorPropagator::wrapperList() const {
  RooArgList result;
  for (int i = 0; i < funcs_.getSize(); ++i) {
    CMSHistFunc const* hf = dynamic_cast<CMSHistFunc const*>(funcs_.at(i));
    if (hf) {
      CMSHistFuncWrapper const* wrapper = hf->wrapper();
      if (wrapper) result.add(*wrapper);
    }
  }
  return result;
}

std::map<std::string, Double_t> CMSHistErrorPropagator::getProcessNorms() const {

      std::map<std::string, Double_t> vals_;
      RooArgList clist(coefList());
      RooArgList plist(funcList());
      /*if (plist.getSize() == 1) {
         CMSHistErrorPropagator *err = dynamic_cast<CMSHistErrorPropagator*>(plist.at(0));
         if (err) {
           clist.removeAll();
           plist.removeAll();
           clist.add(err->coefList());
           plist.add(err->wrapperList());
         }
      }
      */
      for (int i = 0, n = clist.getSize(); i < n; ++i) {
        RooAbsReal *coeff = (RooAbsReal *) clist.at(i);
        std::string coeffName = coeff->GetName();
        RooAbsReal* shape = (RooAbsReal*)plist.at(i);
        std::unique_ptr<RooArgSet> myobs(shape->getObservables(*x_));
        TString normProdName = TString::Format("%s", coeff->GetName());
        RooAbsReal * normProd = nullptr;
        if (coeff->ownedComponents()) {
          normProd = dynamic_cast<RooAbsReal*>(coeff->ownedComponents()->find(normProdName));
        }
        if (!normProd) {
          RooAbsReal* integral = shape->createIntegral(*myobs);
      	  RooArgList normProdInputs;
      	  normProdInputs.add(*integral);
      	  normProdInputs.add(*coeff);
          normProd = new RooProduct(normProdName, "", normProdInputs);
          normProd->addOwnedComponents(normProdInputs);
        }
        vals_[normProdName.Data()] = normProd->getVal();
      }
      return vals_;
}
#undef HFVERBOSE

