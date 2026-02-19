#include "../interface/CMSHistSum.h"
#include "../interface/CMSHistFuncWrapper.h"
#include "../interface/ProcessNormalization.h"
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
#include <charconv>
#include <string>

#define HFVERBOSE 0

bool CMSHistSum::enable_fast_vertical_ = false;

CMSHistSum::CMSHistSum() : initialized_(false), fast_mode_(0) {}

CMSHistSum::CMSHistSum(const char* name,
                                               const char* title,
                                               RooRealVar& x,
                                               RooArgList const& funcs,
                                               RooArgList const& coeffs)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      morphpars_("morphpars", "", this),
      coeffpars_("coeffpars", "", this),
      binpars_("binpars", "", this),
      n_procs_(0),
      n_morphs_(0),
      sentry_(TString(name) + "_sentry", ""),
      binsentry_(TString(name) + "_binsentry", ""),
      initialized_(false),
      analytic_bb_(false),
      fast_mode_(0),
      external_morphs_("external_morphs", "", this) {

  n_procs_ = funcs.getSize();
  assert(n_procs_ == coeffs.getSize());

  coeffpars_.add(coeffs);

  // Resize vectors
  binerrors_.resize(n_procs_);
  vtype_.resize(n_procs_);
  vsmooth_par_.resize(n_procs_);

  // First pass through the list of processes - we don't know the number of
  // vertical morphing parameters yet.
  std::map<std::string, RooRealVar const*> all_vmorphs;
  std::vector<std::map<std::string, int>> vmorph_idx_mapping(n_procs_);
  int n_storage = 0;
  for (int ip = 0; ip < n_procs_; ++ip) {
    CMSHistFunc const* func = dynamic_cast<CMSHistFunc const*>(funcs.at(ip));

    // Verify that each func does not use horizontal morphing or rebinning
    if (func->hmorphs_.getSize() > 0 || func->rebin_) {
      throw std::runtime_error("CMSHistSum does not support horizontal morphing or on-the-fly rebinning");
    }

    // Add this CMSHistFunc to vfuncstmp_ - this is not persisted but will be used to create bbb
    // parameters
    vfuncstmp_.push_back(func);

    // Need to initialize our cache_ with the correct binning
    if (ip == 0) {
      cache_ = func->cache();
      cache_.Clear();
    }

    // Determine the size we'll need for storage_
    n_storage += func->storage_.size();

    // Construct a list of the unique vertical morphing parameters and
    // make a per-process index of the local list positions of each vmorph
    for (int iv = 0; iv < func->vmorphs_.getSize(); ++iv) {
      RooRealVar const* rrv = static_cast<RooRealVar const*>(func->vmorphs_.at(iv));
      all_vmorphs[rrv->GetName()] = rrv;
      vmorph_idx_mapping[ip][rrv->GetName()] = iv;
    }

    // Copy bin errors, vtype and vsmooth_par
    binerrors_[ip] = func->binerrors_;
    vtype_[ip] = func->vtype_;
    vsmooth_par_[ip] = func->vsmooth_par_;
  }

  // Populate a list of the morphing parameter names,
  // and add each parameter to the RooListProxy
  std::vector<std::string> morph_names;
  for (auto const& it : all_vmorphs) {
    morph_names.push_back(it.first);
    morphpars_.add(*it.second);
  }

  // Now we know the number of morphing parameters
  n_morphs_ = morph_names.size();


  // Now we populate storage_, by stealing from each CMSHistFunc in turn
  storage_.clear();
  storage_.resize(n_storage);
  auto copy_it = storage_.begin();
  for (int ip = 0; ip < n_procs_; ++ip) {
    CMSHistFunc const* func = dynamic_cast<CMSHistFunc const*>(funcs.at(ip));
#if HFVERBOSE > 0
    std::cout << "Before copy: size = " << storage_.size() << ", distance = " << std::distance(storage_.begin(), copy_it) << "\n";
#endif
    copy_it = std::copy(func->storage_.begin(), func->storage_.end(), copy_it);
#if HFVERBOSE > 0
    std::cout << "After copy: size = " << storage_.size() << ", distance = " << std::distance(storage_.begin(), copy_it) << "\n";
#endif
  }

  // Next populate the storage_ index positions for fast lookup later
  process_fields_ = std::vector<int>(n_procs_, 0);
  vmorph_fields_ = std::vector<int>(n_procs_ * n_morphs_, -1);
  int current_proc = 0;
  for (int ip = 0; ip < n_procs_; ++ip) {
    process_fields_[ip] = current_proc;
    // current_proc += 1;
    for (int iv = 0; iv < n_morphs_; ++iv) {
      auto it = vmorph_idx_mapping[ip].find(morph_names[iv]);
      if (it != vmorph_idx_mapping[ip].end()) {
        // Position for this vmorph
        morphField(ip, iv) = current_proc + 1 + it->second * 2;
        // Convert to sum/diff
        unsigned idx = process_fields_[ip];
        unsigned idxLo = morphField(ip, iv) + 0;
        unsigned idxHi = morphField(ip, iv) + 1;
        FastTemplate lo = storage_[idxLo];
        FastTemplate hi = storage_[idxHi];
        if (vtype_[ip] == CMSHistFunc::VerticalSetting::QuadLinear) {
          hi.Subtract(storage_[idx]);
          lo.Subtract(storage_[idx]);
        } else if (vtype_[ip] == CMSHistFunc::VerticalSetting::LogQuadLinear) {
          hi.LogRatio(storage_[idx]);
          lo.LogRatio(storage_[idx]);
        }
        FastTemplate::SumDiff(hi, lo, storage_[idxLo], storage_[idxHi]);
      }
    }
    current_proc += (1 +vmorph_idx_mapping[ip].size() * 2);
  }

#if HFVERBOSE > 0
  std::cout << "FUNCTIONS\n";
  for (int ip = 0; ip < n_procs_; ++ip) {
    std::cout << ip << ": " << funcs.at(ip)->GetName() << std::endl;
  }
  std::cout << "MORPHS\n";
  for (int imorph = 0; imorph < n_morphs_; ++imorph) {
    std::cout << imorph << ": " << morph_names[imorph] << std::endl;
  }
  for (int ip = 0; ip < n_procs_; ++ip) {
    std::cout << TString::Format("%4i: %6i |", ip, process_fields_[ip]);
    for (int iv = 0; iv < n_morphs_; ++iv) {
      std::cout << TString::Format("%6i", morphField(ip, iv));
    }
    std::cout << std::endl;
  }
#endif

  // exit(0);
}

CMSHistSum::CMSHistSum(
    CMSHistSum const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      morphpars_("morphpars", this, other.morphpars_),
      coeffpars_("coeffpars", this, other.coeffpars_),
      binpars_("binpars", this, other.binpars_),
      n_procs_(other.n_procs_),
      n_morphs_(other.n_morphs_),
      storage_(other.storage_),
      process_fields_(other.process_fields_),
      vmorph_fields_(other.vmorph_fields_),
      binerrors_(other.binerrors_),
      vtype_(other.vtype_),
      vsmooth_par_(other.vsmooth_par_),
      bintypes_(other.bintypes_),
      cache_(other.cache_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      binsentry_(name ? TString(name) + "_binsentry" : TString(other.GetName())+"_binsentry", ""),
      initialized_(false),
      fast_mode_(0),
      external_morphs_("external_morphs", this, other.external_morphs_),
      external_morph_indices_(other.external_morph_indices_)
{
      initialize();
}

void CMSHistSum::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");

#if HFVERBOSE > 0
  std::cout << "Initializing vectors\n";
#endif

  vmorphpars_.resize(n_morphs_);
  vcoeffpars_.resize(n_procs_);
  compcache_.resize(n_procs_);

  for (int ip = 0; ip < n_procs_; ++ip) {
    vcoeffpars_[ip] = dynamic_cast<RooAbsReal const*>(coeffpars_.at(ip));
    compcache_[ip] = cache_;
    compcache_[ip].Clear();
  }
  for (int iv = 0; iv < n_morphs_; ++iv) {
    vmorphpars_[iv] = dynamic_cast<RooAbsReal const*>(morphpars_.at(iv));
  }

  unsigned nb = cache_.size();
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

  valsum_ = cache_;
  staging_ = cache_;
  valsum_.Clear();
  cache_.Clear();
  staging_.Clear();
  err2sum_.resize(nb, 0.);
  toterr_.resize(nb, 0.);
  binmods_.resize(n_procs_, std::vector<double>(nb, 0.));
  scaledbinmods_.resize(n_procs_, std::vector<double>(nb, 0.));
  coeffvals_.resize(n_procs_, 0.);

  sentry_.addVars(morphpars_);
  sentry_.addVars(coeffpars_);
  binsentry_.addVars(binpars_);

  for (const auto* morph : external_morphs_) {
    RooArgSet* deps = morph->getParameters({*x_});
    sentry_.addVars(*deps);
    delete deps;
  }

  sentry_.setValueDirty();
  binsentry_.setValueDirty();

  initialized_ = true;
}

void CMSHistSum::updateMorphs() const {
  // set up pointers ahead of time for quick loop
  std::vector<CMSExternalMorph*> process_morphs(compcache_.size(), nullptr);
  // if any external morphs are dirty, disable fast_mode_
  for(size_t i=0; i < external_morph_indices_.size(); ++i) {
    auto* morph = static_cast<CMSExternalMorph*>(external_morphs_.at(i));
    process_morphs[external_morph_indices_[i]] = morph;
    if (morph->hasChanged()) {
      fast_mode_ = 0;
    }
  }
  // If we're not in fast mode, need to reset all the compcache_
  #if HFVERBOSE > 0
  std::cout << "fast_mode_ = " << fast_mode_ << std::endl;
  #endif
  for (unsigned ip = 0; ip < compcache_.size(); ++ip) {
    if (fast_mode_ == 0) {
      compcache_[ip].CopyValues(storage_[process_fields_[ip]]);
      if ( process_morphs[ip] != nullptr ) {
        auto& extdata = process_morphs[ip]->batchGetBinValues();
        for(size_t ibin=0; ibin<extdata.size(); ++ibin) {
          compcache_[ip][ibin] *= extdata[ibin];
        }
      }
      if (vtype_[ip] == CMSHistFunc::VerticalSetting::LogQuadLinear) {
        compcache_[ip].Log();
      }
    }
  }
  int n_morphs = vmorphpars_.size();

  if (vertical_prev_vals_.size() == 0) {
    vertical_prev_vals_.resize(n_morphs);
  }
  // Loop through vmorphs
  for (unsigned iv = 0; iv < vmorphpars_.size(); ++iv) {
    double x = vmorphpars_[iv]->getVal();
    if (fast_mode_ == 1 && (x == vertical_prev_vals_[iv])) {
      #if HFVERBOSE > 0
      std::cout << "Skipping " << vmorphpars_[iv]->GetName() << ", prev = now = " << x << std::endl;
      #endif
      continue;
    }
    #if HFVERBOSE > 0
    if (fast_mode_ == 1) {
      std::cout << "Updating " << vmorphpars_[iv]->GetName() << ", prev =  " << vertical_prev_vals_[iv] << ", now = " << x << std::endl;
    }
    #endif


    // If in fast made, check if the value has changed since the last eval,
    // and if it hasn't, skip it

    // For each vmorph need to know the list of processes we apply to
    for (unsigned ip = 0; ip < compcache_.size(); ++ip) {
      int code = vmorph_fields_[ip * n_morphs + iv];
      if (code == -1) continue;
      if (fast_mode_ == 1) {
        double xold = vertical_prev_vals_[iv];
        compcache_[ip].DiffMeld(storage_[code + 1], storage_[code + 0], 0.5*x, smoothStepFunc(x, ip), 0.5*xold, smoothStepFunc(xold, ip));
      } else {
        compcache_[ip].Meld(storage_[code + 1], storage_[code + 0], 0.5*x, smoothStepFunc(x, ip));
      }
    }
    vertical_prev_vals_[iv] = x;
  }

  if (enable_fast_vertical_) fast_mode_ = 1;
}

inline double CMSHistSum::smoothStepFunc(double x, int const& ip) const {
  if (fabs(x) >= vsmooth_par_[ip]) return x > 0 ? +1 : -1;
  double xnorm = x / vsmooth_par_[ip];
  double xnorm2 = xnorm * xnorm;
  return 0.125 * xnorm * (xnorm2 * (3. * xnorm2 - 10.) + 15);
}

void CMSHistSum::updateCache() const {
  initialize();

#if HFVERBOSE > 0
  std::cout << "Sentry: " << sentry_.good() << "\n";
#endif
  if (!sentry_.good()) {
    #if HFVERBOSE > 0
      std::cout << "Calling updateMorphs\n";
    #endif
    updateMorphs();
    for (unsigned i = 0; i < vcoeffpars_.size(); ++i) {
      // vfuncs_[i]->updateCache();
      coeffvals_[i] = vcoeffpars_[i]->getVal();
    }
    #if HFVERBOSE > 0
      std::cout << "Updated coeffs\n";
    #endif

    valsum_.Clear();
    std::fill(err2sum_.begin(), err2sum_.end(), 0.);
    for (unsigned i = 0; i < vcoeffpars_.size(); ++i) {
      staging_ = compcache_[i];
      if (vtype_[i] == CMSHistFunc::VerticalSetting::LogQuadLinear) {
        staging_.Exp();
        staging_.Scale(storage_[process_fields_[i]].Integral() / staging_.Integral());
      }
      staging_.CropUnderflows();
      vectorized::mul_add(valsum_.size(), coeffvals_[i], &(staging_[0]), &valsum_[0]);
      vectorized::mul_add_sqr(valsum_.size(), coeffvals_[i], &(binerrors_[i][0]), &err2sum_[0]);
    }
    vectorized::sqrt(valsum_.size(), &err2sum_[0], &toterr_[0]);
    cache_ = valsum_;
    #if HFVERBOSE > 0
      std::cout << "Updated cache\n";
    #endif
    sentry_.reset();
    binsentry_.setValueDirty();
  }


  if (!binsentry_.good()) {
    // bintypes might have size == 0 if we never ran setupBinPars()
    #if HFVERBOSE > 0
      std::cout << "Assigning bin shifts\n";
    #endif
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      cache_[j] = valsum_[j];
      if (bintypes_[j][0] == 0) {
        continue;
      } else if (bintypes_[j][0] == 1) {
        double x = vbinpars_[j][0]->getVal();
        cache_[j] += toterr_[j] * x;
      } else {
        for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
          if (bintypes_[j][i] == 2) {
            // Poisson: this is a multiplier on the process yield
            scaledbinmods_[i][j] = ((vbinpars_[j][i]->getVal() - 1.) *
                 compcache_[i][j]);
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          } else if (bintypes_[j][i] == 3) {
            // Gaussian This is the addition of the scaled error
            scaledbinmods_[i][j] = vbinpars_[j][i]->getVal() * binerrors_[i][j];
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          }
        }
      }
    }
    #if HFVERBOSE > 0
      std::cout << "Done assigning bin shifts\n";
    #endif
    cache_.CropUnderflows();
    binsentry_.reset();
  }
}

void CMSHistSum::runBarlowBeeston() const {
  updateCache();
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
    double b = bb_.toterr[j] + (bb_.valsum[j] / bb_.toterr[j]) - bb_.gobs[j];
    double c = bb_.valsum[j] - bb_.dat[j] - (bb_.valsum[j] / bb_.toterr[j]) * bb_.gobs[j];
    double tmp = -0.5 * (b + std::copysign(1.0, b) * std::sqrt(b * b - 4. * c));
    double x1 = tmp;
    double x2 = c / tmp;
    bb_.res[j] = std::max(x1, x2);
  }
  for (unsigned j = 0; j < n; ++j) {
    if (toterr_[bb_.use[j]] > 0.) bb_.push_res[j]->setVal(bb_.res[j]);
  }
}

void CMSHistSum::setAnalyticBarlowBeeston(bool flag) const {
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
    bb_.res.clear();
    bb_.gobs.clear();
    bb_.push_res.clear();
    bb_.init = false;
  }
  if (flag && data_.size()) {
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      if (bintypes_[j][0] == 1 && !vbinpars_[j][0]->isConstant()) {
        bb_.use.push_back(j);
        double gobs_val = 0.;
        for (RooAbsArg *arg : vbinpars_[j][0]->valueClients()) {
          if (arg == this || arg == &binsentry_) {
            // std::cout << "Skipping " << this << " " << this->GetName() << "\n";
          } else {
            // std::cout << "Adding " << arg << " " << arg->GetName() << "\n";
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
    bb_.res.resize(n);
    bb_.init = true;
  }
}


RooArgList * CMSHistSum::setupBinPars(double poissonThreshold) {
  RooArgList * res = new RooArgList();
  if (bintypes_.size()) {
    std::cout << "setupBinPars() already called for " << this->GetName() << "\n";
    return res;
  }

  // First initialize all the storage
  initialize();
  // Now fill the bin contents and errors
  updateCache(); // the arg (1) forces updateCache to fill the caches for all bins

  bintypes_.resize(valsum_.size(), std::vector<unsigned>(1, 0));


  std::cout << std::string(60, '=') << "\n";
  std::cout << "Analyzing bin errors for: " << this->GetName() << "\n";
  std::cout << "Poisson cut-off: " << poissonThreshold << "\n";
  std::set<unsigned> skip_idx;
  std::vector<std::string> skipped_procs;
  for (unsigned i = 0; i < vfuncstmp_.size(); ++i) {
    vfuncstmp_[i]->updateCache();
    if (vfuncstmp_[i]->attributes().count("skipForErrorSum")) {
      skipped_procs.push_back(vfuncstmp_[i]->getStringAttribute("combine.process"));
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
    for (unsigned i = 0; i < vfuncstmp_.size(); ++i) {
      if (skip_idx.count(i)) {
        continue;
      }
      sub_sum += vfuncstmp_[i]->cache()[j] * coeffvals_[i];
      sub_err += std::pow(vfuncstmp_[i]->errors()[j] * coeffvals_[i], 2.);;
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
      std::cout << TString::Format("  %-30s\n", "=> Number of weighted events is below Poisson threshold");

      bintypes_[j].resize(vfuncstmp_.size(), 4);

      for (unsigned i = 0; i < vfuncstmp_.size(); ++i) {
        std::string proc =
            vfuncstmp_[i]->stringAttributes().count("combine.process")
                ? vfuncstmp_[i]->getStringAttribute("combine.process")
                : vfuncstmp_[i]->GetName();
        double v_p = vfuncstmp_[i]->cache()[j];
        double e_p = vfuncstmp_[i]->errors()[j];
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


std::unique_ptr<RooArgSet> CMSHistSum::getSentryArgs() const {
  // We can do this without initialising because we're going to hand over
  // the sentry directly
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
  std::unique_ptr<RooArgSet> args(new RooArgSet(sentry_, binsentry_));
  return args;
}

Double_t CMSHistSum::evaluate() const {
  updateCache();
  return cache().GetAt(x_);
}

std::map<std::string, Double_t> CMSHistSum::getProcessNorms() const {
  std::map<std::string, Double_t> vals_;
  initialize();
  if (!sentry_.good()) {
    updateMorphs();

    for (unsigned i = 0; i < vcoeffpars_.size(); ++i) {
      Double_t coeffval = vcoeffpars_[i]->getVal();    

      staging_ = compcache_[i];
      if (vtype_[i] == CMSHistFunc::VerticalSetting::LogQuadLinear) {
        staging_.Exp();
      }
      staging_.CropUnderflows();

      double const * valarray = &(staging_[0]);
      Double_t valsum = 0;
      for(unsigned j=0; j < valsum_.size(); ++j ){
        valsum += coeffval * valarray[j];
      }

      vals_[vcoeffpars_[i]->GetName()] = valsum;
    }
  }
  return vals_;
}

void CMSHistSum::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
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

Int_t CMSHistSum::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (allVars.find(x_.arg().GetName()) == nullptr) allVars.add(x_.arg());
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistSum::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache();
      return cache().IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

void CMSHistSum::setData(RooAbsData const& data) const {
  updateCache();
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

void CMSHistSum::EnableFastVertical() {
  enable_fast_vertical_ = true;
}

void CMSHistSum::injectExternalMorph(int idx, CMSExternalMorph& morph) {
  if ( idx >= coeffpars_.getSize() ) {
    throw std::runtime_error("Process index larger than number of processes in CMSHistSum");
  }
  if ( morph.batchGetBinValues().size() != cache_.size() ) {
    throw std::runtime_error("Mismatched binning between external morph and CMSHistSum");
    // equal edges are user responsibility for now
  }

  for (auto other_idx : external_morph_indices_) {
    if ( idx == other_idx ) {
      external_morphs_.replace(external_morphs_[idx], morph);
      return;
    }
  }
  external_morph_indices_.push_back(idx);
  external_morphs_.add(morph);
}

#undef HFVERBOSE

ClassImp(CMSHistSum)


#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitHS3/JSONIO.h"

using RooFit::Detail::JSONNode;

// Helper: try to pick a reasonable process name.
static std::string cms_hist_sum_process_name(const CMSHistSum &hs, unsigned ip)
{
  auto *coeff = hs.coeffAt(ip);   
  
  if (coeff && coeff->stringAttributes().count("combine.process")) {
    return coeff->getStringAttribute("combine.process");
  }
  
  if(coeff){
    return coeff->GetName();
  }
  
  return "proc_" + std::to_string(ip);
}

bool tryExportHistFactory(RooJSONFactoryWSTool *tool, const CMSHistSum *hs, JSONNode &elem){

    // if we have external morphs injected, you probably want to error out for now.
    if (hs->hasExternalMorphs()) { // optional accessor; otherwise remove and always allow
      RooJSONFactoryWSTool::error("CMSHistSum '" + std::string(hs->GetName()) +
                                 "': exporting with external morphs is not implemented (cannot express as HistFactory modifiers).");
    }

    // Axes: HistFactory exporter uses RooRealVar axis metadata.
    // Use the independent variable name and its range.
    RooAbsReal const &xAbs = hs->getXVar(); 
    auto const *x = dynamic_cast<RooRealVar const *>(&xAbs);
    if (!x) {
      RooJSONFactoryWSTool::error("CMSHistSum '" + std::string(hs->GetName()) +
                                 "': x is not a RooRealVar; cannot export axis metadata.");
    }

    {
      auto &axes = elem["axes"].set_seq();
      auto &ax = axes.append_child().set_map();
      std::string axname = x->GetName();
      RooJSONFactoryWSTool::testValidName(axname, false);
      ax["name"] << axname;
      ax["nbins"] << static_cast<int>(hs->nBins());
      ax["min"] << x->getMin();
      ax["max"] << x->getMax();
    }

    // Samples array
    auto &samples = elem["samples"].set_seq();

    const unsigned nProcs = hs->nProcs();
    const unsigned nMorphs = hs->nMorphs();
    const unsigned nBins  = hs->nBins();

    for (unsigned ip = 0; ip < nProcs; ++ip) {
      const std::string sampleName = cms_hist_sum_process_name(*hs, ip);

      auto &s = samples.append_child().set_map();
      s["name"] << sampleName;

      // --- Nominal histogram (bin yields)
      std::vector<double> nom(nBins, 0.0);
      std::vector<double> err(nBins, 0.0);

      auto const &nomT = hs->nominalTemplate(ip);

      for (unsigned ib = 0; ib < nBins; ++ib) {
        const double w = hs->binWidth(ib);
        nom[ib] = nomT[ib] * w;
        err[ib] = hs->binErrors(ip)[ib] * w;
      }

      auto &dataNode = s["data"].set_map();
      RooJSONFactoryWSTool::exportArray(nBins, nom.data(), dataNode["contents"]);

      // If there are nonzero bin errors, export them and add staterror modifier (HistFactory style).
      bool anyErr = false;
      for (double e : err) { if (e != 0.0) { anyErr = true; break; } }
      if (anyErr) {
        RooJSONFactoryWSTool::exportArray(nBins, err.data(), dataNode["errors"]);
      }

      // --- Modifiers
      auto &mods = s["modifiers"].set_seq();

      auto *coeff = hs->coeffAt(ip);
      if (!coeff) {
	// do nothing
      } else if (auto * rcv = dynamic_cast<RooConstVar const *>(coeff)) {
	// the const'ness of the const parameter is handled elsewhere, we can just pretend like it's a normal NF here
	auto &m = mods.append_child().set_map();
	m["name"] << rcv->GetName();
	m["parameter"] << rcv->GetName();
	m["type"] << "normfactor";
      } else if (auto *rrv = dynamic_cast<RooRealVar const *>(coeff)) {
	auto &m = mods.append_child().set_map();
	m["name"] << rrv->GetName();
	m["parameter"] << rrv->GetName();
	m["type"] << "normfactor";
      } else if (auto* pn = dynamic_cast<const ProcessNormalization*>(coeff)) {
	double nom = pn->nominalValue();
	auto variations = pn->sigmaVariationsAll();
	for(const auto& v:variations){
	  auto &m = mods.append_child().set_map();
	  m["name"] << std::get<0>(v);
	  m["parameter"] << std::get<0>(v);
	  m["type"] << "overallsys";
	  auto& data = m["data"].set_map();
	  data["lo"] << std::get<1>(v)/nom;
	  data["hi"] << std::get<2>(v)/nom;
	}
	if(nom != 1.){
	  auto &m = mods.append_child().set_map();
	  m["name"] << sampleName+"_norm_nominal";
	  m["parameter"] << nom;
	  m["type"] << "normfactor";
	}
      } else {
	auto &m = mods.append_child().set_map();
	m["name"] << coeff->GetName();
	m["type"] << "custom";
	tool->queueExport(*coeff);
      }

      // Vertical morphs -> histosys modifiers
      for (unsigned iv = 0; iv < nMorphs; ++iv) {
        if(!hs->isMorphUsed(ip, iv)) continue;

        RooAbsReal const *par = hs->morphPar(iv);
        if (!par) continue;
        const std::string sysName = par->GetName();

        std::vector<double> lo(nBins, 0.0);
        std::vector<double> hi(nBins, 0.0);

        auto const &upT  = hs->upTemplate(ip, iv); 
        auto const &dnT = hs->downTemplate(ip, iv);
	for(size_t i=0; i<nBins; ++i){
	  lo[i] = dnT[i];
	  hi[i] = upT[i];
	}
	  
        const auto vsetting = hs->vtype(ip);

        auto &m = mods.append_child().set_map();
        m["name"] << sysName;
        m["parameter"] << sysName;
        m["type"] << "histosys";

        auto &d = m["data"].set_map();
        RooJSONFactoryWSTool::exportArray(nBins, lo.data(), d["lo"].set_map()["contents"]);
        RooJSONFactoryWSTool::exportArray(nBins, hi.data(), d["hi"].set_map()["contents"]);

        tool->queueExport(*par);
      }

      // Staterror modifier (if errors were exported)
      if (anyErr) {
        auto &m = mods.append_child().set_map();
        m["name"] << "mcstat";
        m["type"] << "staterror";
      }
    }
    
    return true;
}


#include "RooRealSumPdf.h"
#include "RooProdPdf.h"

const CMSHistSum* getCMSHistSum(const RooRealSumPdf* pdf){
  CMSHistSum* sum = 0;
  if(pdf->funcList().size() == 1 && pdf->coefList().size() == 1 && pdf->coefList().at(0)->GetName() == std::string("ONE") && (sum = dynamic_cast<CMSHistSum*>(pdf->funcList().at(0))) && sum)
    return sum;
  return nullptr;
}

const CMSHistSum* getCMSHistSum(const RooProdPdf* pdf){
  RooRealSumPdf* sum = 0;
  if(pdf->getComponents()->getSize() == 1 && (sum = dynamic_cast<RooRealSumPdf*>(pdf->getComponents()->first())) && sum){
    return getCMSHistSum(sum);
  }
  return nullptr;
}

class CMSHistSumStreamer_CMSHistSum : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto sum = dynamic_cast<const CMSHistSum*>(p);
     if(sum){
       auto retval = tryExportHistFactory(tool,sum,elem);
       if(!retval) return false;
       elem["type"] << key();
       return true;
     }
     return false;
   }
};


class CMSHistSumStreamer_SumPdf : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto sum = getCMSHistSum(dynamic_cast<const RooRealSumPdf*>(p));
     if(sum){
       auto retval = tryExportHistFactory(tool,sum,elem);
       if(retval){
	 elem["type"] << key();
	 return true;
       }
       std::cout << "failure" << std::endl;
     }
     return false;
   }
};


class CMSHistSumStreamer_ProdPdf : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto sum = getCMSHistSum(dynamic_cast<const RooProdPdf*>(p));
     if(sum){
       auto retval = tryExportHistFactory(tool,sum,elem);
       if(!retval) return false;
       elem["type"] << key();
       return true;
     }
     return false;
   }     
};

#include "static_execute.h"

STATIC_EXECUTE([]() {
  using namespace RooFit::JSONIO;
  
  registerExporter<CMSHistSumStreamer_SumPdf>(RooRealSumPdf::Class(), true);
  registerExporter<CMSHistSumStreamer_ProdPdf>(RooProdPdf::Class(), true);
  registerExporter<CMSHistSumStreamer_CMSHistSum>("CMSHistSum", false);    
 });

namespace {

/// Helper: read an array of doubles
static std::vector<double> readDoubleArray(RooFit::Detail::JSONNode const &arr, std::string const &what)
{
   std::vector<double> out;
   out.reserve(arr.num_children());
   for(const auto& c:arr.children()){
     out.push_back(c.val_double());
   }
   return out;
}

FastHisto readFastHisto(const JSONNode &n,
					 const std::string &name,
					 RooArgSet const &vars)
{
  if (!n.has_child("contents"))
    RooJSONFactoryWSTool::error("no contents given");
  
  JSONNode const &contents = n["contents"];
  if (!contents.is_seq())
    RooJSONFactoryWSTool::error("contents are not in list form");
  
  JSONNode const *errors = nullptr;
  if (n.has_child("errors")) {
    errors = &n["errors"];
    if (!errors->is_seq())
      RooJSONFactoryWSTool::error("errors are not in list form");
  }

  if (vars.size() != 1) {
    std::stringstream errMsg;
    errMsg << "readFastHisto only supports 1D, got vars.size()=" << vars.size();
    RooJSONFactoryWSTool::error(errMsg.str());
  }

  auto *x = dynamic_cast<RooRealVar *>(vars.first());
  if (!x)
    RooJSONFactoryWSTool::error("readFastHisto expects a RooRealVar observable");
  
  const RooAbsBinning &binning = x->getBinning();
  const int nbins = binning.numBins();
  
  if (static_cast<int>(contents.num_children()) != nbins) {
    std::stringstream errMsg;
    errMsg << "inconsistent bin numbers: contents=" << contents.num_children() << ", nbins=" << nbins;
    RooJSONFactoryWSTool::error(errMsg.str());
  }
  if (errors && static_cast<int>(errors->num_children()) != nbins) {
    std::stringstream errMsg;
    errMsg << "inconsistent bin numbers: errors=" << errors->num_children() << ", nbins=" << nbins;
    RooJSONFactoryWSTool::error(errMsg.str());
  }
  
  // Read contents
  std::vector<double> vals;
  vals.reserve(contents.num_children());
  for (auto const &cont : contents.children()) {
    vals.push_back(cont.val_double());
  }
  
  // Build bin edges (nbins+1)
  std::vector<double> edges;
  edges.reserve(static_cast<std::size_t>(nbins) + 1u);
  for (int i = 0; i < nbins; ++i) edges.push_back(binning.binLow(i));
  edges.push_back(binning.binHigh(nbins - 1));
  
  // Construct FastHisto directly
  return FastHisto(edges, vals);
}
  

/// Helper: ensure/lookup observable RooRealVar by name
static RooRealVar &getObservable(RooJSONFactoryWSTool *tool, std::string const &obsName,
                                 double xmin, double xmax)
{
   // If it exists already, use it.
   if (auto *arg = tool->workspace()->arg(obsName.c_str())) {
      auto *rrv = dynamic_cast<RooRealVar *>(arg);
      if (!rrv)
         throw std::runtime_error("histfactory_dist importer: observable '" + obsName + "' exists but is not RooRealVar");
      return *rrv;
   }

   // Otherwise create it (range from axis definition)
   auto *rrv = new RooRealVar(obsName.c_str(), obsName.c_str(), 0.0, xmin, xmax);
   tool->workspace()->import(*rrv, RooFit::RecycleConflictNodes());
   return *dynamic_cast<RooRealVar *>(tool->workspace()->arg(obsName.c_str()));
}

/// Helper: ensure/lookup a parameter (RooRealVar) used by normfactor modifiers
static RooRealVar &getOrMakeParameter(RooJSONFactoryWSTool *tool, std::string const &parName,
                                     double init = 1.0, double lo = 0.0, double hi = 10.0)
{
   if (auto *arg = tool->workspace()->arg(parName.c_str())) {
      auto *rrv = dynamic_cast<RooRealVar *>(arg);
      if (!rrv)
         throw std::runtime_error("histfactory_dist importer: parameter '" + parName + "' exists but is not RooRealVar");
      return *rrv;
   }
   auto *rrv = new RooRealVar(parName.c_str(), parName.c_str(), init, lo, hi);
   tool->workspace()->import(*rrv, RooFit::RecycleConflictNodes());
   return *dynamic_cast<RooRealVar *>(tool->workspace()->arg(parName.c_str()));
}

bool parse_double(const std::string& input, double& out)
{
    const char* begin = input.data();
    const char* end   = input.data() + input.size();

    double tmp;
    auto [ptr, ec] = std::from_chars(begin, end, tmp);

    // Successful conversion requires:
    // 1) No error
    // 2) Entire string was consumed
    if((ec == std::errc() && ptr == end)){
      out = tmp;
      return true;
    }
    return false;
}

inline bool approxSymmetricKappas(double kLo, double kHi,
                                 double absLogTol = 1e-10,
                                 double relLogTol = 1e-10)
{
    // kappas must be positive and finite
    if (!(kLo > 0.0) || !(kHi > 0.0)) return false;
    if (!std::isfinite(kLo) || !std::isfinite(kHi)) return false;

    const double a = std::log(kLo);
    const double b = std::log(kHi);

    // symmetric if a ≈ -b  <=>  (a + b) ≈ 0
    const double diff  = std::abs(a + b);
    const double scale = std::max(std::abs(a), std::abs(b)); // 0 if both ~1
    return diff <= std::max(absLogTol, relLogTol * scale);
}

/// The actual importer
class CMSHistSumHistFactoryImporter final : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, RooFit::Detail::JSONNode const &elem) const override
   {
      using RooFit::Detail::JSONNode;

      // Basic validation
      if (!elem.has_child("name"))
         return false;

      const std::string distName = elem["name"].val();

      // Axes: we only support 1D here (which is what Combine uses for CMS_th1x)
      if (!elem.has_child("axes") || elem["axes"].num_children() != 1)
         throw std::runtime_error("CMSHistSum importer: '" + distName + "' must have exactly one axis");

      RooArgSet observables = RooJSONFactoryWSTool::readAxes(elem);
      RooRealVar* x = dynamic_cast<RooRealVar*>(observables.first());
      
      // Samples
      if (!elem.has_child("samples") || elem["samples"].num_children() == 0)
         throw std::runtime_error("CMSHistSum importer: '" + distName + "' has no samples");

      RooArgList funcs;
      RooArgList coeffs;

      // Ensure ONE exists (fixed at 1) for the wrapper RooRealSumPdf
      RooAbsReal *ONE = dynamic_cast<RooAbsReal *>(tool->workspace()->arg("ONE"));
      if (!ONE) {
         auto *one = new RooConstVar("ONE", "ONE", 1.0);
         tool->workspace()->import(*one, RooFit::RecycleConflictNodes());
         ONE = dynamic_cast<RooAbsReal *>(tool->workspace()->arg("ONE"));
      }      

      for(const auto&s:elem["samples"].children()){
	const std::string sampleName = s["name"].val();

         // data.contents (required)
         if (!s.has_child("data") || !s["data"].has_child("contents"))
            throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' missing data.contents");
         std::vector<double> contents = readDoubleArray(s["data"]["contents"], sampleName + ".data.contents");
         if (static_cast<int>(contents.size()) != x->numBins())
            throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' contents size != nbins");

         // data.errors (optional, but we need it for staterror)
         std::vector<double> errors;
         if (s["data"].has_child("errors")) {
            errors = readDoubleArray(s["data"]["errors"], sampleName + ".data.errors");
            if (static_cast<int>(errors.size()) != x->numBins())
               throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' errors size != nbins");
         } else {
	   errors.assign(x->numBins(), 0.0);
         }

         // Build a TH1D with that binning (uniform, defined by axis min/max/nbins)
         const std::string hname = distName + "_" + sampleName + "_nominal";
         auto h = std::make_unique<TH1D>(hname.c_str(), hname.c_str(), x->numBins(), x->getMin(), x->getMax());
         for (int ib = 1; ib <= x->numBins(); ++ib) {
            h->SetBinContent(ib, contents[ib - 1]);
            h->SetBinError(ib, errors[ib - 1]);
         }

	 const std::string funcName = distName + "_histosys_" + sampleName;
         auto *func = new CMSHistFunc(funcName.c_str(), funcName.c_str(), *x, *h);

         RooArgList coeffFactors;
	 struct VMorph {
	   FastHisto low;
	   FastHisto high;
	   VMorph(const FastHisto& low_, const FastHisto& high_):low(std::move(low_)),high(std::move(high_)){}
	 };
	 RooArgList vmorph_pars;
	 std::vector<VMorph> vmorphs;
	 RooArgList overall_pars;
	 std::vector< std::pair<double,double> > overall_vals;
	 std::string norm_nom_name;
	 double norm_nom = 1.;
	 
         if (s.has_child("modifiers")) {
	   for(const auto& m:s["modifiers"].children()){
	     const std::string mtype = m["type"].val();
	     const std::string mname = m["name"].val();

               if (mtype == "normfactor") {
                  // Typical: a floating POI or nuisance scaling a process yield
 		  if(parse_double(m["parameter"].val(),norm_nom)){
		    norm_nom_name = m["parameter"].val();
		  } else {
		    RooRealVar &p = getOrMakeParameter(tool, mname, 1.0, 0.0, 10.0);
		    coeffFactors.add(p);
		  }
	       } else if (mtype == "custom"){
                  // This might be normalization factor that is given by a ProcessNormalization or RooFormulaVar
                  RooAbsReal *obj = tool->workspace()->function(mname);
		  if(!obj){
		    throw std::runtime_error(
					     "CMSHistSum importer: custom modifier '" + mname + "' on sample '" + sampleName +
					     "' not yet imported, there seems to be a priority/ordering issue");
		  }
		  coeffFactors.add(*obj);
               } else if (mtype == "staterror") {
                  // Already encoded via per-bin errors; CMSHistSum will handle BBB if configured.
                  continue;
	       } else if (mtype == "overallsys"){
		 auto *parameter = m.find("parameter");
		 std::string parname(parameter ? parameter->val() : "alpha_" + mname);
		 RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10., 10.0);
		 auto& data = m["data"];
		 double lo = data["lo"].val_double();
		 double hi = data["hi"].val_double();
		 overall_pars.add(p);
		 overall_vals.push_back(std::make_pair(lo,hi));
		 
	       } else if (mtype == "histosys"){
		 // collect all the histosys first, will assign them later when we know how many there are
		 auto *parameter = m.find("parameter");
		 std::string parname(parameter ? parameter->val() : "alpha_" + mname);
		 RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10., 10.0);		 
		 vmorph_pars.add(p);
		 auto &data = m["data"];
		 auto low  = readFastHisto(data["lo"], mname + "Low_" + sampleName, *x);
		 auto high = readFastHisto(data["hi"], mname + "High_" + sampleName, *x);		 
		 vmorphs.emplace_back(std::move(low),std::move(high));
               } else {
                  // Shape and other modifiers: not safely mappable without CMSHistFunc’s public morphing API.
                  throw std::runtime_error(
                      "CMSHistSum importer: modifier type '" + mtype + "' on sample '" + sampleName +
                      "' is not supported by this CMSHistSum importer yet");
               }
	   }
         }

	 assert(vmorphs.size() == vmorph_pars.size());

	 // actually assign the morphs
	 func->setVerticalMorphs(vmorph_pars);
	 func->prepareStorage();
	 for(size_t i=0; i<vmorphs.size(); ++i){
	   func->setShape(0, 0, i+1, 0, vmorphs[i].low);
	   func->setShape(0, 0, i+1, 1, vmorphs[i].high);
	 }
         funcs.add(*func);
	 
	 // actually assign the overall sys
	 if(overall_pars.size() > 0){
	   RooAbsReal* norm = ONE;
	   if(norm_nom != 1.){
	     norm = new RooRealVar(norm_nom_name.c_str(),norm_nom_name.c_str(),norm_nom);
	   }
	   const std::string pnName = distName + "_overallsys_" + sampleName;
	   auto* pn = new ProcessNormalization(pnName.c_str(), pnName.c_str(), *norm);
	   for(size_t i=0; i<overall_pars.size(); ++i){
	     auto par = dynamic_cast<RooAbsReal*>(&(overall_pars[i]));
	     if(approxSymmetricKappas(overall_vals[1].first,overall_vals[i].second)){
	       pn->addLogNormal(overall_vals[i].first * norm_nom, *par);
	     } else {
	       pn->addAsymmLogNormal(overall_vals[i].first * norm_nom, overall_vals[i].second * norm_nom, *par);
	     }
	   }
	   coeffFactors.add(*pn);
	 }
 
         RooAbsReal *coeff = ONE; // default coefficient 1
         if (coeffFactors.getSize() == 1) {
            coeff = dynamic_cast<RooAbsReal *>(coeffFactors.at(0));
         } else if (coeffFactors.getSize() > 1) {
            const std::string coeffName = distName + "_" + sampleName + "_rate";
            auto *prod = new RooProduct(coeffName.c_str(), coeffName.c_str(), coeffFactors);
            tool->workspace()->import(*prod, RooFit::RecycleConflictNodes());
            coeff = dynamic_cast<RooAbsReal *>(tool->workspace()->function(coeffName.c_str()));
            if (!coeff)
               throw std::runtime_error("CMSHistSum importer: failed to import coefficient product '" + coeffName + "'");
         }
         coeffs.add(*coeff);
      }

      // Construct CMSHistSum
      const std::string sumName = distName + "_cmshistsum";
      auto *sum = new CMSHistSum(sumName.c_str(), sumName.c_str(), *x, funcs, coeffs);
      tool->workspace()->import(*sum, RooFit::RecycleConflictNodes());
      sum = dynamic_cast<CMSHistSum *>(tool->workspace()->function(sumName.c_str()));
      if (!sum)
         throw std::runtime_error("CMSHistSum importer: failed to import CMSHistSum '" + sumName + "'");

      // Wrap into RooRealSumPdf with ONE coefficient
      // RooRealSumPdf(name,title, funcs, coeffs) uses:
      //   sum over i coeff_i * func_i ; here: one term: 1 * CMSHistSum
      RooArgList pdfFuncs;
      RooArgList pdfCoeffs;
      pdfFuncs.add(*sum);
      pdfCoeffs.add(*ONE);
      
      auto *pdf = new RooRealSumPdf(distName.c_str(), distName.c_str(), pdfFuncs, pdfCoeffs, true /*extended*/);
      tool->workspace()->import(*pdf, RooFit::RecycleConflictNodes());

      return true;
   }
};

STATIC_EXECUTE([]() {
   using namespace RooFit::JSONIO;

   // priority=true so this wins over the standard HistFactoryImporter
   registerImporter<CMSHistSumHistFactoryImporter>("histfactory_dist", /*priority*/ true);
});

} // namespace




