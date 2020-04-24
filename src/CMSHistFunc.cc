#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFuncWrapper.h"
#include "HiggsAnalysis/CombinedLimit/interface/Accumulators.h"
#include <vector>
#include <ostream>
#include <memory>
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "TH1F.h"
#include "TVector.h"
#include "TMatrix.h"
#include "vectorized.h"
#include "TMath.h"

#define HFVERBOSE 0

bool CMSHistFunc::enable_fast_vertical_ = false;

CMSHistFunc::CMSHistFunc() {
  morph_strategy_ = 0;
  initialized_ = false;
  htype_ = HorizontalType::Integral;
  mtype_ = MomentSetting::NonLinearPosFractions;
  vtype_ = VerticalSetting::QuadLinear;
  divide_by_width_ = true;
  rebin_ = false;
  vsmooth_par_ = 1.0;
  fast_vertical_ = false;
}

CMSHistFunc::CMSHistFunc(const char* name, const char* title, RooRealVar& x,
                         TH1 const& hist, bool divide_by_width)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      vmorphs_("vmorphs", "", this),
      hmorphs_("hmorphs", "", this),
      vmorph_sentry_(TString(name) + "_vmorph_sentry", ""),
      hmorph_sentry_(TString(name) + "_hmorph_sentry", ""),
      cache_(hist),
      binerrors_(FastTemplate(hist.GetNbinsX())),
      morph_strategy_(0),
      initialized_(false),
      rebin_(false),
      htype_(HorizontalType::Integral),
      mtype_(MomentSetting::Linear),
      vtype_(VerticalSetting::QuadLinear),
      divide_by_width_(divide_by_width),
      vsmooth_par_(1.0),
      fast_vertical_(false) {
  if (divide_by_width_) {
    for (unsigned i = 0; i < cache_.size(); ++i) {
      cache_[i] /= cache_.GetWidth(i);
    }
  }
  prepareStorage();  // Prepare storage of size 1 -> the cache_ will be copied in there
  for (unsigned i = 0; i < cache_.size(); ++i) {
    binerrors_[i] = hist.GetBinError(i + 1);
    if (divide_by_width_) {
      binerrors_[i] /= cache_.GetWidth(i);
    }
  }

  if (x.getBins() < int(cache_.size()) &&
      std::fabs(x.getMin() - cache_.GetEdge(0)) < 1E-7 &&
      std::fabs(x.getMax() - cache_.GetEdge(cache_.size())) < 1E-7) {
    std::cout << "[" << this->GetName()
              << "]: histogram binning is finer than observable binning, auto "
                 "rebinning will be enabled\n";
    rebin_ = true;
    TH1F tmp("tmp", "", x.getBins(), x.getBinning().array());
    tmp.SetDirectory(0);
    rebin_cache_ = FastHisto(tmp);
    for (unsigned i = 0; i < cache_.size(); ++i) {
      unsigned target_bin = rebin_cache_.FindBin((cache_.GetEdge(i) + cache_.GetEdge(i+1)) / 2.);
      // std::cout << "Rebin: " << i << " [" << cache_.GetEdge(i) << ","
      //           << cache_.GetEdge(i + 1) << "] --> " << target_bin << " ["
      //           << rebin_cache_.GetEdge(target_bin) << ","
      //           << rebin_cache_.GetEdge(target_bin + 1) << "]\n";
      rebin_scheme_.push_back(target_bin);
    }
    applyRebin();
    binerrors_ = FastTemplate(rebin_cache_.size());
    binerrors_.Clear();
    for (unsigned i = 0; i < cache_.size(); ++i) {
      binerrors_[rebin_scheme_[i]] += std::pow(hist.GetBinError(i + 1), 2.);
    }
    for (unsigned i = 0; i < rebin_cache_.size(); ++i) {
      binerrors_[i] = std::sqrt(binerrors_[i]);
      if (divide_by_width_) {
        binerrors_[i] /= rebin_cache_.GetWidth(i);
      }
    }
  }
}

CMSHistFunc::CMSHistFunc(CMSHistFunc const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      vmorphs_("vmorphs", this, other.vmorphs_),
      hmorphs_("hmorphs", this, other.hmorphs_),
      hpoints_(other.hpoints_),
      vmorph_sentry_(name ? TString(name) + "_vmorph_sentry"
                          : TString(other.vmorph_sentry_.GetName()),
                     ""),
      hmorph_sentry_(name ? TString(name) + "_hmorph_sentry"
                          : TString(other.hmorph_sentry_.GetName()),
                     ""),
      cache_(other.cache_),
      rebin_cache_(other.rebin_cache_),
      binerrors_(other.binerrors_),
      storage_(other.storage_),
      morph_strategy_(other.morph_strategy_),
      initialized_(false),
      rebin_(other.rebin_),
      rebin_scheme_(other.rebin_scheme_),
      htype_(other.htype_),
      mtype_(other.mtype_),
      vtype_(other.vtype_),
      divide_by_width_(other.divide_by_width_),
      vsmooth_par_(other.vsmooth_par_),
      fast_vertical_(false) {
  // initialize();
}

void CMSHistFunc::initialize() const {
  if (initialized_) return;
  vmorph_sentry_.SetName(TString(this->GetName()) + "_vmorph_sentry");
  hmorph_sentry_.SetName(TString(this->GetName()) + "_hmorph_sentry");
  hmorph_sentry_.addVars(hmorphs_);
  vmorph_sentry_.addVars(vmorphs_);
  hmorph_sentry_.setValueDirty();
  vmorph_sentry_.setValueDirty();
  vmorphs_vec_.resize(vmorphs_.getSize());
  for (int i = 0; i < vmorphs_.getSize(); ++i) {
    vmorphs_vec_[i] = (RooAbsReal*)(&vmorphs_[i]);
  }

  setGlobalCache();

  initialized_ = true;
}

void CMSHistFunc::setGlobalCache() const {
  unsigned n_hpoints = 1;
  if (hmorphs_.getSize() == 1) {
    n_hpoints = hpoints_[0].size();

    // This portion of code adapted from:
    // RooMomentMorph.cxx
    TVectorD dm(n_hpoints);
    global_.m.ResizeTo(n_hpoints, n_hpoints);
    TMatrixD M(n_hpoints, n_hpoints);

    for (unsigned i = 0; i < hpoints_[0].size(); ++i) {
      dm[i] = hpoints_[0][i] - hpoints_[0][0];
      M(i, 0) = 1.;
      if (i > 0) M(0, i) = 0.;
    }
    for (unsigned i = 1; i < hpoints_[0].size(); ++i) {
      for (unsigned j = 1; j < hpoints_[0].size(); ++j) {
        M(i, j) = std::pow(dm[i], (double)j);
      }
    }
    global_.m = M.Invert();

    global_.c_scale.resize(hpoints_[0].size(), 0.);
    global_.c_sum.resize(hpoints_[0].size(), 0.);
    // End
    global_.means.resize(hpoints_[0].size(), 0.);
    global_.sigmas.resize(hpoints_[0].size(), 0.);
    global_.slopes.resize(hpoints_[0].size(), 0.);
    global_.offsets.resize(hpoints_[0].size(), 0.);
  }
}

void CMSHistFunc::setActiveBins(unsigned bins) {
  cache_.SetActiveSize(bins);
  for (unsigned i = 0; i < storage_.size(); ++i) {
    storage_[i].SetActiveSize(bins);
  }
  resetCaches();
  setGlobalCache();
}


void CMSHistFunc::resetCaches() {
  mcache_.clear();
  hmorph_sentry_.setValueDirty();
  vmorph_sentry_.setValueDirty();
}

std::unique_ptr<RooArgSet> CMSHistFunc::getSentryArgs() const {
  initialize();
  std::unique_ptr<RooArgSet> vargs(vmorph_sentry_.getComponents());
  std::unique_ptr<RooArgSet> hargs(hmorph_sentry_.getComponents());
  vargs->add(*hargs);
  return vargs;
}

void CMSHistFunc::addHorizontalMorph(RooAbsReal & hvar, TVectorD hpoints) {
  hmorphs_.add(hvar);
  std::vector<double> points(hpoints.GetNrows());
  for (int i = 0; i < hpoints.GetNrows(); ++i) {
    points[i] = hpoints[i];
  }
  hpoints_.push_back(points);
  hmorph_sentry_.addVars(RooArgList(hvar));
  hmorph_sentry_.setValueDirty();
}

void CMSHistFunc::setVerticalMorphs(RooArgList const& vvars) {
  for (int i = 0; i < vvars.getSize(); ++i) {
    vmorphs_.add(*vvars.at(i));
  }
  vmorph_sentry_.addVars(vvars);
  vmorph_sentry_.setValueDirty();
}

void CMSHistFunc::prepareStorage() {
  storage_.clear();
  assert(hmorphs_.getSize() <= 1);
  unsigned n_hpoints = 1;
  if (hmorphs_.getSize() == 1) {
    n_hpoints = hpoints_[0].size();
  }
  unsigned n_vpoints = vmorphs_.getSize();
  storage_.resize(n_hpoints * (1 + n_vpoints * 2));
#if HFVERBOSE > 0
  std::cout << "Storage size set to: " << storage_.size() << "\n";
#endif
  storage_[getIdx(0, 0, 0, 0)] = cache_;
}

void CMSHistFunc::setShape(unsigned hindex, unsigned hpoint, unsigned vindex,
                           unsigned vpoint, TH1 const& hist) {
  unsigned idx = getIdx(hindex, hpoint, vindex, vpoint);
#if HFVERBOSE > 0
  std::cout << "hindex: " << hindex << " hpoint: " << hpoint
            << " vindex: " << vindex << " vpoint: " << vpoint
            << " mapped to idx: " << idx << "\n";
#endif
  storage_.at(idx) = FastHisto(hist);
  if (divide_by_width_) {
    for (unsigned i = 0; i < storage_[idx].size(); ++i) {
      storage_.at(idx)[i] /= cache_.GetWidth(i);
    }
  }
}

void CMSHistFunc::updateCache() const {
  initialize();

  // Quick escape if cache is up-to-date
  if (hmorph_sentry_.good() && vmorph_sentry_.good()) return;

#if HFVERBOSE > 0
  std::cout << "Morphing strategy 0\n";
  std::cout << "Number of horizontal morphs: " << hmorphs_.getSize() << "\n";
  std::cout << "Horizontal morph sentry: " << hmorph_sentry_.good() << "\n";
  std::cout << "Vertical morph sentry: " << vmorph_sentry_.good() << "\n";
  std::cout << "single_point,p1,p2: " << global_.single_point << " " << global_.p1 << " " << global_.p2 << "\n";
#endif

  // Figure out what we're doing:
  //  - singlePoint p1 OR
  //  - h-morphing between p1 and p2
  bool step1 = !hmorph_sentry_.good();
  bool step2 = step1 || !vmorph_sentry_.good();

  if (step1 && hmorphs_.getSize() == 1) {
    double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();
#if HFVERBOSE > 0
    std::cout << "Updating horizontal points...\n";
    std::cout << "Horizontal morph parameter is: " << val << "\n";
#endif
    auto upper = std::lower_bound(hpoints_[0].begin(), hpoints_[0].end(), val);
    auto lower = std::upper_bound(hpoints_[0].begin(), hpoints_[0].end(), val);
    if (upper == hpoints_[0].begin()) {
      global_.p1 = 0;
      global_.single_point = true;
#if HFVERBOSE > 0
      std::cout << "Considered morphing parameter point below the lowest template parameter value. Performing single point morphing using template from lowest parameter value."<< "\n";
#endif
    } else if (lower == hpoints_[0].end()) {
      global_.p1 = hpoints_[0].size() - 1;
      global_.single_point = true;
#if HFVERBOSE > 0
      std::cout << "Considered morphing parameter point above the highest template parameter value. Performing single point morphing using template from highest parameter value."<< "\n";
#endif
    } else {
      lower = upper;
      --lower;
      global_.p1 = lower - hpoints_[0].begin();
      global_.p2 = upper - hpoints_[0].begin();
      global_.single_point = false;
      if (htype_ == HorizontalType::Closest) {
        global_.single_point = true;
        if (fabs(*upper - val) <= fabs(val - *lower)) {
          global_.p1 = upper - hpoints_[0].begin();
        } else {
          global_.p1 = lower - hpoints_[0].begin();
        }
      }
    }
#if HFVERBOSE > 0
    std::cout << "single_point,p1,p2: " << global_.single_point << " " << global_.p1 << " " << global_.p2 << "\n";
#endif
  }

  if (step1 || step2) {
    if (mcache_.size() == 0) mcache_.resize(storage_.size());
  }

  if (step1) {
    fast_vertical_ = false;
  }

  if (morph_strategy_ == 0) {

    if (step1 && !global_.single_point) {
#if HFVERBOSE > 0
      std::cout << "Checking step 1\n";
#endif
      double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();

      if (htype_ == HorizontalType::Moment) {
        updateMomentFractions(val);  // updates fractions in global cache, only depends on hmorph par
      }

      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        for (int vi = 0; vi < (v == 0 ? 1 : 2); ++vi) {
          if (htype_ == HorizontalType::Moment) {
#if HFVERBOSE > 0
            std::cout << "Doing vpoint,vindex = " << v << "\t" << vi << "\n";
#endif
            for (unsigned hi = 0; hi < hpoints_[0].size(); ++hi) {
              // define vec of mean vals
              unsigned idx = getIdx(0, hi, v, vi);
              if (!mcache_[idx].meansig_set) {
                setMeanSig(mcache_[idx], storage_[idx]);
              }
              global_.means[hi] = mcache_[idx].mean;
              global_.sigmas[hi] = mcache_[idx].sigma;
            }
            double M = vectorized::dot_product(global_.means.size(), &global_.means[0], &global_.c_scale[0]);
            double C = vectorized::dot_product(global_.sigmas.size(), &global_.sigmas[0], &global_.c_scale[0]);

            unsigned cidx = getIdx(0, global_.p1, v, vi);
            // The step1 cache might not have been allocated yet...
            mcache_[cidx].step1.Resize(storage_[cidx].size());
            mcache_[cidx].step1.Clear();

            for (unsigned hi = 0; hi < hpoints_[0].size(); ++hi) {
              // We can skip all this if the weight is zero
              if (global_.c_sum[hi] == 0.) continue;

              global_.slopes[hi] = global_.sigmas[hi] / C;
              global_.offsets[hi] = global_.means[hi] - (M * global_.slopes[hi]);

#if HFVERBOSE > 0
                std::cout << hi << "\t" << hpoints_[0][hi]
                          << "\t mean = " << global_.means[hi]
                          << "\t sigma = " << global_.sigmas[hi]
                          << "\t slope = " << global_.slopes[hi]
                          << "\t offset = " << global_.offsets[hi] << "\n";
#endif

              unsigned idx = getIdx(0, hi, v, vi);


              // TODO: Scope for optimisation here if we know binning is regular?
              double xl = cache_.GetEdge(0) * global_.slopes[hi] + global_.offsets[hi];
              int il =  cache_.FindBin(xl);
              double xh = 0.;
              int ih = 0;
              int n = storage_[idx].size();

#if HFVERBOSE > 0
              storage_[idx].Dump();
#endif

              for (unsigned ib = 0; ib < storage_[idx].size(); ++ib) {

                double sum = 0.;

                xh = cache_.GetEdge(ib+1) * global_.slopes[hi] + global_.offsets[hi];
                ih =  cache_.FindBin(xh);

#if HFVERBOSE > 1
                std::cout << "Calculating bin " << ib << "\txl = " << xl
                          << "\til = " << il << "\txh = " << xh
                          << "\tih = " << ih << "\n";
#endif

                if (il != -1 && il != n && il != ih) {
                  sum += (cache_.GetEdge(il + 1) - xl) * storage_[idx][il];
#if HFVERBOSE > 1
                  std::cout << "Adding from lower edge: to boundary = "
                            << cache_.GetEdge(il + 1)
                            << "\tcontent = " << storage_[idx][il] << "\n";
#endif
                }

                for (int step = il + 1; step < ih; ++step) {
                  sum += cache_.GetWidth(step) * storage_[idx][step];
#if HFVERBOSE > 1
                  std::cout << "Adding whole bin: bin = " << step
                            << "\tcontent = " << storage_[idx][step] << "\n";
#endif
                }
                // Add the fraction of the last bin
                if (ih != -1 && ih != n && il != ih) {
                  sum += (xh - cache_.GetEdge(ih)) * storage_[idx][ih];
#if HFVERBOSE > 1
                  std::cout
                      << "Adding to upper edge: from boundary = "
                      << cache_.GetEdge(ih)
                      << "\tcontent = " << storage_[idx][ih] << "\n";
#endif
                }

                if (il == ih && il != -1 && il != n) {
                  sum += (xh - xl) * storage_[idx][il];
#if HFVERBOSE > 1
                  std::cout << "Adding partial bin: bin = " << il
                            << "\tcontent = " << storage_[idx][il] << "\n";
#endif
                }

                mcache_[cidx].step1[ib] += (global_.c_sum[hi] * sum / cache_.GetWidth(ib));
                il = ih;
                xl = xh;
              }
            }
#if HFVERBOSE > 0
            mcache_[cidx].step1.Dump();
#endif
          }
          if (htype_ == HorizontalType::Integral) {
            unsigned idx1 = getIdx(0, global_.p1, v, vi);
            unsigned idx2 = getIdx(0, global_.p2, v, vi);
            if (!mcache_[idx1].cdf_set) {
#if HFVERBOSE > 0
              std::cout << "Setting cdf for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
#endif
              setCdf(mcache_[idx1], storage_[idx1]);
            }
            if (!mcache_[idx2].cdf_set) {
#if HFVERBOSE > 0
              std::cout << "Setting cdf for " << 0 << " " << global_.p2 << " " << v << " " << vi << "\n";
#endif
              setCdf(mcache_[idx2], storage_[idx2]);
            }
            if (!mcache_[idx1].interp_set) {
#if HFVERBOSE > 0
              std::cout << "Setting interp for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
#endif
              prepareInterpCache(mcache_[idx1], mcache_[idx2]);
            }
#if HFVERBOSE > 0
            std::cout << "Doing horizontal morph for  " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
#endif
            double x1 = hpoints_[0][global_.p1];
            double x2 = hpoints_[0][global_.p2];
            double y1 = mcache_[idx1].integral;
            double y2 = mcache_[idx2].integral;
            if(y1 <= 0.0 || y2 <= 0.0)
            {
              FastTemplate emptyhist(cache_.size());
              mcache_[idx1].step1 = emptyhist;
            }
            else
            {
              mcache_[idx1].step1 = cdfMorph(idx1, x1, x2, val);
              mcache_[idx1].step1.CropUnderflows();
              double ym = y1 + ((y2 - y1) / (x2 - x1)) * (val - x1);
              mcache_[idx1].step1.Scale(ym / integrateTemplate(mcache_[idx1].step1));
            }
          }
        }

        // This next block is identical to the Integral morphing, so should make it common code
        if (v >= 1) {
#if HFVERBOSE > 0
          std::cout << "Setting sumdiff for vmorph " << v << "\n";
#endif
          unsigned idx = getIdx(0, global_.p1, 0, 0);
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = mcache_[idxLo].step1;
          FastTemplate hi = mcache_[idxHi].step1;
          if (vtype_ == VerticalSetting::QuadLinear) {
            hi.Subtract(mcache_[idx].step1);
            lo.Subtract(mcache_[idx].step1);
          } else if (vtype_ == VerticalSetting::LogQuadLinear) {
            hi.LogRatio(mcache_[idx].step1);
            lo.LogRatio(mcache_[idx].step1);
          }
          // TODO: could skip the next two lines if .sum and .diff have been set before
          mcache_[idxLo].sum = mcache_[idx].step1;
          mcache_[idxLo].diff = mcache_[idx].step1;
          FastTemplate::SumDiff(hi, lo, mcache_[idxLo].sum, mcache_[idxLo].diff);
        }
      }
      hmorph_sentry_.reset();
    }

    if (step1 && global_.single_point) {
#if HFVERBOSE > 0
      std::cout << "Checking step 1 (single point)\n";
#endif
      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        unsigned idx = getIdx(0, global_.p1, 0, 0);
        if (v == 0) {
          mcache_[idx].step1 = storage_[idx];
        }
        if (v >= 1) {
#if HFVERBOSE > 0
          std::cout << "Setting sumdiff for vmorph " << v << "\n";
#endif
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = storage_[idxLo];
          FastTemplate hi = storage_[idxHi];
          if (vtype_ == VerticalSetting::QuadLinear) {
            hi.Subtract(storage_[idx]);
            lo.Subtract(storage_[idx]);
          } else if (vtype_ == VerticalSetting::LogQuadLinear) {
            hi.LogRatio(storage_[idx]);
            lo.LogRatio(storage_[idx]);
          }
          // TODO: could skip the next two lines if .sum and .diff have been set before
          mcache_[idxLo].sum = storage_[idx];
          mcache_[idxLo].diff = storage_[idx];
          FastTemplate::SumDiff(hi, lo, mcache_[idxLo].sum, mcache_[idxLo].diff);
        }
      }
      hmorph_sentry_.reset();
    }

    if (step2) {
#if HFVERBOSE > 0
      std::cout << "Checking step 2\n";
#endif

      unsigned idx = getIdx(0, global_.p1, 0, 0);
      if (vertical_prev_vals_.size() == 0) {
        vertical_prev_vals_.resize(vmorphs_.getSize());
      }
      if (!fast_vertical_) {
        mcache_[idx].step2 = mcache_[idx].step1;
        if (vtype_ == VerticalSetting::LogQuadLinear) {
          mcache_[idx].step2.Log();
        }
      }

#if HFVERBOSE > 0
      std::cout << "Template before vmorph: " << mcache_[idx].step2.Integral() << "\n";
      mcache_[idx].step2.Dump();
#endif

      for (int v = 0; v < vmorphs_.getSize(); ++v) {
        double x = vmorphs_vec_[v]->getVal();
        // if we're in fast_vertical then need to check if this vmorph value has changed.
        // if it hasn't then we skip immediately.
        if (fast_vertical_ && (x == vertical_prev_vals_[v])) continue;

        unsigned vidx = getIdx(0, global_.p1, v+1, 0);

        if (fast_vertical_) {
          double xold = vertical_prev_vals_[v];
          mcache_[idx].step2.DiffMeld(mcache_[vidx].diff, mcache_[vidx].sum, 0.5*x, smoothStepFunc(x), 0.5*xold, smoothStepFunc(xold));
        } else {
          mcache_[idx].step2.Meld(mcache_[vidx].diff, mcache_[vidx].sum, 0.5*x, smoothStepFunc(x));
        }
        vertical_prev_vals_[v] = x;

#if HFVERBOSE > 1
        std::cout << "Morphing for " << vmorphs_[v].GetName() << " with value: " << x << "\n";
        std::cout << "Template after vmorph " << v+1 << ": " << mcache_[idx].step2.Integral() << "\n";
        mcache_[idx].step2.Dump();
#endif
      }
      cache_.CopyValues(mcache_[idx].step2);
      if (vtype_ == VerticalSetting::LogQuadLinear) {
        cache_.Exp();
      }
      cache_.CropUnderflows();
      if (enable_fast_vertical_) fast_vertical_ = true;

#if HFVERBOSE > 0
      std::cout << "Final cache: " << cache_.Integral() << "\n";
      cache_.Dump();
#endif

      vmorph_sentry_.reset();
    }
  }

  if (rebin_ && (step1 || step2)) applyRebin();
}


void CMSHistFunc::applyRebin() const {
  rebin_cache_.Clear();
  for (unsigned i = 0; i < cache_.size(); ++i) {
    rebin_cache_[rebin_scheme_[i]] += (cache_[i] * cache_.GetWidth(i));
  }
  for (unsigned i = 0; i < rebin_cache_.size(); ++i) {
    rebin_cache_[i] /= rebin_cache_.GetWidth(i);
  }
}


Double_t CMSHistFunc::evaluate() const {
  updateCache();
  return cache().GetAt(x_);
}

unsigned CMSHistFunc::getIdx(unsigned hindex, unsigned hpoint, unsigned vindex,
                             unsigned vpoint) const {
  unsigned n_vpoints = vmorphs_.getSize();
  return hindex + hpoint * (1 + n_vpoints * 2) +
         (vindex * 2 - (vindex > 0 ? 1 : 0)) + vpoint;
}

inline double CMSHistFunc::smoothStepFunc(double x) const {
  if (fabs(x) >= vsmooth_par_) return x > 0 ? +1 : -1;
  double xnorm = x / vsmooth_par_;
  double xnorm2 = xnorm * xnorm;
  return 0.125 * xnorm * (xnorm2 * (3. * xnorm2 - 10.) + 15);
}

void CMSHistFunc::setCdf(Cache& c, FastTemplate const& h) const {
  c.cdf = FastTemplate(h.size() + 1);
  c.integral = integrateTemplate(h);
  for (unsigned i = 1; i < c.cdf.size(); ++i) {
    c.cdf[i] = (h[i - 1] * cache_.GetWidth(i-1)) / c.integral + c.cdf[i - 1];
  }
  c.cdf_set = true;
}

void CMSHistFunc::setMeanSig(Cache& c, FastTemplate const& h) const {
  double iw = 0.;
  double ixw = 0.;
  double ixxw = 0.;
  for (unsigned  i = 0; i < h.size(); ++i) {
    double x = (cache_.GetEdge(i + 1) + cache_.GetEdge(i)) / 2.;
    iw += cache_.GetWidth(i) * h[i];
    ixw += cache_.GetWidth(i) * h[i] * x;
    ixxw += cache_.GetWidth(i) * h[i] * x * x;
  }
  c.mean = ixw / iw;
  c.sigma = (ixxw / iw) - c.mean * c.mean;
  c.meansig_set = true;
}

void CMSHistFunc::updateMomentFractions(double m) const {
  int nPdf = hpoints_[0].size();

  double dm = m - hpoints_[0][0];

#if HFVERBOSE > 0
  std::cout << "dm =  " << dm << "\n";
  std::cout << "mtype =  " << mtype_ << "\n";
#endif

  // fully non-linear
  double sumposfrac = 0.;
  for (int i = 0; i < nPdf; ++i) {
    double ffrac = 0.;
    for (int j = 0; j < nPdf; ++j) {
      ffrac += global_.m(j, i) * (j == 0 ? 1. : std::pow(dm, (double)j));
    }
    if (ffrac >= 0) sumposfrac += ffrac;
    // fractions for pdf
    global_.c_sum[i] = ffrac;
    global_.c_scale[i] = ffrac;
  }

  // various mode settings
  int imin = global_.p1;
  int imax = global_.single_point ? global_.p1 : global_.p2;
  double mfrac =
      (m - hpoints_[0][imin]) / (hpoints_[0][imax] - hpoints_[0][imin]);
  switch (mtype_) {
    case NonLinear:
      // default already set above
      break;

    case SineLinear:
      mfrac = TMath::Sin(TMath::PiOver2() * mfrac);
    // now fall through to Linear case

    case Linear:
      for (int i = 0; i < nPdf; ++i) {
        global_.c_sum[i] = 0.;
        global_.c_scale[i] = 0.;
      }
      if (imax > imin) {  // m in between mmin and mmax
        global_.c_sum[imin] = 1. - mfrac;
        global_.c_scale[imin] = 1. - mfrac;
        global_.c_sum[imax] = mfrac;
        global_.c_scale[imax] = mfrac;
      } else if (imax == imin) {  // m outside mmin and mmax
        global_.c_sum[imin] = 1.;
        global_.c_scale[imin] = 1.;
      }
      break;
    case NonLinearLinFractions:
      for (int i = 0; i < nPdf; ++i) {
        global_.c_sum[i] = 0.;
      }
      if (imax > imin) {  // m in between mmin and mmax
        global_.c_sum[imin] = 1. - mfrac;
        global_.c_sum[imax] = mfrac;
      } else if (imax == imin) {  // m outside mmin and mmax
        global_.c_sum[imin] = 1.;
      }
      break;
    case NonLinearPosFractions:
      for (Int_t i = 0; i < nPdf; ++i) {
        if (global_.c_sum[i] < 0) {
          global_.c_sum[i] = 0;
        }
        global_.c_sum[i] = global_.c_sum[i] / sumposfrac;
      }
      break;
  }
#if HFVERBOSE > 0
  for (unsigned i = 0; i < global_.c_sum.size(); ++i) {
    std::cout << hpoints_[0][i] << "\t c_sum = " << global_.c_sum[i] << "\t c_scale = " << global_.c_scale[i] << "\n";
  }
#endif
}

void CMSHistFunc::prepareInterpCache(Cache& c1,
                                     Cache const& c2) const {
  // *
  // *......We are going to step through all the edges of both input
  // *      cdf's ordered by increasing value of y. We start at the
  // *      lower edge, but first we should identify the upper ends of the
  // *      curves. These (ixl1, ixl2) are the first point in each cdf from
  // *      above that has the same integral as the last edge.
  // *

  Int_t ix1l = cache_.size();
  Int_t ix2l = cache_.size();
  while (c1.cdf[ix1l - 1] >= c1.cdf[ix1l]) {
    ix1l = ix1l - 1;
  }
  while (c2.cdf[ix2l - 1] >= c2.cdf[ix2l]) {
    ix2l = ix2l - 1;
  }

  // *
  // *......Step up to the beginnings of the curves. These (ix1, ix2) are the
  // *      first non-zero points from below.

  Int_t ix1 = -1;
  do {
    ix1 = ix1 + 1;
  } while (c1.cdf[ix1 + 1] <= c1.cdf[0]);

  Int_t ix2 = -1;
  do {
    ix2 = ix2 + 1;
  } while (c2.cdf[ix2 + 1] <= c2.cdf[0]);

#if HFVERBOSE > 2
  std::cout << "First and last edge of hist1: " << ix1 << " " << ix1l
            << std::endl;
  std::cout << "Relevant bins of cdf1: (x, y) = " << std::endl;
  for (int ind=ix1; ind<=ix1l; ind++)
  {
    std::cout << "\t(" << cache_.GetEdge(ind) << "," << c1.cdf[ind] << ")" << std::endl;
  }
  std::cout << "First and last edge of hist2: " << ix2 << " " << ix2l
            << std::endl;
  std::cout << "Relevant bins of cdf2: (x, y) = " << std::endl;
  for (int ind=ix2; ind<=ix2l; ind++)
  {
    std::cout << "\t(" << cache_.GetEdge(ind) << "," << c2.cdf[ind] << ")" << std::endl;
  }
#endif


  // ....The first interpolated point should be computed now.

  Int_t nx3 = 0;
  // Double_t x1, x2, x;
  Double_t x1, x2;
  x1 = cache_.GetEdge(ix1);
  x2 = cache_.GetEdge(ix2);

  c1.x1.push_back(x1);
  c1.x2.push_back(x2);
  c1.y.push_back(0.);

  // x = wt1 * x1 + wt2 * x2;
  // xdisn[nx3] = x;
  // sigdisn[nx3] = 0;
  // if (idebug >= 1) {
  //   std::cout << "First interpolated point: " << xdisn[nx3] << " "
  //             << sigdisn[nx3] << std::endl;
  //   std::cout << "                          " << x1 << " <= " << x
  //             << " <= " << x2 << std::endl;
  // }

  // .....Loop over the remaining point in both curves. Getting the last
  //      points may be a bit tricky due to limited floating point
  //      precision.

#if HFVERBOSE > 2
  std::cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l
       << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;
#endif

  Double_t yprev = -1;  // The probability y of the previous point, it will
  // get updated and used in the loop.
  Double_t y = 0;
  while ((ix1 < ix1l) | (ix2 < ix2l)) {

#if HFVERBOSE > 2
    std::cout << "----Top of while with ix1=" << ix1 << ", ix1l=" << ix1l
         << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;
#endif
    // .....Increment to the next lowest point. Step up to the next
    //      kink in case there are several empty (flat in the integral)
    //      bins.

    Int_t i12type = -1;  // Tells which input distribution we need to
                         // see next point of.
    if ((c1.cdf[ix1 + 1] <= c2.cdf[ix2 + 1] || ix2 == ix2l) && ix1 < ix1l) {
      ix1 = ix1 + 1;
      while (c1.cdf[ix1 + 1] < c1.cdf[ix1] && ix1 < ix1l) {
        ix1 = ix1 + 1;
      }
      i12type = 1;
    } else if (ix2 < ix2l) {
      ix2 = ix2 + 1;
      while (c2.cdf[ix2 + 1] <= c2.cdf[ix2] && ix2 < ix2l) {
        ix2 = ix2 + 1;
      }
      i12type = 2;
    }
    if (i12type == 1) {
#if HFVERBOSE > 2
      std::cout << "Pair for i12type=1: " << c2.cdf[ix2] << " "
                << c1.cdf[ix1] << " " << c2.cdf[ix2 + 1] << std::endl;
#endif
      x1 = cache_.GetEdge(ix1);
      y = c1.cdf[ix1];
      Double_t x20 = cache_.GetEdge(ix2);
      Double_t y20 = c2.cdf[ix2];
      Double_t x21 = x20;
      Double_t y21 = y20;
      if (ix2 < ix2l) {
        x21 = cache_.GetEdge(ix2 + 1);
        y21 = c2.cdf[ix2 + 1];
      }
      // .....Calculate where the cummulative probability y in distribution 1
      //      intersects between the 2 points from distribution 2 which
      //      bracket it.

      if (y21 > y20) {
        x2 = x20 + (x21 - x20) * (y - y20) / (y21 - y20);
      } else {
        x2 = x20;
      }
    } else {
#if HFVERBOSE > 2
      std::cout << "Pair for i12type=2: " << c1.cdf[ix1] << " "
                << c2.cdf[ix2] << " " << c1.cdf[ix1 + 1] << std::endl;
#endif
      x2 = cache_.GetEdge(ix2);
      y = c2.cdf[ix2];
      Double_t x10 = cache_.GetEdge(ix1);
      Double_t y10 = c1.cdf[ix1];
      Double_t x11 = x10;
      Double_t y11 = y10;
      if (ix1 < ix1l) {
        x11 = cache_.GetEdge(ix1 + 1);
        y11 = c1.cdf[ix1 + 1];
      }

      // .....Calculate where the cummulative probability y in distribution 2
      //      intersects between the 2 points from distribution 1 which
      //      brackets it.

      if (y11 > y10) {
        x1 = x10 + (x11 - x10) * (y - y10) / (y11 - y10);
      } else {
        x1 = x10;
      }
    }

    // .....Interpolate between the x's in the 2 distributions at the
    //      cummulative probability y. Store the (x,y) for provisional
    //      edge nx3 in (xdisn[nx3],sigdisn[nx3]). nx3 grows for each point
    //      we add the the arrays. Note: Should probably turn the pair into
    //      a structure to make the code more object-oriented and readable.

    // x = wt1 * x1 + wt2 * x2;
    if (y > yprev) {
      nx3 = nx3 + 1;
#if HFVERBOSE > 2
        std::cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" << nx3
             << ", y=" << y << ", yprev=" << yprev << std::endl;
             // << ", x= " << x << ", y=" << y << ", yprev=" << yprev << std::endl;
#endif
      yprev = y;
      c1.x1.push_back(x1);
      c1.x2.push_back(x2);
      c1.y.push_back(y);
      // xdisn[nx3] = x;
      // sigdisn[nx3] = y;
#if HFVERBOSE > 2
      std::cout << "    ix1=" << ix1 << ", ix2= " << ix2
           << ", i12type= " << i12type << ", sigdis1[ix1]=" << c1.cdf[ix1]
           << std::endl;
      std::cout << "        "
           << ", nx3=" << nx3 << ", y= " << c1.y[nx3]
           // << ", nx3=" << nx3 << ", x=" << x << ", y= " << c1.y[nx3]
           << std::endl;
#endif
    }
  }
#if HFVERBOSE > 2
    for (unsigned i = 0; i < c1.y.size(); i++) {
      std::cout << " nx " << i << " " << c1.x1[i] << " " << c1.x2[i] << " " << c1.y[i]
                << std::endl;
    }
#endif
  c1.interp_set = true;
}

FastTemplate CMSHistFunc::cdfMorph(unsigned idx, double par1, double par2,
                                 double parinterp) const {
  double wt1;
  double wt2;
  if (par2 != par1) {
    wt1 = 1. - (parinterp - par1) / (par2 - par1);
    wt2 = 1. + (parinterp - par2) / (par2 - par1);
  } else {
    wt1 = 0.5;
    wt2 = 0.5;
  }

  // ......Give a warning if this is an extrapolation.

  if (wt1 < 0 || wt1 > 1. || wt2 < 0. || wt2 > 1. ||
      fabs(1 - (wt1 + wt2)) > 1.0e-4) {
    std::cout << "Warning! th1fmorph: This is an extrapolation!! Weights are "
              << wt1 << " and " << wt2 << " (sum=" << wt1 + wt2 << ")"
              << std::endl;
  }
#if HFVERBOSE > 2
    std::cout << "th1morph - Weights: " << wt1 << " " << wt2 << std::endl;
#endif

  // Treatment for empty histograms: Return an empty histogram
  // with interpolated bins.

  // FIXME should restore zero integral behaviour
  // if (hist1.Integral() <= 0 || hist2.Integral() <= 0) {
  //   std::cout << "Warning! th1morph detects an empty input histogram. Empty "
  //           "interpolated histogram returned: " << std::endl;
  //   return (TH1F("morphed", "morphed", mc_.nbn, mc_.xminn, mc_.xmaxn));
  // }
  // if (idebug >= 1)
  //   std::cout << "Input histogram content sums: " << hist1.Integral() << " "
  //        << hist2.Integral() << std::endl;

  // *
  // *......Extract the single precision histograms into double precision arrays
  // *      for the interpolation computation. The offset is because sigdis(i)
  // *      describes edge i (there are nbins+1 of them) while dist1/2
  // *      describe bin i. Be careful, ROOT does not use C++ convention to
  // *      number bins: dist1[ibin] is content of bin ibin where ibin runs from
  // *      1 to nbins. We allocate some extra space for the derived
  // distributions
  // *      because there may be as many as nb1+nb2+2 edges in the intermediate
  // *      interpolated cdf described by xdisn[i] (position of edge i) and
  // *      sigdisn[i] (cummulative probability up this edge) before we project
  // *      into the final binning.
  // Float_t const* dist1 = hist1.GetArray();
  // Float_t const* dist2 = hist2.GetArray();

  // *......Now we loop over the edges of the bins of the interpolated
  // *      histogram and find out where the interpolated cdf 3
  // *      crosses them. This projection defines the result and will
  // *      be stored (after differention and renormalization) in the
  // *      output histogram.
  // *
  // *......We set all the bins following the final edge to the value
  // *      of the final edge.

  Cache const& c1 = mcache_[idx];

  int nbn = cache_.size();
  int nbe = cache_.size() + 1;

  double x  = cache_.GetEdge(nbn);
  int ix = nbn;

  int nx3 = c1.y.size();
  std::vector<double> xdisn(nx3);
  for (int i = 0; i < nx3; ++i) {
    xdisn[i] = wt1 * c1.x1[i] + wt2 * c1.x2[i];
  }
#if HFVERBOSE > 2
    std::cout << "relevant x,y for the morphed cdf: "  << std::endl;
    for (unsigned int i=0; i < c1.y.size(); ++i){ std::cout << "\t (x,y) = (" << xdisn[i] << "," << c1.y[i] << ")" << std::endl; }
#endif
  std::vector<double> sigdisf(nbe, 0.);


  nx3 = nx3 - 1;

#if HFVERBOSE > 2
    std::cout << "------> Any final bins to set? " << x << " " << xdisn[nx3]
              << std::endl;
#endif

  while (x >= xdisn[nx3]) {
    sigdisf[ix] = c1.y[nx3];

#if HFVERBOSE > 2
      std::cout << "   Setting final bins " << ix << " " << x << " "
                << sigdisf[ix] << std::endl;
#endif

    ix = ix - 1;
    x = cache_.GetEdge(ix);
  }
  Int_t ixl = ix + 1;

#if HFVERBOSE > 2
  std::cout << " Now ixl=" << ixl << " ix=" << ix << std::endl;
#endif

  // *
  // *......The beginning may be empty, so we have to step up to the first
  // *      edge where the result is nonzero. We zero the bins which have
  // *      and upper (!) edge which is below the first point of the
  // *      cummulative distribution we are going to project to this
  // *      output histogram binning.
  // *

  ix = 0;
  x = cache_.GetEdge(ix + 1);

#if HFVERBOSE > 2
  std::cout << "Start setting initial bins at x=" << x << std::endl;
#endif

  while (x <= xdisn[0]) {
    sigdisf[ix] = c1.y[0];
#if HFVERBOSE > 2
    std::cout << "   Setting initial bins " << ix << " " << x << " "
              << xdisn[1] << " " << sigdisf[ix] << std::endl;
#endif
    ix = ix + 1;
    x = cache_.GetEdge(ix + 1);
  }
  Int_t ixf = ix;

#if HFVERBOSE > 2
  std::cout << "Bins left to loop over:" << ixf << "-" << ixl << std::endl;
#endif
  // *......Also the end (from y to 1.0) often comes before the last edge
  // *      so we have to set the following to 1.0 as well.

  Int_t ix3 = 0;  // Problems with initial edge!!!
  double y = 0;
  for (ix = ixf; ix < ixl; ix++) {
    x = cache_.GetEdge(ix);
    if (x < xdisn[0]) {
      y = 0;
    } else if (x > xdisn[nx3]) {
      y = 1.;
    } else {
      while (xdisn[ix3 + 1] <= x && ix3 < 2 * nbn) {
        ix3 = ix3 + 1;
      }
      int bin = cache_.FindBin(x);
      if (bin == -1) bin = 0;
      // std::cout << "x=" << x << ", bin=" << bin << "\n";
      Double_t dx2 = cache_.GetWidth(bin);
      if (xdisn[ix3 + 1] - x >= 1.0 * dx2) {  // Empty bin treatment
#if HFVERBOSE > 2
        std::cout << "Warning - th1fmorph: encountered empty bin." << std::endl;
#endif
        y = c1.y[ix3];
      } else if (xdisn[ix3 + 1] > xdisn[ix3]) {  // Normal bins
        y = c1.y[ix3] + (c1.y[ix3 + 1] - c1.y[ix3]) *
                               (x - xdisn[ix3]) / (xdisn[ix3 + 1] - xdisn[ix3]);
      } else {  // Is this ever used?
        y = 0;
        std::cout << "Warning - th1fmorph: This probably shoudn't happen! "
                  << std::endl;
        std::cout << "Warning - th1fmorph: Zero slope solving x(y)"
                  << std::endl;
        // for (unsigned z = 0; z < c1.y.size(); ++z) {
        //   printf("x %.6f x1 %.6f x2 %.6f y %.6f\n", xdisn[z], c1.x1[z], c1.x2[z], c1.y[z]);
        // }
      }
    }
    sigdisf[ix] = y;
#if HFVERBOSE > 2
      std::cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3]
                << ", x=" << x << ", next xdisn=" << xdisn[ix3 + 1]
                << std::endl;
      std::cout << "   cdf n=" << c1.y[ix3] << ", y=" << y
           << ", next point=" << c1.y[ix3 + 1] << std::endl;
#endif
  }
#if HFVERBOSE > 2
  std::cout << "CDF mapped into the histogram binning:" << std::endl;
  for (int ind = std::max(0, ixf - 1); ind < ixl+1; ind++) {
    std::cout << "\t(bin,x,y) = (" << ind << "," << cache_.GetEdge(ind) << "," << sigdisf[ind] << ")" << std::endl;
  }
#endif

  // .....Differentiate interpolated cdf and return renormalized result in
  //      new histogram.

  FastTemplate morphedhist(nbn);

  for (ix = nbn - 1; ix > -1; ix--) {
    y = sigdisf[ix + 1] - sigdisf[ix];
    morphedhist[ix] = y / cache_.GetWidth(ix);
  }

  // ......All done, return the result.

  return morphedhist;
}

double CMSHistFunc::integrateTemplate(FastTemplate const& t) const {
  DefaultAccumulator<double> total = 0;
  for (unsigned int i = 0; i < t.size(); ++i) total += t[i] * cache_.GetWidth(i);
  return total.sum();
}


void CMSHistFunc::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  cache().Dump();
  std::cout << ">> Errors:\n";
  binerrors_.Dump();
  std::cout << ">> Horizontal morphing sentry: " << hmorph_sentry_.good() << "\n";
  hmorph_sentry_.Print("v");
  std::cout << ">> Vertical morphing sentry: " << vmorph_sentry_.good() << "\n";
  vmorph_sentry_.Print("v");
}

Int_t CMSHistFunc::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistFunc::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  switch (code) {
    case 1: {
      updateCache();
      return cache().IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

CMSHistFuncWrapper const* CMSHistFunc::wrapper() const {
  RooFIter iter = this->valueClientMIterator();
  RooAbsArg *arg = nullptr;
  while((arg = iter.next())) {
    CMSHistFuncWrapper const* wrapper = dynamic_cast<CMSHistFuncWrapper const*>(arg);
    if (wrapper) return wrapper;
  }
  return nullptr;
}

RooAbsReal const& CMSHistFunc::getXVar() const {
  return x_.arg();
}

void CMSHistFunc::EnableFastVertical() {
  enable_fast_vertical_ = true;
}

#undef HFVERBOSE
