#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
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

CMSHistFunc::CMSHistFunc() {
  morph_strategy_ = 0;
  veval = 0;
  initialized_ = false;
  htype_ = HorizontalType::Closest;
  mtype_ = MomentSetting::NonLinearPosFractions;
}

CMSHistFunc::CMSHistFunc(const char* name, const char* title, RooRealVar& x,
                         TH1 const& hist)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      vmorphs_("vmorphs", "", this),
      hmorphs_("hmorphs", "", this),
      vmorph_sentry_(TString(name) + "_vmorph_sentry", ""),
      hmorph_sentry_(TString(name) + "_hmorph_sentry", ""),
      cache_(hist),
      binerrors_(FastTemplate(hist.GetNbinsX())),
      morph_strategy_(0),
      veval(0),
      initialized_(false),
      htype_(HorizontalType::Moment),
      mtype_(MomentSetting::Linear) {
  prepareStorage();  // Prepare storage of size 1 -> the cache_ will be copied in there
  for (unsigned i = 0; i < cache_.size(); ++i) {
    binerrors_[i] = hist.GetBinError(i + 1);
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
      binerrors_(other.binerrors_),
      storage_(other.storage_),
      morph_strategy_(other.morph_strategy_),
      veval(other.veval),
      initialized_(false),
      htype_(other.htype_),
      mtype_(other.mtype_) {
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

  unsigned n_hpoints = 1;
  if (hmorphs_.getSize() == 1) {
    n_hpoints = hpoints_[0].size();
    global_.bedgesn.resize(cache_.size() + 1);
    for (unsigned i = 0; i < global_.bedgesn.size(); ++i) {
     global_.bedgesn[i] = cache_.GetEdge(i);
    }

    // This portion of code adapted from:
    // RooMomentMorph.cxx
    TVectorD dm(n_hpoints);
    global_.m.ResizeTo(n_hpoints, n_hpoints);
    // global_.m = TMatrixD(n_hpoints, n_hpoints);
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
  }
  initialized_ = true;
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
  FNLOGC(std::cout, veval) << "Storage size set to: " << storage_.size() << "\n";
  storage_[getIdx(0, 0, 0, 0)] = cache_;
}

void CMSHistFunc::setShape(unsigned hindex, unsigned hpoint, unsigned vindex,
                           unsigned vpoint, TH1 const& hist) {
  unsigned idx = getIdx(hindex, hpoint, vindex, vpoint);
  FNLOGC(std::cout, veval) << "hindex: " << hindex << " hpoint: " << hpoint
                           << " vindex: " << vindex << " vpoint: " << vpoint
                           << " mapped to idx: " << idx << "\n";
  storage_.at(idx) = FastHisto(hist);
}

void CMSHistFunc::updateCache() const {
  initialize();

  // Quick escape if cache is up-to-date
  if (hmorph_sentry_.good() && vmorph_sentry_.good()) return;

  if (morph_strategy_ == 0) {
    // FNLOGC(std::cout, veval) << "Morphing strategy 0\n";
    // FNLOGC(std::cout, veval) << "Number of horizontal morphs: " << hmorphs_.getSize() << "\n";
    // FNLOGC(std::cout, veval) << "Horizontal morph sentry: " << hmorph_sentry_.good() << "\n";
    // FNLOGC(std::cout, veval) << "Vertical morph sentry: " << vmorph_sentry_.good() << "\n";
    // FNLOGC(std::cout, veval) << "single_point,p1,p2: " << global_.single_point << " " << global_.p1 << " " << global_.p2 << "\n";
    // Figure out what we're doing:
    //  - singlePoint p1 OR
    //  - h-morphing between p1 and p2
    bool step1 = !hmorph_sentry_.good();
    bool step2 = !vmorph_sentry_.good();

    if (step1 && hmorphs_.getSize() == 1) {
      double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();
      FNLOGC(std::cout, veval) << "Updating horizontal points...\n";
      FNLOGC(std::cout, veval) << "Horizontal morph parameter is: " << val << "\n";
      auto upper = std::lower_bound(hpoints_[0].begin(), hpoints_[0].end(), val);
      if (upper == hpoints_[0].begin()) {
        global_.p1 = 0;
        global_.single_point = true;
      } else if (upper == hpoints_[0].end()) {
        global_.p1 = hpoints_[0].size() - 1;
        global_.single_point = true;
      } else {
        auto lower = upper;
        --lower;
        global_.p1 = lower - hpoints_[0].begin();
        global_.p2 = upper - hpoints_[0].begin();
        global_.single_point = false;
        // mh_lo_ = lower->first;
        // mh_hi_ = upper->first;
        // if (!can_morph_) {
        //   single_point_ = true;
        //   if (fabs(upper->first - mh_) <= fabs(mh_ - lower->first)) {
        //     p1_ = upper->second;
        //   } else {
        //     p1_ = lower->second;
        //   }
        // } else {
        //   single_point_ = false;
        // }
      }
      FNLOGC(std::cout, veval) << "single_point,p1,p2: " << global_.single_point << " " << global_.p1 << " " << global_.p2 << "\n";
    }

    if (step1 || step2) {
      if (mcache_.size() == 0) mcache_.resize(storage_.size());
    }

    if (step1 && htype_ == HorizontalType::Moment && !global_.single_point) {
      FNLOGC(std::cout, veval) << "Checking step 1\n";
      double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();
      updateMomentFractions(val);  // updates fractions in global cache, only depends on hmorph par

      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        for (int vi = 0; vi < (v == 0 ? 1 : 2); ++vi) {
          // FNLOGC(std::cout, veval) << "Doing vpoint,vindex = " << v << "\t" << vi << "\n";
          std::vector<double> means(hpoints_[0].size(), 0.);
          std::vector<double> sigmas(hpoints_[0].size(), 0.);
          std::vector<double> slopes(hpoints_[0].size(), 0.);
          std::vector<double> offsets(hpoints_[0].size(), 0.);
          for (unsigned hi = 0; hi < hpoints_[0].size(); ++hi) {
            // define vec of mean vals
            unsigned idx = getIdx(0, hi, v, vi);
            if (!mcache_[idx].meansig_set) {
              setMeanSig(mcache_[idx], storage_[idx]);
            }
            means[hi] = mcache_[idx].mean;
            sigmas[hi] = mcache_[idx].sigma;
          }
          double M = vectorized::dot_product(means.size(), &means[0], &global_.c_scale[0]);
          double C = vectorized::dot_product(sigmas.size(), &sigmas[0], &global_.c_scale[0]);

          unsigned cidx = getIdx(0, global_.p1, v, vi);
          // The step1 cache might not have been allocated yet...
          mcache_[cidx].step1.Resize(storage_[cidx].size());
          mcache_[cidx].step1.Clear();

          for (unsigned hi = 0; hi < hpoints_[0].size(); ++hi) {
            // We can skip all this if the weight is zero
            if (global_.c_sum[hi] == 0.) continue;

            slopes[hi] = sigmas[hi] / C;
            offsets[hi] = means[hi] - (M * slopes[hi]);
            // if (veval) {
            //   std::cout << hi << "\t" << hpoints_[0][hi]
            //             << "\t mean = " << means[hi]
            //             << "\t sigma = " << sigmas[hi]
            //             << "\t slope = " << slopes[hi]
            //             << "\t offset = " << offsets[hi] << "\n";
            // }
            unsigned idx = getIdx(0, hi, v, vi);


            // TODO: Scope for optimisation here if we know binning is regular?
            double xl = storage_[idx].GetEdge(0) * slopes[hi] + offsets[hi];
            int il =  storage_[idx].FindBin(xl);
            double xh = 0.;
            int ih = 0;
            int n = storage_[idx].size();

            if (veval) storage_[idx].Dump();

            for (unsigned ib = 0; ib < storage_[idx].size(); ++ib) {

              double sum = 0.;

              xh = storage_[idx].GetEdge(ib+1) * slopes[hi] + offsets[hi];
              ih =  storage_[idx].FindBin(xh);
              // FNLOGC(std::cout, veval) << "Calculating bin " << ib
              //                          << "\txl = " << xl << "\til = " << il
              //                          << "\txh = " << xh << "\tih = " << ih
              //                          << "\n";

              if (il != -1 && il != n && il != ih) {
                sum += (storage_[idx].GetEdge(il + 1) - xl) * storage_[idx][il];
                // FNLOGC(std::cout, veval)
                //     << "Adding from lower edge: to boundary = "
                //     << storage_[idx].GetEdge(il + 1)
                //     << "\tcontent = " << storage_[idx][il] << "\n";
              }

              for (int step = il + 1; step < ih; ++step) {
                sum += storage_[idx].GetWidth(step) * storage_[idx][step];
                // FNLOGC(std::cout, veval)
                //     << "Adding whole bin: bin = "
                //     << step
                //     << "\tcontent = " << storage_[idx][step] << "\n";
              }
              // Add the fraction of the last bin
              if (ih != -1 && ih != n && il != ih) {
                sum += (xh - storage_[idx].GetEdge(ih)) * storage_[idx][ih];
                // FNLOGC(std::cout, veval)
                //     << "Adding to upper edge: from boundary = "
                //     << storage_[idx].GetEdge(ih)
                //     << "\tcontent = " << storage_[idx][ih] << "\n";
              }

              if (il == ih && il != -1 && il != n) {
                sum += (xh - xl) * storage_[idx][il];
                // FNLOGC(std::cout, veval)
                //     << "Adding partial bin: bin = "
                //     << il
                //     << "\tcontent = " << storage_[idx][il] << "\n";
              }

              mcache_[cidx].step1[ib] += (global_.c_sum[hi] * sum / storage_[idx].GetWidth(ib));
              il = ih;
              xl = xh;
            }
          }

          if (veval) mcache_[cidx].step1.Dump();
        }

        // This next block is identical to the Integral morphing, so should make it common code
        if (v >= 1) {
          FNLOGC(std::cout, veval) << "Setting sumdiff for vmorph " << v << "\n";
          unsigned idx = getIdx(0, global_.p1, 0, 0);
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = mcache_[idxLo].step1;
          FastTemplate hi = mcache_[idxHi].step1;
          hi.Subtract(mcache_[idx].step1);
          lo.Subtract(mcache_[idx].step1);
          // TODO: could skip the next two lines if .sum and .diff have been set before
          mcache_[idxLo].sum = mcache_[idx].step1;
          mcache_[idxLo].diff = mcache_[idx].step1;
          FastTemplate::SumDiff(hi, lo, mcache_[idxLo].sum, mcache_[idxLo].diff);
        }
      }
      hmorph_sentry_.reset();
    }

    if (step1 && htype_ == HorizontalType::Integral && !global_.single_point) {
      FNLOGC(std::cout, veval) << "Checking step 1\n";
      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        for (int vi = 0; vi < (v == 0 ? 1 : 2); ++vi) {
          unsigned idx1 = getIdx(0, global_.p1, v, vi);
          unsigned idx2 = getIdx(0, global_.p2, v, vi);
          if (!mcache_[idx1].cdf_set) {
            FNLOGC(std::cout, veval) << "Setting cdf for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
            setCdf(mcache_[idx1], storage_[idx1]);
          }
          if (!mcache_[idx2].cdf_set) {
            FNLOGC(std::cout, veval) << "Setting cdf for " << 0 << " " << global_.p2 << " " << v << " " << vi << "\n";
            setCdf(mcache_[idx2], storage_[idx2]);
          }
          if (!mcache_[idx1].interp_set) {
            FNLOGC(std::cout, veval) << "Setting interp for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
            prepareInterpCache(mcache_[idx1], mcache_[idx2]);
          }
          FNLOGC(std::cout, veval) << "Doing horizontal morph for  " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
          double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();
          double x1 = hpoints_[0][global_.p1];
          double x2 = hpoints_[0][global_.p2];
          double y1 = mcache_[idx1].integral;
          double y2 = mcache_[idx2].integral;
          mcache_[idx1].step1 = cdfMorph(idx1, x1, x2, val);
          mcache_[idx1].step1.CropUnderflows();
          double ym = y1 + ((y2 - y1) / (x2 - x1)) * (val - x1);
          mcache_[idx1].step1.Scale(ym/mcache_[idx1].step1.Integral());
          if (veval >= 2) {
            std::cout << "Template left: " << storage_[idx1].Integral() << "\n";
            storage_[idx1].Dump();
            std::cout << "Template right: " << storage_[idx2].Integral() << "\n";
            storage_[idx2].Dump();
            std::cout << "Template morphed: " << mcache_[idx1].step1.Integral() << "\n";
            mcache_[idx1].step1.Dump();
          }
        }
        if (v >= 1) {
          FNLOGC(std::cout, veval) << "Setting sumdiff for vmorph " << v << "\n";
          unsigned idx = getIdx(0, global_.p1, 0, 0);
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = mcache_[idxLo].step1;
          FastTemplate hi = mcache_[idxHi].step1;
          hi.Subtract(mcache_[idx].step1);
          lo.Subtract(mcache_[idx].step1);
          // TODO: could skip the next two lines if .sum and .diff have been set before
          mcache_[idxLo].sum = mcache_[idx].step1;
          mcache_[idxLo].diff = mcache_[idx].step1;
          FastTemplate::SumDiff(hi, lo, mcache_[idxLo].sum, mcache_[idxLo].diff);
        }
      }
      hmorph_sentry_.reset();
    }

    if (step1 && global_.single_point) {
      FNLOGC(std::cout, veval) << "Checking step 1 (single point)\n";
      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        unsigned idx = getIdx(0, global_.p1, 0, 0);
        if (v == 0) {
          mcache_[idx].step1 = storage_[idx];
        }
        if (v >= 1) {
          FNLOGC(std::cout, veval) << "Setting sumdiff for vmorph " << v << "\n";
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = storage_[idxLo];
          FastTemplate hi = storage_[idxHi];
          hi.Subtract(storage_[idx]);
          lo.Subtract(storage_[idx]);
          // TODO: could skip the next two lines if .sum and .diff have been set before
          mcache_[idxLo].sum = storage_[idx];
          mcache_[idxLo].diff = storage_[idx];
          FastTemplate::SumDiff(hi, lo, mcache_[idxLo].sum, mcache_[idxLo].diff);
          if (veval >= 2) {
            std::cout << "Template nominal: " << storage_[idx].Integral() << "\n";
            storage_[idx].Dump();
          }
        }
      }
      hmorph_sentry_.reset();
    }

    if (step2) {
      // FNLOGC(std::cout, veval) << "Checking step 2\n";
      unsigned idx = getIdx(0, global_.p1, 0, 0);
      mcache_[idx].step2 = mcache_[idx].step1;
      // if (veval >= 2) {
      //   std::cout << "Template before vmorph: " << mcache_[idx].step2.Integral() << "\n";
      //   mcache_[idx].step2.Dump();
      // }
      for (int v = 0; v < vmorphs_.getSize(); ++v) {
        unsigned vidx = getIdx(0, global_.p1, v+1, 0);

        double x = ((RooRealVar&)vmorphs_[v]).getVal();
        //     // std::cout << "Morphing for " << vmorphs_[i].GetName() << " with value: " << x << "\n";
        // double a = 0.5*x;
        // double b = smoothStepFunc(x);
        mcache_[idx].step2.Meld(mcache_[vidx].diff, mcache_[vidx].sum, 0.5*x, smoothStepFunc(x));
        // if (veval >= 2) {
        //   std::cout << "Template after vmorph " << v+1 << ": " << mcache_[idx].step2.Integral() << "\n";
        //   mcache_[idx].step2.Dump();
        // }
      }
      mcache_[idx].step2.CropUnderflows();
      cache_.CopyValues(mcache_[idx].step2);
      // if (veval >= 1) {
      //   std::cout << "Final cache: " << cache_.Integral() << "\n";
      //   cache_.Dump();
      // }
      vmorph_sentry_.reset();
    }
  }
}

Double_t CMSHistFunc::evaluate() const {
  // LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  updateCache();

  return cache_.GetAt(x_);
}

unsigned CMSHistFunc::getIdx(unsigned hindex, unsigned hpoint, unsigned vindex,
                             unsigned vpoint) const {
  unsigned n_vpoints = vmorphs_.getSize();
  return hindex + hpoint * (1 + n_vpoints * 2) +
         (vindex * 2 - (vindex > 0 ? 1 : 0)) + vpoint;
}

inline double CMSHistFunc::smoothStepFunc(double x) const {
  double _smoothRegion = 1.0;
  if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
  double xnorm = x / _smoothRegion, xnorm2 = xnorm * xnorm;
  return 0.125 * xnorm * (xnorm2 * (3. * xnorm2 - 10.) + 15);
}

void CMSHistFunc::setCdf(Cache& c, FastHisto const& h) const {
  c.cdf = FastTemplate(h.size() + 1);
  c.integral = h.Integral();
  for (unsigned i = 1; i < c.cdf.size(); ++i) {
    c.cdf[i] = h[i - 1] / c.integral + c.cdf[i - 1];
  }
  c.cdf_set = true;
}

void CMSHistFunc::setMeanSig(Cache& c, FastHisto const& h) const {
  double iw = 0.;
  double ixw = 0.;
  double ixxw = 0.;
  for (unsigned  i = 0; i < h.size(); ++i) {
    double x = (h.GetEdge(i + 1) + h.GetEdge(i)) / 2.;
    iw += h.GetWidth(i) * h[i];
    ixw += h.GetWidth(i) * h[i] * x;
    ixxw += h.GetWidth(i) * h[i] * x * x;
  }
  c.mean = ixw / iw;
  c.sigma = (ixxw / iw) - c.mean * c.mean;
  c.meansig_set = true;
}

void CMSHistFunc::updateMomentFractions(double m) const {
  int nPdf = hpoints_[0].size();

  double dm = m - hpoints_[0][0];

  FNLOGC(std::cout, veval) << "dm =  " << dm << "\n";
  FNLOGC(std::cout, veval) << "mtype =  " << mtype_ << "\n";

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
    // if (verbose) { cout << ffrac << endl; }
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
  if (veval) {
    for (unsigned i = 0; i < global_.c_sum.size(); ++i) {
      std::cout << hpoints_[0][i] << "\t c_sum = " << global_.c_sum[i] << "\t c_scale = " << global_.c_scale[i] << "\n";
    }
  }
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
  int idebug = 0;

  Int_t ix1l = global_.bedgesn.size() - 1;
  Int_t ix2l = global_.bedgesn.size() - 1;
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

  if (idebug >= 1) {
    std::cout << "First and last edge of hist1: " << ix1 << " " << ix1l
              << std::endl;
    std::cout << "   " << c1.cdf[ix1] << " " << c1.cdf[ix1 + 1] << std::endl;
    std::cout << "First and last edge of hist2: " << ix2 << " " << ix2l
              << std::endl;
    std::cout << "   " << c2.cdf[ix2] << " " << c2.cdf[ix2 + 1] << std::endl;
  }

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

  if (idebug >= 1) {
    std::cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l
         << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;
  }

  Double_t yprev = -1;  // The probability y of the previous point, it will
  // get updated and used in the loop.
  Double_t y = 0;
  while ((ix1 < ix1l) | (ix2 < ix2l)) {
    if (idebug >= 1)
      std::cout << "----Top of while with ix1=" << ix1 << ", ix1l=" << ix1l
           << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;

    // .....Increment to the next lowest point. Step up to the next
    //      kink in case there are several empty (flat in the integral)
    //      bins.

    Int_t i12type = -1;  // Tells which input distribution we need to
                         // see next point of.
    if ((c1.cdf[ix1 + 1] <= c2.cdf[ix2 + 1] || ix2 == ix2l) && ix1 < ix1l) {
      ix1 = ix1 + 1;
      while (c1.cdf[ix1 + 1] <= c1.cdf[ix1] && ix1 < ix1l) {
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
      if (idebug >= 3) {
        std::cout << "Pair for i12type=1: " << c2.cdf[ix2] << " "
                  << c1.cdf[ix1] << " " << c2.cdf[ix2 + 1] << std::endl;
      }
      x1 = cache_.GetEdge(ix1);
      y = c1.cdf[ix1];
      Double_t x20 = cache_.GetEdge(ix2);
      Double_t x21 = cache_.GetEdge(ix2 + 1);
      Double_t y20 = c2.cdf[ix2];
      Double_t y21 = c2.cdf[ix2 + 1];

      // .....Calculate where the cummulative probability y in distribution 1
      //      intersects between the 2 points from distribution 2 which
      //      bracket it.

      if (y21 > y20) {
        x2 = x20 + (x21 - x20) * (y - y20) / (y21 - y20);
      } else {
        x2 = x20;
      }
    } else {
      if (idebug >= 3) {
        std::cout << "Pair for i12type=2: " << c1.cdf[ix1] << " "
                  << c2.cdf[ix2] << " " << c1.cdf[ix1 + 1] << std::endl;
      }
      x2 = cache_.GetEdge(ix2);
      y = c2.cdf[ix2];
      Double_t x10 = cache_.GetEdge(ix1);
      Double_t x11 = cache_.GetEdge(ix1 + 1);
      Double_t y10 = c1.cdf[ix1];
      Double_t y11 = c1.cdf[ix1 + 1];

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
      if (idebug >= 1) {
        std::cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" << nx3
             << ", y=" << y << ", yprev=" << yprev << std::endl;
             // << ", x= " << x << ", y=" << y << ", yprev=" << yprev << std::endl;
      }
      yprev = y;
      c1.x1.push_back(x1);
      c1.x2.push_back(x2);
      c1.y.push_back(y);
      // xdisn[nx3] = x;
      // sigdisn[nx3] = y;
      if (idebug >= 1) {
        std::cout << "    ix1=" << ix1 << ", ix2= " << ix2
             << ", i12type= " << i12type << ", sigdis1[ix1]=" << c1.cdf[ix1]
             << std::endl;
        std::cout << "        "
             << ", nx3=" << nx3 << ", y= " << c1.y[nx3]
             // << ", nx3=" << nx3 << ", x=" << x << ", y= " << c1.y[nx3]
             << std::endl;
      }
    }
  }
  if (idebug >= 3)
    for (Int_t i = 0; i < nx3; i++) {
      std::cout << " nx " << i << " " << c1.x1[i] << " " << c1.x2[i] << " " << c1.y[i]
                << std::endl;
    }
  c1.interp_set = true;
}

FastTemplate CMSHistFunc::cdfMorph(unsigned idx, double par1, double par2,
                                 double parinterp) const {
  unsigned idebug = 0;

  // Extract bin parameters of input histograms 1 and 2.
  // Supports the cases of non-equidistant as well as equidistant binning
  // and also the case that binning of histograms 1 and 2 is different.
 // TArrayD *bedgesn = new TArrayD(cache_.size() + 1);
 // for (int i = 0; i < bedgesn->GetSize(); ++i) {
 //  (*bedgesn)[i] = cache_.GetEdge(i);
 // }

  // ......The weights (wt1,wt2) are the complements of the "distances" between
  //       the values of the parameters at the histograms and the desired
  //       interpolation point. For example, wt1=0, wt2=1 means that the
  //       interpolated histogram should be identical to input histogram 2.
  //       Check that they make sense. If par1=par2 then we can choose any
  //       valid set of wt1,wt2 so why not take the average?

  Double_t wt1, wt2;
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
  if (idebug >= 1)
    std::cout << "th1morph - Weights: " << wt1 << " " << wt2 << std::endl;


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
  // HMorphNewCDFCache const& c2 = mcache_[globalIdx2];
  std::vector<double> const& bedgesn = global_.bedgesn;

  int nbn = bedgesn.size() - 1;
  int nbe = bedgesn.size();

  double x  = bedgesn.back();
  int ix = nbn;

  int nx3 = c1.y.size();
  std::vector<double> xdisn(nx3);
  for (int i = 0; i < nx3; ++i) {
    xdisn[i] = wt1 * c1.x1[i] + wt2 * c1.x2[i];
  }
  std::vector<double> sigdisf(nbe, 0.);


  nx3 = nx3 - 1;
  if (idebug >= 1)
    std::cout << "------> Any final bins to set? " << x << " " << xdisn[nx3]
              << std::endl;
  while (x >= xdisn[nx3]) {
    sigdisf[ix] = c1.y[nx3];
    if (idebug >= 2)
      std::cout << "   Setting final bins" << ix << " " << x << " "
                << sigdisf[ix] << std::endl;
    ix = ix - 1;
    x = bedgesn[ix];
  }
  Int_t ixl = ix + 1;
  if (idebug >= 1) std::cout << " Now ixl=" << ixl << " ix=" << ix << std::endl;

  // *
  // *......The beginning may be empty, so we have to step up to the first
  // *      edge where the result is nonzero. We zero the bins which have
  // *      and upper (!) edge which is below the first point of the
  // *      cummulative distribution we are going to project to this
  // *      output histogram binning.
  // *

  ix = 0;
  x = bedgesn[ix + 1];
  if (idebug >= 1)
    std::cout << "Start setting initial bins at x=" << x << std::endl;
  while (x <= xdisn[0]) {
    sigdisf[ix] = c1.y[0];
    if (idebug >= 1)
      std::cout << "   Setting initial bins " << ix << " " << x << " "
                << xdisn[1] << " " << sigdisf[ix] << std::endl;
    ix = ix + 1;
    x = bedgesn[ix + 1];
  }
  Int_t ixf = ix;

  if (idebug >= 1)
    std::cout << "Bins left to loop over:" << ixf << "-" << ixl << std::endl;

  // *......Also the end (from y to 1.0) often comes before the last edge
  // *      so we have to set the following to 1.0 as well.

  Int_t ix3 = 0;  // Problems with initial edge!!!
  double y = 0;
  for (ix = ixf; ix < ixl; ix++) {
    x = bedgesn[ix];
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
      if (xdisn[ix3 + 1] - x > 1.1 * dx2) {  // Empty bin treatment
        y = c1.y[ix3 + 1];
      } else if (xdisn[ix3 + 1] > xdisn[ix3]) {  // Normal bins
        y = c1.y[ix3] + (c1.y[ix3 + 1] - c1.y[ix3]) *
                               (x - xdisn[ix3]) / (xdisn[ix3 + 1] - xdisn[ix3]);
      } else {  // Is this ever used?
        y = 0;
        std::cout << "Warning - th1fmorph: This probably shoudn't happen! "
                  << std::endl;
        std::cout << "Warning - th1fmorph: Zero slope solving x(y)"
                  << std::endl;
      }
    }
    sigdisf[ix] = y;
    if (idebug >= 3) {
      std::cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3]
                << ", x=" << x << ", next xdisn=" << xdisn[ix3 + 1]
                << std::endl;
      std::cout << "   cdf n=" << c1.y[ix3] << ", y=" << y
           << ", next point=" << c1.y[ix3 + 1] << std::endl;
    }
  }

  // .....Differentiate interpolated cdf and return renormalized result in
  //      new histogram.

  // TH1F morphedhist("morphed", "morphed", mc_.nbn, 0, static_cast<float>(mc_.nbn));

  FastTemplate morphedhist(nbn);

  for (ix = nbn - 1; ix > -1; ix--) {
    y = sigdisf[ix + 1] - sigdisf[ix];
    morphedhist[ix] = y;
  }

  // ......All done, return the result.

  return morphedhist;
}

void CMSHistFunc::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  cache_.Dump();
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
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache();
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

CMSHistFuncV::CMSHistFuncV(const CMSHistFunc& hpdf, const RooAbsData& data,
                           bool includeZeroWeights)
    : hpdf_(hpdf), begin_(0), end_(0) {
  hpdf.updateCache();
  std::vector<int> bins;
  RooArgSet obs(hpdf.x_.arg());
  const RooRealVar& x = static_cast<const RooRealVar&>(*obs.first());
  bool aligned = true;
  for (int i = 0, n = data.numEntries(); i < n; ++i) {
    obs = *data.get(i);
    if (data.weight() == 0 && !includeZeroWeights) continue;
    int idx = hpdf.cache_.FindBin(x.getVal());
    if (!bins.empty() && idx != bins.back() + 1) aligned = false;
    bins.push_back(idx);
  }
  if (bins.empty()) {
    // nothing to do.
  } else if (aligned) {
    begin_ = bins.front();
    end_ = bins.back() + 1;
    // std::cout << "Created CMSHistFuncV from " << hpdf.GetName() << ",
    // aligned, " << (end_-begin_) << " bins." << std::endl;
  } else {
    nbins_ = bins.size();
    bins_.swap(bins);
    blocks_.clear();
    int start = bins_[0], istart = 0;
    for (int i = 1, n = bins_.size(); i < n; ++i) {
      if (bins_[i] != bins_[i - 1] + 1) {
        blocks_.push_back(Block(istart, start, bins_[i - 1] + 1));
        start = bins_[i];
        istart = i;
      }
    }
    blocks_.push_back(Block(istart, start, bins_.back() + 1));
    if (blocks_.size() < 4 * bins_.size()) {
      // std::cout << "Created CMSHistFuncV from " << hpdf.GetName() << ",
      // block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      bins_.clear();
    } else {
      // std::cout << "Created CMSHistFuncV from " << hpdf.GetName() << ",
      // non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      blocks_.clear();
    }
  }
}

void CMSHistFuncV::fill(std::vector<Double_t>& out) const {
  hpdf_.updateCache();
  if (begin_ != end_) {
    out.resize(end_ - begin_);
    std::copy(&hpdf_.cache_.GetBinContent(begin_),
              &hpdf_.cache_.GetBinContent(end_), out.begin());
  } else if (!blocks_.empty()) {
    out.resize(nbins_);
    for (auto b : blocks_)
      std::copy(&hpdf_.cache_.GetBinContent(b.begin),
                &hpdf_.cache_.GetBinContent(b.end), out.begin() + b.index);
  } else {
    out.resize(bins_.size());
    for (int i = 0, n = bins_.size(); i < n; ++i) {
      out[i] = hpdf_.cache_.GetBinContent(bins_[i]);
    }
  }
}
