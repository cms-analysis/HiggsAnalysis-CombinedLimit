#include "HiggsAnalysis/CombinedLimit/interface/RooMorphingPdf2.h"
#include <stdexcept>
#include <vector>
#include <ostream>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"
#include "TVector.h"

CMSHistFunc::CMSHistFunc() {
  morph_strategy_ = 0;
}

CMSHistFunc::CMSHistFunc(const char* name, const char* title, RooRealVar& x,
                         TH1 const& hist)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      vmorphs_("vmorphs", "", this),
      hmorphs_("hmorphs", "", this),
      cache_(hist) {}

CMSHistFunc::CMSHistFunc(CMSHistFunc const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      vmorphs_("vmorphs", this, other.vmorphs_),
      hmorphs_("hmorphs", this, other.hmorphs_),
      cache_(other.cache_) {}


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
    // mc_.nbn = cache_.size();
    // mc_.xminn = cache_.GetEdge(0);
    // mc_.xmaxn = cache_.GetEdge(mc_.nbn);
    // mc_.sigdis1.clear();
    // mc_.sigdis2.clear();
    // mc_.sigdisn.clear();
    // mc_.xdisn.clear();
    // mc_.sigdisf.clear();
    // mc_.sigdis1.resize(1 + mc_.nbn);
    // mc_.sigdis2.resize(1 + mc_.nbn);
    // mc_.sigdisn.resize(2 * (1 + mc_.nbn));
    // mc_.xdisn.resize(2 * (1 + mc_.nbn));
    // mc_.sigdisf.resize(mc_.nbn + 1);
    // FNLOG(std::cout) << "Prepared HMorphCache with nbn = " << mc_.nbn
    //                  << ", xminn=" << mc_.xminn << ", xmaxn=" << mc_.xmaxn
    //                  << "\n";
    global_.bedgesn.resize(cache_.size() + 1);
    for (unsigned i = 0; i < global_.bedgesn.size(); ++i) {
     global_.bedgesn[i] = cache_.GetEdge(i);
    }

  }
  unsigned n_vpoints = vmorphs_.getSize();
  storage_.resize(n_hpoints * (1 + n_vpoints * 2));
  std::cout << "Storage size set to: " << storage_.size() << "\n";

  hmorph_cache_.resize(n_hpoints * (1 + n_vpoints * 2));
}

unsigned CMSHistFunc::getIdx(unsigned hindex, unsigned hpoint, unsigned vindex, unsigned vpoint) const {
  unsigned n_vpoints = vmorphs_.getSize();
  return hindex + hpoint * (1 + n_vpoints * 2) + (vindex * 2 - (vindex > 0 ? 1 : 0)) + vpoint;
}


void CMSHistFunc::setShape(unsigned hindex, unsigned hpoint, unsigned vindex, unsigned vpoint, TH1 const& hist) {
  unsigned idx = getIdx(hindex, hpoint, vindex, vpoint);
  std::cout << "hindex: " << hindex << " hpoint: " << hpoint
            << " vindex: " << vindex << " vpoint: " << vpoint
            << " mapped to idx: " << idx << "\n";
  storage_.at(idx) = FastHisto(hist);
}

void CMSHistFunc::prepareInterpCaches() const {
  for (unsigned h = 0; h < hpoints_.at(0).size() - 1; ++h) {
    std::cout << "Setting interpolation point caches for hpoints " << h << ":" << h+1 << "\n";
    prepareInterpCache(hmorph_cache_.at(getIdx(0, h, 0, 0)), hmorph_cache_.at(getIdx(0, h+1, 0, 0)));
    for (int v = 0; v < vmorphs_.getSize(); ++v) {

    }
  }
}

void CMSHistFunc::prepareInterpCache(HMorphCache& c1,
                                     HMorphCache const& c2) const {
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

void CMSHistFunc::setCdf(HMorphCache& c, FastHisto const& h) const {
  c.cdf = FastTemplate(h.size() + 1);
  c.integral = h.Integral();
  for (unsigned i = 1; i < c.cdf.size(); ++i) {
    c.cdf[i] = h[i - 1] / c.integral + c.cdf[i - 1];
  }
  c.cdf_set = true;
}

Double_t CMSHistFunc::evaluate() const {
  veval = 1;
  LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  // FNLOG(std::cout) << x_ << "\n";

  if (morph_strategy_ == 0) {
    FNLOGC(std::cout, veval) << "Morphing strategy 0\n";
    FNLOGC(std::cout, veval) << "Number of horizontal morphs: " << hmorphs_.getSize() << "\n";
    FNLOGC(std::cout, veval) << "Horizontal morph sentry: " << hmorph_sentry_.good() << "\n";
    FNLOGC(std::cout, veval) << "Vertical morph sentry: " << vmorph_sentry_.good() << "\n";
    FNLOGC(std::cout, veval) << "single_point,p1,p2: " << global_.single_point << " " << global_.p1 << " " << global_.p2 << "\n";
    // Figure out what we're doing:
    //  - singlePoint p1 OR
    //  - h-morphing between p1 and p2
    bool step1 = !hmorph_sentry_.good();
    bool step2 = !vmorph_sentry_.good();
    // step1 = true;

    if (!hmorph_sentry_.good() && hmorphs_.getSize() == 1) {
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
    hmorph_sentry_.reset();

    if (step1 && !global_.single_point) {
      FNLOGC(std::cout, veval) << "Checking step 1\n";
      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        for (int vi = 0; vi < (v == 0 ? 1 : 2); ++vi) {
          unsigned idx1 = getIdx(0, global_.p1, v, vi);
          unsigned idx2 = getIdx(0, global_.p2, v, vi);
          if (!hmorph_cache_[idx1].cdf_set) {
            FNLOGC(std::cout, veval) << "Setting cdf for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
            setCdf(hmorph_cache_[idx1], storage_[idx1]);
          }
          if (!hmorph_cache_[idx2].cdf_set) {
            FNLOGC(std::cout, veval) << "Setting cdf for " << 0 << " " << global_.p2 << " " << v << " " << vi << "\n";
            setCdf(hmorph_cache_[idx2], storage_[idx2]);
          }
          if (!hmorph_cache_[idx1].interp_set) {
            FNLOGC(std::cout, veval) << "Setting interp for " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
            prepareInterpCache(hmorph_cache_[idx1], hmorph_cache_[idx2]);
          }
          FNLOGC(std::cout, veval) << "Doing horizontal morph for  " << 0 << " " << global_.p1 << " " << v << " " << vi << "\n";
          double val = ((RooRealVar*)(hmorphs_.at(0)))->getVal();
          double x1 = hpoints_[0][global_.p1];
          double x2 = hpoints_[0][global_.p2];
          double y1 = hmorph_cache_[idx1].integral;
          double y2 = hmorph_cache_[idx2].integral;
          hmorph_cache_[idx1].step1 = morph2(idx1, idx2, x1, x2, val);
          hmorph_cache_[idx1].step1.CropUnderflows();
          double ym = y1 + ((y2 - y1) / (x2 - x1)) * (val - x1);
          hmorph_cache_[idx1].step1.Scale(ym/hmorph_cache_[idx1].step1.Integral());
        }
        if (v >= 1) {
          FNLOGC(std::cout, veval) << "Setting sumdiff for vmorph " << v << "\n";
          unsigned idx = getIdx(0, global_.p1, 0, 0);
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = hmorph_cache_[idxLo].step1;
          FastTemplate hi = hmorph_cache_[idxHi].step1;
          hi.Subtract(hmorph_cache_[idx].step1);
          lo.Subtract(hmorph_cache_[idx].step1);
          // TODO: could skip the next two lines if .sum and .diff have been set before
          hmorph_cache_[idxLo].sum = hmorph_cache_[idx].step1;
          hmorph_cache_[idxLo].diff = hmorph_cache_[idx].step1;
          FastTemplate::SumDiff(hi, lo, hmorph_cache_[idxLo].sum, hmorph_cache_[idxLo].diff);
        }
      }
    }

    if (step1 && global_.single_point) {
      FNLOGC(std::cout, veval) << "Checking step 1 (single point)\n";
      for (int v = 0; v < vmorphs_.getSize() + 1; ++v) {
        if (v >= 1) {
          FNLOGC(std::cout, veval) << "Setting sumdiff for vmorph " << v << "\n";
          unsigned idx = getIdx(0, global_.p1, 0, 0);
          unsigned idxLo = getIdx(0, global_.p1, v, 0);
          unsigned idxHi = getIdx(0, global_.p1, v, 1);
          FastTemplate lo = storage_[idxLo];
          FastTemplate hi = storage_[idxHi];
          hi.Subtract(storage_[idx]);
          lo.Subtract(storage_[idx]);
          // TODO: could skip the next two lines if .sum and .diff have been set before
          hmorph_cache_[idxLo].sum = storage_[idx];
          hmorph_cache_[idxLo].diff = storage_[idx];
          FastTemplate::SumDiff(hi, lo, hmorph_cache_[idxLo].sum, hmorph_cache_[idxLo].diff);
        }
      }
    }

    if (step2) {
      FNLOGC(std::cout, veval) << "Checking step 2\n";
      unsigned idx = getIdx(0, global_.p1, 0, 0);
      hmorph_cache_[idx].step2 = storage_[idx];
      for (int v = 0; v < vmorphs_.getSize(); ++v) {
        unsigned vidx = getIdx(0, global_.p1, v+1, 0);

        double x = ((RooRealVar&)vmorphs_[v]).getVal();
        //     // std::cout << "Morphing for " << vmorphs_[i].GetName() << " with value: " << x << "\n";
        double a = 0.5*x;
        double b = smoothStepFunc(x);
        hmorph_cache_[idx].step2.Meld(hmorph_cache_[vidx].diff, hmorph_cache_[vidx].sum, a, b);

      }
      hmorph_cache_[idx].step2.CropUnderflows();
      cache_ = hmorph_cache_[idx].step2;
    }
    vmorph_sentry_.reset();
  }
  
  return cache_.GetAt(x_);
}

void CMSHistFunc::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  cache_.Dump();
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
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

FastTemplate CMSHistFunc::morph2(unsigned globalIdx1, unsigned globalIdx2,
                                 double par1, double par2,
                                 double parinterp) const {

  // FIXME not a very elegant way to check if the cache is ready
  // if (hmorph_cache_.at(0).y.size() == 0) {
  //   prepareInterpCaches();
  // }
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

  HMorphCache const& c1 = hmorph_cache_[globalIdx1];
  // HMorphNewCDFCache const& c2 = hmorph_cache_[globalIdx2];
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


// FastTemplate CMSHistFunc::morph(FastTemplate const& hist1,
//                                    FastTemplate const& hist2, double par1,
//                                    double par2, double parinterp) const {
//  unsigned idebug = 3;

//   // Extract bin parameters of input histograms 1 and 2.
//   // Supports the cases of non-equidistant as well as equidistant binning
//   // and also the case that binning of histograms 1 and 2 is different.
//  TArrayD *bedgesn = new TArrayD(cache_.size() + 1);
//  for (int i = 0; i < bedgesn->GetSize(); ++i) {
//   (*bedgesn)[i] = cache_.GetEdge(i);
//  }

//   // ......The weights (wt1,wt2) are the complements of the "distances" between
//   //       the values of the parameters at the histograms and the desired
//   //       interpolation point. For example, wt1=0, wt2=1 means that the
//   //       interpolated histogram should be identical to input histogram 2.
//   //       Check that they make sense. If par1=par2 then we can choose any
//   //       valid set of wt1,wt2 so why not take the average?

//   Double_t wt1, wt2;
//   if (par2 != par1) {
//     wt1 = 1. - (parinterp - par1) / (par2 - par1);
//     wt2 = 1. + (parinterp - par2) / (par2 - par1);
//   } else {
//     wt1 = 0.5;
//     wt2 = 0.5;
//   }

//   // ......Give a warning if this is an extrapolation.

//   if (wt1 < 0 || wt1 > 1. || wt2 < 0. || wt2 > 1. ||
//       fabs(1 - (wt1 + wt2)) > 1.0e-4) {
//     std::cout << "Warning! th1fmorph: This is an extrapolation!! Weights are "
//               << wt1 << " and " << wt2 << " (sum=" << wt1 + wt2 << ")"
//               << std::endl;
//   }
//   if (idebug >= 1)
//     std::cout << "th1morph - Weights: " << wt1 << " " << wt2 << std::endl;

//   if (idebug >= 1)
//     std::cout << "New hist: " << mc_.nbn << " " << mc_.xminn << " " << mc_.xmaxn
//               << std::endl;

//   // Treatment for empty histograms: Return an empty histogram
//   // with interpolated bins.

//   if (hist1.Integral() <= 0 || hist2.Integral() <= 0) {
//     std::cout << "Warning! th1morph detects an empty input histogram. Empty "
//             "interpolated histogram returned: " << std::endl;
//     return (TH1F("morphed", "morphed", mc_.nbn, mc_.xminn, mc_.xmaxn));
//   }
//   if (idebug >= 1)
//     std::cout << "Input histogram content sums: " << hist1.Integral() << " "
//          << hist2.Integral() << std::endl;
//   // *
//   // *......Extract the single precision histograms into double precision arrays
//   // *      for the interpolation computation. The offset is because sigdis(i)
//   // *      describes edge i (there are nbins+1 of them) while dist1/2
//   // *      describe bin i. Be careful, ROOT does not use C++ convention to
//   // *      number bins: dist1[ibin] is content of bin ibin where ibin runs from
//   // *      1 to nbins. We allocate some extra space for the derived
//   // distributions
//   // *      because there may be as many as nb1+nb2+2 edges in the intermediate
//   // *      interpolated cdf described by xdisn[i] (position of edge i) and
//   // *      sigdisn[i] (cummulative probability up this edge) before we project
//   // *      into the final binning.
//   // Float_t const* dist1 = hist1.GetArray();
//   // Float_t const* dist2 = hist2.GetArray();

//   // AG: Use already-allocated arrays
//   std::fill(mc_.sigdis1.begin(), mc_.sigdis1.end(), 0.);
//   std::fill(mc_.sigdis2.begin(), mc_.sigdis2.end(), 0.);
//   std::fill(mc_.sigdisn.begin(), mc_.sigdisn.end(), 0.);
//   std::fill(mc_.xdisn.begin(), mc_.xdisn.end(), 0.);
//   std::fill(mc_.sigdisf.begin(), mc_.sigdisf.end(), 0.);
//   Double_t *sigdis1 = &(mc_.sigdis1[0]);
//   Double_t *sigdis2 = &(mc_.sigdis2[0]);
//   Double_t *sigdisn = &(mc_.sigdisn[0]);
//   Double_t *xdisn = &(mc_.xdisn[0]);
//   Double_t *sigdisf = &(mc_.sigdisf[0]);

//   sigdis1[0] = 0;
//   sigdis2[0] = 0;  // Start with cdf=0 at left edge

//   for (Int_t i = 0; i < mc_.nbn; i++) {  // Remember, bin i has edges at i-1 and
//     sigdis1[i+1] = hist1[i];               // i and i runs from 1 to nb.
//     sigdis2[i+1] = hist2[i];
//   }

//   if (idebug >= 3) {
//     for (Int_t i = 0; i < mc_.nbn; i++) {
//       std::cout << i << " dist1[" << i << "] " << hist1[i] << std::endl;
//     }
//     for (Int_t i = 0; i < mc_.nbn; i++) {
//       std::cout << i << " dist2[" << i << "] " << hist2[i] << std::endl;
//     }
//   }

//   // ......Normalize the distributions to 1 to obtain pdf's and integrate
//   //      (sum) to obtain cdf's.
//   Double_t total = 0;
//   for (Int_t i = 0; i < mc_.nbn + 1; i++) {
//     total += sigdis1[i];
//   }
//   if (idebug >= 1) std::cout << "Total histogram 1: " << total << std::endl;
//   for (Int_t i = 1; i < mc_.nbn + 1; i++) {
//     sigdis1[i] = sigdis1[i] / total + sigdis1[i - 1];
//   }

//   total = 0.;
//   for (Int_t i = 0; i < mc_.nbn + 1; i++) {
//     total += sigdis2[i];
//   }
//   if (idebug >= 1) std::cout << "Total histogram 2: " << total << std::endl;
//   for (Int_t i = 1; i < mc_.nbn + 1; i++) {
//     sigdis2[i] = sigdis2[i] / total + sigdis2[i - 1];
//   }

//   if (idebug >= 3) {
//     for (Int_t i = 0; i < mc_.nbn + 1; i++) {
//       std::cout << i << " sigdis1[" << i << "] " << sigdis1[i] << std::endl;
//     }
//     for (Int_t i = 0; i < mc_.nbn + 1; i++) {
//       std::cout << i << " sigdis2[" << i << "] " << sigdis2[i] << std::endl;
//     }
//   }

//   // *
//   // *......We are going to step through all the edges of both input
//   // *      cdf's ordered by increasing value of y. We start at the
//   // *      lower edge, but first we should identify the upper ends of the
//   // *      curves. These (ixl1, ixl2) are the first point in each cdf from
//   // *      above that has the same integral as the last edge.
//   // *

//   Int_t ix1l = mc_.nbn;
//   Int_t ix2l = mc_.nbn;
//   while (sigdis1[ix1l - 1] >= sigdis1[ix1l]) {
//     ix1l = ix1l - 1;
//   }
//   while (sigdis2[ix2l - 1] >= sigdis2[ix2l]) {
//     ix2l = ix2l - 1;
//   }

//   // *
//   // *......Step up to the beginnings of the curves. These (ix1, ix2) are the
//   // *      first non-zero points from below.

//   Int_t ix1 = -1;
//   do {
//     ix1 = ix1 + 1;
//   } while (sigdis1[ix1 + 1] <= sigdis1[0]);

//   Int_t ix2 = -1;
//   do {
//     ix2 = ix2 + 1;
//   } while (sigdis2[ix2 + 1] <= sigdis2[0]);

//   if (idebug >= 1) {
//     std::cout << "First and last edge of hist1: " << ix1 << " " << ix1l
//               << std::endl;
//     std::cout << "   " << sigdis1[ix1] << " " << sigdis1[ix1 + 1] << std::endl;
//     std::cout << "First and last edge of hist2: " << ix2 << " " << ix2l
//               << std::endl;
//     std::cout << "   " << sigdis2[ix2] << " " << sigdis2[ix2 + 1] << std::endl;
//   }

//   // ....The first interpolated point should be computed now.

//   Int_t nx3 = 0;
//   Double_t x1, x2, x;
//   // x1 = morph_axis_.GetBinLowEdge(ix1 + 1);
//   // x2 = morph_axis_.GetBinLowEdge(ix2 + 1);
//   x1 = cache_.GetEdge(ix1);
//   x2 = cache_.GetEdge(ix2);
//   x = wt1 * x1 + wt2 * x2;
//   xdisn[nx3] = x;
//   sigdisn[nx3] = 0;
//   if (idebug >= 1) {
//     std::cout << "First interpolated point: " << xdisn[nx3] << " "
//               << sigdisn[nx3] << std::endl;
//     std::cout << "                          " << x1 << " <= " << x
//               << " <= " << x2 << std::endl;
//   }

//   // .....Loop over the remaining point in both curves. Getting the last
//   //      points may be a bit tricky due to limited floating point
//   //      precision.

//   if (idebug >= 1) {
//     std::cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l
//          << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;
//   }

//   Double_t yprev = -1;  // The probability y of the previous point, it will
//   // get updated and used in the loop.
//   Double_t y = 0;
//   while ((ix1 < ix1l) | (ix2 < ix2l)) {
//     if (idebug >= 1)
//       std::cout << "----Top of while with ix1=" << ix1 << ", ix1l=" << ix1l
//            << ", ix2=" << ix2 << ", ix2l=" << ix2l << std::endl;

//     // .....Increment to the next lowest point. Step up to the next
//     //      kink in case there are several empty (flat in the integral)
//     //      bins.

//     Int_t i12type = -1;  // Tells which input distribution we need to
//                          // see next point of.
//     if ((sigdis1[ix1 + 1] <= sigdis2[ix2 + 1] || ix2 == ix2l) && ix1 < ix1l) {
//       ix1 = ix1 + 1;
//       while (sigdis1[ix1 + 1] <= sigdis1[ix1] && ix1 < ix1l) {
//         ix1 = ix1 + 1;
//       }
//       i12type = 1;
//     } else if (ix2 < ix2l) {
//       ix2 = ix2 + 1;
//       while (sigdis2[ix2 + 1] <= sigdis2[ix2] && ix2 < ix2l) {
//         ix2 = ix2 + 1;
//       }
//       i12type = 2;
//     }
//     if (i12type == 1) {
//       if (idebug >= 3) {
//         std::cout << "Pair for i12type=1: " << sigdis2[ix2] << " "
//                   << sigdis1[ix1] << " " << sigdis2[ix2 + 1] << std::endl;
//       }
//       x1 = cache_.GetEdge(ix1);
//       y = sigdis1[ix1];
//       Double_t x20 = cache_.GetEdge(ix2);
//       Double_t x21 = cache_.GetEdge(ix2 + 1);
//       Double_t y20 = sigdis2[ix2];
//       Double_t y21 = sigdis2[ix2 + 1];

//       // .....Calculate where the cummulative probability y in distribution 1
//       //      intersects between the 2 points from distribution 2 which
//       //      bracket it.

//       if (y21 > y20) {
//         x2 = x20 + (x21 - x20) * (y - y20) / (y21 - y20);
//       } else {
//         x2 = x20;
//       }
//     } else {
//       if (idebug >= 3) {
//         std::cout << "Pair for i12type=2: " << sigdis1[ix1] << " "
//                   << sigdis2[ix2] << " " << sigdis1[ix1 + 1] << std::endl;
//       }
//       x2 = cache_.GetEdge(ix2);
//       y = sigdis2[ix2];
//       Double_t x10 = cache_.GetEdge(ix1);
//       Double_t x11 = cache_.GetEdge(ix1 + 1);
//       Double_t y10 = sigdis1[ix1];
//       Double_t y11 = sigdis1[ix1 + 1];

//       // .....Calculate where the cummulative probability y in distribution 2
//       //      intersects between the 2 points from distribution 1 which
//       //      brackets it.

//       if (y11 > y10) {
//         x1 = x10 + (x11 - x10) * (y - y10) / (y11 - y10);
//       } else {
//         x1 = x10;
//       }
//     }

//     // .....Interpolate between the x's in the 2 distributions at the
//     //      cummulative probability y. Store the (x,y) for provisional
//     //      edge nx3 in (xdisn[nx3],sigdisn[nx3]). nx3 grows for each point
//     //      we add the the arrays. Note: Should probably turn the pair into
//     //      a structure to make the code more object-oriented and readable.

//     x = wt1 * x1 + wt2 * x2;
//     if (y > yprev) {
//       nx3 = nx3 + 1;
//       if (idebug >= 1) {
//         std::cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" << nx3
//              << ", x= " << x << ", y=" << y << ", yprev=" << yprev << std::endl;
//       }
//       yprev = y;
//       xdisn[nx3] = x;
//       sigdisn[nx3] = y;
//       if (idebug >= 1) {
//         std::cout << "    ix1=" << ix1 << ", ix2= " << ix2
//              << ", i12type= " << i12type << ", sigdis1[ix1]=" << sigdis1[ix1]
//              << std::endl;
//         std::cout << "        "
//              << ", nx3=" << nx3 << ", x=" << x << ", y= " << sigdisn[nx3]
//              << std::endl;
//       }
//     }
//   }
//   if (idebug >= 3)
//     for (Int_t i = 0; i < nx3; i++) {
//       std::cout << " nx " << i << " " << xdisn[i] << " " << sigdisn[i]
//                 << std::endl;
//     }

//   // *......Now we loop over the edges of the bins of the interpolated
//   // *      histogram and find out where the interpolated cdf 3
//   // *      crosses them. This projection defines the result and will
//   // *      be stored (after differention and renormalization) in the
//   // *      output histogram.
//   // *
//   // *......We set all the bins following the final edge to the value
//   // *      of the final edge.

//   x = mc_.xmaxn;
//   Int_t ix = mc_.nbn;

//   if (idebug >= 1)
//     std::cout << "------> Any final bins to set? " << x << " " << xdisn[nx3]
//               << std::endl;
//   while (x >= xdisn[nx3]) {
//     sigdisf[ix] = sigdisn[nx3];
//     if (idebug >= 2)
//       std::cout << "   Setting final bins" << ix << " " << x << " "
//                 << sigdisf[ix] << std::endl;
//     ix = ix - 1;
//     x = (*bedgesn)[ix];
//   }
//   Int_t ixl = ix + 1;
//   if (idebug >= 1) std::cout << " Now ixl=" << ixl << " ix=" << ix << std::endl;

//   // *
//   // *......The beginning may be empty, so we have to step up to the first
//   // *      edge where the result is nonzero. We zero the bins which have
//   // *      and upper (!) edge which is below the first point of the
//   // *      cummulative distribution we are going to project to this
//   // *      output histogram binning.
//   // *

//   ix = 0;
//   x = (*bedgesn)[ix + 1];
//   if (idebug >= 1)
//     std::cout << "Start setting initial bins at x=" << x << std::endl;
//   while (x <= xdisn[0]) {
//     sigdisf[ix] = sigdisn[0];
//     if (idebug >= 1)
//       std::cout << "   Setting initial bins " << ix << " " << x << " "
//                 << xdisn[1] << " " << sigdisf[ix] << std::endl;
//     ix = ix + 1;
//     x = (*bedgesn)[ix + 1];
//   }
//   Int_t ixf = ix;

//   if (idebug >= 1)
//     std::cout << "Bins left to loop over:" << ixf << "-" << ixl << std::endl;

//   // *......Also the end (from y to 1.0) often comes before the last edge
//   // *      so we have to set the following to 1.0 as well.

//   Int_t ix3 = 0;  // Problems with initial edge!!!
//   for (ix = ixf; ix < ixl; ix++) {
//     x = (*bedgesn)[ix];
//     if (x < xdisn[0]) {
//       y = 0;
//     } else if (x > xdisn[nx3]) {
//       y = 1.;
//     } else {
//       while (xdisn[ix3 + 1] <= x && ix3 < 2 * mc_.nbn) {
//         ix3 = ix3 + 1;
//       }
//       int bin = cache_.FindBin(x);
//       if (bin == -1) bin = 0;
//       std::cout << "x=" << x << ", bin=" << bin << "\n";
//       Double_t dx2 = cache_.GetWidth(bin);
//       if (xdisn[ix3 + 1] - x > 1.1 * dx2) {  // Empty bin treatment
//         y = sigdisn[ix3 + 1];
//       } else if (xdisn[ix3 + 1] > xdisn[ix3]) {  // Normal bins
//         y = sigdisn[ix3] + (sigdisn[ix3 + 1] - sigdisn[ix3]) *
//                                (x - xdisn[ix3]) / (xdisn[ix3 + 1] - xdisn[ix3]);
//       } else {  // Is this ever used?
//         y = 0;
//         std::cout << "Warning - th1fmorph: This probably shoudn't happen! "
//                   << std::endl;
//         std::cout << "Warning - th1fmorph: Zero slope solving x(y)"
//                   << std::endl;
//       }
//     }
//     sigdisf[ix] = y;
//     if (idebug >= 3) {
//       std::cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3]
//                 << ", x=" << x << ", next xdisn=" << xdisn[ix3 + 1]
//                 << std::endl;
//       std::cout << "   cdf n=" << sigdisn[ix3] << ", y=" << y
//            << ", next point=" << sigdisn[ix3 + 1] << std::endl;
//     }
//   }

//   // .....Differentiate interpolated cdf and return renormalized result in
//   //      new histogram.

//   // TH1F morphedhist("morphed", "morphed", mc_.nbn, 0, static_cast<float>(mc_.nbn));

//   FastTemplate morphedhist(mc_.nbn);

//   for (ix = mc_.nbn - 1; ix > -1; ix--) {
//     y = sigdisf[ix + 1] - sigdisf[ix];
//     morphedhist[ix] = y;
//   }

//   // ......All done, return the result.

//   return morphedhist;
// }