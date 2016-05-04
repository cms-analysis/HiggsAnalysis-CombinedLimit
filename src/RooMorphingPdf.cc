#include "HiggsAnalysis/CombinedLimit/interface/RooMorphingPdf.h"
#include <stdexcept>
#include <vector>
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"

RooMorphingPdf::RooMorphingPdf()
    : x_(RooRealProxy()),
      mh_(RooRealProxy()),
      pdfs_(RooListProxy()),
      masses_(std::vector<double>()),
      can_morph_(false),
      rebin_(TArrayI()),
      target_axis_(TAxis()),
      morph_axis_(TAxis()),
      init_(false),
      cache_(FastHisto()),
      single_point_(false),
      p1_(nullptr),
      p2_(nullptr),
      mh_lo_(0.),
      mh_hi_(0.) {}

RooMorphingPdf::RooMorphingPdf(const char* name, const char* title,
                               RooRealVar& x, RooAbsReal& mh,
                               RooArgList const& pdfs,
                               std::vector<double> const& masses,
                               bool const& can_morph, TAxis const& target_axis,
                               TAxis const& morph_axis)
    : RooAbsPdf(name, title),
      x_(x.GetName(), x.GetTitle(), this, x),
      mh_(mh.GetName(), mh.GetTitle(), this, mh),
      pdfs_(pdfs.GetName(), pdfs.GetTitle(), this),
      masses_(masses),
      can_morph_(can_morph),
      rebin_(TArrayI()),
      target_axis_(target_axis),
      morph_axis_(morph_axis),
      init_(false),
      cache_(FastHisto()),
      single_point_(false),
      p1_(nullptr),
      p2_(nullptr),
      mh_lo_(0.),
      mh_hi_(0.) {
  SetAxisInfo();
  TIterator* pdf_iter = pdfs.createIterator();
  RooAbsArg* pdf;
  while ((pdf = reinterpret_cast<RooAbsArg*>(pdf_iter->Next())))
    pdfs_.add(*pdf);
  delete pdf_iter;
}

RooMorphingPdf::RooMorphingPdf(const RooMorphingPdf& other, const char* name)
    : RooAbsPdf(other, name),
      x_(other.x_.GetName(), this, other.x_),
      mh_(other.mh_.GetName(), this, other.mh_),
      pdfs_(other.pdfs_.GetName(), this, other.pdfs_),
      masses_(other.masses_),
      can_morph_(other.can_morph_),
      rebin_(other.rebin_),
      target_axis_(other.target_axis_),
      morph_axis_(other.morph_axis_),
      init_(false),  // force a re-init if we copy
      cache_(other.cache_),
      single_point_(false),
      p1_(nullptr),
      p2_(nullptr),
      mh_lo_(0.),
      mh_hi_(0.) {}

void RooMorphingPdf::SetAxisInfo() {
  if (!morph_axis_.IsVariableBinSize()) {
    TArrayD new_bins(morph_axis_.GetNbins() + 1);
    new_bins[0] = morph_axis_.GetBinLowEdge(1);
    for (int i = 1; i <= morph_axis_.GetNbins(); ++i) {
      new_bins[i] = morph_axis_.GetBinUpEdge(i);
    }
    morph_axis_.Set(morph_axis_.GetNbins(), &(new_bins[0]));
  }

  rebin_ = TArrayI(morph_axis_.GetNbins());
  for (int i = 0; i < rebin_.GetSize(); ++i) {
    rebin_[i] = target_axis_.FindFixBin(morph_axis_.GetBinCenter(i+1)) - 1;
    // std::cout << "Morph bin " << i << " goes to target bin " << rebin_[i]
    //           << "\n";
  }
}

void RooMorphingPdf::Init() const {
  hmap_.clear();
  if (masses_.size() != unsigned(pdfs_.getSize())) {
    throw std::runtime_error(
        "RooMorphingPdf: Number of mass points differs from the number of "
        "stored pdfs!");
  }
  for (unsigned i = 0; i < masses_.size(); ++i) {
    FastVerticalInterpHistPdf2* tpdf =
        dynamic_cast<FastVerticalInterpHistPdf2*>(pdfs_.at(i));
    if (!tpdf)
      throw std::runtime_error(
          "RooMorphingPdf: Supplied pdf is not of type "
          "FastVerticalInterpHistPdf2");
    tpdf->evaluate();
    if (tpdf->cache().size() != unsigned(morph_axis_.GetNbins())) {
      throw std::runtime_error(
          "RooMorphingPdf: Component pdfs must have the same number of bins as "
          "the TAxis supplied for morphing");
    }
    hmap_[masses_[i]] = tpdf;
  }
  cache_ = FastHisto(TH1F("tmp", "tmp", target_axis_.GetNbins(), 0,
                          static_cast<float>(target_axis_.GetNbins())));


  mc_.nbn = morph_axis_.GetNbins();
  mc_.xminn = (*morph_axis_.GetXbins())[0];
  mc_.xmaxn = (*morph_axis_.GetXbins())[mc_.nbn];

  mc_.sigdis1.clear();
  mc_.sigdis2.clear();
  mc_.sigdisn.clear();
  mc_.xdisn.clear();
  mc_.sigdisf.clear();
  mc_.sigdis1.resize(1 + mc_.nbn);
  mc_.sigdis2.resize(1 + mc_.nbn);
  mc_.sigdisn.resize(2 * (1 + mc_.nbn));
  mc_.xdisn.resize(2 * (1 + mc_.nbn));
  mc_.sigdisf.resize(mc_.nbn + 1);

  // We don't know if mh_ is a value or a function,
  // but luckily the addVars() method will figure it
  // out for us
  sentry_.addVars(RooArgList(mh_.arg()));
  sentry_.setValueDirty();

  init_ = true;
}

Double_t RooMorphingPdf::evaluate() const {
  if (!init_ || hmap_.empty()) Init();

  // cache is empty: throw exception
  // (something went very wrong...)
  if (hmap_.empty()) {
    std::cout << "name:   " << this->GetName() << "\n";
    std::cout << "init?:  " << init_ << "\n";
    std::cout << "masses: " << masses_.size() << "\n";
    std::cout << "pdfs:   " << pdfs_.getSize() << "\n";
    throw std::runtime_error("RooMorphingPdf: Cache is empty!");
  }

  if (!sentry_.good()) {
    // Mass value has changed so we need to figure out what
    // new masspoints to take
    MassMapIter upper = hmap_.lower_bound(mh_);
    if (upper == hmap_.begin()) {
      single_point_ = true;
      p1_ = upper->second;
    } else if (upper == hmap_.end()) {
      single_point_ = true;
      --upper;
      p1_ = upper->second;
    } else {
      MassMapIter lower = upper;
      --lower;
      p1_ = lower->second;
      p2_ = upper->second;
      mh_lo_ = lower->first;
      mh_hi_ = upper->first;
      if (!can_morph_) {
        single_point_ = true;
        if (fabs(upper->first - mh_) <= fabs(mh_ - lower->first)) {
          p1_ = upper->second;
        } else {
          p1_ = lower->second;
        }
      } else {
        single_point_ = false;
      }
    }
  }

  if (single_point_) {
    if (!(p1_->cacheIsGood() && sentry_.good())) {
      p1_->evaluate();
      cache_.Clear();
      for (unsigned i = 0; i < p1_->cache().size(); ++i) {
        cache_[rebin_[i]] += p1_->cache()[i];
      }
      cache_.CropUnderflows();
      cache_.Normalize();
    }
  } else {
    if (!(p1_->cacheIsGood() && p2_->cacheIsGood() && sentry_.good())) {
      p1_->evaluate();
      p2_->evaluate();
      FastTemplate result =
          morph(p1_->cache(), p2_->cache(), mh_lo_, mh_hi_, mh_);
      cache_.Clear();
      for (unsigned i = 0; i < result.size(); ++i) {
        cache_[rebin_[i]] += result[i];
      }
      cache_.CropUnderflows();
      cache_.Normalize();
    }
  }
  sentry_.reset();
  return cache_.GetAt(x_);
}

FastTemplate RooMorphingPdf::morph(FastTemplate const& hist1,
                                   FastTemplate const& hist2, double par1,
                                   double par2, double parinterp) const {
 unsigned idebug = 0;

  // Extract bin parameters of input histograms 1 and 2.
  // Supports the cases of non-equidistant as well as equidistant binning
  // and also the case that binning of histograms 1 and 2 is different.

  TArrayD const* bedgesn = morph_axis_.GetXbins();

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

  if (idebug >= 1)
    std::cout << "New hist: " << mc_.nbn << " " << mc_.xminn << " " << mc_.xmaxn
              << std::endl;

  // Treatment for empty histograms: Return an empty histogram
  // with interpolated bins.

  if (hist1.Integral() <= 0 || hist2.Integral() <= 0) {
    std::cout << "Warning! th1morph detects an empty input histogram. Empty "
            "interpolated histogram returned: " << std::endl;
    return (TH1F("morphed", "morphed", mc_.nbn, mc_.xminn, mc_.xmaxn));
  }
  if (idebug >= 1)
    std::cout << "Input histogram content sums: " << hist1.Integral() << " "
         << hist2.Integral() << std::endl;
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

  // AG: Use already-allocated arrays
  std::fill(mc_.sigdis1.begin(), mc_.sigdis1.end(), 0.);
  std::fill(mc_.sigdis2.begin(), mc_.sigdis2.end(), 0.);
  std::fill(mc_.sigdisn.begin(), mc_.sigdisn.end(), 0.);
  std::fill(mc_.xdisn.begin(), mc_.xdisn.end(), 0.);
  std::fill(mc_.sigdisf.begin(), mc_.sigdisf.end(), 0.);
  Double_t *sigdis1 = &(mc_.sigdis1[0]);
  Double_t *sigdis2 = &(mc_.sigdis2[0]);
  Double_t *sigdisn = &(mc_.sigdisn[0]);
  Double_t *xdisn = &(mc_.xdisn[0]);
  Double_t *sigdisf = &(mc_.sigdisf[0]);

  sigdis1[0] = 0;
  sigdis2[0] = 0;  // Start with cdf=0 at left edge

  for (Int_t i = 0; i < mc_.nbn; i++) {  // Remember, bin i has edges at i-1 and
    sigdis1[i+1] = hist1[i];               // i and i runs from 1 to nb.
    sigdis2[i+1] = hist2[i];
  }

  if (idebug >= 3) {
    for (Int_t i = 0; i < mc_.nbn; i++) {
      std::cout << i << " dist1" << hist1[i] << std::endl;
    }
    for (Int_t i = 0; i < mc_.nbn; i++) {
      std::cout << i << " dist2" << hist2[i] << std::endl;
    }
  }

  // ......Normalize the distributions to 1 to obtain pdf's and integrate
  //      (sum) to obtain cdf's.
  Double_t total = 0;
  for (Int_t i = 0; i < mc_.nbn + 1; i++) {
    total += sigdis1[i];
  }
  if (idebug >= 1) std::cout << "Total histogram 1: " << total << std::endl;
  for (Int_t i = 1; i < mc_.nbn + 1; i++) {
    sigdis1[i] = sigdis1[i] / total + sigdis1[i - 1];
  }

  total = 0.;
  for (Int_t i = 0; i < mc_.nbn + 1; i++) {
    total += sigdis2[i];
  }
  if (idebug >= 1) std::cout << "Total histogram 22: " << total << std::endl;
  for (Int_t i = 1; i < mc_.nbn + 1; i++) {
    sigdis2[i] = sigdis2[i] / total + sigdis2[i - 1];
  }

  // *
  // *......We are going to step through all the edges of both input
  // *      cdf's ordered by increasing value of y. We start at the
  // *      lower edge, but first we should identify the upper ends of the
  // *      curves. These (ixl1, ixl2) are the first point in each cdf from
  // *      above that has the same integral as the last edge.
  // *

  Int_t ix1l = mc_.nbn;
  Int_t ix2l = mc_.nbn;
  while (sigdis1[ix1l - 1] >= sigdis1[ix1l]) {
    ix1l = ix1l - 1;
  }
  while (sigdis2[ix2l - 1] >= sigdis2[ix2l]) {
    ix2l = ix2l - 1;
  }

  // *
  // *......Step up to the beginnings of the curves. These (ix1, ix2) are the
  // *      first non-zero points from below.

  Int_t ix1 = -1;
  do {
    ix1 = ix1 + 1;
  } while (sigdis1[ix1 + 1] <= sigdis1[0]);

  Int_t ix2 = -1;
  do {
    ix2 = ix2 + 1;
  } while (sigdis2[ix2 + 1] <= sigdis2[0]);

  if (idebug >= 1) {
    std::cout << "First and last edge of hist1: " << ix1 << " " << ix1l
              << std::endl;
    std::cout << "   " << sigdis1[ix1] << " " << sigdis1[ix1 + 1] << std::endl;
    std::cout << "First and last edge of hist2: " << ix2 << " " << ix2l
              << std::endl;
    std::cout << "   " << sigdis2[ix2] << " " << sigdis2[ix2 + 1] << std::endl;
  }

  // ....The first interpolated point should be computed now.

  Int_t nx3 = 0;
  Double_t x1, x2, x;
  x1 = morph_axis_.GetBinLowEdge(ix1 + 1);
  x2 = morph_axis_.GetBinLowEdge(ix2 + 1);
  x = wt1 * x1 + wt2 * x2;
  xdisn[nx3] = x;
  sigdisn[nx3] = 0;
  if (idebug >= 1) {
    std::cout << "First interpolated point: " << xdisn[nx3] << " "
              << sigdisn[nx3] << std::endl;
    std::cout << "                          " << x1 << " <= " << x
              << " <= " << x2 << std::endl;
  }

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
    if ((sigdis1[ix1 + 1] <= sigdis2[ix2 + 1] || ix2 == ix2l) && ix1 < ix1l) {
      ix1 = ix1 + 1;
      while (sigdis1[ix1 + 1] <= sigdis1[ix1] && ix1 < ix1l) {
        ix1 = ix1 + 1;
      }
      i12type = 1;
    } else if (ix2 < ix2l) {
      ix2 = ix2 + 1;
      while (sigdis2[ix2 + 1] <= sigdis2[ix2] && ix2 < ix2l) {
        ix2 = ix2 + 1;
      }
      i12type = 2;
    }
    if (i12type == 1) {
      if (idebug >= 3) {
        std::cout << "Pair for i12type=1: " << sigdis2[ix2] << " "
                  << sigdis1[ix1] << " " << sigdis2[ix2 + 1] << std::endl;
      }
      x1 = morph_axis_.GetBinLowEdge(ix1 + 1);
      y = sigdis1[ix1];
      Double_t x20 = morph_axis_.GetBinLowEdge(ix2 + 1);
      Double_t x21 = morph_axis_.GetBinUpEdge(ix2 + 1);
      Double_t y20 = sigdis2[ix2];
      Double_t y21 = sigdis2[ix2 + 1];

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
        std::cout << "Pair for i12type=2: " << sigdis1[ix1] << " "
                  << sigdis2[ix2] << " " << sigdis1[ix1 + 1] << std::endl;
      }
      x2 = morph_axis_.GetBinLowEdge(ix2 + 1);
      y = sigdis2[ix2];
      Double_t x10 = morph_axis_.GetBinLowEdge(ix1 + 1);
      Double_t x11 = morph_axis_.GetBinUpEdge(ix1 + 1);
      Double_t y10 = sigdis1[ix1];
      Double_t y11 = sigdis1[ix1 + 1];

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

    x = wt1 * x1 + wt2 * x2;
    if (y > yprev) {
      nx3 = nx3 + 1;
      if (idebug >= 1) {
        std::cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" << nx3
             << ", x= " << x << ", y=" << y << ", yprev=" << yprev << std::endl;
      }
      yprev = y;
      xdisn[nx3] = x;
      sigdisn[nx3] = y;
      if (idebug >= 1) {
        std::cout << "    ix1=" << ix1 << ", ix2= " << ix2
             << ", i12type= " << i12type << ", sigdis1[ix1]=" << sigdis1[ix1]
             << std::endl;
        std::cout << "        "
             << ", nx3=" << nx3 << ", x=" << x << ", y= " << sigdisn[nx3]
             << std::endl;
      }
    }
  }
  if (idebug >= 3)
    for (Int_t i = 0; i < nx3; i++) {
      std::cout << " nx " << i << " " << xdisn[i] << " " << sigdisn[i]
                << std::endl;
    }

  // *......Now we loop over the edges of the bins of the interpolated
  // *      histogram and find out where the interpolated cdf 3
  // *      crosses them. This projection defines the result and will
  // *      be stored (after differention and renormalization) in the
  // *      output histogram.
  // *
  // *......We set all the bins following the final edge to the value
  // *      of the final edge.

  x = mc_.xmaxn;
  Int_t ix = mc_.nbn;

  if (idebug >= 1)
    std::cout << "------> Any final bins to set? " << x << " " << xdisn[nx3]
              << std::endl;
  while (x >= xdisn[nx3]) {
    sigdisf[ix] = sigdisn[nx3];
    if (idebug >= 2)
      std::cout << "   Setting final bins" << ix << " " << x << " "
                << sigdisf[ix] << std::endl;
    ix = ix - 1;
    x = (*bedgesn)[ix];
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
  x = (*bedgesn)[ix + 1];
  if (idebug >= 1)
    std::cout << "Start setting initial bins at x=" << x << std::endl;
  while (x <= xdisn[0]) {
    sigdisf[ix] = sigdisn[0];
    if (idebug >= 1)
      std::cout << "   Setting initial bins " << ix << " " << x << " "
                << xdisn[1] << " " << sigdisf[ix] << std::endl;
    ix = ix + 1;
    x = (*bedgesn)[ix + 1];
  }
  Int_t ixf = ix;

  if (idebug >= 1)
    std::cout << "Bins left to loop over:" << ixf << "-" << ixl << std::endl;

  // *......Also the end (from y to 1.0) often comes before the last edge
  // *      so we have to set the following to 1.0 as well.

  Int_t ix3 = 0;  // Problems with initial edge!!!
  for (ix = ixf; ix < ixl; ix++) {
    x = (*bedgesn)[ix];
    if (x < xdisn[0]) {
      y = 0;
    } else if (x > xdisn[nx3]) {
      y = 1.;
    } else {
      while (xdisn[ix3 + 1] <= x && ix3 < 2 * mc_.nbn) {
        ix3 = ix3 + 1;
      }
      Double_t dx2 = morph_axis_.GetBinWidth(morph_axis_.FindFixBin(x));
      if (xdisn[ix3 + 1] - x > 1.1 * dx2) {  // Empty bin treatment
        y = sigdisn[ix3 + 1];
      } else if (xdisn[ix3 + 1] > xdisn[ix3]) {  // Normal bins
        y = sigdisn[ix3] + (sigdisn[ix3 + 1] - sigdisn[ix3]) *
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
      std::cout << "   cdf n=" << sigdisn[ix3] << ", y=" << y
           << ", next point=" << sigdisn[ix3 + 1] << std::endl;
    }
  }

  // .....Differentiate interpolated cdf and return renormalized result in
  //      new histogram.

  // TH1F morphedhist("morphed", "morphed", mc_.nbn, 0, static_cast<float>(mc_.nbn));

  FastTemplate morphedhist(mc_.nbn);

  for (ix = mc_.nbn - 1; ix > -1; ix--) {
    y = sigdisf[ix + 1] - sigdisf[ix];
    morphedhist[ix] = y;
  }

  // ......All done, return the result.

  return morphedhist;
}
