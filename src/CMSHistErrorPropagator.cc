#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include <stdexcept>
#include <vector>
#include <ostream>
#include <memory>
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
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
  vbinpars_.resize(nb, nullptr);
  for (unsigned j = 0, r = 0; j < nb; ++j) {
    if (bintypes_.size() && bintypes_[j] == 1) {
      vbinpars_[j] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
      ++r;
    }
  }
  // for (unsigned i = 0; i < nb; ++i) {
  //   vbinpars_[i] = dynamic_cast<RooAbsReal const*>(binpars_.at(i));
  // }
  valsum_ = vfuncs_[0]->getCacheHisto();
  valsum_.Clear();
  scaledvalsum_ = vfuncs_[0]->getCacheHisto();
  scaledvalsum_.Clear();
  err2sum_.resize(nb, 0.);
  toterr_.resize(nb, 0.);
  valvec_.resize(nf, std::vector<double>(nb, 0.));
  err2vec_.resize(nf, std::vector<double>(nb, 0.));
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

// eval 0 = default
// eval 1 = force
// eval 2 = optimised
void CMSHistErrorPropagator::fillSumAndErr(int eval) const {
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
    if (eval == 1) {
      valsum_.Clear();
      std::fill(err2sum_.begin(), err2sum_.end(), 0.);
      for (unsigned i = 0; i < vfuncs_.size(); ++i) {
        vectorized::mul_add(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->getCacheHisto()[0]), &valsum_[0]);
        // This next one probably won't actually vectorize because of the sqrt. Check if gcc can do this (
        // may raise an error flag which prevents it)
        vectorized::mul_add_sqr(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->getBinErrors()[0]), &err2sum_[0]);
      }
      vectorized::sqrt(valsum_.size(), &err2sum_[0], &toterr_[0]);
      scaledvalsum_ = valsum_;
    } else {
      for (unsigned j = 0; j < valsum_.size(); ++j) {
        if (bintypes_.size() && bintypes_[j] == 1) {
          valsum_[j] = 0.;
          err2sum_[j] = 0.;
          // if (v > 1) std::cout << "Bin " << j << "\n";
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            valvec_[i][j] = vfuncs_[i]->getCacheHisto()[j] * coeffvals_[i];
            valsum_[j] += valvec_[i][j];
            double e =  vfuncs_[i]->getBinErrors()[j] * coeffvals_[i];
            err2vec_[i][j] = e * e;
            // if (v > 1) printf("%.6f/%.6f   ", valvec_[i][j], err2vec_[i][j]);
            err2sum_[j] += err2vec_[i][j];
          }
          scaledvalsum_[j] = valsum_[j];
          toterr_[j] = std::sqrt(err2sum_[j]);
          // if (v > 1) printf(" | %.6f/%.6f/%.6f\n", valsum_[j], err2sum_[j], toterr_[j]);
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            if (err2sum_[j] > 0. && coeffvals_[i] > 0.) {
              binmods_[i][j] = (toterr_[j] * err2vec_[i][j]) / (err2sum_[j] * coeffvals_[i]);
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
      if (bintypes_[j] == 1) {
        if (data_.size() && toterr_[j] > 0. && vbinpars_[j]->isConstant()) {
          //// Original:
          // double val = valsum_[j];
          // double err = toterr_[j] / valsum_[j];
          // double a = 1;
          // double b = val * err * err - 1.;
          // double c = - data_[j] * err * err;
          // double tmp = -0.5 * (b + copysign(1.0, b) * sqrt(b * b - 4. * a * c));
          // double x1 = tmp / a;
          // double x2 = c /tmp;
          //// Opt:

          // double val = valsum_[j];
          double err = toterr_[j] / valsum_[j];
          // double a = 1;
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
          ((RooRealVar*)vbinpars_[j])->setVal((std::max(x1, x2) - 1.) / err);
          RooAbsArg::setDirtyInhibit(false);
          // Rely on the first value client being the constraint
          vbinpars_[j]->valueClientMIterator().next()->setValueDirty();
        }
        double x = vbinpars_[j]->getVal();
        scaledvalsum_[j] = valsum_[j] + toterr_[j] * x;
        // if (valsum_[j] > 0.) {
        //   scaledvalsum_[j] = valsum_[j] * std::pow(((valsum_[j] + toterr_[j]) / valsum_[j]), x);
        // }
        if (eval == 0) {
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            scaledbinmods_[i][j] = binmods_[i][j] * x;
          }
        }
      }
    }
    RooAbsArg::setDirtyInhibit(false);
    scaledvalsum_.CropUnderflows();
    binsentry_.reset();
  }
}

RooArgList * CMSHistErrorPropagator::setupBinPars() {
  RooArgList * res = new RooArgList();
  if (bintypes_.size()) {
    std::cout << "setupBinPars() already called for " << this->GetName() << "\n";
    return res;
  }

  // First initialize all the storage
  initialize();
  // Now fill the bin contents and errors
  fillSumAndErr(1); // the arg (1) forces fillSumAndError to fill the caches for all bins

  bintypes_.resize(valsum_.size(), 0);



  std::cout << "Analysing bin errors for: " << this->GetName() << "\n";

  for (unsigned j = 0; j < valsum_.size(); ++j) {
    std::cout << "Bin " << j << ": contents=" << valsum_[j] << " error=" << toterr_[j] << "\n";
    if (toterr_[j] > 0.) {
      double alpha = (toterr_[j] * toterr_[j]) / valsum_[j];
      double n = valsum_[j] / alpha;
      double nround = int(n + 0.5);
      double alpharound = valsum_[j] / nround;
      std::cout << " -- equivalent events:         " << n << " +/- " << std::sqrt(n) << " alpha = " << alpha << "\n";
      std::cout << " -- equivalent events (round): " << nround << " +/- " << std::sqrt(nround) << " alpha = " << alpharound << "\n";
      bintypes_[j] = 1;
      RooRealVar var(TString::Format("%s_bin%i", this->GetName(), j), "", 0, -7, 7);
      res->addClone(var);
    }
  }

  binpars_.add(*res);
  binsentry_.addVars(binpars_);
  binsentry_.setValueDirty();

  for (unsigned j = 0, r = 0; j < valsum_.size(); ++j) {
    if (bintypes_[j] == 1) {
      vbinpars_[j] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
      ++r;
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
  fillSumAndErr();
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
  fillSumAndErr(1);
  return scaledvalsum_.GetAt(x_);
}


void CMSHistErrorPropagator::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                    TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  if (bintypes_.size()) {
    std::cout << ">> Bin parameters are initialized:\n";
    for (unsigned i = 0; i < bintypes_.size(); ++i) {
      std::cout << ">> " << i << "\t" << bintypes_[i]  << "\t" << vbinpars_[i] << "\n";
    }
  }
  std::cout << ">> Current cache:\n";
  valsum_.Dump();
  std::cout << ">> Current cache (bin scaled):\n";
  scaledvalsum_.Dump();
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
      fillSumAndErr(1);
      return scaledvalsum_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

CMSHistErrorPropagatorV::CMSHistErrorPropagatorV(const CMSHistErrorPropagator& hpdf, const RooAbsData& data,
                           bool includeZeroWeights)
    : hpdf_(hpdf), begin_(0), end_(0) {
  hpdf.fillSumAndErr(1);
  hpdf.data_.resize(hpdf.scaledvalsum_.size(), 0.);
  std::vector<int> bins;
  RooArgSet obs(hpdf.x_.arg());
  const RooRealVar& x = static_cast<const RooRealVar&>(*obs.first());
  bool aligned = true;
  for (int i = 0, n = data.numEntries(); i < n; ++i) {
    obs = *data.get(i);
    if (data.weight() == 0 && !includeZeroWeights) continue;
    int idx = hpdf.scaledvalsum_.FindBin(x.getVal());
    hpdf.data_[idx] = data.weight();
    if (!bins.empty() && idx != bins.back() + 1) aligned = false;
    bins.push_back(idx);
  }
  if (bins.empty()) {
    // nothing to do.
  } else if (aligned) {
    begin_ = bins.front();
    end_ = bins.back() + 1;
    // std::cout << "Created CMSHistErrorPropagatorV from " << hpdf.GetName() << ",
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
      // std::cout << "Created CMSHistErrorPropagatorV from " << hpdf.GetName() << ",
      // block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      bins_.clear();
    } else {
      // std::cout << "Created CMSHistErrorPropagatorV from " << hpdf.GetName() << ",
      // non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      blocks_.clear();
    }
  }
}

void CMSHistErrorPropagatorV::fill(std::vector<Double_t>& out) const {
  hpdf_.fillSumAndErr(1);
  if (begin_ != end_) {
    out.resize(end_ - begin_);
    std::copy(&hpdf_.scaledvalsum_.GetBinContent(begin_),
              &hpdf_.scaledvalsum_.GetBinContent(end_), out.begin());
  } else if (!blocks_.empty()) {
    out.resize(nbins_);
    for (auto b : blocks_)
      std::copy(&hpdf_.scaledvalsum_.GetBinContent(b.begin),
                &hpdf_.scaledvalsum_.GetBinContent(b.end), out.begin() + b.index);
  } else {
    out.resize(bins_.size());
    for (int i = 0, n = bins_.size(); i < n; ++i) {
      out[i] = hpdf_.scaledvalsum_.GetBinContent(bins_[i]);
    }
  }
}

