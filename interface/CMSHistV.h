#ifndef CMSHistV_h
#define CMSHistV_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooArgProxy.h"
#include "RooAbsData.h"
#include "Rtypes.h"
#include "RooRealVar.h"

template <class T>
class CMSHistV {
 public:
  CMSHistV(const T&, const RooAbsData& data, bool includeZeroWeights = false);
  void fill(std::vector<Double_t>& out) const;

 private:
  const T& hpdf_;
  int begin_, end_, nbins_;
  struct Block {
    int index, begin, end;
    Block(int i, int begin_, int end_) : index(i), begin(begin_), end(end_) {}
  };
  std::vector<Block> blocks_;
  std::vector<int> bins_;
};

template <class T>
CMSHistV<T>::CMSHistV(const T& hpdf, const RooAbsData& data,
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
    int idx = hpdf.cache().FindBin(x.getVal());
    if (!bins.empty() && idx != bins.back() + 1) aligned = false;
    bins.push_back(idx);
  }
  if (bins.empty()) {
    // nothing to do.
  } else if (aligned) {
    begin_ = bins.front();
    end_ = bins.back() + 1;
    // std::cout << "Created CMSHistV from " << hpdf.GetName() << ",
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
      // std::cout << "Created CMSHistV from " << hpdf.GetName() << ",
      // block-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      bins_.clear();
    } else {
      // std::cout << "Created CMSHistV from " << hpdf.GetName() << ",
      // non-aligned, " << bins_.size() << " bins, " <<  blocks_.size() << "
      // blocks." << std::endl;
      blocks_.clear();
    }
  }
}

template <class T>
void CMSHistV<T>::fill(std::vector<Double_t>& out) const {
  hpdf_.updateCache();
  if (begin_ != end_) {
    out.resize(end_ - begin_);
    std::copy(&hpdf_.cache().GetBinContent(begin_),
              (&hpdf_.cache().GetBinContent(end_-1))+1, out.begin());
  } else if (!blocks_.empty()) {
    out.resize(nbins_);
    for (auto b : blocks_)
      std::copy(&hpdf_.cache().GetBinContent(b.begin),
                (&hpdf_.cache().GetBinContent(b.end-1))+1, out.begin() + b.index);
  } else {
    out.resize(bins_.size());
    for (int i = 0, n = bins_.size(); i < n; ++i) {
      out[i] = hpdf_.cache().GetBinContent(bins_[i]);
      // if ((int)hpdf_.cache().GetNbinsX()>bins_[i]) out[i] = hpdf_.cache().GetBinContent(bins_[i]);
      // else out[i]=0;
    }
  }
}

#endif
