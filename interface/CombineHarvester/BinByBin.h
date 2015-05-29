#ifndef CombineTools_BinByBin_h
#define CombineTools_BinByBin_h
#include <string>
#include "CombineTools/interface/CombineHarvester.h"

namespace ch {

class BinByBinFactory {
 public:
  BinByBinFactory();

  void MergeBinErrors(CombineHarvester &cb);

  void AddBinByBin(CombineHarvester &src, CombineHarvester &dest);

  void MergeAndAdd(CombineHarvester &src, CombineHarvester &dest);

  inline BinByBinFactory& SetVerbosity(unsigned verbosity) {
    v_ = verbosity;
    return *this;
  }

  inline BinByBinFactory& SetAddThreshold(double val) {
    bbb_threshold_ = val;
    return *this;
  }
  inline BinByBinFactory& SetMergeThreshold(double val) {
    merge_threshold_ = val;
    return *this;
  }

  inline BinByBinFactory& SetPattern(std::string const& pattern) {
    pattern_ = pattern;
    return *this;
  }

  inline BinByBinFactory& SetFixNorm(bool fix) {
    fix_norm_ = fix;
    return *this;
  }

 private:
  std::string pattern_;
  unsigned v_;
  double bbb_threshold_;
  double merge_threshold_;
  bool fix_norm_;
};
}

#endif
