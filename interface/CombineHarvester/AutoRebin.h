#ifndef CombineTools_AutoRebin_h
#define CombineTools_AutoRebin_h
#include <string>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
/**
 * Merges bin uncertainties and creates bin-by-bin statistical uncertainties
 *
 * Typical usage:
 * 
 *     auto bbb = ch::AutoRebin()
 *         .SetAddThreshold(0.1)
 *         .SetMergeThreshold(0.5)
 *         .SetFixNorm(true);
 *     bbb.MergeBinErrors(cb.cp().backgrounds());
 *     bbb.AddBinByBin(cb.cp().backgrounds(), cb);
 * 
 * See below for details on each class method.
 */
class AutoRebin {
 public:
  AutoRebin();

  /**
   * Create bin-by-bin shape uncertainties for every process in **src**, and
   * add these to **dest**
   *
   * The behaviour of this function is controlled by three parameters:
   * 
   *   * The first is the uncertainty threshold which determines whether a 
   *     systematic will be created for a given histogram bin
   *     (\ref SetAddThreshold). For a bin with content \f$x\f$ and bin error
   *     \f$e\f$, this threshold is thefractional error \f$e/x\f$.
   *   * The second parameters (\ref SetFixNorm) is a flag that when set to
   *     `true` will re-normalize the up and down shape templates to have the
   *     same integral as the nominal. This means the nuisance parameter for
   *     this bin is only able to change the shape of the distribution: if a
   *     bin is shifted up, all other bins are shifted down slightly to
   *     preserve the normalisation. This feature is useful when the
   *     statistical uncertainty in the process yield is considered as an
   *     independent nuisance parameter.
   *   * The third parameters is a pattern-string which prescribes how the
   *     bin-by-bin uncertainties should be named. By default this is
   *     `CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bin_$#`, but can be
   *     changed with \ref SetPattern. Note that the special term `$#`
   *     should always be included as this is replaced with the bin index.
   */
  void Rebin(CombineHarvester &src, CombineHarvester &dest);

  /**
   * By default this class only produces output on the screen when an error
   * occurs, set to a value greater than zero for more verbose output
   */
  inline AutoRebin& SetVerbosity(unsigned verbosity) {
    v_ = verbosity;
    return *this;
  }

  /**
   * The threshold for the merging algorithm
   */
  inline AutoRebin& SetEmptyBinThreshold(double val) {
    empty_bin_threshold_ = val;
    return *this;
  }


 private:
  unsigned v_;
  double empty_bin_threshold_;
};
}

#endif
