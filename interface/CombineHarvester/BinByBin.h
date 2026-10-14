#ifndef CombineTools_BinByBin_h
#define CombineTools_BinByBin_h
#include <string>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
/**
 * Merges bin uncertainties and creates bin-by-bin statistical uncertainties
 *
 * Typical usage:
 * 
 *     auto bbb = ch::BinByBinFactory()
 *         .SetAddThreshold(0.1)
 *         .SetMergeThreshold(0.5)
 *         .SetFixNorm(true);
 *     bbb.MergeBinErrors(cb.cp().backgrounds());
 *     bbb.AddBinByBin(cb.cp().backgrounds(), cb);
 * 
 * See below for details on each class method.
 */
class BinByBinFactory {
 public:
  BinByBinFactory();

  /**
   * Merges histogram bin errors between processes
   *
   * The method will loop through each analysis category and apply the merging
   * algorithm to the histograms of all processes within that category:
   *
   *  * In a given histogram bin consider contributions from the \f$j\f$
   *    processes with bin contents \f$x_{i}\f$, bin errors \f$e_{i}\f$
   *    and with \f$e_{i}/x_{i}\f$ larger than the bin-by-bin addition
   *    threshold (\ref SetAddThreshold). Note that by setting this threshold
   *    to zero the contributions from all processes will be considered.
   *  * Evaluate the total squared uncertainty from these processes: 
   *    \f$e_{TOT}^{2} = e_{1}^{2} + e_{2}^{2} + ... + e_{j}^{2} \f$
   *  * Rank bins from lowest to highest fraction error: \f$e_{i}^{2}/e_{TOT}^{2}\f$
   *  * Start removing bin errors from the lowest ranked bins as long as the total
   *    sum of squared error removed (\f$r_{sub}\f$) is less than the merge threshold
   *    (\ref SetMergeThreshold).
   *  * Once as many bin errors have been removed as possible, scale up the
   *    remaining errors by a factor \f$\sqrt{1/(1-r_{sub})}\f$, which ensure the
   *    total bin error from all processes remains the same.
   *  * Note that the BinByBinFactory default is a merging threshold of zero,
   *    in which no merging takes place. The most extreme threshold is 1, in
   *    which all errors will be merged into a single process.
   *
   * @note For a valid result all analysis event categories must have a unique
   * `bin()` value.
   */
  void MergeBinErrors(CombineHarvester &cb);

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
  void AddBinByBin(CombineHarvester &src, CombineHarvester &dest);

  /**
   * A convenience function which calls \ref MergeBinErrors and \ref
   * AddBinByBin in turn
   */
  void MergeAndAdd(CombineHarvester &src, CombineHarvester &dest);

  /**
   * By default this class only produces output on the screen when an error
   * occurs, set to a value greater than zero for more verbose output
   */
  inline BinByBinFactory& SetVerbosity(unsigned verbosity) {
    v_ = verbosity;
    return *this;
  }

  /**
   * Set the fractional bin error threshold for bin-by-bin creation and
   * for participation in the merging algorithm
   */
  inline BinByBinFactory& SetAddThreshold(double val) {
    bbb_threshold_ = val;
    return *this;
  }

  /**
   * The threshold for the merging algorithm
   */
  inline BinByBinFactory& SetMergeThreshold(double val) {
    merge_threshold_ = val;
    return *this;
  }

  /**
   * The pattern-string for the systematic naming convention
   */
  inline BinByBinFactory& SetPattern(std::string const& pattern) {
    pattern_ = pattern;
    return *this;
  }

  /**
   * Whether or not the bin-by-bin systematics are allowed to vary the process
   * normalisation
   */
  inline BinByBinFactory& SetFixNorm(bool fix) {
    fix_norm_ = fix;
    return *this;
  }

  /**
   * Construct approximate Poisson uncertainties instead of default Gaussian
   */
  inline BinByBinFactory& SetPoissonErrors(bool poisson_errors) {
    poisson_errors_ = poisson_errors;
    return *this;
  }

  /**
   * Set whether bins with zero content can participate in the merging procedure
   */
  inline BinByBinFactory& SetMergeZeroBins(bool merge) {
    merge_zero_bins_ = merge;
    return *this;
  }

  /**
   * Set whether bins with error >= content participate in the merging procedure
   */
  inline BinByBinFactory& SetMergeSaturatedBins(bool merge) {
    merge_saturated_bins_ = merge;
    return *this;
  }

  /**
   * Getter functions for class attributes
   */
  inline unsigned GetVerbosity() {return v_;}
  inline double GetAddThreshold() {return bbb_threshold_;}
  inline double GetMergeThreshold() {return merge_threshold_;}
  inline std::string GetPattern() {return pattern_;}
  inline bool GetFixNorm() {return fix_norm_;}
  inline bool GetPoissonErrors() {return poisson_errors_;}
  inline bool GetMergeZeroBins() {return merge_zero_bins_;}
  inline bool GetMergeSaturatedBins() {return merge_saturated_bins_;}

 private:
  std::string pattern_;
  unsigned v_;
  double bbb_threshold_;
  double merge_threshold_;
  bool fix_norm_;
  bool poisson_errors_;
  bool merge_zero_bins_;
  bool merge_saturated_bins_;
};
}

#endif
