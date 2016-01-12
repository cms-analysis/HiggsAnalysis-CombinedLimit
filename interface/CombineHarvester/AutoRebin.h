#ifndef CombineTools_AutoRebin_h
#define CombineTools_AutoRebin_h
#include <string>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
/**
 * Tests for any bins below a certain threshold and if they exist merges them
 * with neighborouring bins
 *
 * Typical usage:
 * 
 *     auto rebin = ch::AutoRebin()
 *         .SetEmptyBinThreshold(0.);
 *     rebin.Rebin(cb.cp().channel({"tt"}), cb);
 * 
 * See below for details on each class method.
 */
class AutoRebin {
 public:
  AutoRebin();

  /**
   * Work out optimal binning using the total background histogram built from 
   * **src**, and apply the binning to **dest**
   * 
   *   * Builds the total background histogram and creates vectors for the
   *   original binning and the new binning 
   *
   *   * Calls the FindNewBinning function passing it the empty new_bins vector,
   *   the total background histogram, and the condition for what we consider to
   *   be an empty bin
   */
  void Rebin(CombineHarvester &src, CombineHarvester &dest);
  /**
   * Pass through the total background histogram to find bins failing the
   * required condition ("empty" bins)
   * 
   *   * Starts from the bin with maximum bin content and moves left first. For
   *   any binnings passing the condition their value of BinLowEdge is added to
   *   new_bins vector.
   *
   *   * Then moves from maximum bin content to end of distribution moving right
   *   and does the same.
   *
   *   * Creates new total background hist with the new binning and checks if
   *   all bins pass the condition. Continues recursively until they all do.
   */
  void FindNewBinning(TH1F &total_bkg, std::vector<double> &new_bins, double bin_condition);

  /**
   * Set to a value greater than zero for more verbose output
   */
  inline AutoRebin& SetVerbosity(unsigned verbosity) {
    v_ = verbosity;
    return *this;
  }

  /**
   * The threshold for what we consider an "empty" bin (usually 0)
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
