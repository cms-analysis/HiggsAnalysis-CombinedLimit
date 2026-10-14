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
   * required condition ("empty" bins). For RebinMode==0:
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
   *
   *   For RebinMode==1:
   *
   *   * Starts from bin with lowest content which fails the condition. Tries
   *   moving left and right merging bins until threshold is met.
   *
   *   * Chooses from left and right to minimise number of bins lost
   *
   *   * Repeats with new lowest bin until all bins pass threshold
   */
  void FindNewBinning(TH1F &total_bkg, std::vector<double> &new_bins, 
     double bin_condition, double bin_uncert_fraction, int mode);

  /**
   * Return bin with maximum value of fractional error
   */
  int GetMaximumFracUncertBin(TH1F &total_bkg);
  
  /**
   * Set to a value greater than zero for more verbose output
   */
  inline AutoRebin& SetVerbosity(unsigned verbosity) {
    v_ = verbosity;
    return *this;
  }

  /**
   * The threshold for which we consider merging bins containing less than this
   * value
   */
  inline AutoRebin& SetBinThreshold(double val) {
    bin_threshold_ = val;
    return *this;
  }
  
  /**
   * The threshold on the bin uncertainty fraction for which we consider merging
   * bins containing less than this value
   */
  inline AutoRebin& SetBinUncertFraction(double val) {
    bin_uncert_fraction_ = val;
    return *this;
  }
  
  /**
   * Flag to perform different versions of rebinning algorithm
   */
  inline AutoRebin& SetRebinMode(int val) {
    rebin_mode_ = val;
    return *this;
  }
  
  /**
   * Whether to actually perform the rebinning or just print it out
   */
  inline AutoRebin& SetPerformRebin(bool val) {
    perform_rebin_ = val;
    return *this;
  }


 private:
  unsigned v_;
  int rebin_mode_;
  bool perform_rebin_;
  double bin_threshold_;
  double bin_uncert_fraction_;
};
}

#endif
