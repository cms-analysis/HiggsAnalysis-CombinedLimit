#include "CombineHarvester/CombineTools/interface/AutoRebin.h" 
#include <iostream> 
#include <string> 
#include <vector> 
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

namespace ch {

AutoRebin::AutoRebin()
    : v_(0),
      rebin_mode_(0), 
      perform_rebin_(true),
      bin_threshold_(0.), 
      bin_uncert_fraction_(100.) {}


void AutoRebin::Rebin(CombineHarvester &src, CombineHarvester &dest) {

  auto bins = src.bin_set();
  //Perform rebinning separately for each channel-category of the CH instance
    
  for (auto bin : bins) {
    //Build histogram containing total of all backgrounds, avoid using
    //GetShapeWithUncertainty so that possible negative bins are retained
    TH1F data_obs = src.cp().bin({bin}).GetObservedShape();
    TH1F total_bkg;
    if (data_obs.GetXaxis()->GetXbins()->GetArray()){
      total_bkg = TH1F("","",data_obs.GetNbinsX(), 
                  data_obs.GetXaxis()->GetXbins()->GetArray()) ;
    } else {
      total_bkg = TH1F("","",data_obs.GetNbinsX(),
                  data_obs.GetXaxis()->GetBinLowEdge(1),data_obs.GetXaxis()->GetBinLowEdge(data_obs.GetNbinsX()+1));
    }
    src.cp().bin({bin}).backgrounds().ForEachProc([&](ch::Process *proc) {
        total_bkg.Add((proc->ClonedScaledShape()).get());
    });
      
    //Create a vector to store the original and new binning
    int nbins = total_bkg.GetNbinsX();
    std::vector<double> init_bins; 
    for(int i=1; i<=nbins+1; ++i) init_bins.push_back(total_bkg.GetBinLowEdge(i));
    std::cout << "[AutoRebin] Searching for bins failing conditions "
    "for analysis bin id " << bin << std::endl; 
    std::vector<double> new_bins = init_bins;
        
    //Call to function to fill new_bins vector with the recommended binning
    AutoRebin::FindNewBinning(total_bkg, new_bins, bin_threshold_, 
                bin_uncert_fraction_, rebin_mode_);
        
    //Inform user if binning has been modified 
    if(new_bins.size() != init_bins.size()) { 
      std::cout << "[AutoRebin] Some bins not satisfying requested condition "
      "found for analysis bin " << bin << ", merging with neighbouring bins"
      << std::endl; 
      if(v_ > 0) {
        unsigned i_old = 0;
        unsigned i_new = 0;
        std::cout << (
          boost::format("%-21s %-10s %-21s\n") % "Init Edges/Widths" % "Content" % "New Edges/Widths");

        for (; i_old < init_bins.size(); ++i_old) {
          double i_old_width = (i_old < (init_bins.size() - 1))
                                   ? (init_bins[i_old + 1] - init_bins[i_old])
                                   : 0.;
          std::cout << (
            boost::format("%-10.0f %-10.0f %-10.2g") % init_bins[i_old] % i_old_width % total_bkg.GetBinContent(i_old+1));
          bool new_aligned = (i_new < new_bins.size()) ? std::fabs(init_bins[i_old] - new_bins[i_new]) < 1E-8 : false;
          if (new_aligned) {
            double i_new_width = (i_new < (new_bins.size() - 1))
                                     ? (new_bins[i_new + 1] - new_bins[i_new])
                                     : 0.;

            std::cout << (
              boost::format("%-10.0f %-10.0f\n") % new_bins[i_new] % i_new_width);
            ++i_new;
          } else {
            std::cout << (
              boost::format("%-10s %-10s\n") % "-" % "-");
          }
        }
        // std::cout << "Original binning: " << std::endl; 
        // for (std::vector<double>::const_iterator i =
        // init_bins.begin(); i != init_bins.end(); ++i) std::cout << *i << ", ";
        // std::cout << std::endl; 
        // std::cout << "New binning: " << std::endl; 
        // for (std::vector<double>::const_iterator i = new_bins.begin(); 
        //         i != new_bins.end(); ++i) 
        //     std::cout << *i << ", "; 
        // std::cout << std::endl; 
      }
      //Altering binning in CH instance for all distributions if requested
      if(perform_rebin_) {
        std::cout << "[AutoRebin] Applying binning to all relevant distributions "
          "for analysis bin id " << bin << std::endl; 
      dest.cp().bin({bin}).VariableRebin(new_bins); 
      }
    } else std::cout << "[AutoRebin] Did not find any bins to merge for analysis "
        "bin id: " << bin << std::endl;
  }
}


void AutoRebin::FindNewBinning(TH1F &total_bkg, std::vector<double> &new_bins, 
                double bin_condition, double bin_uncert_fraction, int mode) {

  bool all_bins = true;
  //Find the maximum bin
  int hbin_idx = total_bkg.GetMaximumBin();
  int nbins = total_bkg.GetNbinsX();
  std::cout << "[AutoRebin::FindNewBinning] Searching for bins failing condition "
  "using algorithm mode: " << mode << std::endl; 
   
  //Mode 0 is the simplest version of merging algorithm. It loops from the
  //peak bin left to bin 1 then right to final bin. When the "next" bin in
  //either direction fails
  //the threshold, it removes the boundary between the current bin and the
  //next bin. This will be overly conservative in situations where a better
  //solution is to sum a subset of subsequent bins to reach the threshold and
  //will lead to more rebinning than necessary.
  if(mode==0) {
    new_bins.clear();
    //Loop from highest bin, first left and then right, merging bins where
    //necessary to make new list of BinLowEdges
    if(v_>0) std::cout << "[AutoRebin::FindNewBinning] Testing bins of " 
       "total_bkg hist to find those failing condition, starting from "
       "maximum bin: " << hbin_idx << std::endl; 
    for(int idx=hbin_idx; idx>1; --idx) { 
      if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
          total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
          total_bkg.GetBinContent(idx) << ", Bin error fraction: " <<
          total_bkg.GetBinError(idx-1)/total_bkg.GetBinContent(idx-1)
          << std::endl;
      if(total_bkg.GetBinContent(idx-1) > bin_condition 
          && (total_bkg.GetBinError(idx-1)/total_bkg.GetBinContent(idx-1)
          < bin_uncert_fraction) 
          && idx!=1)
          new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
    } 
    for(int idx=hbin_idx+1; idx<=nbins; ++idx) { 
      if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
          total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
          total_bkg.GetBinContent(idx) << ", Bin error fraction: " <<
          total_bkg.GetBinError(idx)/total_bkg.GetBinContent(idx)
          << std::endl;
      if(total_bkg.GetBinContent(idx) > bin_condition
          && (total_bkg.GetBinError(idx)/total_bkg.GetBinContent(idx)
          < bin_uncert_fraction))
          new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
    }
    //Add to vector the extremal bins, which are untreated in the above 
    new_bins.push_back(total_bkg.GetBinLowEdge(1));
    new_bins.push_back(total_bkg.GetBinLowEdge(nbins+1));
    std::sort(new_bins.begin(), new_bins.end()); 
        
  //Mode 1 finds bins below threshold, then tries merging with bins left or
  //right until threshold is met, the final choice being the one which
  //minimises the number of required bin merges.     
  } else if(mode == 1 || mode == 2) {
    //Start from finding the minimum bin and the bin with largest fractional
    //error
    int lbin_idx = 0;
    if (mode == 1) {
      lbin_idx = total_bkg.GetMinimumBin();
    } else if (mode == 2) {
      double maxval = std::numeric_limits<double>::max();
      // loop through bins from left to right - we take the rightmost
      // bin in the case multiple bins are tied on the smallest value
      for (int idx = 1; idx <= total_bkg.GetNbinsX(); ++idx) {
        if (total_bkg.GetBinContent(idx) <= maxval) {
          lbin_idx = idx;
          maxval = total_bkg.GetBinContent(idx);
        }
      }
    }

    int herrbin_idx = GetMaximumFracUncertBin(total_bkg);
    bool bin_tot_flag = total_bkg.GetBinContent(lbin_idx) <= bin_condition;
    bool bin_err_flag = total_bkg.GetBinError(herrbin_idx)/
            total_bkg.GetBinContent(herrbin_idx) >= bin_uncert_fraction;
    if(bin_tot_flag || bin_err_flag) {
      //Fix issues with minimal bins in first pass and bin errors after 
      lbin_idx = (bin_tot_flag ? lbin_idx : herrbin_idx);
      if(v_>0 && bin_tot_flag) std::cout << "[AutoRebin::FindNewBinning] Testing "
          "bins of total_bkg hist to find those with entry <= " << bin_condition 
          << ", starting from bin: " << lbin_idx << std::endl; 
      if(v_>0 && !bin_tot_flag) std::cout << "[AutoRebin::FindNewBinning] Testing "
          "bins of total_bkg hist to find those with fractional error >= " <<
          bin_uncert_fraction << ", starting from bin: " << lbin_idx << std::endl; 
      //loop left first
      double left_tot = total_bkg.GetBinContent(lbin_idx);
      double left_err_tot = total_bkg.GetBinError(lbin_idx);
      int idx = lbin_idx - 1;
      //Note: this logic works so long as bins with 0 entries are merged
      //first, to give sensible values for the fractional error
      while (((bin_tot_flag && left_tot <= bin_condition) || (!bin_tot_flag && 
          left_err_tot/left_tot >= bin_uncert_fraction)) && idx > 0) {
          left_err_tot = sqrt(pow(left_err_tot,2) + pow(total_bkg.GetBinError(idx),2));
          left_tot = left_tot + total_bkg.GetBinContent(idx);
          if(v_>0) std::cout << "Moving left, bin index: " << idx << ", BinLowEdge: " << 
          total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
          total_bkg.GetBinContent(idx) << ", Running total from combining bins: "
          << left_tot << ", Current fractional error: " << left_err_tot/left_tot 
          << std::endl;
          --idx;
      }
      int left_bins = lbin_idx - idx;
      double right_tot = total_bkg.GetBinContent(lbin_idx);
      double right_err_tot = total_bkg.GetBinError(lbin_idx);
      idx = lbin_idx + 1;
      //now try right
      while (((bin_tot_flag && right_tot <= bin_condition) || (!bin_tot_flag && 
          right_err_tot/right_tot >= bin_uncert_fraction)) && idx <= nbins) {
          right_err_tot = sqrt(pow(right_err_tot,2) + pow(total_bkg.GetBinError(idx),2));
          right_tot = right_tot + total_bkg.GetBinContent(idx);
          if(v_>0) std::cout << "Moving right, bin index: " << idx << ", BinLowEdge: " << 
          total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
          total_bkg.GetBinContent(idx) << ", Running total from combining bins: "
          << right_tot << ", Current fractional error: " << right_err_tot/right_tot 
          << std::endl;
          ++idx;
      }
      int right_bins = idx - lbin_idx;
      bool left_pass = bin_tot_flag ? left_tot > bin_condition : 
            left_err_tot/left_tot < bin_uncert_fraction;
      bool right_pass = bin_tot_flag ? right_tot > bin_condition : 
            right_err_tot/right_tot < bin_uncert_fraction;
      if (mode == 2) right_pass = false;  // never merge right in mode 2
      //Decision on which way to merge - remove relevant entries from
      //binning vector as applicable
      if(left_pass && !right_pass) {
          if(v_>0) std::cout << "Merging left using " << left_bins - 1 << 
              " bins" << std::endl;
          for(int i = lbin_idx; i > lbin_idx - left_bins + 1; --i) {
              if(v_>0) std::cout << "Erasing bin low edge for bin " << i << 
                  ", BinLowEdge: " << total_bkg.GetBinLowEdge(i) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i)),new_bins.end() );
          }
      } else if(right_pass && !left_pass) {
          if(v_>0) std::cout << "Merging right using " << right_bins - 1 << 
              " bins" << std::endl;   
          for(int i = lbin_idx; i < lbin_idx + right_bins - 1; ++i) {
              if(v_>0) std::cout << "Erasing bin upper edge for bin " << i << 
                  ", BinUpperEdge: " << total_bkg.GetBinLowEdge(i+1) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i+1)),new_bins.end() );
          }
      } else if(left_pass && right_pass && left_bins < right_bins) {
          if(v_>0) std::cout << "Merging left using " << left_bins - 1 << 
              " bins" << std::endl;
          for(int i = lbin_idx; i > lbin_idx - left_bins + 1; --i) {
              if(v_>0) std::cout << "Erasing bin low edge for bin " << i << 
                  ", BinLowEdge: " << total_bkg.GetBinLowEdge(i) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i)),new_bins.end() );
          }
      } else if(left_pass && right_pass && left_bins > right_bins) {
          if(v_>0) std::cout << "Merging right using " <<  right_bins - 1 << 
              " bins" << std::endl;
          for(int i = lbin_idx; i < lbin_idx + right_bins - 1 ; ++i) {
              if(v_>0) std::cout << "Erasing bin upper edge for bin " << i << 
                  ", BinUpperEdge: " << total_bkg.GetBinLowEdge(i+1) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i+1)),new_bins.end() );
          }
      //In case where left and right are equally good, merge towards the
      //peak (unclear to me if this is really always better but it generally more
      //conservative in leading to bins with higher content)
      } else if(left_pass && right_pass && left_bins == right_bins &&
          lbin_idx < hbin_idx) {
          if(v_>0) std::cout << "Merging right using " <<  right_bins - 1 << 
              " bins" << std::endl;
          for(int i = lbin_idx; i < lbin_idx + right_bins - 1; ++i) {
              if(v_>0) std::cout << "Erasing bin upper edge for bin " << i << 
                  ", BinUpperEdge: " << total_bkg.GetBinLowEdge(i+1) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i+1)),new_bins.end() );
          }
      } else if(left_pass && right_pass && left_bins == right_bins && 
          lbin_idx > hbin_idx) {
          if(v_>0) std::cout << "Merging left using " << left_bins - 1 << 
              " bins" << std::endl;
          for(int i = lbin_idx; i > lbin_idx - left_bins + 1; --i) {
              if(v_>0) std::cout << "Erasing bin low edge for bin " << i << 
                  ", BinLowEdge: " << total_bkg.GetBinLowEdge(i) << std::endl;
              new_bins.erase( std::remove( new_bins.begin(), new_bins.end(), 
                total_bkg.GetBinLowEdge(i)),new_bins.end() );
          }
      } else if(!left_pass && !right_pass) {
          std::cout << "[AutoRebin::FindNewBinning] WARNING: No solution found "
          "to satisfy condition, try merging all bins" << std::endl;
          for(int i = 2; i<=nbins; ++i) new_bins.erase( std::remove( new_bins.begin(), 
              new_bins.end(), total_bkg.GetBinLowEdge(i)),new_bins.end() );
      }
    }
  } else {
      std::cout << "[AutoRebin::FindNewBinning] Chosen mode " << mode << " not"
      " currently supported, exiting without altering binning" << std::endl;
      return;
  }
    
  TH1F* total_bkg_new;
  total_bkg_new = (TH1F*)total_bkg.Rebin(new_bins.size()-1, "total_bkg_new",
    &new_bins[0]); 
    
  //Flag checks if all bins now satify the condition after pass through
  int nbins_new = total_bkg_new->GetNbinsX();
  if(nbins_new!=nbins) {
    if(v_>0) std::cout << std::endl;
    if(v_>0) std::cout << "[AutoRebin::FindNewBinning] New binning found: " 
        << std::endl;
    for(int i=1; i<=nbins_new+1; ++i) {
      if(v_>0)  std::cout << "Bin index: " << i << ", BinLowEdge: " << 
            total_bkg_new->GetBinLowEdge(i) <<  ", Bin content: " << 
            total_bkg_new->GetBinContent(i) << ", Bin error fraction: " <<
            total_bkg_new->GetBinError(i)/total_bkg_new->GetBinContent(i) 
            << std::endl;
       //Last bin is allowed to be failing condition, just there to enclose the
       //penultimate bin
       if((total_bkg_new->GetBinContent(i) <= bin_condition 
           || total_bkg_new->GetBinError(i)/total_bkg_new->GetBinContent(i) 
           >= bin_uncert_fraction)
           && i!=nbins_new+1) all_bins=false;
    }
  } 
     
  if(all_bins) return;
  //Run function again if all_bins is not true, i.e. if it's possible that not
  //all bins pass condition after one pass 
  else FindNewBinning(*total_bkg_new, new_bins, bin_condition, bin_uncert_fraction, mode);

}


int AutoRebin::GetMaximumFracUncertBin(TH1F &total_bkg) {
   int idx = 1;
   double maximum = 0;
   for(int i = 1; i <= total_bkg.GetNbinsX(); ++i){
       if(total_bkg.GetBinError(i)/total_bkg.GetBinContent(i) > maximum) {
           maximum = total_bkg.GetBinError(i)/total_bkg.GetBinContent(i);
           idx = i;
       }
   }

   return idx; 
}




}
