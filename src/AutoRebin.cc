#include "CombineHarvester/CombineTools/interface/AutoRebin.h" 
#include <iostream> 
#include <string> 
#include <vector> 
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

namespace ch {

AutoRebin::AutoRebin()
    : v_(0), 
      empty_bin_threshold_(0.) {}


void AutoRebin::Rebin(CombineHarvester &src, CombineHarvester &dest) {

    auto bins = src.bin_set();
    //Perform rebinning separately for each channel-category of the CH instance
    for (auto bin : bins) {
        //Build histogram containing total of all backgrounds, avoid using
        //GetShapeWithUncertainty so that possible negative bins are retained
        //NB need to replace this with data_obs when figure out how
        TH1F data_obs = src.cp().bin({bin}).GetObservedShape();
        TH1F total_bkg("","",data_obs.GetNbinsX(), data_obs.GetXaxis()->GetXbins()->GetArray()) ;
        src.cp().bin({bin}).backgrounds().ForEachProc([&](ch::Process *proc) {
            total_bkg.Add((proc->ClonedScaledShape()).get());
        });
      
        //Find the maximum bin
        int hbin_idx = total_bkg.GetMaximumBin(); 
        int nbins = total_bkg.GetNbinsX();
        //Create a vector to store the original binning
        std::vector<double> init_bins; 
        for(int i=1; i<=nbins+1; ++i) init_bins.push_back(total_bkg.GetBinLowEdge(i));
        //Loop from highest bin, first left and then right, merging bins where
        //necessary to make new list of BinLowEdges
        std::cout << "[AutoRebin] Searching for bins failing condition starting from "
        "maximum bin at position " << hbin_idx << ", for analysis bin id " << bin << std::endl; 
        std::vector<double> new_bins; 
        
        for(int idx=hbin_idx; idx>1; --idx) { 
            if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
                total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
                total_bkg.GetBinContent(idx) << std::endl;
            if(total_bkg.GetBinContent(idx-1) > empty_bin_threshold_ && idx!=1)
                new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
        } 
        for(int idx=hbin_idx+1; idx<=nbins; ++idx) { 
            if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
                total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
                total_bkg.GetBinContent(idx) << std::endl;
            if(total_bkg.GetBinContent(idx) > empty_bin_threshold_)
                new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
        }
        
        //Retain enclosing bin edges in the new list 
        new_bins.push_back(total_bkg.GetBinLowEdge(1));
        new_bins.push_back(total_bkg.GetBinLowEdge(nbins+1));
        std::sort(new_bins.begin(), new_bins.end()); 
        if(new_bins.size() != init_bins.size()) { 
            std::cout << "[AutoRebin] Some bins not satisfying requested condition "
            "found for analysis bin " << bin << ", merging with neighbouring bins" << std::endl; 
            if(v_>0) { 
                std::cout << "Original binning: " << std::endl; 
                for (std::vector<double>::const_iterator i =
                init_bins.begin(); i != init_bins.end(); ++i) 
                    std::cout << *i << ", ";
                std::cout << std::endl; 
                std::cout << "New binning: " << std::endl; 
                for (std::vector<double>::const_iterator i = new_bins.begin(); i != new_bins.end(); ++i) 
                    std::cout << *i << ", "; 
                    std::cout << std::endl; 
            }
        } else std::cout << "[AutoRebin] Did not find any bins to merge for analysis "
        "bin id: " << bin << std::endl;
        
        dest.cp().bin({bin}).VariableRebin(new_bins); 
        
        /*TH1F data_obs_new = dest.cp().bin({bin}).GetObservedShape(); 
        TH1F total_bkg_new("","",data_obs_new.GetNbinsX(), data_obs_new.GetXaxis()->GetXbins()->GetArray()) ;
        dest.cp().bin({bin}).backgrounds().ForEachProc([&](ch::Process *proc) {
            total_bkg_new.Add((proc->ClonedScaledShape()).get());
        });
        int nbins_new = total_bkg_new.GetNbinsX();
        for(int i=1; i<=nbins_new+1; ++i) {
           std::cout << "Bin index: " << i << ", BinLowEdge: " << 
                total_bkg_new.GetBinLowEdge(i) <<  ", Bin content: " << 
                total_bkg_new.GetBinContent(i) << std::endl;
        }*/
    
    }
   

}


}
