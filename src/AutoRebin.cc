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
        TH1F data_obs = src.cp().bin({bin}).GetObservedShape();
        TH1F total_bkg("","",data_obs.GetNbinsX(), data_obs.GetXaxis()->GetXbins()->GetArray()) ;
        src.cp().bin({bin}).backgrounds().ForEachProc([&](ch::Process *proc) {
            total_bkg.Add((proc->ClonedScaledShape()).get());
        });
      
        //Create a vector to store the original and new binning
        int nbins = total_bkg.GetNbinsX();
        std::vector<double> init_bins; 
        for(int i=1; i<=nbins+1; ++i) init_bins.push_back(total_bkg.GetBinLowEdge(i));
        std::cout << "[AutoRebin] Searching for bins failing condition "
        "for analysis bin id " << bin << std::endl; 
        std::vector<double> new_bins;
        
        //Call to function to fill new_bins with the updated binning
        AutoRebin::FindNewBinning(total_bkg, new_bins, empty_bin_threshold_);
        
        //Inform user if binning has been modified 
        if(new_bins.size() != init_bins.size()) { 
            std::cout << "[AutoRebin] Some bins not satisfying requested condition "
            "found for analysis bin " << bin << ", merging with neighbouring bins" << std::endl; 
            if(v_>0) { 
                std::cout << "Original binning: " << std::endl; 
                for (std::vector<double>::const_iterator i =
                init_bins.begin(); i != init_bins.end(); ++i) std::cout << *i << ", ";
                std::cout << std::endl; 
                std::cout << "New binning: " << std::endl; 
                for (std::vector<double>::const_iterator i = new_bins.begin(); i != new_bins.end(); ++i) 
                    std::cout << *i << ", "; 
                std::cout << std::endl; 
            }
        } else std::cout << "[AutoRebin] Did not find any bins to merge for analysis "
        "bin id: " << bin << std::endl;
        
        //Altering binning in CH instance for all distributions
        dest.cp().bin({bin}).VariableRebin(new_bins); 
        
    
    }
}


void AutoRebin::FindNewBinning(TH1F &total_bkg, std::vector<double> &new_bins, double bin_condition) {
    
    int nbins = total_bkg.GetNbinsX();
    new_bins.clear();
    //Find the maximum bin
    int hbin_idx = total_bkg.GetMaximumBin(); 
    std::cout << "[AutoRebin::FindNewBinning] Searching for bins failing condition starting from "
    "maximum bin at position " << hbin_idx << std::endl; 
    //Loop from highest bin, first left and then right, merging bins where
    //necessary to make new list of BinLowEdges
    if(v_>0) std::cout << "[AutoRebin::FindNewBinning] Testing bins of total_bkg hist to find"
        " those failing condition: " << std::endl; 
    for(int idx=hbin_idx; idx>1; --idx) { 
        if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
            total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
            total_bkg.GetBinContent(idx) << std::endl;
        if(total_bkg.GetBinContent(idx-1) > bin_condition && idx!=1)
            new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
    } 
    for(int idx=hbin_idx+1; idx<=nbins; ++idx) { 
        if(v_>0) std::cout << "Bin index: " << idx << ", BinLowEdge: " << 
            total_bkg.GetBinLowEdge(idx) <<  ", Bin content: " << 
            total_bkg.GetBinContent(idx) << std::endl;
        if(total_bkg.GetBinContent(idx) > bin_condition)
            new_bins.push_back(total_bkg.GetBinLowEdge(idx)); 
    }
    //Add to vector the extremal bins, which are untreated in the above 
    new_bins.push_back(total_bkg.GetBinLowEdge(1));
    new_bins.push_back(total_bkg.GetBinLowEdge(nbins+1));
    std::sort(new_bins.begin(), new_bins.end()); 
    
    TH1F* total_bkg_new = (TH1F*)total_bkg.Rebin(new_bins.size()-1, "total_bkg_new", &new_bins[0]); 
    
    //Flag checks if all bins now satify the condition. This would be true by
    //construction for a condition of bin content > 0 if the input histogram 
    //had no negative bins, but to support negative bins and possibility 
    //of a different choice of bin threshold we make the function recursive.
    bool all_bins = true;
    int nbins_new = total_bkg_new->GetNbinsX();
    if(v_>0) std::cout << std::endl;
    if(v_>0) std::cout << "[AutoRebin::FindNewBinning] Producing new histogram with following binning: " << std::endl;
    for(int i=1; i<=nbins_new+1; ++i) {
    if(v_>0)  std::cout << "Bin index: " << i << ", BinLowEdge: " << 
            total_bkg_new->GetBinLowEdge(i) <<  ", Bin content: " << 
            total_bkg_new->GetBinContent(i) << std::endl;
       //Last bin is allowed to be failing condition, just there to enclose the
       //penultimate bin
       if(total_bkg_new->GetBinContent(i) <= bin_condition && i!=nbins_new+1) all_bins=false;
    }
     
    if(all_bins) return;
    else FindNewBinning(*total_bkg_new, new_bins, bin_condition);

}






}
