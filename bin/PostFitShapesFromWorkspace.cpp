#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "CombineTools/interface/CombineHarvester.h"
#include "CombineTools/interface/TFileIO.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]){
  string datacard   = "";
  string fitresult  = "";
  string mass       = "";
  bool postfit      = false;
  bool sampling     = false;
  string output     = "";
  bool factors      = false;

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");
  po::options_description config("Configuration");
  config.add_options()
    ("datacard,d",      po::value<string>(&datacard)->required(),
        "The datacard .txt file [REQUIRED]")
    ("output,o ",       po::value<string>(&output)->required(),
        "Name of the output root file to create [REQUIRED]")
    ("fitresult,f",     po::value<string>(&fitresult)->default_value(fitresult),
        "Path to a RooFitResult, only needed for postfit")
    ("mass,m",          po::value<string>(&mass)->default_value(""),
        "Signal mass point of the input datacard")
    ("postfit",  po::value<bool>(&postfit)->default_value(postfit)->implicit_value(true),
        "Create post-fit histograms in addition to pre-fit")
    ("sampling", po::value<bool>(&sampling)->default_value(sampling)->implicit_value(true),
        "Use the cov. matrix sampling method for the post-fit uncertainty")
    ("print",    po::value<bool>(&factors)->default_value(factors)->implicit_value(true),
        "Print tables of background shifts and relative uncertainties");

  if (sampling && !postfit)
    throw std::logic_error(
        "Can't sample the fit covariance matrix for pre-fit!");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
                .options(help_config)
                .allow_unregistered()
                .run(),
            vm);
  po::notify(vm);
  if (vm.count("help")) {
    cout << config << "\n";
    cout << "Example usage: " << endl;
    cout << "PostFitShapes -d htt_mt_125.txt -o htt_mt_125_shapes.root -m 125 "
            "-f mlfit.root:fit_s --postfit --sampling --print \n";
    return 1;
  }
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  ch::CombineHarvester cmb;
  cmb.ParseDatacard(datacard, "", "", "", 0, mass);

  RooFitResult *res = nullptr;
  auto bins = cmb.bin_set();

  TFile outfile(output.c_str(), "RECREATE");
  TH1::AddDirectory(false);

  map<string, map<string, TH1F>> pre_shapes;
  // We can always do the prefit version,
  // Loop through the bins writing the shapes to the output file
  for (auto bin : bins) {
    ch::CombineHarvester cmb_bin = move(cmb.cp().bin({bin}));
    pre_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
    auto procs = cmb.process_set();
    for (auto proc : procs) {
      pre_shapes[bin][proc] =
          cmb_bin.cp().process({proc}).GetShapeWithUncertainty();
    }
    pre_shapes[bin]["TotalBkg"] =
        cmb_bin.cp().backgrounds().GetShapeWithUncertainty();
    pre_shapes[bin]["TotalSig"] =
        cmb_bin.cp().signals().GetShapeWithUncertainty();
    outfile.cd();
    for (auto& iter : pre_shapes[bin]) {
      ch::WriteToTFile(&(iter.second), &outfile, bin + "_prefit/" + iter.first);
    }
  }

  if (factors) {
    std::cout << boost::format("%-25s %-32s\n") % "Bin" % "Total relative uncert. (prefit)";
    std::cout << string(58, '-') << "\n";
  }
  for (auto bin : bins) {
    ch::CombineHarvester cmb_bin = move(cmb.cp().bin({bin}));
    double rate = cmb_bin.cp().backgrounds().GetRate();
    double err  = cmb_bin.cp().backgrounds().GetUncertainty();
    if (factors)
      std::cout << boost::format("%-25s %-10.5f\n") % bin % (err / rate);
  }

  // Here we'll do some horrendous copy-n-paste coding to do the same again
  // for the postfit
  if (postfit) {
    res = new RooFitResult(ch::OpenFromTFile<RooFitResult>(fitresult));
    cmb.UpdateParameters(ch::ExtractFitParameters(*res));

    if (factors) {
      std::cout << boost::format("\n%-25s %-32s\n") % "Bin" % "Total relative uncert. (postfit)";
      std::cout << string(58, '-') << "\n";
    }
    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = move(cmb.cp().bin({bin}));
      double rate = cmb_bin.cp().backgrounds().GetRate();
      double err  = sampling ? cmb_bin.cp().backgrounds().GetUncertainty(res, 500)
                             : cmb_bin.cp().backgrounds().GetUncertainty();
      if (factors)
        std::cout << boost::format("%-25s %-10.5f\n") % bin % (err / rate);
    }

    map<string, map<string, TH1F>> post_shapes;
    if (factors) {
      std::cout << boost::format("\n%-25s %-20s %-10s\n") % "Bin" % "Process" %
                       "Scale factor";
      std::cout << string(58, '-') << "\n";
    }

    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = move(cmb.cp().bin({bin}));
      post_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
      auto procs = cmb.process_set();
      for (auto proc : procs) {
        post_shapes[bin][proc] =
            sampling
                ? cmb_bin.cp().process({proc}).GetShapeWithUncertainty(res, 500)
                : cmb_bin.cp().process({proc}).GetShapeWithUncertainty();
        if (factors) std::cout << boost::format("%-25s %-20s %-10.5f\n") % bin % proc %
                         (post_shapes[bin][proc].Integral() /
                          pre_shapes[bin][proc].Integral());
      }
      post_shapes[bin]["TotalBkg"] =
          sampling ? cmb_bin.cp().backgrounds().GetShapeWithUncertainty(res, 500)
                   : cmb_bin.cp().backgrounds().GetShapeWithUncertainty();
      post_shapes[bin]["TotalSig"] =
          sampling ? cmb_bin.cp().signals().GetShapeWithUncertainty(res, 500)
                   : cmb_bin.cp().signals().GetShapeWithUncertainty();
      outfile.cd();
      for (auto & iter : post_shapes[bin]) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         bin + "_postfit/" + iter.first);
      }
    }
  }
  outfile.Close();
  return 0;
}

