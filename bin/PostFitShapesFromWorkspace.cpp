#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "CombineTools/interface/CombineHarvester.h"
#include "CombineTools/interface/TFileIO.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]) {
  /*
  In this example we will build a fully-functional program to serve a few common
  post-fit plotting tasks in one go:

   - Create pre- and post-fit histograms for each process in each bin
   - Create pre- and post-fit histograms for the total signal and total
  background contributions. The latter can be used to draw the "Bkg.
  uncertainty" band
   - Calculate a post-fit/pre-fit scale factor for each background to be used as
  input for making post-fit control plots of other distributions

  We will create the user configuration using the boost program_options library.
  We first declare the variables we need, then connect them to named
  command-line options. As well as providing the essential information like the
  input datacard and fit result, we also let the user choose whether or not to
  do the post-fit calculation, whether to use the covariance matrix sampling
  method (to get the effect of nuisance parameter correlations in the
  uncertainty) and whether or not to calculate and print the post/pre scale
  factors.
  */
  string datacard   = "";
  string fitresult  = "";
  string mass       = "";
  bool postfit      = false;
  bool sampling     = false;
  string output     = "";
  bool factors      = false;
  unsigned samples  = 500;

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  config.add_options()
    ("datacard,d",
      po::value<string>(&datacard)->required(),
      "The datacard .txt file [REQUIRED]")
    ("output,o ",
      po::value<string>(&output)->required(),
      "Name of the output root file to create [REQUIRED]")
    ("fitresult,f",
      po::value<string>(&fitresult)->default_value(fitresult),
      "Path to a RooFitResult, only needed for postfit")
    ("mass,m",
      po::value<string>(&mass)->default_value(""),
      "Signal mass point of the input datacard")
    ("postfit",
      po::value<bool>(&postfit)
      ->default_value(postfit)->implicit_value(true),
      "Create post-fit histograms in addition to pre-fit")
    ("sampling",
      po::value<bool>(&sampling)->default_value(sampling)->implicit_value(true),
      "Use the cov. matrix sampling method for the post-fit uncertainty")
    ("samples",
      po::value<unsigned>(&samples)->default_value(samples),
      "Number of samples to make in each evaluate call")
    ("print",
      po::value<bool>(&factors)->default_value(factors)->implicit_value(true),
      "Print tables of background shifts and relative uncertainties");

  if (sampling && !postfit) {
    throw logic_error(
        "Can't sample the fit covariance matrix for pre-fit!");
  }

  po::variables_map vm;

  // First check if the user has set the "--help" or "-h" option, and if so
  // just prin the usage information and quit
  po::store(po::command_line_parser(argc, argv)
    .options(help_config).allow_unregistered().run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    cout << config << "\nExample usage:\n";
    cout << "PostFitShapes -d htt_mt_125.txt -o htt_mt_125_shapes.root -m 125 "
            "-f mlfit.root:fit_s --postfit --sampling --print\n";
    return 1;
  }

  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Create CH instance and parse the datacard - the only metadata we care
  // about here is the mass
  ch::CombineHarvester cmb;
  cmb.ParseDatacard(datacard, "", "", "", 0, mass);

  // Drop any process that has no hist/data/pdf
  cmb.FilterProcs([&](ch::Process * proc) {
    bool no_shape = !proc->shape() && !proc->data() && !proc->pdf();
    if (no_shape) {
      cout << "Filtering process with no shape:\n";
      cout << ch::Process::PrintHeader << *proc << "\n";
    }
    return no_shape;
  });

  auto bins = cmb.bin_set();

  TFile outfile(output.c_str(), "RECREATE");
  TH1::AddDirectory(false);

  // Create a map of maps for storing histograms in the form:
  //   pre_shapes[<bin>][<process>]
  map<string, map<string, TH1F>> pre_shapes;
  // We can always do the prefit version,
  // Loop through the bins writing the shapes to the output file
  for (auto bin : bins) {
    ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});
    // This next line is a temporary fix for models with parameteric RooFit pdfs
    // - we try and set the number of bins to evaluate the pdf to be the same as
    // the number of bins in data
    cmb_bin.SetPdfBins(cmb_bin.GetObservedShape().GetNbinsX());

    // Fill the data and process histograms
    pre_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
    for (auto proc : cmb_bin.process_set()) {
      pre_shapes[bin][proc] =
          cmb_bin.cp().process({proc}).GetShapeWithUncertainty();
    }

    // The fill total signal and total bkg hists
    pre_shapes[bin]["TotalBkg"] =
        cmb_bin.cp().backgrounds().GetShapeWithUncertainty();
    pre_shapes[bin]["TotalSig"] =
        cmb_bin.cp().signals().GetShapeWithUncertainty();

    // Can write these straight into the output file
    outfile.cd();
    for (auto& iter : pre_shapes[bin]) {
      ch::WriteToTFile(&(iter.second), &outfile, bin + "_prefit/" + iter.first);
    }
  }

  // Print out the relative uncert. on the bkg
  if (factors) {
    cout << boost::format("%-25s %-32s\n") % "Bin" %
                "Total relative bkg uncert. (prefit)";
    cout << string(58, '-') << "\n";
    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});
      double rate = cmb_bin.cp().backgrounds().GetRate();
      double err = cmb_bin.cp().backgrounds().GetUncertainty();
      cout << boost::format("%-25s %-10.5f\n") % bin %
                  (rate > 0. ? (err / rate) : 0.);
    }
  }

  // Now we can do the same again but for the post-fit model
  if (postfit) {
    // Get the fit result and update the parameters to the post-fit model
    RooFitResult res = ch::OpenFromTFile<RooFitResult>(fitresult);
    cmb.UpdateParameters(res);

    // Calculate the post-fit fractional background uncertainty in each bin
    if (factors) {
      cout << boost::format("\n%-25s %-32s\n") % "Bin" %
                  "Total relative bkg uncert. (postfit)";
      cout << string(58, '-') << "\n";
      for (auto bin : bins) {
        ch::CombineHarvester cmb_bkgs = cmb.cp().bin({bin}).backgrounds();
        double rate = cmb_bkgs.GetRate();
        double err = sampling ? cmb_bkgs.GetUncertainty(res, samples)
                              : cmb_bkgs.GetUncertainty();
        cout << boost::format("%-25s %-10.5f\n") % bin %
                    (rate > 0. ? (err / rate) : 0.);
      }
    }

    map<string, map<string, TH1F>> post_shapes;
    // As we calculate the post-fit yields can also print out the post/pre scale
    // factors
    if (factors) {
      cout << boost::format("\n%-25s %-20s %-10s\n") % "Bin" % "Process" %
                       "Scale factor";
      cout << string(58, '-') << "\n";
    }

    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});
      post_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
      for (auto proc : cmb_bin.process_set()) {
        auto cmb_proc = cmb_bin.cp().process({proc});
        // Method to get the shape uncertainty depends on whether we are using
        // the sampling method or the "wrong" method (assumes no correlations)
        post_shapes[bin][proc] =
            sampling ? cmb_proc.GetShapeWithUncertainty(res, samples)
                     : cmb_proc.GetShapeWithUncertainty();
        // Print out the post/pre scale factors
        if (factors) {
          TH1 const& pre = pre_shapes[bin][proc];
          TH1 const& post = post_shapes[bin][proc];
          cout << boost::format("%-25s %-20s %-10.5f\n") % bin % proc %
                      (pre.Integral() > 0. ? (post.Integral() / pre.Integral())
                                           : 1.0);
        }
      }
      // Fill the total sig. and total bkg. hists
      auto cmb_bkgs = cmb_bin.cp().backgrounds();
      auto cmb_sigs = cmb_bin.cp().signals();
      post_shapes[bin]["TotalBkg"] =
          sampling ? cmb_bkgs.GetShapeWithUncertainty(res, samples)
                   : cmb_bkgs.GetShapeWithUncertainty();
      post_shapes[bin]["TotalSig"] =
          sampling ? cmb_sigs.GetShapeWithUncertainty(res, samples)
                   : cmb_sigs.GetShapeWithUncertainty();
      outfile.cd();
      // Write the post-fit histograms
      for (auto & iter : post_shapes[bin]) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         bin + "_postfit/" + iter.first);
      }
    }
  }
  // And we're done!
  outfile.Close();
  return 0;
}

