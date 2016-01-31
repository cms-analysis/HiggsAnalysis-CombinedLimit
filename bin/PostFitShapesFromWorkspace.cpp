#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]) {
  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  string datacard   = "";
  string workspace  = "";
  string fitresult  = "";
  string mass       = "";
  bool postfit      = false;
  bool sampling     = false;
  string output     = "";
  bool factors      = false;
  unsigned samples  = 500;
  std::string freeze_arg = "";

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  config.add_options()
    ("workspace,w",
      po::value<string>(&workspace)->required(),
      "The input workspace-containing file [REQUIRED]")
    ("datacard,d",
      po::value<string>(&datacard),
      "The input datacard, only used for rebinning")
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
      "Print tables of background shifts and relative uncertainties")
    ("freeze",
      po::value<string>(&freeze_arg)->default_value(freeze_arg),
      "Format PARAM1,PARAM2=X,PARAM3=Y where the values X and Y are optional");

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
    cout << "PostFitShapesFromWorkspace.root -d htt_mt_125.txt -w htt_mt_125.root -o htt_mt_125_shapes.root -m 125 "
            "-f mlfit.root:fit_s --postfit --sampling --print\n";
    return 1;
  }

  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  TFile infile(workspace.c_str());

  RooWorkspace *ws = dynamic_cast<RooWorkspace*>(gDirectory->Get("w"));

  if (!ws) {
    throw std::runtime_error(
        FNERROR("Could not locate workspace in input file"));
  }

  // Create CH instance and parse the workspace
  ch::CombineHarvester cmb;
  cmb.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(cmb, *ws, "ModelConfig", "data_obs", false);

  vector<string> freeze_vec;
  boost::split(freeze_vec, freeze_arg, boost::is_any_of(","));
  for (auto const& item : freeze_vec) {
    vector<string> parts;
    boost::split(parts, item, boost::is_any_of("="));
    if (parts.size() == 1) {
      ch::Parameter *par = cmb.GetParameter(parts[0]);
      if (par) par->set_frozen(true);
    } else {
      if (parts.size() == 2) {
        ch::Parameter *par = cmb.GetParameter(parts[0]);
        if (par) {
          par->set_val(boost::lexical_cast<double>(parts[1]));
          par->set_frozen(true);
        }
      }
    }
  }
  // cmb.GetParameter("r")->set_frozen(true);

  ch::CombineHarvester cmb_card;
  if (datacard != "") {
    cmb_card.ParseDatacard(datacard, "", "", "", 0, mass);
  }

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
          pre_shapes[bin][proc].SetName(proc.c_str());
    }

    // The fill total signal and total bkg hists
    pre_shapes[bin]["TotalBkg"] =
        cmb_bin.cp().backgrounds().GetShapeWithUncertainty();
    pre_shapes[bin]["TotalSig"] =
        cmb_bin.cp().signals().GetShapeWithUncertainty();
    pre_shapes[bin]["TotalProcs"] =
        cmb_bin.cp().GetShapeWithUncertainty();

    if (datacard != "") {
      TH1F ref = cmb_card.cp().bin({bin}).GetObservedShape();
      for (auto & it : pre_shapes[bin]) {
        it.second = ch::RestoreBinning(it.second, ref);
      }
    }

    // Can write these straight into the output file
    outfile.cd();
    for (auto& iter : pre_shapes[bin]) {
      ch::WriteHistToTFileAndSetName(&(iter.second), &outfile, bin + "_prefit/" + iter.first);
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
      post_shapes[bin]["TotalProcs"] =
          sampling ? cmb_bin.cp().GetShapeWithUncertainty(res, samples)
                   : cmb_bin.cp().GetShapeWithUncertainty();

      if (datacard != "") {
        TH1F ref = cmb_card.cp().bin({bin}).GetObservedShape();
        for (auto & it : post_shapes[bin]) {
          it.second = ch::RestoreBinning(it.second, ref);
        }
      }

      outfile.cd();
      // Write the post-fit histograms
      for (auto & iter : post_shapes[bin]) {
        ch::WriteHistToTFileAndSetName(&(iter.second), &outfile,
                         bin + "_postfit/" + iter.first);
      }
    }
  }
  // And we're done!
  outfile.Close();
  return 0;
}

