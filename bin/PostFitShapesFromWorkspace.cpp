#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH2F.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace po = boost::program_options;

using namespace std;

void ReverseBins(TH1F & h) {
  std::vector<float> contents(h.GetNbinsX());
  std::vector<float> errors(h.GetNbinsX());
  for (int i = 0; i < h.GetNbinsX(); ++i) {
    contents[i] = h.GetBinContent(i + 1);
    errors[i] = h.GetBinError(i + 1);
  }
  for (int i = 0; i < h.GetNbinsX(); ++i) {
    h.SetBinContent(h.GetNbinsX() - i, contents[i]);
    h.SetBinError(h.GetNbinsX() - i, errors[i]);
  }
  // return h;
}

int main(int argc, char* argv[]) {
  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  string datacard   = "";
  string workspace  = "";
  string fitresult  = "";
  string mass       = "";
  bool postfit      = false;
  bool sampling     = false;
  bool no_sampling  = false;
  string output     = "";
  bool factors      = false;
  unsigned samples  = 500;
  std::string freeze_arg = "";
  bool covariance   = false;
  string data       = "data_obs";
  bool skip_prefit  = false;
  bool skip_proc_errs = false;
  bool total_shapes = false;
  std::vector<std::string> reverse_bins_;

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  config.add_options()
    ("workspace,w",
      po::value<string>(&workspace)->required(),
      "The input workspace-containing file [REQUIRED]")
    ("dataset",
      po::value<string>(&data)->default_value(data),
      "The input dataset name")
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
      "Use the cov. matrix sampling method for the post-fit uncertainty (deprecated, this is the default)")
    ("no-sampling",
      po::value<bool>(&no_sampling)->default_value(no_sampling)->implicit_value(true),
      "Do not use the cov. matrix sampling method for the post-fit uncertainty")
    ("samples",
      po::value<unsigned>(&samples)->default_value(samples),
      "Number of samples to make in each evaluate call")
    ("print",
      po::value<bool>(&factors)->default_value(factors)->implicit_value(true),
      "Print tables of background shifts and relative uncertainties")
    ("freeze",
      po::value<string>(&freeze_arg)->default_value(freeze_arg),
      "Format PARAM1,PARAM2=X,PARAM3=Y where the values X and Y are optional")
    ("covariance",
      po::value<bool>(&covariance)->default_value(covariance)->implicit_value(true),
      "Save the covariance and correlation matrices of the process yields")
    ("skip-prefit",
      po::value<bool>(&skip_prefit)->default_value(skip_prefit)->implicit_value(true),
      "Skip the pre-fit evaluation")
    ("skip-proc-errs",
      po::value<bool>(&skip_proc_errs)->default_value(skip_proc_errs)->implicit_value(true),
      "Skip evaluation of errors on individual processes")
    ("total-shapes",
      po::value<bool>(&total_shapes)->default_value(total_shapes)->implicit_value(true),
      "Save signal- and background shapes added for all channels/categories")
    ("reverse-bins", po::value<vector<string>>(&reverse_bins_)->multitoken(), "List of bins to reverse the order for");


  po::variables_map vm;

  // First check if the user has set the "--help" or "-h" option, and if so
  // just prin the usage information and quit
  po::store(po::command_line_parser(argc, argv)
    .options(help_config).allow_unregistered().run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    cout << config << "\nExample usage:\n";
    cout << "PostFitShapesFromWorkspace.root -d htt_mt_125.txt -w htt_mt_125.root -o htt_mt_125_shapes.root -m 125 "
            "-f mlfit.root:fit_s --postfit --print\n";
    return 1;
  }

  // Parse the main config options
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  if (sampling) {
    std::cout<<"WARNING: the default behaviour of PostFitShapesFromWorkspace is to use the covariance matrix sampling method for the post-fit uncertainty. The option --sampling is deprecated and will be removed in future versions of CombineHarvester"<<std::endl;
  }

  TFile infile(workspace.c_str());

  RooWorkspace *ws = dynamic_cast<RooWorkspace*>(gDirectory->Get("w"));

  if (!ws) {
    throw std::runtime_error(
        FNERROR("Could not locate workspace in input file"));
  }

  // Create CH instance and parse the workspace
  ch::CombineHarvester cmb;
  cmb.SetFlag("workspaces-use-clone", true);
  ch::ParseCombineWorkspace(cmb, *ws, "ModelConfig", data, false);

  // Only evaluate in case parameters to freeze are provided
  if(! freeze_arg.empty())
  {
    vector<string> freeze_vec;
    boost::split(freeze_vec, freeze_arg, boost::is_any_of(","));
    for (auto const& item : freeze_vec) {
      vector<string> parts;
      boost::split(parts, item, boost::is_any_of("="));
      if (parts.size() == 1) {
        ch::Parameter *par = cmb.GetParameter(parts[0]);
        if (par) par->set_frozen(true);
        else throw std::runtime_error(
          FNERROR("Requested variable to freeze does not exist in workspace"));
      } else {
        if (parts.size() == 2) {
          ch::Parameter *par = cmb.GetParameter(parts[0]);
          if (par) {
            par->set_val(boost::lexical_cast<double>(parts[1]));
            par->set_frozen(true);
          }
          else throw std::runtime_error(
            FNERROR("Requested variable to freeze does not exist in workspace"));
        }
      }
    }
  }
  // cmb.GetParameter("r")->set_frozen(true);

  ch::CombineHarvester cmb_card;
  cmb_card.SetFlag("workspaces-use-clone",true);
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

  auto bins = cmb.cp().bin_set();

  TFile outfile(output.c_str(), "RECREATE");
  TH1::AddDirectory(false);

  // Create a map of maps for storing histograms in the form:
  //   pre_shapes[<bin>][<process>]
  map<string, map<string, TH1F>> pre_shapes;

  // Also create a simple map for storing total histograms, summed 
  // over all bins, in the form:
  //   pre_shapes_tot[<process>]
  map<string, TH1F> pre_shapes_tot;

  // We can always do the prefit version,
  // Loop through the bins writing the shapes to the output file
  if (!skip_prefit) {
    if(total_shapes){
      pre_shapes_tot["data_obs"] = cmb.GetObservedShape();
      // Then fill total signal and total bkg hists
      std::cout << ">> Doing prefit: TotalBkg" << std::endl;
      pre_shapes_tot["TotalBkg"] =
          cmb.cp().backgrounds().GetShapeWithUncertainty();
      std::cout << ">> Doing prefit: TotalSig" << std::endl;
      pre_shapes_tot["TotalSig"] =
          cmb.cp().signals().GetShapeWithUncertainty();
      std::cout << ">> Doing prefit: TotalProcs" << std::endl;
      pre_shapes_tot["TotalProcs"] =
          cmb.cp().GetShapeWithUncertainty();

      if (datacard != "") {
        TH1F ref = cmb_card.cp().GetObservedShape();
        for (auto & it : pre_shapes_tot) {
          it.second = ch::RestoreBinning(it.second, ref);
        }
      }

      // Can write these straight into the output file
      outfile.cd();
      for (auto& iter : pre_shapes_tot) {
        ch::WriteToTFile(&(iter.second), &outfile, "prefit/" + iter.first);
      }
    }
    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});
      // This next line is a temporary fix for models with parameteric RooFit pdfs
      // - we try and set the number of bins to evaluate the pdf to be the same as
      // the number of bins in data
      // cmb_bin.SetPdfBins(cmb_bin.GetObservedShape().GetNbinsX());

      // Fill the data and process histograms
      pre_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
      for (auto proc : cmb_bin.process_set()) {
        std::cout << ">> Doing prefit: " << bin << "," << proc << std::endl;
        if (skip_proc_errs) {
          pre_shapes[bin][proc] =
              cmb_bin.cp().process({proc}).GetShape();
        } else {
          pre_shapes[bin][proc] =
              cmb_bin.cp().process({proc}).GetShapeWithUncertainty();
        }
      }

      // The fill total signal and total bkg hists
      std::cout << ">> Doing prefit: " << bin << "," << "TotalBkg" << std::endl;
      pre_shapes[bin]["TotalBkg"] =
          cmb_bin.cp().backgrounds().GetShapeWithUncertainty();
      std::cout << ">> Doing prefit: " << bin << "," << "TotalSig" << std::endl;
      pre_shapes[bin]["TotalSig"] =
          cmb_bin.cp().signals().GetShapeWithUncertainty();
      std::cout << ">> Doing prefit: " << bin << "," << "TotalProcs" << std::endl;
      pre_shapes[bin]["TotalProcs"] =
          cmb_bin.cp().GetShapeWithUncertainty();


      if (datacard != "") {
        TH1F ref = cmb_card.cp().bin({bin}).GetObservedShape();
        for (auto & it : pre_shapes[bin]) {
          it.second = ch::RestoreBinning(it.second, ref);
        }
      }

      for (auto const& rbin : reverse_bins_) {
        if (rbin != bin) continue;
        auto & hists = pre_shapes[bin];
        for (auto it = hists.begin(); it != hists.end(); ++it) {
          ReverseBins(it->second);
        }
      }
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
        cout << boost::format("%-25s %-10.5f") % bin %
                    (rate > 0. ? (err / rate) : 0.) << std::endl;
      }
    }
  }


  // Now we can do the same again but for the post-fit model
  if (postfit) {
    // Get the fit result and update the parameters to the post-fit model
    RooFitResult res = ch::OpenFromTFile<RooFitResult>(fitresult);
    cmb.UpdateParameters(res);

    // Calculate the post-fit fractional background uncertainty in each bin

    map<string, map<string, TH1F>> post_shapes;
    map<string, TH2F> post_yield_cov;
    map<string, TH2F> post_yield_cor;

    map<string, TH1F> post_shapes_tot;

    if(total_shapes){
      post_shapes_tot["data_obs"] = cmb.GetObservedShape();
      // Fill the total sig. and total bkg. hists
      auto cmb_bkgs = cmb.cp().backgrounds();
      auto cmb_sigs = cmb.cp().signals();
      std::cout << ">> Doing postfit: TotalBkg" << std::endl;
      post_shapes_tot["TotalBkg"] =
          no_sampling ? cmb_bkgs.GetShapeWithUncertainty()
                   : cmb_bkgs.GetShapeWithUncertainty(res, samples);
      std::cout << ">> Doing postfit: TotalSig" << std::endl;
      post_shapes_tot["TotalSig"] =
          no_sampling ? cmb_sigs.GetShapeWithUncertainty()
                   : cmb_sigs.GetShapeWithUncertainty(res, samples);
      std::cout << ">> Doing postfit: TotalProcs" << std::endl;
      post_shapes_tot["TotalProcs"] =
          no_sampling ? cmb.cp().GetShapeWithUncertainty()
                   : cmb.cp().GetShapeWithUncertainty(res, samples);

      if (datacard != "") {
        TH1F ref = cmb_card.cp().GetObservedShape();
        for (auto & it : post_shapes_tot) {
          it.second = ch::RestoreBinning(it.second, ref);
        }
      }

      outfile.cd();
      // Write the post-fit histograms
      for (auto & iter : post_shapes_tot) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         "postfit/" + iter.first);
      }
    }


    for (auto bin : bins) {
      ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});
      post_shapes[bin]["data_obs"] = cmb_bin.GetObservedShape();
      for (auto proc : cmb_bin.process_set()) {
        auto cmb_proc = cmb_bin.cp().process({proc});
        // Method to get the shape uncertainty depends on whether we are using
        // the sampling method or the "wrong" method (assumes no correlations)
        std::cout << ">> Doing postfit: " << bin << "," << proc << std::endl;
        if (skip_proc_errs) {
          post_shapes[bin][proc] = cmb_proc.GetShape();
        } else {
          post_shapes[bin][proc] =
              no_sampling ? cmb_proc.GetShapeWithUncertainty()
                       : cmb_proc.GetShapeWithUncertainty(res, samples);
        }
      }
      if (!no_sampling && covariance) {
        post_yield_cov[bin] = cmb_bin.GetRateCovariance(res, samples);
        post_yield_cor[bin] = cmb_bin.GetRateCorrelation(res, samples);
      }
      // Fill the total sig. and total bkg. hists
      auto cmb_bkgs = cmb_bin.cp().backgrounds();
      auto cmb_sigs = cmb_bin.cp().signals();
      std::cout << ">> Doing postfit: " << bin << "," << "TotalBkg" << std::endl;
      post_shapes[bin]["TotalBkg"] =
          no_sampling ? cmb_bkgs.GetShapeWithUncertainty()
                   : cmb_bkgs.GetShapeWithUncertainty(res, samples);
      std::cout << ">> Doing postfit: " << bin << "," << "TotalSig" << std::endl;
      post_shapes[bin]["TotalSig"] =
          no_sampling ? cmb_sigs.GetShapeWithUncertainty()
                   : cmb_sigs.GetShapeWithUncertainty(res, samples);
      std::cout << ">> Doing postfit: " << bin << "," << "TotalProcs" << std::endl;
      post_shapes[bin]["TotalProcs"] =
          no_sampling ? cmb_bin.cp().GetShapeWithUncertainty()
                   : cmb_bin.cp().GetShapeWithUncertainty(res, samples);

      if (datacard != "") {
        TH1F ref = cmb_card.cp().bin({bin}).GetObservedShape();
        for (auto & it : post_shapes[bin]) {
          it.second = ch::RestoreBinning(it.second, ref);
        }
      }

      outfile.cd();
      // Write the post-fit histograms

      for (auto const& rbin : reverse_bins_) {
        if (rbin != bin) continue;
        std::cout << ">> reversing hists in bin " << bin << "\n";
        auto & hists = post_shapes[bin];
        for (auto it = hists.begin(); it != hists.end(); ++it) {
          ReverseBins(it->second);
        }
      }

      for (auto & iter : post_shapes[bin]) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         bin + "_postfit/" + iter.first);
      }
      for (auto & iter : post_yield_cov) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         iter.first+"_cov");
      }
      for (auto & iter : post_yield_cor) {
        ch::WriteToTFile(&(iter.second), &outfile,
                         iter.first+"_cor");
      }

    }

    if (factors) {
      cout << boost::format("\n%-25s %-32s\n") % "Bin" %
                  "Total relative bkg uncert. (postfit)";
      cout << string(58, '-') << "\n";
      for (auto bin : bins) {
        ch::CombineHarvester cmb_bkgs = cmb.cp().bin({bin}).backgrounds();
        double rate = cmb_bkgs.GetRate();
        double err = no_sampling ? cmb_bkgs.GetUncertainty()
                              : cmb_bkgs.GetUncertainty(res, samples);
        cout << boost::format("%-25s %-10.5f") % bin %
                    (rate > 0. ? (err / rate) : 0.) << std::endl;
      }
    }

    // As we calculate the post-fit yields can also print out the post/pre scale
    // factors
    if (factors && postfit) {
      cout << boost::format("\n%-25s %-20s %-10s\n") % "Bin" % "Process" %
                  "Scale factor";
      cout << string(58, '-') << "\n";
      for (auto bin : bins) {
        ch::CombineHarvester cmb_bin = cmb.cp().bin({bin});

        for (auto proc : cmb_bin.process_set()) {
          // Print out the post/pre scale factors
          TH1 const& pre = pre_shapes[bin][proc];
          TH1 const& post = post_shapes[bin][proc];
          cout << boost::format("%-25s %-20s %-10.5f\n") % bin % proc %
                      (pre.Integral() > 0. ? (post.Integral() / pre.Integral())
                                           : 1.0);
        }
      }
    }
  }
  // And we're done!
  outfile.Close();
  return 0;
}

