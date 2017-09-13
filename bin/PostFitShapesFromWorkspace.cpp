#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "TSystem.h"
#include "TH2F.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
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
  bool covariance   = false;
  string data       = "data_obs";
  bool skip_prefit  = false;
  bool skip_proc_errs = false;
  bool make_yield_tables=false;

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
      "Use the cov. matrix sampling method for the post-fit uncertainty")
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
    ("make_yield_tables",
      po::value<bool>(&make_yield_tables)->default_value(make_yield_tables)->implicit_value(true),
      "Print yield tables in Latex format");

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

  cmb.ForEachObj(boost::bind(ch::SetFromBinName,_1,"$ANALYSIS_$CHANNEL_$BINID_$ERA"));
  auto bins = cmb.cp().bin_set();

  TFile outfile(output.c_str(), "RECREATE");
  TH1::AddDirectory(false);

  // Create a map of maps for storing histograms in the form:
  //   pre_shapes[<bin>][<process>]
  map<string, map<string, TH1F>> pre_shapes;
  map<string, map<string, std::pair<double,double>>> pre_yields;
  // We can always do the prefit version,
  // Loop through the bins writing the shapes to the output file
  if (!skip_prefit) {
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

   //Setting up some maps of lables we'll need later
     std::vector<string> chns = {"mt","et","tt","em"};
     std::map<string,string> bkg_proc_labels;
     std::map<string,string> header_labels;
     std::map<string,std::vector<int>> channel_binids;
     std::map<string,std::map<int,string>> channel_binlabels;
     bkg_proc_labels["ZTT"] = "$\\PZ\\rightarrow\\Pgt\\Pgt$";
     bkg_proc_labels["ZL"] = "$\\PZ$+jets (l faking $\\Pgt$)";
     bkg_proc_labels["ZLL"] = "$\\PZ\\rightarrow\\ell\\ell$";
     bkg_proc_labels["W"] = "$\\PW$+jets";
     bkg_proc_labels["QCD"] = "QCD";
     bkg_proc_labels["jetFakes"] = "jet$\\rightarrow\\Pgt$ fakes";
     bkg_proc_labels["TT"] = "$\\cPqt\\cPaqt$";
     bkg_proc_labels["TTT"] = "$\\cPqt\\cPaqt$ (genuine $\\Pgt$)";
     bkg_proc_labels["VV"] = "Di-bosons+single top";
     bkg_proc_labels["VVT"] = "Di-bosons+single top (genuine $\\Pgt$)";
     header_labels["et"] = "$e\\tau_h$-channel";
     header_labels["mt"] = "$\\mu\\tau_h$-channel";
     header_labels["tt"] = "$\\tau_h\\tau_h$-channel";
     header_labels["em"] = "$e\\mu$-channel";
     channel_binids["et"] = {8,9,10,11};
     channel_binids["mt"] = {8,9,10,11};
     channel_binids["tt"] = {8,9};
     channel_binids["em"] = {8,9,10,11,12,13};
     channel_binlabels["tt"][8] = "No b-tag";
     channel_binlabels["tt"][9] = "B-tag";
     channel_binlabels["mt"][8] = "No b-tag tight";
     channel_binlabels["mt"][9] = "B-tag tight";
     channel_binlabels["mt"][10] = "No b-tag loose $m_{\\text{T}}$";
     channel_binlabels["mt"][11] = "B-tag loose $m_{\\text{T}}$";
     channel_binlabels["et"][8] = "No b-tag tight";
     channel_binlabels["et"][9] = "B-tag tight";
     channel_binlabels["et"][10] = "No b-tag loose $m_{\\text{T}}$";
     channel_binlabels["et"][11] = "B-tag loose $m_{\\text{T}}$";
     channel_binlabels["em"][8] = "No b-tag low $D_{\\zeta}$";
     channel_binlabels["em"][9] = "B-tag low $D_{\\zeta}$";
     channel_binlabels["em"][10] = "No b-tag medium $D_{\\zeta}$";
     channel_binlabels["em"][11] = "B-tag medium $D_{\\zeta}$";
     channel_binlabels["em"][12] = "No b-tag high $D_{\\zeta}$";
     channel_binlabels["em"][13] = "B-tag high $D_{\\zeta}$";
     if(make_yield_tables&&!postfit){
       for (string chn: chns){
         auto chn_bins = cmb.cp().channel({chn}).bin_id(channel_binids[chn]).bin_set();
         unsigned ncols = chn_bins.size();
         cout<< "\\begin{tabular}{l";
         for (unsigned i=0; i<ncols;++i) cout<< "r@{$ \\,\\,\\pm\\,\\,$}l";
         cout <<"}\n\\hline\n";
         cout<< boost::format("\\multicolumn{%i}{c}{%s} \\\\\n") % (ncols*2+1) % header_labels[chn];
         cout<< "\\hline\n";
         std::string cm_energy = "13~\\TeV";
         std::cout <<boost::format("& \\multicolumn{%i}{c}{$\\sqrt{s}$=%s} ")% (ncols * 2) %cm_energy;
         cout <<" \\\\\n";
         cout << "Process";
         for (auto bin: chn_bins){
           std::string bin_label = "";
           if (bin.find("_8_") != std::string::npos) bin_label = channel_binlabels[chn][8];
           else if (bin.find("_9_") != std::string::npos) bin_label = channel_binlabels[chn][9];
           else if (bin.find("_10_") != std::string::npos) bin_label = channel_binlabels[chn][10];
           else if (bin.find("_11_") != std::string::npos) bin_label = channel_binlabels[chn][11];
           else if (bin.find("_12_") != std::string::npos) bin_label = channel_binlabels[chn][12];
           else bin_label = channel_binlabels[chn][13];
           cout << "& \\multicolumn{2}{c}{" + bin_label+ "}";
         }
         cout<< "\\\\\n";
         cout<<"\\hline\n";
         string fmt = "& %-10.0f & %-10.0f";
         string fmts = "& %-10.0f & %-10.1f";
         ch::CombineHarvester chn_bin;
         for (auto bin :chn_bins){
           chn_bin = cmb.cp().bin({bin});
           auto chn_bkgs = chn_bin.cp().backgrounds();
           auto chn_sig = chn_bin.cp().signals();
           pre_yields[bin]["data_obs"] = std::make_pair(0.,0.);
           double tot_unc = chn_bkgs.GetUncertainty();
           pre_yields[bin]["TotalBkg"] = std::make_pair(chn_bkgs.GetRate(),tot_unc);
           double sig_unc = chn_sig.GetUncertainty();
           pre_yields[bin]["TotalSignal"] = std::make_pair(chn_sig.GetRate(),sig_unc);
           for( auto proc: chn_bkgs.process_set()){
              auto chn_proc = chn_bin.cp().process({proc});
              double proc_unc =chn_proc.GetUncertainty();
              pre_yields[bin][proc] = std::make_pair(chn_proc.GetRate(), proc_unc);
           }
         }
         auto chn_bkgs = chn_bin.cp().backgrounds();
         for (auto proc :chn_bkgs.process_set()){
           cout<< boost::format("%-38s") %bkg_proc_labels[proc];
           for (auto bin :chn_bins){
           cout <<boost::format(pre_yields[bin][proc].second <1. ? fmts:fmt) %
                  pre_yields[bin][proc].first % pre_yields[bin][proc].second;
           }
           cout << "\\\\\n";
        }
        cout << "\\hline\n";
        cout << boost::format("%-38s") % "Total Background";
        for (auto bin: chn_bins)
           cout <<boost::format(fmt) %pre_yields[bin]["TotalBkg"].first % pre_yields[bin]["TotalBkg"].second;
        cout<< "\\\\\n";
        cout<< "\\hline\n";
        cout << boost::format("%-38s") % "A/H$\\rightarrow\\Pgt\\Pgt$";
        for (auto bin : chn_bins)
            cout << boost::format(pre_yields[bin]["TotalSignal"].second<1.? fmts:fmt) % pre_yields[bin]["TotalSignal"].first % pre_yields[bin]["TotalSignal"].second;
          cout << "\\\\\n";
          cout << "\\hline\n";
          cout << boost::format("%-38s") % "Data";
          for (auto bin : chn_bins)
            cout << boost::format("& \\multicolumn{2}{c}{%-10.0f}") % (pre_yields[bin]["data_obs"]).first;
          cout << "\\\\\n";
          cout << "\\hline\n";
          cout << "\\end{tabular}\n";
         }
       }

  // Now we can do the same again but for the post-fit model
  if (postfit) {
    // Get the fit result and update the parameters to the post-fit model
    RooFitResult res = ch::OpenFromTFile<RooFitResult>(fitresult);
    cmb.UpdateParameters(res);

    // Calculate the post-fit fractional background uncertainty in each bin

    map<string, map<string, TH1F>> post_shapes;
    map<string, map<string, std::pair<double,double>>> post_yields;
    map<string, TH2F> post_yield_cov;
    map<string, TH2F> post_yield_cor;

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
              sampling ? cmb_proc.GetShapeWithUncertainty(res, samples)
                       : cmb_proc.GetShapeWithUncertainty();
        }
      }
      if (sampling && covariance) {
        post_yield_cov[bin] = cmb_bin.GetRateCovariance(res, samples);
        post_yield_cor[bin] = cmb_bin.GetRateCorrelation(res, samples);
      }
      // Fill the total sig. and total bkg. hists
      auto cmb_bkgs = cmb_bin.cp().backgrounds();
      auto cmb_sigs = cmb_bin.cp().signals();
      std::cout << ">> Doing postfit: " << bin << "," << "TotalBkg" << std::endl;
      post_shapes[bin]["TotalBkg"] =
          sampling ? cmb_bkgs.GetShapeWithUncertainty(res, samples)
                   : cmb_bkgs.GetShapeWithUncertainty();
      std::cout << ">> Doing postfit: " << bin << "," << "TotalSig" << std::endl;
      post_shapes[bin]["TotalSig"] =
          sampling ? cmb_sigs.GetShapeWithUncertainty(res, samples)
                   : cmb_sigs.GetShapeWithUncertainty();
      std::cout << ">> Doing postfit: " << bin << "," << "TotalProcs" << std::endl;
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
        double err = sampling ? cmb_bkgs.GetUncertainty(res, samples)
                              : cmb_bkgs.GetUncertainty();
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
    if(make_yield_tables){
     for (string chn: chns){
       auto chn_bins = cmb.cp().channel({chn}).bin_id(channel_binids[chn]).bin_set();
       unsigned ncols = chn_bins.size();
       cout<< "\\begin{tabular}{l";
       for (unsigned i=0; i<ncols;++i) cout<< "r@{$ \\,\\,\\pm\\,\\,$}l";
       cout <<"}\n\\hline\n";
       cout<< boost::format("\\multicolumn{%i}{c}{%s} \\\\\n") % (ncols*2+1) % header_labels[chn];
       cout<< "\\hline\n";
       std::string cm_energy = "13~\\TeV";
       std::cout <<boost::format("& \\multicolumn{%i}{c}{$\\sqrt{s}$=%s} ")% (ncols * 2) %cm_energy;
       cout <<" \\\\\n";
       cout << "Process";
       for (auto bin: chn_bins){
         std::string bin_label = "";
         if (bin.find("_8_") != std::string::npos) bin_label = channel_binlabels[chn][8];
         else if (bin.find("_9_") != std::string::npos) bin_label = channel_binlabels[chn][9];
         else if (bin.find("_10_") != std::string::npos) bin_label = channel_binlabels[chn][10];
         else if (bin.find("_11_") != std::string::npos) bin_label = channel_binlabels[chn][11];
         else if (bin.find("_12_") != std::string::npos) bin_label = channel_binlabels[chn][12];
         else bin_label = channel_binlabels[chn][13];
         cout << "& \\multicolumn{2}{c}{" + bin_label+ "}";
       }
       cout<< "\\\\\n";
       cout<<"\\hline\n";
       string fmt = "& %-10.0f & %-10.0f";
       string fmts = "& %-10.0f & %-10.1f";
       ch::CombineHarvester chn_bin;
       for (auto bin :chn_bins){
         chn_bin = cmb.cp().bin({bin});
         auto chn_bkgs = chn_bin.cp().backgrounds();
         auto chn_sig = chn_bin.cp().signals();
         post_yields[bin]["data_obs"] = std::make_pair(chn_bin.GetObservedRate(),0.);
         double tot_unc = sampling ? chn_bkgs.GetUncertainty(res,samples) : chn_bkgs.GetUncertainty();
         post_yields[bin]["TotalBkg"] = std::make_pair(chn_bkgs.GetRate(),tot_unc);
         double sig_unc = sampling ? chn_sig.GetUncertainty(res,samples) : chn_sig.GetUncertainty();
         post_yields[bin]["TotalSignal"] = std::make_pair(chn_sig.GetRate(),sig_unc);
         for( auto proc: chn_bkgs.process_set()){
            auto chn_proc = chn_bin.cp().process({proc});
            double proc_unc = sampling ? chn_proc.GetUncertainty(res,samples) : chn_proc.GetUncertainty();
            post_yields[bin][proc] = std::make_pair(chn_proc.GetRate(), proc_unc);
         }
       }
       auto chn_bkgs = chn_bin.cp().backgrounds();
       for (auto proc :chn_bkgs.process_set()){
         cout<< boost::format("%-38s") %bkg_proc_labels[proc];
         for (auto bin :chn_bins){
         cout <<boost::format(post_yields[bin][proc].second <1. ? fmts:fmt) %
                 post_yields[bin][proc].first % post_yields[bin][proc].second;
         }
         cout << "\\\\\n";
      }
      cout << "\\hline\n";
      cout << boost::format("%-38s") % "Total Background";
      for (auto bin: chn_bins)
         cout <<boost::format(fmt) %post_yields[bin]["TotalBkg"].first % post_yields[bin]["TotalBkg"].second;
      cout<< "\\\\\n";
      cout<< "\\hline\n";
      cout << boost::format("%-38s") % "A/H$\\rightarrow\\Pgt\\Pgt$";
      for (auto bin : chn_bins)
          cout << boost::format(post_yields[bin]["TotalSignal"].second<1.? fmts:fmt) % post_yields[bin]["TotalSignal"].first % post_yields[bin]["TotalSignal"].second;
        cout << "\\\\\n";
        cout << "\\hline\n";
        cout << boost::format("%-38s") % "Data";
        for (auto bin : chn_bins)
          cout << boost::format("& \\multicolumn{2}{c}{%-10.0f}") % (post_yields[bin]["data_obs"]).first;
        cout << "\\\\\n";
        cout << "\\hline\n";
        cout << "\\end{tabular}\n";
       }
     }
  }
  // And we're done!
  outfile.Close();
  return 0;
}

