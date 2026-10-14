#include <iostream>
#include <vector>
#include <map>
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"

#include "TLatex.h"
#include "TH2F.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLine.h"
#include "TColor.h"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"

// #include "Utilities/interface/FnRootTools.h"

#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/SOverBTools.h"

namespace po = boost::program_options;

using namespace std;
// using namespace ic;


int main(int argc, char* argv[]){
  string cfg;                                   // The configuration file

  vector<string> datacards;
  string fitresult_file = "";
  string parse_rule = "";
  string texname          = "";
  string shortname        = "";
  bool postfit            = true;

  po::variables_map vm;
  po::notify(vm);
  po::options_description config("configuration");
  config.add_options()
    ("help,h",  "print the help message")
    ("datacard,d",
      po::value<vector<string>> (&datacards)->multitoken()->required())
    ("fitresult,f",
      po::value<string> (&fitresult_file))
    ("parse_rule,p",
      po::value<string> (&parse_rule)
      ->default_value("$MASS/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt"))
    ("texname",
      po::value<string>(&texname)->required(),
      "[REQUIRED] output name (no extension)")
    ("shortname",
      po::value<string>(&shortname)->default_value(""),
      "[REQUIRED] output name (no extension)")
    ("postfit",
      po::value<bool>(&postfit)->required(),
      "[REQUIRED] use the pulls file to make a post-fit plot");
  po::store(po::command_line_parser(argc, argv)
    .options(config).allow_unregistered().run(), vm);
  po::notify(vm);

  ch::CombineHarvester cmb;
  for (auto const& d : datacards) {
    cmb.ParseDatacard(d, parse_rule);
  }
  ch::SetStandardBinNames(cmb);


  RooFitResult *fitresult = nullptr;
  if (fitresult_file.length() && postfit) {
    fitresult = new RooFitResult(ch::OpenFromTFile<RooFitResult>(fitresult_file));
    auto fitparams = ch::ExtractFitParameters(*fitresult);
    cmb.UpdateParameters(fitparams);
  }

  TH1F sig = cmb.cp().signals().GetShape();
  TH1F bkg = cmb.cp().backgrounds().GetShape();
  ch::SOverBInfo info(ch::SOverBInfo(&sig, &bkg, 3500, 0.682));
  double s_over_sb = info.s/(info.s + info.b);
  double s_over_root_sb = info.s/std::sqrt(info.s + info.b);
  double width = (info.x_hi-info.x_lo)/2.;

  double tot_sig     = cmb.cp()
      .process({"ggH", "qqH", "WH", "ZH"}).GetRate();
  double tot_sig_err = 0.0;
  if (fitresult) {
    tot_sig_err = cmb.cp()
        .process({"ggH", "qqH", "WH", "ZH"})
        .GetUncertainty(fitresult, 500);
  } else {
    tot_sig_err = cmb.cp()
        .process({"ggH", "qqH", "WH", "ZH"})
        .GetUncertainty();
  }

  double tot_bkg     = cmb.cp()
      .backgrounds().GetRate();
  double tot_bkg_err = 0.0;
  if (fitresult) {
    tot_bkg_err = cmb.cp()
        .backgrounds()
        .GetUncertainty(fitresult, 500);
  } else {
    tot_bkg_err = cmb.cp()
        .backgrounds()
        .GetUncertainty();
  }

  if (tot_bkg_err > 100.) {
    tot_bkg = std::floor((tot_bkg/10.) + 0.5) * 10.;
    tot_bkg_err = std::floor((tot_bkg_err/10.) + 0.5) * 10.;
  }

  double tot_dat = cmb.GetObservedRate();

  double frac_ggh = (cmb.cp()
      .process({"ggH"}).GetRate()) / tot_sig;
  // double frac_qqh = (cmb.cp()
  //     .process(true, {"qqH"}).GetRate()) / tot_sig;
  double frac_vh = (cmb.cp()
      .process({"WH", "ZH"}).GetRate()) / tot_sig;

  // std::cout << frac_ggh << std::endl;
  // std::cout << frac_qqh << std::endl;
  // std::cout << frac_vh << std::endl;

  string def_sb = "%.3f";
  if (s_over_sb >= 0.1) def_sb = "%.2f";
  string def_rsb = "%.3f";
  if (s_over_root_sb >= 0.1) def_rsb = "%.2f";

  std::cout << boost::format("%s & "
    "\\multicolumn{3}{c}{\\includegraphics[height=3mm]{figs/htt-boxes/%s.pdf}} &"
    " %.1f,%.1f & %.0f,%.0f & %.0f & "+def_sb+" & "+def_rsb+" & %.1f \\\\ \n")
    % texname % shortname % tot_sig % tot_sig_err % tot_bkg % tot_bkg_err % tot_dat % s_over_sb % s_over_root_sb % width;

    TCanvas canv("canv","canv",350,70);
    TBox x;
    x.SetFillColor(kViolet+10); // VH
    x.DrawBox(1.-frac_vh,0,1.0,1);
    x.SetFillColor(kAzure+7); // ggH
    x.DrawBox(0.,0,frac_ggh,1);
    x.SetFillColor(kTeal+5); // qqH
    x.DrawBox(frac_ggh,0,1.-frac_vh,1);
    // canv.SaveAs((shortname+".pdf").c_str());

  return 0;
}

