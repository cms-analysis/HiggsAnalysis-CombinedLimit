#include <iostream>
#include <vector>
#include <map>

#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"

#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/SOverBTools.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include "CombineHarvester/CombineTools/interface/Plotting.h"
#include "CombineHarvester/CombineTools/interface/Plotting_Style.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]){
  TH1::AddDirectory(0);
  ModTDRStyle();
  string cfg;                                   // The configuration file

  bool do_ratio = false;
  bool do_logy = true;

  vector<string> datacards;
  string fitresult_file = "";
  string parse_rule = "";

  string output           = "";
  string text1            = "";
  string text2            = "";
  bool postfit            = true;

  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
    ("help,h",  "print the help message")
    ("datacards,d",
      po::value<vector<string>> (&datacards)->multitoken()->required())
    ("fitresult,f",
      po::value<string> (&fitresult_file))
    ("parse_rule,p",
      po::value<string> (&parse_rule)
      ->default_value("$MASS/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt"))
    ("output",
      po::value<string>(&output)->required(),
      "[REQUIRED] output name (no extension)")
    ("text1",
      po::value<string>(&text1)->default_value(""),
      "[REQUIRED] output name (no extension)")
    ("text2",
      po::value<string>(&text2)->default_value(""),
      "[REQUIRED] output name (no extension)")
    ("postfit",
      po::value<bool>(&postfit)->required(),
      "[REQUIRED] use the pulls file to make a post-fit plot");
  po::store(po::command_line_parser(argc, argv)
    .options(config).allow_unregistered().run(), vm);
  po::notify(vm);

  ch::CombineHarvester cmb;
  // cmb.SetVerbosity(2);
  for (auto const& d : datacards) {
    cmb.ParseDatacard(d, parse_rule);
  }
  ch::SetStandardBinNames(cmb);

  RooFitResult* fitresult = nullptr;
  if (fitresult_file.length() && postfit) {
    fitresult =
        new RooFitResult(ch::OpenFromTFile<RooFitResult>(fitresult_file));
    auto fitparams = ch::ExtractFitParameters(*fitresult);
    cmb.UpdateParameters(fitparams);
  }

  auto bins = cmb.bin_set();



  TH1F proto("proto", "proto", 17, -3, 0.4);

  for (auto const& bin : bins) {
    ch::CombineHarvester cmb_bin = std::move(cmb.cp().bin({bin}));
    TH1F bkg = cmb_bin.cp().backgrounds().GetShape();
    TH1F sig = cmb_bin.cp().signals().GetShape();
    std::map<int, int> bin_mapping;
    for (int b = 1; b <= bkg.GetNbinsX(); ++b) {
      float i_sig = sig.GetBinContent(b);
      float i_bkg = bkg.GetBinContent(b);
      int bin_result = 1;
      if (i_sig > 0.) {
        float s_sb = i_sig / (i_sig + i_bkg);
        bin_result = proto.FindFixBin(std::log10(s_sb));
        if (bin_result == 0) {
          bin_result = 1;
        } else if (bin_result == 16) {
          std::cout << "Bin " << bin << " has an s/s+b == 1\n";
        } else if (bin_result == 18) {
          std::cout << "Bin " << bin << " has an s/s+b: " << s_sb << "\n";
          std::cout << "Bin " << bin << " has s: " << i_sig << "\n";
          std::cout << "Bin " << bin << " has b: " << i_bkg << "\n";
        }
      }
      bin_mapping[b] = bin_result;
    }
    cmb_bin.ForEachObs([&](ch::Observation *e) {
      TH1 const* old_h = e->shape();
      std::unique_ptr<TH1> new_h = ch::make_unique<TH1F>(TH1F(proto));
      for (int b = 1; b <= old_h->GetNbinsX(); ++b) {
        int new_bin = (bin_mapping.find(b))->second;
        new_h->SetBinContent(
            new_bin,
            new_h->GetBinContent(new_bin) + old_h->GetBinContent(b));
      }
      e->set_shape(std::move(new_h), false);
    });

    cmb_bin.ForEachProc([&](ch::Process *e) {
      TH1 const* old_h = e->shape();
      std::unique_ptr<TH1> new_h = ch::make_unique<TH1F>(TH1F(proto));
      for (int b = 1; b <= old_h->GetNbinsX(); ++b) {
        int new_bin = (bin_mapping.find(b))->second;
        new_h->SetBinContent(
            new_bin,
            new_h->GetBinContent(new_bin) + old_h->GetBinContent(b));
      }
      e->set_shape(std::move(new_h), false);
    });
    cmb_bin.ForEachSyst([&](ch::Systematic *e) {
      if (e->type() != "shape") return;
      TH1 const* old_hd = e->shape_d();
      TH1 const* old_hu = e->shape_u();
      std::unique_ptr<TH1> new_hd = ch::make_unique<TH1F>(TH1F(proto));
      std::unique_ptr<TH1> new_hu = ch::make_unique<TH1F>(TH1F(proto));
      for (int b = 1; b <= old_hd->GetNbinsX(); ++b) {
        int new_bin = (bin_mapping.find(b))->second;
        new_hd->SetBinContent(
            new_bin,
            new_hd->GetBinContent(new_bin) + old_hd->GetBinContent(b));
        new_hu->SetBinContent(
            new_bin,
            new_hu->GetBinContent(new_bin) + old_hu->GetBinContent(b));
      }
      e->set_shapes(std::move(new_hu), std::move(new_hd), nullptr);
    });
  }

  std::vector<TH1F *> by_chn;

  auto channels = cmb.channel_set();

  THStack stack_bkg;

  for (auto chn : channels) {
    by_chn.push_back(
        new TH1F(cmb.cp().channel({chn}).backgrounds().GetShape()));
    if (chn == "et") by_chn.back()->SetFillColor(TColor::GetColor(84, 166, 226));
    if (chn == "mt") by_chn.back()->SetFillColor(TColor::GetColor(212, 66,  88));
    if (chn == "em") by_chn.back()->SetFillColor(TColor::GetColor(201, 198, 116));
    if (chn == "tt") by_chn.back()->SetFillColor(TColor::GetColor(245, 197, 85));
    if (chn == "ee") by_chn.back()->SetFillColor(TColor::GetColor(247, 186, 254));
    if (chn == "mm") by_chn.back()->SetFillColor(TColor::GetColor(112, 142, 122));
    if (chn == "vhtt") by_chn.back()->SetFillColor(TColor::GetColor(137, 132, 192));
    if (chn == "et") by_chn.back()->SetTitle("e_{}#tau_{h}");
    if (chn == "mt") by_chn.back()->SetTitle("#mu_{}#tau_{h}");
    if (chn == "em") by_chn.back()->SetTitle("e#mu");
    if (chn == "ee") by_chn.back()->SetTitle("ee");
    if (chn == "mm") by_chn.back()->SetTitle("#mu#mu");
    if (chn == "tt") by_chn.back()->SetTitle("#tau_{h}#tau_{h}");
    if (chn == "vhtt") by_chn.back()->SetTitle("VH");
    by_chn.back()->SetBinContent(17, 0.);
    stack_bkg.Add(by_chn.back());

  }

  TH1F sb_sig = cmb.cp().signals().GetShape();
  TH1F sb_bkg = cmb.cp().backgrounds().GetShapeWithUncertainty();
  TH1F sb_err = sb_bkg;
  TH1F sb_obs = cmb.cp().GetObservedShape();

  sb_sig.SetBinContent(17, 0.);
  sb_bkg.SetBinContent(17, 0.);
  sb_obs.SetBinContent(17, 0.);

  sb_bkg.SetFillColor(18);
  sb_sig.SetFillColor(kRed);
  sb_sig.SetFillStyle(3004);
  sb_sig.SetLineColor(kRed);

  int new_idx = CreateTransparentColor(12, 0.4);
  sb_err.SetFillColor(new_idx);
  sb_err.SetMarkerSize(0);

  THStack *stack = new THStack("stack", "stack");
  stack->Add((TH1F*)sb_bkg.Clone(), "hist");
  stack->Add((TH1F*)sb_sig.Clone(), "hist");


  TCanvas* canv = new TCanvas(output.c_str(), output.c_str());
  canv->cd();

  std::vector<TPad*> pads =
      do_ratio ? TwoPadSplit(0.29, 0.005, 0.005) : OnePad();

  pads[0]->SetLogy(do_logy);
  std::vector<TH1*> h = CreateAxisHists(2, &sb_obs);
  if (do_ratio) {
    SetupTwoPadSplitAsRatio(pads, h[0], h[1], "Obs - Bkg", true, -10., 10.);
    h[1]->GetXaxis()->SetTitle("log(S/S+B)");
    h[0]->GetYaxis()->SetTitle("Events");
  } else {
    h[0]->GetXaxis()->SetTitle("log(S/S+B)");
    h[0]->GetYaxis()->SetTitle("Events");
  }

  h[0]->Draw("");

  if (pads[0]->GetLogy()) h[0]->SetMinimum(0.1);

  stack->Draw("histsame");
  stack_bkg.Draw("histsame");
  sb_err.Draw("e2same");

  sb_obs.Draw("esamex0");

  FixTopRange(pads[0], GetPadYMax(pads[0]), 0.15);

  TLegend *legend = PositionedLegend(0.30, 0.30, 3, 0.03);
  legend->AddEntry(&sb_obs, "Observed", "pe");
  // legend->AddEntry(&sb_bkg, "Background", "f");
  for (auto chn: by_chn) {
    legend->AddEntry(chn, chn->GetTitle(), "f");
  }
  legend->AddEntry(&sb_sig, "Signal", "f");
  // legend->AddEntry(&sb_err, "Bkg. Uncertainty", "f");
  legend->Draw();

  for (auto pad : pads) {
    pad->cd();
    FixOverlay();
  }

  DrawTitle(pads[0], "CombineHarvester", 1);



  canv->Print(".pdf");



  // TFile outfile((output+".root").c_str(), "RECREATE");
  // outfile.cd();
  // ch::WriteToTFile(&sb_sig, &outfile, "TotalSig");
  // ch::WriteToTFile(&sb_bkg, &outfile, "TotalBkg");
  // ch::WriteToTFile(&sb_obs, &outfile, "TotalObs");
  // outfile.Close();



  return 0;
}

