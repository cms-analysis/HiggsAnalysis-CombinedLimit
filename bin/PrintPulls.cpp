#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "TTree.h"

#include "boost/regex.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"

#include "CombineHarvester/CombineTools/interface/Plotting_Style.h"

namespace po = boost::program_options;
using namespace std;


struct Vals {
  std::string parameter;
  double par_best;
  double par_lo;
  double par_hi;
  double par_best_post;
  double par_lo_post;
  double par_hi_post;
  double par_best_pre;
  double par_lo_pre;
  double par_hi_pre;
  double impact;
  double impact_post;
  double impact_pre;
  int rank_impact;
  int rank_impact_post;
  int rank_impact_pre;
};


// A struct and some functions extracted from framework and reproduced here
struct Pull {
  std::string name;
  double      prefit;
  double      prefit_err;
  double      bonly;
  double      bonly_err;
  double      splusb;
  double      splusb_err;
  double      rho;
  void Print() const;
};
void PullsFromFile(std::string const& filename, std::vector<Pull> & pullvec, bool verbose);
std::vector<std::string> ParseFileLines(std::string const& file_name);
bool BvsSBComparator(Pull const& pull1, Pull const& pull2);
bool ImpactComparator(Pull const& pull1, Pull const& pull2);

int main(int argc, char* argv[]) {

  TFile f("impact.root");
  TTree *t = (TTree*)f.Get("impact");

  Vals reader;
  TString *param_str = nullptr;
  t->SetBranchAddress("parameter",        &param_str);
  t->SetBranchAddress("par_best",         &reader.par_best);
  t->SetBranchAddress("par_lo",           &reader.par_lo);
  t->SetBranchAddress("par_hi",           &reader.par_hi);
  t->SetBranchAddress("par_best_post",    &reader.par_best_post);
  t->SetBranchAddress("par_lo_post",      &reader.par_lo_post);
  t->SetBranchAddress("par_hi_post",      &reader.par_hi_post);
  t->SetBranchAddress("par_best_pre",     &reader.par_best_pre);
  t->SetBranchAddress("par_lo_pre",       &reader.par_lo_pre);
  t->SetBranchAddress("par_hi_pre",       &reader.par_hi_pre);
  t->SetBranchAddress("impact",           &reader.impact);
  t->SetBranchAddress("impact_post",      &reader.impact_post);
  t->SetBranchAddress("impact_pre",       &reader.impact_pre);
  t->SetBranchAddress("rank_impact",      &reader.rank_impact);
  t->SetBranchAddress("rank_impact_post", &reader.rank_impact_post);
  t->SetBranchAddress("rank_impact_pre",  &reader.rank_impact_pre);

  std::map<std::string, Vals> impact_map;
  for (unsigned i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    reader.parameter = *param_str;
    impact_map[reader.parameter] = reader;
  }
  // Define and parse arguments
  string input1, name1;
  bool splusb_1;
  string input2, name2;
  bool splusb_2;
  bool draw_difference;
  bool draw_first;
  string output;
  vector<string> filter_regex_str;
  int max = -1;

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");
  po::options_description config("Configuration");
  config.add_options()
    ("input1",                po::value<string>(&input1)->required(), "The first mlfit file [REQUIRED]")
    ("input2",                po::value<string>(&input2)->required(), "The second mlfit file [REQUIRED]")
    ("label1",                po::value<string>(&name1)->default_value("input1"), "A label for the first input")
    ("label2",                po::value<string>(&name2)->default_value("input2"), "A label for the second input")
    ("output",                po::value<string>(&output)->default_value(""), "output filename")
    ("sb1",                   po::value<bool>(&splusb_1)->required(), "Use s+b pulls from the first file? [REQUIRED]")
    ("sb2",                   po::value<bool>(&splusb_2)->required(), "Use s+b pulls from the second file? [REQUIRED]")
    ("draw_difference",       po::value<bool>(&draw_difference)->default_value(true), "Draw the difference between inputs")
    ("draw_first",            po::value<bool>(&draw_first)->default_value(false), "Draw only first")
    ("max",                   po::value<int>(&max)->default_value(max), "Max number to print")
    ("filter_regex",          po::value<vector<string>>(&filter_regex_str), "A regular expression to filter pulls");
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(help_config).allow_unregistered().run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    cout << config << "\n";
    cout << "Example usage: " << endl;
    cout << "print-pulls --input1=mlfit_a.root --input2=mlfit_b.root --label1=\"PullsA\" --label2=\"PullsB\" --sb1=true --sb2=false --filter_regex=\".*_bin_.*\"" << endl;
    return 1;
  }
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::notify(vm);

  // Build a list of regular expressions to filter nuisances
  std::vector<boost::regex> filter_regex;
  for (unsigned i = 0; i < filter_regex_str.size(); ++i) {
    filter_regex.push_back(boost::regex(filter_regex_str[i]));
    std::cout << "Filter nuisances with regex: " << filter_regex_str[i]
              << std::endl;
  }

  // Set a nice drawing style
  ModTDRStyle();

  // Extract pulls from the two input files
  std::vector<Pull> pulls1;
  std::vector<Pull> pulls2;
  PullsFromFile(input1, pulls1, false);
  PullsFromFile(input2, pulls2, false);

  // Build new lists of the pulls common to both inputs, and in the same order
  std::vector<Pull> pulls1sorted;
  std::vector<Pull> pulls2sorted;
  for (unsigned i = 0; i < pulls1.size(); ++i) {
    bool found = false;
    for (unsigned j = 0; j < pulls2.size(); ++j) {
      if (pulls1[i].name == pulls2[j].name) {
        found = true;
        pulls2sorted.push_back(pulls2[j]);
        continue;
      }
    }
    if (found == true) pulls1sorted.push_back(pulls1[i]);
  }

  // Build a final vector of pulls taking b-only or s+b for each input
  // respectively - skip pulls that match one of the filter regex
  std::vector<Pull> final;
  for (unsigned i = 0; i < pulls1sorted.size(); ++i) {
    bool skip = false;
    for (unsigned j = 0; j < filter_regex.size(); ++j) {
      if (boost::regex_match(pulls1sorted[i].name, filter_regex[j]))
        skip = true;
    }
    if (skip) continue;
    Pull pull;
    pull.name = pulls1sorted[i].name;
    pull.bonly = splusb_1 ? pulls1sorted[i].splusb : pulls1sorted[i].bonly ;
    pull.bonly_err = splusb_1 ? pulls1sorted[i].splusb_err : pulls1sorted[i].bonly_err ;
    pull.splusb = splusb_2 ? pulls2sorted[i].splusb : pulls2sorted[i].bonly ;
    pull.splusb_err = splusb_2 ? pulls2sorted[i].splusb_err : pulls2sorted[i].bonly_err ;
    if (!impact_map.count(pull.name)) {
      std::cout << pull.name << "\n";
      pull.rho = 0.;
    } else {
      pull.rho = impact_map.at(pull.name).impact;
      std::cout << pull.rho << "\n";
    }
    final.push_back(pull);
  }
  // std::sort(final.begin(), final.end(), ImpactComparator);
  std::sort(final.begin(), final.end(), BvsSBComparator);
  unsigned npulls = final.size();
  if (max > 0 && max < int(npulls)) npulls = max;


  // Prepare the canvas
  TCanvas canv("canvas", "canvas", draw_difference ? 1200 : 800 , 1200);
  canv.cd();
  canv.Update();
  TPad* pad1 = NULL;
  TPad* pad2 = NULL;
  if (draw_difference) {
    pad1 = new TPad("pad1","pad1",0, 0, 0.85, 1);
    pad2 = new TPad("pad2","pad2",0.85, 0, 1, 1);
    pad1->Draw();
    pad2->Draw();
  } else {
    pad1 = new TPad("pad1","pad1",0, 0, 1, 1);
  }
  pad1->SetLeftMargin(0.45);
  pad1->SetRightMargin(0.03);
  pad1->SetGrid(1,0);
  pad1->Draw();
  pad1->cd();

  // Create a dummy TH2F to attach the nuisance labels
  // and TGraphErrors for each set of pulls
  TH2F *hpulls = new TH2F("pulls","pulls", 6, -3, 3, npulls, 0, npulls);
  TGraphErrors vals(npulls);
  TGraphErrors vals2(npulls);
  float ts1 = 30./(pad1->GetWw()*pad1->GetAbsWNDC());
  for (unsigned i = 0; i < npulls; ++i) {
    vals.SetPoint(i, final[i].bonly, double(i) + 0.5);
    vals.SetPointError(i, final[i].bonly_err, 0);
    vals2.SetPoint(i, final[i].splusb, double(i) + 0.5);
    vals2.SetPointError(i, final[i].splusb_err, 0);
    hpulls->GetYaxis()->SetBinLabel(i+1, final[i].name.c_str());
  }
  hpulls->GetYaxis()->LabelsOption("v");
  hpulls->SetStats(false);
  hpulls->GetYaxis()->SetLabelSize(0.025);
  gStyle->SetEndErrorSize(5);
  hpulls->GetXaxis()->SetTitle("Fit Value #hat{#theta}");
  hpulls->GetXaxis()->CenterTitle();
  // hpulls->GetXaxis()->SetTitleSize(0.04);
  hpulls->Draw("");
  hpulls->GetXaxis()->SetLabelSize(ts1);
  vals.Draw("pSAME");
  vals.SetLineWidth(2);
  vals2.SetLineWidth(2);
  vals.SetMarkerStyle(20);
  vals2.SetMarkerStyle(24);
  vals2.SetLineColor(4);
  vals2.SetMarkerColor(4);
  if(!draw_first) vals2.Draw("pSAME");
  canv.Update();

  TGraphErrors valsd;
  if (pad2) {
    pad2->cd();
    // pad2->SetBottomMargin(0.07);
    pad2->SetLeftMargin(0.1);
    pad2->SetRightMargin(0.2);
    // pad2->SetTopMargin(0.05);
    pad2->SetGrid(1,0);
    TH2F *hpulldiff = new TH2F("pulldiff","pulldiff", 2, 0, 0.08, npulls, 0, npulls);
    hpulldiff->GetYaxis()->SetLabelSize(0.0);
    hpulldiff->GetYaxis()->SetTitleSize(0.0);
    hpulldiff->GetXaxis()->SetTitle("#splitline{Impact}{(rank)}");

    valsd = TGraphErrors(npulls);
    for (unsigned i = 0; i < npulls; ++i) {
      valsd.SetPoint(i, (final[i].rho)/2., double(i) + 0.5);
      valsd.SetPointError(i, (final[i].rho)/2., 0.3);
    }
    hpulldiff->SetStats(false);
    if (name1 == "") name1 = input1;
    if (name2 == "") name2 = input2;
    float ts2 = 30./(pad2->GetWw()*pad2->GetAbsWNDC());
    hpulldiff->GetXaxis()->SetLabelSize(ts2);
    hpulldiff->GetXaxis()->SetTitleSize(0.20);
    hpulldiff->GetXaxis()->SetTitleOffset(0.25);
    hpulldiff->GetXaxis()->SetLabelOffset(-0.105);
    hpulldiff->GetXaxis()->SetNdivisions(502);
    hpulldiff->Draw("");
    valsd.SetFillColor(15);
    valsd.SetLineColor(15);
    valsd.SetMarkerColor(15);
    valsd.Draw("E2SAME");
    valsd.SetMarkerStyle(20);
    valsd.SetMarkerColor(4);
    canv.Update();
    TLatex latex;
    latex.SetTextAlign(32);
    latex.SetTextSize(0.11);
    for (unsigned i = 0; i < npulls; ++i) {
      TString txt = Form("%i", -1);
      if (impact_map.count(final[i].name)) {
        txt = Form("%i", impact_map.at(final[i].name).rank_impact);
      }
      latex.DrawLatex(0.075, double(i)+0.5, txt);
    }
    pad2->RedrawAxis();
  }
  canv.cd();
  TLatex *title_latex = new TLatex();
  title_latex->SetNDC();
  title_latex->SetTextSize(0.03);
  title_latex->SetTextAlign(31);
  if(!draw_first) title_latex->DrawLatex(0.95,0.965, (std::string(name1) + (splusb_1?"(s+b)":"(b-only)") + " vs. #color[4]{" + name2 + "}" + (splusb_2?"(s+b)":"(b-only)")).c_str());
  else title_latex->DrawLatex(0.95,0.965, (std::string(name1) + (splusb_1?"(s+b)":"(b-only)")).c_str());
  canv.Update();
  if (output == "") output = ("compare_pulls_"+name1+(splusb_1?"_sb_":"_b_")+name2+(splusb_2?"_sb":"_b"));
  canv.SaveAs((output+".pdf").c_str());
  canv.SaveAs((output+".png").c_str());
  return 0;
}

void PullsFromFile(std::string const& filename, std::vector<Pull> & pullvec, bool verbose) {
  // std::vector<std::string> lines = ParseFileLines(filename);
  TFile f(filename.c_str());
  RooFitResult *fit_b = static_cast<RooFitResult*>(gDirectory->Get("fit_b"));
  RooFitResult *fit_s = static_cast<RooFitResult*>(gDirectory->Get("fit_s"));
  if (!fit_s) return;
  RooArgList const& list_s = fit_s->floatParsFinal();
  std::set<std::string> names;
  unsigned n = list_s.getSize();
  for (unsigned i = 0; i < n; ++i) {
    RooRealVar const* var = dynamic_cast<RooRealVar const*>(list_s.at(i));
    if (!var) continue;
    if (std::string(var->GetName()) == "r") continue;
    pullvec.push_back(Pull());
    Pull & new_pull = pullvec.back();
    new_pull.name = var->GetName();
    new_pull.prefit = 0.;
    new_pull.prefit_err = 1.;
    if (fit_b) {
      RooRealVar const* bvar = dynamic_cast<RooRealVar const*>(
          fit_b->floatParsFinal().find(var->GetName()));
      new_pull.bonly = bvar->getVal();
      new_pull.bonly_err = bvar->getError();
    }
    new_pull.splusb = var->getVal();
    new_pull.splusb_err = var->getError();
    new_pull.rho = 0.;
    if (verbose) new_pull.Print();
  }
}

// bool BvsSBComparator(Pull const& pull1, Pull const& pull2) {
//   return (fabs(pull1.bonly - pull1.splusb) > fabs(pull2.bonly - pull2.splusb));
// }

bool BvsSBComparator(Pull const& pull1, Pull const& pull2) {
  double err1 = std::min(pull1.bonly_err, pull1.splusb_err);
  double err2 = std::min(pull2.bonly_err, pull2.splusb_err);

  return (fabs(pull1.bonly - pull1.splusb)/err1 > fabs(pull2.bonly - pull2.splusb)/err2);
}

bool ImpactComparator(Pull const& pull1, Pull const& pull2) {
  return pull1.rho > pull2.rho;
}

std::vector<std::string> ParseFileLines(std::string const& file_name) {
  // Build a vector of input files
  std::vector<std::string> files;
  std::ifstream file;
  file.open(file_name.c_str());
  if (!file.is_open()) {
    std::cerr << "Warning: File " << file_name << " cannot be opened." << std::endl;
    return files;
  }
  std::string line = "";
  while(std::getline(file, line)) {//while loop through lines
    files.push_back(line);
  }
  file.close();
  return files;
}

void Pull::Print() const {
  std::cout << boost::format("%-60s %-4.2f +/- %-4.2f   %+-4.2f +/- %-4.2f   %+-4.2f +/- %-4.2f   %+-4.2f \n")
    % name % prefit %prefit_err % bonly % bonly_err % splusb % splusb_err % rho;
}




