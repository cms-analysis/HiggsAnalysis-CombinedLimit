#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooWorkspace.h"
#include "RooMomentMorph.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFuncWrapper.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"
#include "TRandom3.h"
#include "TVector.h"
#include "RooBinning.h"

namespace po = boost::program_options;

using namespace std;

RooRealVar VarFromHist(TString name, TString title, TH1 const& hist) {
  RooBinning * binning;
  if (hist.GetXaxis()->GetXbins()) {
    binning = new RooBinning(hist.GetNbinsX(), hist.GetXaxis()->GetXbins()->GetArray());
  } else {
    binning = new RooBinning(hist.GetNbinsX(), hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
  }
  RooRealVar x(name, title, hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
  // x.getBinning().Print();
  x.setBinning(*binning);
  delete binning;
  return x;
}


int main(int argc, char* argv[]) {
  RooWorkspace wsp("test");
  TH1::AddDirectory(false);
  // Need this to read combine workspaces
  gSystem->Load("libHiggsAnalysisCombinedLimit");

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");

  po::options_description config("Configuration");
  // config.add_options()
  //   ("workspace,w",
  //     po::value<string>(&workspace)->required(),
  //     "The input workspace-containing file [REQUIRED]")
  //   ("datacard,d",
  //     po::value<string>(&datacard),
  //     "The input datacard, only used for rebinning");

  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);


  // TH1F x1("x1", "", 10, 0, 10);
  // TH1F x2("x1", "", 10, 0, 10);

  // x1.SetBinContent(1, 5.);
  // x1.SetBinContent(2, 10.);
  // x1.SetBinContent(3, 5.);

  // x2.SetBinContent(8, 5.);
  // x2.SetBinContent(9, 10.);
  // x2.SetBinContent(10, 5.);

  // RooRealVar x("x", "", 0.);

  // RooDataHist rdh1("rdh1", "", RooArgSet(x), &x1);
  // RooDataHist rdh2("rdh2", "", RooArgSet(x), &x2);

  // RooHistPdf h1("h1", "", RooArgSet(x), rdh1);
  // RooHistPdf h2("h2", "", RooArgSet(x), rdh2);

  // RooRealVar m("m", "", 0, 0, 1);

  // TVectorD mvec(2);
  // mvec[0] = 0.;
  // mvec[1] = 1.;

  // RooMomentMorph morph("morph", "", m, RooArgList(x), RooArgList(h1, h2), mvec);

  // CMSHistFunc c_test("test", "test", x, x1);
  // c_test.addHorizontalMorph(m, mvec);
  // // c_test.setVerticalMorphs(RooArgList(tes));
  // c_test.prepareStorage();
  // c_test.setShape(0, 0, 0, 0, x1);
  // c_test.setShape(0, 1, 0, 0, x2);
  // c_test.setEvalVerbose(0);

  // m.setVal(0.5);

  // for (unsigned i = 0; i < 100; ++i) {
  //   x.setVal(double(i) / 10.);
  //   std::cout << x.getVal() << "\t" << morph.getVal() << "\t" << c_test.getVal() << "\n";
  // }


  // c_test.evaluate();


  TFile fshapes("htt_mt.inputs-mssm-13TeV-mttot.root");
  TH1F *h_ggH140 = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140");
  // TH1F *h_ggH140_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVDown");
  // TH1F *h_ggH140_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVUp");
  TH1F *h_ggH160 = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160");
  // TH1F *h_ggH160_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVDown");
  // TH1F *h_ggH160_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVUp");
  TH1F *h_ggH180 = (TH1F*)gDirectory->Get("/mt_nobtag/ggH180");
  // TH1F *h_ggH180_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ggH180_CMS_scale_t_mt_13TeVDown");
  // TH1F *h_ggH180_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ggH180_CMS_scale_t_mt_13TeVUp");
  // TH1F *h_dat = (TH1F*)h_ggH140->Clone();
  // h_dat->Reset();
  // h_dat->FillRandom(h_ggH140, int(h_ggH140->Integral()));
  // TH1F *h_dat = (TH1F*)gDirectory->Get("/mt_nobtag/data_obs");
  fshapes.Close();
  h_ggH140->Print("range");
  // h_dat->Print("range");

  RooRealVar x = VarFromHist("x", "x", *h_ggH140);
  RooRealVar tes("CMS_scale_t_mt_13TeV", "CMS_scale_t_mt_13TeV", 0, -7, 7);
  CMSHistFunc c_ggH("ggH", "ggH", x, *h_ggH140);

  RooRealVar mH("mH", "mH", 150, 140, 180);
  c_ggH.addHorizontalMorph(mH, TVectorD(3, std::vector<double>{140., 160., 180.}.data()));
  c_ggH.setHorizontalType(CMSHistFunc::HorizontalType::Moment);
  c_ggH.setMomentType(CMSHistFunc::MomentSetting::NonLinearPosFractions);
  // c_ggH.setVerticalMorphs(RooArgList(tes));
  c_ggH.prepareStorage();
  // c_ggH.setEvalVerbose(0);

  c_ggH.setShape(0, 0, 0, 0, *h_ggH140);
  // c_ggH.setShape(0, 0, 1, 0, *h_ggH140_lo);
  // c_ggH.setShape(0, 0, 1, 1, *h_ggH140_hi);
  c_ggH.setShape(0, 1, 0, 0, *h_ggH160);
  // c_ggH.setShape(0, 1, 1, 0, *h_ggH160_lo);
  // c_ggH.setShape(0, 1, 1, 1, *h_ggH160_hi);
  c_ggH.setShape(0, 2, 0, 0, *h_ggH180);
  // c_ggH.setShape(0, 2, 1, 0, *h_ggH180_lo);
  // c_ggH.setShape(0, 2, 1, 1, *h_ggH180_hi);


  RooDataHist rdh1("rdh1", "", RooArgSet(x), h_ggH140);
  RooDataHist rdh2("rdh2", "", RooArgSet(x), h_ggH160);
  RooDataHist rdh3("rdh3", "", RooArgSet(x), h_ggH180);

  RooHistPdf h1("h1", "", RooArgSet(x), rdh1);
  RooHistPdf h2("h2", "", RooArgSet(x), rdh2);
  RooHistPdf h3("h3", "", RooArgSet(x), rdh3);

  RooMomentMorph morph(
      "morph", "", mH, RooArgList(x), RooArgList(h1, h2, h3),
      TVectorD(3, std::vector<double>{140., 160., 180.}.data()));
  morph.setMode(RooMomentMorph::Setting::NonLinearPosFractions);
  // c_ggH.evaluate();
  // mH.setVal(170);
  c_ggH.evaluate();
  RooWorkspace w("w", "");

  // c_ggH.setEvalVerbose(0);

  w.import(c_ggH);

  TRandom3 rng;
  for (unsigned r = 0; r < 1E7; ++r) {
    mH.setVal(rng.Uniform(mH.getMin(), mH.getMax()));
    x.setVal(rng.Uniform(x.getMin(), x.getMax()));
    // mH.Print();
    // tes.setVal(rng.Gaus(0, 1));
    // tes.Print();
    // w.function("ggH")->getVal();
    // c_ggH.getVal();
    morph.getVal();
  }


  // w.writeToFile("workspace.root");


  // h_ggH140->Print("range");
  // h_ggH160->Print("range");

  return 0;
}
