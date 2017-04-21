#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooBinning.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFuncWrapper.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"
#include "TRandom3.h"
#include "TVector.h"
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

TH1F * RebinHist(TH1F * hist) {
  // TH1::AddDirectory(0);
  // TH1F* shape = new TH1F("tmp", "tmp", hist->GetNbinsX(), 0.,
  //            static_cast<float>(hist->GetNbinsX()));
  // for (int i = 1; i <= hist->GetNbinsX(); ++i) {
  //   shape->SetBinContent(i, hist->GetBinContent(i));
  // }
  // return shape;
  return hist;
}

int main(int argc, char* argv[]) {
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

  TFile fshapes("htt_mt.inputs-mssm-13TeV-mttot.root");
  TH1F *h_ggH140 = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH140"));
  TH1F *h_ggH140_lo = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVDown"));
  TH1F *h_ggH140_hi = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVUp"));
  TH1F *h_ggH160 = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH160"));
  TH1F *h_ggH160_lo = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVDown"));
  TH1F *h_ggH160_hi = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVUp"));
  TH1F *h_ggH180 = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH180"));
  TH1F *h_ggH180_lo = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH180_CMS_scale_t_mt_13TeVDown"));
  TH1F *h_ggH180_hi = RebinHist((TH1F*)gDirectory->Get("/mt_nobtag/ggH180_CMS_scale_t_mt_13TeVUp"));
  TH1F *h_dat = (TH1F*)h_ggH140->Clone();

  h_dat->Reset();
  h_dat->FillRandom(h_ggH140, int(h_ggH140->Integral()));
  // TH1F *h_dat = (TH1F*)gDirectory->Get("/mt_nobtag/data_obs");
  std::vector<double> rebins = {0., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 250., 300., 350., 400., 500., 700., 1100., 1500., 1900., 2300., 2700., 3100., 3500., 3900.};
  h_dat->Print("range");
  h_dat = (TH1F*)h_dat->Rebin(rebins.size() - 1, "", &(rebins[0]));
  h_dat->Print("range");
  fshapes.Close();
  // h_ggH140->Print("range");
  // h_dat->Print("range");

  RooRealVar x = VarFromHist("x", "x", *h_dat);
  x.Print();
  RooRealVar tes("CMS_scale_t_mt_13TeV", "CMS_scale_t_mt_13TeV", 0, -7, 7);
  CMSHistFunc c_ggH("ggH", "ggH", x, *h_ggH140);

  c_ggH.setHorizontalType(CMSHistFunc::HorizontalType::Integral);

  RooRealVar mH("mH", "mH", 150, 140, 180);
  c_ggH.addHorizontalMorph(mH, TVectorD(3, std::vector<double>{140., 160., 180.}.data()));
  c_ggH.setVerticalMorphs(RooArgList(tes));
  c_ggH.prepareStorage();
  c_ggH.setEvalVerbose(1);

  c_ggH.setShape(0, 0, 0, 0, *h_ggH140);
  c_ggH.setShape(0, 0, 1, 0, *h_ggH140_lo);
  c_ggH.setShape(0, 0, 1, 1, *h_ggH140_hi);
  c_ggH.setShape(0, 1, 0, 0, *h_ggH160);
  c_ggH.setShape(0, 1, 1, 0, *h_ggH160_lo);
  c_ggH.setShape(0, 1, 1, 1, *h_ggH160_hi);
  c_ggH.setShape(0, 2, 0, 0, *h_ggH180);
  c_ggH.setShape(0, 2, 1, 0, *h_ggH180_lo);
  c_ggH.setShape(0, 2, 1, 1, *h_ggH180_hi);

  c_ggH.evaluate();
  // c_ggH.evaluate();
  c_ggH.Print("v");

  // TRandom3 rng;
  // for (unsigned r = 0; r < 3; ++r) {
  //   mH.setVal(rng.Uniform(mH.getMin(), mH.getMax()));
  //   mH.Print();
  //   // tes.setVal(rng.Gaus(0, 1));
  //   // tes.Print();
  //   c_ggH.evaluate();
  //   c_ggH.Print("v");
  // }


  return 0;
}

