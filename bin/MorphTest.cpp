#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMorphingPdf2.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"
#include "TRandom3.h"
#include "TVector.h"
namespace po = boost::program_options;

using namespace std;

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
  TH1F *h_ggH140 = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140");
  TH1F *h_ggH140_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVDown");
  TH1F *h_ggH140_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ggH140_CMS_scale_t_mt_13TeVUp");
  TH1F *h_ggH160 = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160");
  TH1F *h_ggH160_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVDown");
  TH1F *h_ggH160_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ggH160_CMS_scale_t_mt_13TeVUp");
  TH1F *h_dat = (TH1F*)h_ggH140->Clone();
  h_dat->Reset();
  h_dat->FillRandom(h_ggH140, int(h_ggH140->Integral()));
  // TH1F *h_dat = (TH1F*)gDirectory->Get("/mt_nobtag/data_obs");
  fshapes.Close();
  h_ggH140->Print("range");
  h_dat->Print("range");

  RooRealVar x("x", "x", h_ggH140->GetXaxis()->GetXmin(), h_ggH140->GetXaxis()->GetXmax());
  RooRealVar tes("CMS_scale_t_mt_13TeV", "CMS_scale_t_mt_13TeV", 0, -7, 7);
  CMSHistFunc c_ggH("ggH", "ggH", x, *h_ggH140);

  RooRealVar mH("mH", "mH", 150, 140, 160);
  c_ggH.addHorizontalMorph(mH, TVectorD(2, std::vector<double>{140., 160.}.data()));
  c_ggH.setVerticalMorphs(RooArgList(tes));
  c_ggH.prepareStorage();

  c_ggH.setShape(0, 0, 0, 0, *h_ggH140);
  c_ggH.setShape(0, 0, 1, 0, *h_ggH140_lo);
  c_ggH.setShape(0, 0, 1, 1, *h_ggH140_hi);
  c_ggH.setShape(0, 1, 0, 0, *h_ggH160);
  c_ggH.setShape(0, 1, 1, 0, *h_ggH160_lo);
  c_ggH.setShape(0, 1, 1, 1, *h_ggH160_hi);

  c_ggH.evaluate();

  // TH1F* validate = th1fmorph("test", "test", h_ggH140, h_ggH160, 140., 160., mH.getVal(), 1.0);
  // validate->Print("range");
  // TRandom3 rng;
  // for (unsigned i = 0; i < 1E6; ++i) {
  //   tes.setVal(rng.Gaus(0, 1));
  //   c_ggH.evaluate();
  // }

  // RooRealVar y_ztt("Y_ZTT", "", 1, 0.5, 1.5);

  // RooDataHist c_dat("data", "", x, h_dat);

  // RooRealSumPdf pdf("pdf", "", RooArgList(c_ztt), RooArgList(y_ztt), true);
  // // pdf.Print("v");
  // // pdf.Print("tree");
  // RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG);
  // RooAbsArg::verboseDirty(true);
  // x.setVal(5.);
  // std::cout << pdf.getVal() << "\t" << pdf.getVal(x) << "\n";
  // x.setVal(15.);
  // std::cout << pdf.getVal() << "\t" << pdf.getVal(x) << "\n";


  // auto *nll = pdf.createNLL(c_dat, RooFit::Extended(true));
  // nll->Print();
  // y_ztt.setVal(1.1);
  // nll->Print();

  // std::cout << "PDF integral: " << pdf.createIntegral(RooArgSet(x))->getVal() << "\n";


  // c_ztt.Print("v");
  // c_ztt.Print("");
  // c_ztt.Print("v");

  return 0;
}

