#include <map>
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooWorkspace.h"
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

  TFile fshapes("htt_mt.inputs-mssm-13TeV-mttot.root");
  TH1F *h_ZTT = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT");
  TH1F *h_W = (TH1F*)gDirectory->Get("/mt_nobtag/W");
  TH1F *h_ZTT_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT_CMS_scale_t_mt_13TeVDown");
  TH1F *h_ZTT_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT_CMS_scale_t_mt_13TeVUp");

  TH1F *h_TTT = (TH1F*)gDirectory->Get("/mt_nobtag/TTT");
  TH1F *h_TTJ = (TH1F*)gDirectory->Get("/mt_nobtag/TTJ");
  TH1F *h_QCD = (TH1F*)gDirectory->Get("/mt_nobtag/QCD");
  TH1F *h_VVT = (TH1F*)gDirectory->Get("/mt_nobtag/VVT");
  TH1F *h_VVJ = (TH1F*)gDirectory->Get("/mt_nobtag/VVJ");
  TH1F *h_ZL = (TH1F*)gDirectory->Get("/mt_nobtag/ZL");
  TH1F *h_ZJ = (TH1F*)gDirectory->Get("/mt_nobtag/ZJ");


  // TH1F *h_ZTT160 = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT160");
  // TH1F *h_ZTT160_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT160_CMS_scale_t_mt_13TeVDown");
  // TH1F *h_ZTT160_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT160_CMS_scale_t_mt_13TeVUp");
  // TH1F *h_ZTT180 = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT180");
  // TH1F *h_ZTT180_lo = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT180_CMS_scale_t_mt_13TeVDown");
  // TH1F *h_ZTT180_hi = (TH1F*)gDirectory->Get("/mt_nobtag/ZTT180_CMS_scale_t_mt_13TeVUp");
  // TH1F *h_dat = (TH1F*)h_ZTT->Clone();
  // h_dat->Reset();
  // h_dat->FillRandom(h_ZTT, int(h_ZTT->Integral()));
  // TH1F *h_dat = (TH1F*)gDirectory->Get("/mt_nobtag/data_obs");
  fshapes.Close();
  h_ZTT->Print("range");
  h_W->Print("range");
  // h_dat->Print("range");

  // RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG);
  // RooAbsArg::verboseDirty(true);
  {
    RooRealVar x = VarFromHist("x", "x", *h_ZTT);

    RooRealVar tes("CMS_scale_t_mt_13TeV", "CMS_scale_t_mt_13TeV", 0, -7, 7);
    CMSHistFunc c_ZTT("ZTT", "ZTT", x, *h_ZTT);
    RooArgList binparsx;
    for (int i = 0; i < h_ZTT->GetNbinsX(); ++i) {
      RooRealVar binvar(TString::Format("oldbin%i", i), "", 0, -7 , 7);
      binparsx.addClone(binvar);
    }
    RooArgList vvars(tes);
    // vvars.add(binparsx);
    // c_ZTT.setVerticalMorphs(RooArgList(binparsx));
    c_ZTT.setVerticalMorphs(vvars);
    c_ZTT.prepareStorage();
    c_ZTT.setEvalVerbose(0);
    c_ZTT.setShape(0, 0, 0, 0, *h_ZTT);
    c_ZTT.setShape(0, 0, 1, 0, *h_ZTT_lo);
    c_ZTT.setShape(0, 0, 1, 1, *h_ZTT_hi);

    // for (int i = 0; i < h_ZTT->GetNbinsX(); ++i) {
    //   TH1F* hlo = (TH1F*)h_ZTT->Clone();
    //   hlo->SetBinContent(i+1, h_ZTT->GetBinContent(i+1) - h_ZTT->GetBinError(i+1));
    //   TH1F* hhi = (TH1F*)h_ZTT->Clone();
    //   hhi->SetBinContent(i+1, h_ZTT->GetBinContent(i+1) + h_ZTT->GetBinError(i+1));
    //   c_ZTT.setShape(0, 0, i+2, 0, *hlo);
    //   c_ZTT.setShape(0, 0, i+2, 1, *hhi);
    //   delete hlo;
    //   delete hhi;
    // }

    CMSHistFunc c_W("W", "W", x, *h_W); c_W.prepareStorage(); c_W.setShape(0, 0, 0, 0, *h_W); RooRealVar coeff_W("coeff_W", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_TTJ("TTJ", "TTJ", x, *h_TTJ); c_TTJ.prepareStorage(); c_TTJ.setShape(0, 0, 0, 0, *h_TTJ); RooRealVar coeff_TTJ("coeff_TTJ", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_TTT("TTT", "TTT", x, *h_TTT); c_TTT.prepareStorage(); c_TTT.setShape(0, 0, 0, 0, *h_TTT); RooRealVar coeff_TTT("coeff_TTT", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_VVT("VVT", "VVT", x, *h_VVT); c_VVT.prepareStorage(); c_VVT.setShape(0, 0, 0, 0, *h_VVT); RooRealVar coeff_VVT("coeff_VVT", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_VVJ("VVJ", "VVJ", x, *h_VVJ); c_VVJ.prepareStorage(); c_VVJ.setShape(0, 0, 0, 0, *h_VVJ); RooRealVar coeff_VVJ("coeff_VVJ", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_QCD("QCD", "QCD", x, *h_QCD); c_QCD.prepareStorage(); c_QCD.setShape(0, 0, 0, 0, *h_QCD); RooRealVar coeff_QCD("coeff_QCD", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_ZL("ZL", "ZL", x, *h_ZL); c_ZL.prepareStorage(); c_ZL.setShape(0, 0, 0, 0, *h_ZL); RooRealVar coeff_ZL("coeff_ZL", "", 1.0, 0.5, 1.5);
    CMSHistFunc c_ZJ("ZJ", "ZJ", x, *h_ZJ); c_ZJ.prepareStorage(); c_ZJ.setShape(0, 0, 0, 0, *h_ZJ); RooRealVar coeff_ZJ("coeff_ZJ", "", 1.0, 0.5, 1.5);

    RooRealVar coeff_ZTT("coeff_ZTT", "", 1.0, 0.5, 1.5);
    RooArgList binpars;
    for (int i = 0; i < h_ZTT->GetNbinsX(); ++i) {
      RooRealVar binvar(TString::Format("bin%i", i), "", 0, -7 , 7);
      binpars.addClone(binvar);
    }
    // binpars.Print("v");

    // CMSHistErrorPropagator c_prop(
    //     "prop", "prop", x,
    //     RooArgList(c_ZTT, c_W, c_TTJ, c_TTT, c_VVJ, c_VVT, c_QCD, c_ZL, c_ZJ),
    //     RooArgList(coeff_ZTT, coeff_W, coeff_TTJ, coeff_TTT, coeff_VVJ,
    //                coeff_VVT, coeff_QCD, coeff_ZL, coeff_ZJ));


    CMSHistErrorPropagator c_prop(
        "prop", "prop", x,
        RooArgList(c_ZTT, c_W),
        RooArgList(coeff_ZTT, coeff_W));

    auto newbinargs = c_prop.setupBinPars(0.);
    newbinargs->Print();

    CMSHistFuncWrapper c_ZTT_wrapper("ZTT_wrapper", "ZTT_wrapper", x, c_ZTT, c_prop, 0);
    CMSHistFuncWrapper c_W_wrapper("W_wrapper", "W_wrapper", x, c_W, c_prop, 1);


    c_ZTT.Print("v");
    c_W.Print("v");
    ((RooRealVar*)newbinargs->at(2))->setVal(+2.0);
    c_prop.evaluate();
    c_prop.Print("v");


    // CMSHistFuncWrapper c_TTJ_wrapper("TTJ_wrapper", "TTJ_wrapper", x, c_TTJ, c_prop, 1);
    // CMSHistFuncWrapper c_TTT_wrapper("TTT_wrapper", "TTT_wrapper", x, c_TTT, c_prop, 1);
    // CMSHistFuncWrapper c_VVJ_wrapper("VVJ_wrapper", "VVJ_wrapper", x, c_VVJ, c_prop, 1);
    // CMSHistFuncWrapper c_VVT_wrapper("VVT_wrapper", "VVT_wrapper", x, c_VVT, c_prop, 1);
    // CMSHistFuncWrapper c_QCD_wrapper("QCD_wrapper", "QCD_wrapper", x, c_QCD, c_prop, 1);
    // CMSHistFuncWrapper c_ZL_wrapper("ZL_wrapper", "ZL_wrapper", x, c_ZL, c_prop, 1);
    // CMSHistFuncWrapper c_ZJ_wrapper("ZJ_wrapper", "ZJ_wrapper", x, c_ZJ, c_prop, 1);

    // c_ZTT_wrapper.setEvalVerbose(1);
    // c_ZTT_wrapper.evaluate();
    // x.setVal(3000.);
    // ((RooRealVar &)binpars[3]).setVal(1.2);
    // ((RooRealVar &)binpars[4]).setVal(-1.5);
    // c_ZTT_wrapper.evaluate();
    // c_W_wrapper.evaluate();
    // c_W_wrapper.evaluate();
    // c_W_wrapper.evaluate();
    // c_ZTT_wrapper.evaluate();
    // c_ZTT_wrapper.evaluate();

    // TRandom3 rng;
    // for (unsigned r = 0; r < 1E1; ++r) {
    //   RooRealVar &var = (RooRealVar&)binpars[rng.Integer(binpars.getSize())];
    //   // RooRealVar &var = (RooRealVar&)binparsx[rng.Integer(binparsx.getSize())];
    //   var.setVal(rng.Uniform(-3, +3));
    //   var.Print();

    //   // mH.setVal(rng.Uniform(mH.getMin(), mH.getMax()));
    //   // mH.Print();
    //   // tes.setVal(rng.Gaus(0, 1));
    //   // coeff_ZTT.setVal(rng.Gaus(1, 0.05));
    //   // tes.Print();
    //   c_ZTT_wrapper.evaluate();
    //   c_W_wrapper.evaluate();
    //   // c_ZTT.evaluate();
    //   // c_W.evaluate();
    // }

    // c_ZTT.setEvalVerbose(1);
    wsp.import(c_prop);
    // wsp.import(c_W_wrapper);
  }
  // wsp.Print();
  // ((CMSHistFunc *)wsp.function("ZTT"))->setEvalVerbose(1);
  // wsp.function("ZTT_wrapper")->Print("");
  // wsp.function("ZTT_wrapper")->Print("v");
  // wsp.function("prop")->Print("v");

  wsp.writeToFile("test_ws.root");

  return 0;
}
