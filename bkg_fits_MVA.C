// Fits to background samples
//  to extract parameters
// Aron Soha
// April 5, 2018
// Sep 19, 2018 : This version now uses the DeepESM-based samples
// Nov 15, 2018 : Updated to use RooParametricHist and the latest MVA

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH2.h"

using namespace RooFit;

void construct_formula(string procName, RooArgList& binlist, const RooArgList& paramlist) {

  int max_bin = 20; // 14 means just njets=14, 20 means last bin is inclusive up through njets=20

  for (int i=1; i<=8; i++) {

    stringstream form;
    RooArgList formArgList;

    form << "(@0";
    formArgList.add(paramlist[0]); // N7_tt for this MVA bin

    if (i>=2) {
      form << "*@1";
      formArgList.add(paramlist[1]); // p0_tt
    }

    if (i>=3) { // for bin 3 and up
      for (int j=3; j<=i; j++) {
	form << "*(@2+(@1-@2)*exp((" << j << "-2)*@3))";
      }
      formArgList.add(paramlist[2]); // p1_tt
      formArgList.add(paramlist[3]); // p2_tt
    } // end bin 3 and up

    // The last bin covers from njet=14 through njet=max_bin
    if (i==8) {
      for (int k=9; k<=max_bin-6; k++) {
	form << " + @0*@1";
	for (int j=3; j<=k; j++) {
	  form << "*(@2+(@1-@2)*exp((" << j << "-2)*@3))";
	}
      }
    }
    form << ")";

    // Create RooFormulaVar for this bin
    stringstream binName;
    binName << procName << "_b" << i;
    RooFormulaVar* binvar = new RooFormulaVar(binName.str().c_str(), "", form.str().c_str(), RooArgList(formArgList));
    binlist.add(*binvar);

    cout << "bin i = " << i << " , njets = " << i+6 << endl;
    cout << "process bin name : " << binName.str().c_str() << endl;
    cout << "Formula : " << form.str().c_str() << endl;
    formArgList.Print();
    cout << endl;

  }

}


void bkg_fits_MVA() {

  RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
  //wspace->factory("nj[6.5,14.5]");
  //wspace->var("nj")->setBins(8);
  //RooArgSet vars(*wspace->var("nj"));

  // njet is our variable, 8 bins, from 7 up through 14,
  //   Note that njet=14 is inclusive as >=14
  // D1, D2, D3, D4 are the MVA bins
  wspace->factory("nj_D1[6.5,14.5]");
  //wspace->factory("nj_D1[0,8]");
  wspace->var("nj_D1")->setBins(8);
  RooArgSet vars_D1(*wspace->var("nj_D1"));

  wspace->factory("nj_D2[6.5,14.5]");
  //wspace->factory("nj_D2[0,8]");
  wspace->var("nj_D2")->setBins(8);
  RooArgSet vars_D2(*wspace->var("nj_D2"));

  wspace->factory("nj_D3[6.5,14.5]");
  //wspace->factory("nj_D3[0,8]");
  wspace->var("nj_D3")->setBins(8);
  RooArgSet vars_D3(*wspace->var("nj_D3"));

  wspace->factory("nj_D4[6.5,14.5]");
  //wspace->factory("nj_D4[0,8]");
  wspace->var("nj_D4")->setBins(8);
  RooArgSet vars_D4(*wspace->var("nj_D4"));

  wspace->var("nj_D1")->setRange("low",6.5,11.5);
  wspace->var("nj_D2")->setRange("low",6.5,11.5);
  wspace->var("nj_D3")->setRange("low",6.5,11.5);
  wspace->var("nj_D4")->setRange("low",6.5,11.5);

  TFile* file = TFile::Open("njets_for_Aron_V1.2.3_Nov19.root");

  // tt ---------------------------------------------------------------------------------------

  // D1
  TH1* ttMC_th1_D1 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin1_TT",ttMC_th1_D1);

  RooRealVar p0_tt_D1("p0_tt_D1","p0 of tt bkg shape",0.35,0.00,1.00);
  RooRealVar p1_tt_D1("p1_tt_D1","p1 of tt bkg shape",0.21,-1.00,1.00);
  //RooRealVar p1_tt_D1("p1_tt_D1","p1 of tt bkg shape",0.21,0.00,1.00);
  RooRealVar p2_tt_D1("p2_tt_D1","p2 of tt bkg shape",-0.25,-1.00,0.0);

  //RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",ttMC_th1_D1->GetBinContent(1),ttMC_th1_D1->GetBinContent(1)-5000,ttMC_th1_D1->GetBinContent(1)+5000);
  //RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",ttMC_th1_D1->GetBinContent(1),ttMC_th1_D1->GetBinContent(1),ttMC_th1_D1->GetBinContent(1));
  RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",1);
  
  RooDataHist ttMC_hist_D1("ttMC_obs_D1","tt MC observed in signal region D1",vars_D1,ttMC_th1_D1);
  wspace->import(ttMC_hist_D1);
  // shape for tt bkg MC
  RooArgList parlist_D1(N7_tt_D1,p0_tt_D1,p1_tt_D1,p2_tt_D1);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  string procName_D1 = "background_tt_D1";
  construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1);
  RooParametricHist background_tt_D1(procName_D1.c_str(),"",*wspace->var("nj_D1"),*bkg_tt_bins_D1,*ttMC_th1_D1);
  wspace->import(background_tt_D1,RooFit::RecycleConflictNodes());

  RooAbsPdf* model_tt_D1 = wspace->pdf("background_tt_D1");
  RooFitResult* r_tt_D1 = model_tt_D1->fitTo(ttMC_hist_D1,Save(),SumW2Error(kFALSE),Minos(kTRUE),Range("low"));
  //RooFitResult* r_tt_D1 = model_tt_D1->fitTo(ttMC_hist_D1,Save(),SumW2Error(kFALSE),Range("low"));
  //RooFitResult* r_tt_D1 = model_tt_D1->fitTo(ttMC_hist_D1,Save(),SumW2Error(kFALSE),Minos(kTRUE));
  //RooFitResult* r_tt_D1 = model_tt_D1->fitTo(ttMC_hist_D1,Save(),SumW2Error(kFALSE));
  r_tt_D1->Print("V");

  RooPlot* ttframe_D1 = wspace->var("nj_D1")->frame(Title("njets tt D1 PDF"));
  ttMC_hist_D1.plotOn(ttframe_D1);
  model_tt_D1->plotOn(ttframe_D1,Precision(0.00001));
  model_tt_D1->paramOn(ttframe_D1,Layout(0.35,0.622912,0.773497));
  ttframe_D1->getAttLine()->SetLineWidth(0);

  //new TCanvas("tt D1 fit result","tt D1 fit results",600,600);
  //gPad->SetLeftMargin(0.15);
  //ttframe_D1->GetYaxis()->SetTitleOffset(2.0);
  //ttframe_D1->Draw();

  // D2
  TH1* ttMC_th1_D2 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin2_TT",ttMC_th1_D2);

  RooRealVar p0_tt_D2("p0_tt_D2","p0 of tt bkg shape",0.35,0.0,1.00);
  RooRealVar p1_tt_D2("p1_tt_D2","p1 of tt bkg shape",0.21,-1.00,1.00);
  //RooRealVar p1_tt_D2("p1_tt_D2","p1 of tt bkg shape",0.21,0.00,1.00);
  RooRealVar p2_tt_D2("p2_tt_D2","p2 of tt bkg shape",-0.25,-1.00,0.0);

  //RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",ttMC_th1_D2->GetBinContent(1),ttMC_th1_D2->GetBinContent(1)-5000,ttMC_th1_D2->GetBinContent(1)+5000);
  RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",1);
  
  RooDataHist ttMC_hist_D2("ttMC_obs_D2","tt MC observed in signal region D2",vars_D2,ttMC_th1_D2);
  wspace->import(ttMC_hist_D2);
  // shape for tt bkg MC
  RooArgList parlist_D2(N7_tt_D2,p0_tt_D2,p1_tt_D2,p2_tt_D2);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  string procName_D2 = "background_tt_D2";
  construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2);
  RooParametricHist background_tt_D2(procName_D2.c_str(),"",*wspace->var("nj_D2"),*bkg_tt_bins_D2,*ttMC_th1_D2);
  wspace->import(background_tt_D2,RooFit::RecycleConflictNodes());

  RooAbsPdf* model_tt_D2 = wspace->pdf("background_tt_D2");
  RooFitResult* r_tt_D2 = model_tt_D2->fitTo(ttMC_hist_D2,Save(),SumW2Error(kFALSE),Minos(kTRUE),Range("low"));
  //RooFitResult* r_tt_D2 = model_tt_D2->fitTo(ttMC_hist_D2,Save(),SumW2Error(kFALSE),Range("low"));
  //RooFitResult* r_tt_D2 = model_tt_D2->fitTo(ttMC_hist_D2,Save(),SumW2Error(kFALSE),Minos(kTRUE));
  //RooFitResult* r_tt_D2 = model_tt_D2->fitTo(ttMC_hist_D2,Save(),SumW2Error(kFALSE));
  r_tt_D2->Print("V");

  RooPlot* ttframe_D2 = wspace->var("nj_D2")->frame(Title("njets tt D2 PDF"));
  ttMC_hist_D2.plotOn(ttframe_D2);
  model_tt_D2->plotOn(ttframe_D2,Precision(0.00001));
  model_tt_D2->paramOn(ttframe_D2,Layout(0.35,0.622912,0.773497));
  ttframe_D2->getAttLine()->SetLineWidth(0);

  //new TCanvas("tt D2 fit result","tt D2 fit results",600,600);
  //gPad->SetLeftMargin(0.15);
  //ttframe_D2->GetYaxis()->SetTitleOffset(2.0);
  //ttframe_D2->Draw();

  // D3
  TH1* ttMC_th1_D3 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin3_TT",ttMC_th1_D3);

  RooRealVar p0_tt_D3("p0_tt_D3","p0 of tt bkg shape",0.35,0.0,1.00);
  RooRealVar p1_tt_D3("p1_tt_D3","p1 of tt bkg shape",0.21,-1.00,1.00);
  // RooRealVar p1_tt_D3("p1_tt_D3","p1 of tt bkg shape",0.21,0.00,1.00);
  RooRealVar p2_tt_D3("p2_tt_D3","p2 of tt bkg shape",-0.25,-1.00,0.0);

  //RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",ttMC_th1_D3->GetBinContent(1),ttMC_th1_D3->GetBinContent(1)-5000,ttMC_th1_D3->GetBinContent(1)+5000);
  RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",1);
  
  RooDataHist ttMC_hist_D3("ttMC_obs_D3","tt MC observed in signal region D3",vars_D3,ttMC_th1_D3);
  wspace->import(ttMC_hist_D3);
  // shape for tt bkg MC
  RooArgList parlist_D3(N7_tt_D3,p0_tt_D3,p1_tt_D3,p2_tt_D3);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  string procName_D3 = "background_tt_D3";
  construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3);
  RooParametricHist background_tt_D3(procName_D3.c_str(),"",*wspace->var("nj_D3"),*bkg_tt_bins_D3,*ttMC_th1_D3);
  wspace->import(background_tt_D3,RooFit::RecycleConflictNodes());

  RooAbsPdf* model_tt_D3 = wspace->pdf("background_tt_D3");
  RooFitResult* r_tt_D3 = model_tt_D3->fitTo(ttMC_hist_D3,Save(),SumW2Error(kFALSE),Minos(kTRUE),Range("low"));
  //RooFitResult* r_tt_D3 = model_tt_D3->fitTo(ttMC_hist_D3,Save(),SumW2Error(kFALSE),Range("low"));
  //RooFitResult* r_tt_D3 = model_tt_D3->fitTo(ttMC_hist_D3,Save(),SumW2Error(kFALSE),Minos(kTRUE));
  //RooFitResult* r_tt_D3 = model_tt_D3->fitTo(ttMC_hist_D3,Save(),SumW2Error(kFALSE));
  r_tt_D3->Print("V");

  RooPlot* ttframe_D3 = wspace->var("nj_D3")->frame(Title("njets tt D3 PDF"));
  ttMC_hist_D3.plotOn(ttframe_D3);
  model_tt_D3->plotOn(ttframe_D3,Precision(0.00001));
  model_tt_D3->paramOn(ttframe_D3,Layout(0.35,0.622912,0.773497));
  ttframe_D3->getAttLine()->SetLineWidth(0);

  //new TCanvas("tt D3 fit result","tt D3 fit results",600,600);
  //gPad->SetLeftMargin(0.15);
  //ttframe_D3->GetYaxis()->SetTitleOffset(2.0);
  //ttframe_D3->Draw();

  // D4
  TH1* ttMC_th1_D4 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin4_TT",ttMC_th1_D4);

  RooRealVar p0_tt_D4("p0_tt_D4","p0 of tt bkg shape",0.35,0.0,1.00);
  RooRealVar p1_tt_D4("p1_tt_D4","p1 of tt bkg shape",0.21,-1.00,1.00);
  //RooRealVar p1_tt_D4("p1_tt_D4","p1 of tt bkg shape",0.21,0.00,1.00);
  RooRealVar p2_tt_D4("p2_tt_D4","p2 of tt bkg shape",-0.25,-1.00,0.0);

  //RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",ttMC_th1_D4->GetBinContent(1),ttMC_th1_D4->GetBinContent(1)-5000,ttMC_th1_D4->GetBinContent(1)+5000);
  RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",1);
  
  RooDataHist ttMC_hist_D4("ttMC_obs_D4","tt MC observed in signal region D4",vars_D4,ttMC_th1_D4);
  wspace->import(ttMC_hist_D4);
  // shape for tt bkg MC
  RooArgList parlist_D4(N7_tt_D4,p0_tt_D4,p1_tt_D4,p2_tt_D4);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  string procName_D4 = "background_tt_D4";
  construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4);
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("nj_D4"),*bkg_tt_bins_D4,*ttMC_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());

  RooAbsPdf* model_tt_D4 = wspace->pdf("background_tt_D4");
  RooFitResult* r_tt_D4 = model_tt_D4->fitTo(ttMC_hist_D4,Save(),SumW2Error(kFALSE),Minos(kTRUE),Range("low"));
  //RooFitResult* r_tt_D4 = model_tt_D4->fitTo(ttMC_hist_D4,Save(),SumW2Error(kFALSE),Range("low"));
  //RooFitResult* r_tt_D4 = model_tt_D4->fitTo(ttMC_hist_D4,Save(),SumW2Error(kFALSE),Minos(kTRUE));
  //RooFitResult* r_tt_D4 = model_tt_D4->fitTo(ttMC_hist_D4,Save(),SumW2Error(kFALSE));
  r_tt_D4->Print("V");

  RooPlot* ttframe_D4 = wspace->var("nj_D4")->frame(Title("njets tt D4 PDF"));
  ttMC_hist_D4.plotOn(ttframe_D4);
  model_tt_D4->plotOn(ttframe_D4,Precision(0.00001));
  model_tt_D4->paramOn(ttframe_D4,Layout(0.35,0.622912,0.773497));
  ttframe_D4->getAttLine()->SetLineWidth(0);

  //new TCanvas("tt D4 fit result","tt D4 fit results",600,600);
  //gPad->SetLeftMargin(0.15);
  //ttframe_D4->GetYaxis()->SetTitleOffset(2.0);
  //ttframe_D4->Draw();


  // other ---------------------------------------------------------------------------------------

  // // D1
  // TH1* otherMC_th1_D1 = 0;
  // file->GetObject("h_njets_pt30_1l_deepESMbin1_other",otherMC_th1_D1);

  // wspace->factory("p0_other_D1[0.24,0.0,1.5]");
  // wspace->factory("p1_other_D1[0.27,0.0,0.42]");
  // wspace->factory("p2_other_D1[-0.64,-2.5,0.0]");

  // wspace->factory(TString::Format("N7_other_D1[%f,%f,%f]",otherMC_th1_D1->GetBinContent(1),otherMC_th1_D1->GetBinContent(1)-5000,otherMC_th1_D1->GetBinContent(1)+5000));

  // RooDataHist otherMC_hist_D1("otherMC_obs_D1","other MC observed in signal region 1",vars_D1,otherMC_th1_D1);
  // wspace->import(otherMC_hist_D1);
  // // shape for other bkg MC
  // wspace->factory("RooNjetsPdf::bkgMC_other_D1(nj_D1,p0_other_D1,p1_other_D1,p2_other_D1,N7_other_D1)");
  // wspace->import(*wspace->pdf("bkgMC_other_D1"));

  // RooAbsPdf* model_other_D1 = wspace->pdf("bkgMC_other_D1");
  // RooFitResult* r_other_D1 = model_other_D1->fitTo(otherMC_hist_D1,Save(),SumW2Error(kFALSE));
  // r_other_D1->Print("V");

  // RooPlot* otherframe_D1 = wspace->var("nj_D1")->frame(Title("njets other D1 PDF"));
  // otherMC_hist_D1.plotOn(otherframe_D1);
  // model_other_D1->plotOn(otherframe_D1);
  // model_other_D1->paramOn(otherframe_D1,Layout(0.35,0.622912,0.773497));
  // otherframe_D1->getAttLine()->SetLineWidth(0);

  // //new TCanvas("other D1 fit result","other D1 fit results",600,600);
  // //gPad->SetLeftMargin(0.15);
  // //otherframe_D1->GetYaxis()->SetTitleOffset(2.0);
  // //otherframe_D1->Draw();

  // // D2
  // TH1* otherMC_th1_D2 = 0;
  // file->GetObject("h_njets_pt30_1l_deepESMbin2_other",otherMC_th1_D2);

  // wspace->factory("p0_other_D2[0.24,0.0,1.5]");
  // wspace->factory("p1_other_D2[0.27,0.0,0.42]");
  // wspace->factory("p2_other_D2[-0.64,-2.5,0.0]");

  // wspace->factory(TString::Format("N7_other_D2[%f,%f,%f]",otherMC_th1_D2->GetBinContent(1),otherMC_th1_D2->GetBinContent(1)-5000,otherMC_th1_D2->GetBinContent(1)+5000));
  
  // RooDataHist otherMC_hist_D2("otherMC_obs_D2","other MC observed in signal region 2",vars_D2,otherMC_th1_D2);
  // wspace->import(otherMC_hist_D2);
  // // shape for other bkg MC
  // wspace->factory("RooNjetsPdf::bkgMC_other_D2(nj_D2,p0_other_D2,p1_other_D2,p2_other_D2,N7_other_D2)");
  // wspace->import(*wspace->pdf("bkgMC_other_D2"));

  // RooAbsPdf* model_other_D2 = wspace->pdf("bkgMC_other_D2");
  // RooFitResult* r_other_D2 = model_other_D2->fitTo(otherMC_hist_D2,Save(),SumW2Error(kFALSE));
  // r_other_D2->Print("V");

  // RooPlot* otherframe_D2 = wspace->var("nj_D2")->frame(Title("njets other D2 PDF"));
  // otherMC_hist_D2.plotOn(otherframe_D2);
  // model_other_D2->plotOn(otherframe_D2);
  // model_other_D2->paramOn(otherframe_D2,Layout(0.35,0.622912,0.773497));
  // otherframe_D2->getAttLine()->SetLineWidth(0);

  // //new TCanvas("other D2 fit result","other D2 fit results",600,600);
  // //gPad->SetLeftMargin(0.15);
  // //otherframe_D2->GetYaxis()->SetTitleOffset(2.0);
  // //otherframe_D2->Draw();

  // // D3
  // TH1* otherMC_th1_D3 = 0;
  // file->GetObject("h_njets_pt30_1l_deepESMbin3_other",otherMC_th1_D3);

  // wspace->factory("p0_other_D3[0.24,0.0,1.5]");
  // wspace->factory("p1_other_D3[0.27,0.0,0.42]");
  // wspace->factory("p2_other_D3[-0.64,-2.5,0.0]");

  // wspace->factory(TString::Format("N7_other_D3[%f,%f,%f]",otherMC_th1_D3->GetBinContent(1),otherMC_th1_D3->GetBinContent(1)-5000,otherMC_th1_D3->GetBinContent(1)+5000));

  // RooDataHist otherMC_hist_D3("otherMC_obs_D3","other MC observed in signal region 3",vars_D3,otherMC_th1_D3);
  // wspace->import(otherMC_hist_D3);
  // // shape for other bkg MC
  // wspace->factory("RooNjetsPdf::bkgMC_other_D3(nj_D3,p0_other_D3,p1_other_D3,p2_other_D3,N7_other_D3)");
  // wspace->import(*wspace->pdf("bkgMC_other_D3"));

  // RooAbsPdf* model_other_D3 = wspace->pdf("bkgMC_other_D3");
  // RooFitResult* r_other_D3 = model_other_D3->fitTo(otherMC_hist_D3,Save(),SumW2Error(kFALSE));
  // r_other_D3->Print("V");

  // RooPlot* otherframe_D3 = wspace->var("nj_D3")->frame(Title("njets other D3 PDF"));
  // otherMC_hist_D3.plotOn(otherframe_D3);
  // model_other_D3->plotOn(otherframe_D3);
  // model_other_D3->paramOn(otherframe_D3,Layout(0.35,0.622912,0.773497));
  // otherframe_D3->getAttLine()->SetLineWidth(0);

  // //new TCanvas("other D3 fit result","other D3 fit results",600,600);
  // //gPad->SetLeftMargin(0.15);
  // //otherframe_D3->GetYaxis()->SetTitleOffset(2.0);
  // //otherframe_D3->Draw();

  // // D4
  // TH1* otherMC_th1_D4 = 0;
  // file->GetObject("h_njets_pt30_1l_deepESMbin4_other",otherMC_th1_D4);

  // wspace->factory("p0_other_D4[0.24,0.0,1.5]");
  // wspace->factory("p1_other_D4[0.27,0.0,0.42]");
  // wspace->factory("p2_other_D4[-0.64,-2.5,0.0]");

  // wspace->factory(TString::Format("N7_other_D4[%f,%f,%f]",otherMC_th1_D4->GetBinContent(1),otherMC_th1_D4->GetBinContent(1)-5000,otherMC_th1_D4->GetBinContent(1)+5000));
  
  // RooDataHist otherMC_hist_D4("otherMC_obs_D4","other MC observed in signal region 4",vars_D4,otherMC_th1_D4);
  // wspace->import(otherMC_hist_D4);
  // // shape for other bkg MC
  // wspace->factory("RooNjetsPdf::bkgMC_other_D4(nj_D4,p0_other_D4,p1_other_D4,p2_other_D4,N7_other_D4)");
  // wspace->import(*wspace->pdf("bkgMC_other_D4"));

  // RooAbsPdf* model_other_D4 = wspace->pdf("bkgMC_other_D4");
  // RooFitResult* r_other_D4 = model_other_D4->fitTo(otherMC_hist_D4,Save(),SumW2Error(kFALSE));
  // r_other_D4->Print("V");

  // RooPlot* otherframe_D4 = wspace->var("nj_D4")->frame(Title("njets other D4 PDF"));
  // otherMC_hist_D4.plotOn(otherframe_D4);
  // model_other_D4->plotOn(otherframe_D4);
  // model_other_D4->paramOn(otherframe_D4,Layout(0.35,0.622912,0.773497));
  // otherframe_D4->getAttLine()->SetLineWidth(0);

  // //new TCanvas("other D4 fit result","other D4 fit results",600,600);
  // //gPad->SetLeftMargin(0.15);
  // //otherframe_D4->GetYaxis()->SetTitleOffset(2.0);
  // //otherframe_D4->Draw();

  // ===================================================================================

  // Clone the fit result frames, so that I can plot a log version without the fit params
  RooPlot* clone_ttframe_D1 = (RooPlot*)ttframe_D1->Clone("ttframe_D1");
  RooPlot* clone_ttframe_D2 = (RooPlot*)ttframe_D2->Clone("ttframe_D2");
  RooPlot* clone_ttframe_D3 = (RooPlot*)ttframe_D3->Clone("ttframe_D3");
  RooPlot* clone_ttframe_D4 = (RooPlot*)ttframe_D4->Clone("ttframe_D4");

  // RooPlot* clone_otherframe_D1 = (RooPlot*)otherframe_D1->Clone("otherframe_D1");
  // RooPlot* clone_otherframe_D2 = (RooPlot*)otherframe_D2->Clone("otherframe_D2");
  // RooPlot* clone_otherframe_D3 = (RooPlot*)otherframe_D3->Clone("otherframe_D3");
  // RooPlot* clone_otherframe_D4 = (RooPlot*)otherframe_D4->Clone("otherframe_D4");

  // ===================================================================================

  // Get the N7 for each fit

  TPaveText *text_N7_tt_D1 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  text_N7_tt_D1->AddText(TString::Format("N7 = %.1f",ttMC_th1_D1->GetBinContent(1)));
  text_N7_tt_D1->SetFillStyle(0);
  text_N7_tt_D1->SetBorderSize(1);
  text_N7_tt_D1->SetLineColor(0);
  TPaveText *text_N7_tt_D2 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  text_N7_tt_D2->AddText(TString::Format("N7 = %.1f",ttMC_th1_D2->GetBinContent(1)));
  text_N7_tt_D2->SetFillStyle(0);
  text_N7_tt_D2->SetBorderSize(1);
  text_N7_tt_D2->SetLineColor(0);
  TPaveText *text_N7_tt_D3 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  text_N7_tt_D3->AddText(TString::Format("N7 = %.1f",ttMC_th1_D3->GetBinContent(1)));
  text_N7_tt_D3->SetFillStyle(0);
  text_N7_tt_D3->SetBorderSize(1);
  text_N7_tt_D3->SetLineColor(0);
  TPaveText *text_N7_tt_D4 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  text_N7_tt_D4->AddText(TString::Format("N7 = %.1f",ttMC_th1_D4->GetBinContent(1)));
  text_N7_tt_D4->SetFillStyle(0);
  text_N7_tt_D4->SetBorderSize(1);
  text_N7_tt_D4->SetLineColor(0);

  // TPaveText *text_N7_other_D1 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  // text_N7_other_D1->AddText(TString::Format("N7 = %.1f",otherMC_th1_D1->GetBinContent(1)));
  // text_N7_other_D1->SetFillStyle(0);
  // text_N7_other_D1->SetBorderSize(1);
  // text_N7_other_D1->SetLineColor(0);
  // TPaveText *text_N7_other_D2 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  // text_N7_other_D2->AddText(TString::Format("N7 = %.1f",otherMC_th1_D2->GetBinContent(1)));
  // text_N7_other_D2->SetFillStyle(0);
  // text_N7_other_D2->SetBorderSize(1);
  // text_N7_other_D2->SetLineColor(0);
  // TPaveText *text_N7_other_D3 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  // text_N7_other_D3->AddText(TString::Format("N7 = %.1f",otherMC_th1_D3->GetBinContent(1)));
  // text_N7_other_D3->SetFillStyle(0);
  // text_N7_other_D3->SetBorderSize(1);
  // text_N7_other_D3->SetLineColor(0);
  // TPaveText *text_N7_other_D4 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  // text_N7_other_D4->AddText(TString::Format("N7 = %.1f",otherMC_th1_D4->GetBinContent(1)));
  // text_N7_other_D4->SetFillStyle(0);
  // text_N7_other_D4->SetBorderSize(1);
  // text_N7_other_D4->SetLineColor(0);



  // ===================================================================================

  // Get the chi^2 for each fit

  TPaveText *chiSquare_tt_D1 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_tt_D1->AddText(TString::Format("Chi Square = %0.2f",ttframe_D1->chiSquare(4)));
  chiSquare_tt_D1->SetFillStyle(0);
  chiSquare_tt_D1->SetBorderSize(1);
  chiSquare_tt_D1->SetLineColor(0);
  TPaveText *chiSquare_tt_D2 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_tt_D2->AddText(TString::Format("Chi Square = %0.2f",ttframe_D2->chiSquare(4)));
  chiSquare_tt_D2->SetFillStyle(0);
  chiSquare_tt_D2->SetBorderSize(1);
  chiSquare_tt_D2->SetLineColor(0);
  TPaveText *chiSquare_tt_D3 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_tt_D3->AddText(TString::Format("Chi Square = %0.2f",ttframe_D3->chiSquare(4)));
  chiSquare_tt_D3->SetFillStyle(0);
  chiSquare_tt_D3->SetBorderSize(1);
  chiSquare_tt_D3->SetLineColor(0);
  TPaveText *chiSquare_tt_D4 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_tt_D4->AddText(TString::Format("Chi Square = %0.2f",ttframe_D4->chiSquare(4)));
  chiSquare_tt_D4->SetFillStyle(0);
  chiSquare_tt_D4->SetBorderSize(1);
  chiSquare_tt_D4->SetLineColor(0);

  // TPaveText *chiSquare_other_D1 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  // chiSquare_other_D1->AddText(TString::Format("Chi Square = %0.2f",otherframe_D1->chiSquare(4)));
  // chiSquare_other_D1->SetFillStyle(0);
  // chiSquare_other_D1->SetBorderSize(1);
  // chiSquare_other_D1->SetLineColor(0);
  // TPaveText *chiSquare_other_D2 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  // chiSquare_other_D2->AddText(TString::Format("Chi Square = %0.2f",otherframe_D2->chiSquare(4)));
  // chiSquare_other_D2->SetFillStyle(0);
  // chiSquare_other_D2->SetBorderSize(1);
  // chiSquare_other_D2->SetLineColor(0);
  // TPaveText *chiSquare_other_D3 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  // chiSquare_other_D3->AddText(TString::Format("Chi Square = %0.2f",otherframe_D3->chiSquare(4)));
  // chiSquare_other_D3->SetFillStyle(0);
  // chiSquare_other_D3->SetBorderSize(1);
  // chiSquare_other_D3->SetLineColor(0);
  // TPaveText *chiSquare_other_D4 = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  // chiSquare_other_D4->AddText(TString::Format("Chi Square = %0.2f",otherframe_D4->chiSquare(4)));
  // chiSquare_other_D4->SetFillStyle(0);
  // chiSquare_other_D4->SetBorderSize(1);
  // chiSquare_other_D4->SetLineColor(0);


  // ===================================================================================

  // Make pull plots of each fit
  // The pull is defined as (curve-histogram +/- error_on_histogram) / error_on_histogram
  //   so it is expected that for the case of symmetric errors the error bars on the pull are +/- 1. 

  RooHist* hpull_tt_D1 = ttframe_D1->pullHist();
  RooPlot* framePull_tt_D1 = wspace->var("nj_D1")->frame(Title("Pull"));
  framePull_tt_D1->addPlotable(hpull_tt_D1,"P");

  RooHist* hpull_tt_D2 = ttframe_D2->pullHist();
  RooPlot* framePull_tt_D2 = wspace->var("nj_D2")->frame(Title("Pull"));
  framePull_tt_D2->addPlotable(hpull_tt_D2,"P");

  RooHist* hpull_tt_D3 = ttframe_D3->pullHist();
  RooPlot* framePull_tt_D3 = wspace->var("nj_D3")->frame(Title("Pull"));
  framePull_tt_D3->addPlotable(hpull_tt_D3,"P");

  RooHist* hpull_tt_D4 = ttframe_D4->pullHist();
  RooPlot* framePull_tt_D4 = wspace->var("nj_D4")->frame(Title("Pull"));
  framePull_tt_D4->addPlotable(hpull_tt_D4,"P");


  // RooHist* hpull_other_D1 = otherframe_D1->pullHist();
  // RooPlot* framePull_other_D1 = wspace->var("nj_D1")->frame(Title("Pull"));
  // framePull_other_D1->addPlotable(hpull_other_D1,"P");

  // RooHist* hpull_other_D2 = otherframe_D2->pullHist();
  // RooPlot* framePull_other_D2 = wspace->var("nj_D2")->frame(Title("Pull"));
  // framePull_other_D2->addPlotable(hpull_other_D2,"P");

  // RooHist* hpull_other_D3 = otherframe_D3->pullHist();
  // RooPlot* framePull_other_D3 = wspace->var("nj_D3")->frame(Title("Pull"));
  // framePull_other_D3->addPlotable(hpull_other_D3,"P");

  // RooHist* hpull_other_D4 = otherframe_D4->pullHist();
  // RooPlot* framePull_other_D4 = wspace->var("nj_D4")->frame(Title("Pull"));
  // framePull_other_D4->addPlotable(hpull_other_D4,"P");


  // ===================================================================================

  // Make 2D hists of the correlation matrices
  TH2* hcorr_tt_D1 = r_tt_D1->correlationHist();
  TH2* hcorr_tt_D2 = r_tt_D2->correlationHist();
  TH2* hcorr_tt_D3 = r_tt_D3->correlationHist();
  TH2* hcorr_tt_D4 = r_tt_D4->correlationHist();

  // TH2* hcorr_other_D1 = r_other_D1->correlationHist();
  // TH2* hcorr_other_D2 = r_other_D2->correlationHist();
  // TH2* hcorr_other_D3 = r_other_D3->correlationHist();
  // TH2* hcorr_other_D4 = r_other_D4->correlationHist();


  // ===================================================================================

  // Create the summary canvases

  gStyle->SetPaintTextFormat("2.3f");

  TCanvas* c_tt_D1 = new TCanvas("c_tt_D1","tt D1",900,900);
  c_tt_D1->Divide(2,2);
  c_tt_D1->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  ttframe_D1->GetYaxis()->SetTitleOffset(2.0);
  ttframe_D1->Draw();
  //chiSquare_tt_D1->Draw("same");  
  text_N7_tt_D1->Draw("same");
  c_tt_D1->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_ttframe_D1->GetYaxis()->SetTitleOffset(1.4);
  clone_ttframe_D1->SetMinimum(0.001);
  clone_ttframe_D1->SetMaximum(100000);
  clone_ttframe_D1->remove("background_tt_D1_paramBox");
  clone_ttframe_D1->Draw();
  c_tt_D1->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  framePull_tt_D1->SetYTitle("(Fit-Data)/errData");
  framePull_tt_D1->GetYaxis()->SetTitleOffset(1.4);
  framePull_tt_D1->Draw();
  c_tt_D1->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_tt_D1->SetStats(0);
  hcorr_tt_D1->SetMarkerSize(2.4);
  hcorr_tt_D1->Draw("COLZ,TEXT");
  c_tt_D1->SaveAs("BkgFit_tt_D1.png");

  TCanvas* c_tt_D2 = new TCanvas("c_tt_D2","tt D2",900,900);
  c_tt_D2->Divide(2,2);
  c_tt_D2->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  ttframe_D2->GetYaxis()->SetTitleOffset(2.0);
  ttframe_D2->Draw();
  //chiSquare_tt_D2->Draw("same");
  text_N7_tt_D2->Draw("same");
  c_tt_D2->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_ttframe_D2->GetYaxis()->SetTitleOffset(1.4);
  clone_ttframe_D2->SetMinimum(0.001);
  clone_ttframe_D2->SetMaximum(100000);
  clone_ttframe_D2->remove("background_tt_D2_paramBox");
  clone_ttframe_D2->Draw();
  c_tt_D2->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  framePull_tt_D2->SetYTitle("(Fit-Data)/errData");
  framePull_tt_D2->GetYaxis()->SetTitleOffset(1.4);
  framePull_tt_D2->Draw();
  c_tt_D2->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_tt_D2->SetStats(0);
  hcorr_tt_D2->SetMarkerSize(2.4);
  hcorr_tt_D2->Draw("COLZ,TEXT");
  c_tt_D2->SaveAs("BkgFit_tt_D2.png");

  TCanvas* c_tt_D3 = new TCanvas("c_tt_D3","tt D3",900,900);
  c_tt_D3->Divide(2,2);
  c_tt_D3->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  ttframe_D3->GetYaxis()->SetTitleOffset(2.0);
  ttframe_D3->Draw();
  //chiSquare_tt_D3->Draw("same");
  text_N7_tt_D3->Draw("same");
  c_tt_D3->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_ttframe_D3->GetYaxis()->SetTitleOffset(1.4);
  clone_ttframe_D3->SetMinimum(0.001);
  clone_ttframe_D3->SetMaximum(100000);
  clone_ttframe_D3->remove("background_tt_D3_paramBox");
  clone_ttframe_D3->Draw();
  c_tt_D3->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  framePull_tt_D3->SetYTitle("(Fit-Data)/errData");
  framePull_tt_D3->GetYaxis()->SetTitleOffset(1.4);
  framePull_tt_D3->Draw();
  c_tt_D3->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_tt_D3->SetStats(0);
  hcorr_tt_D3->SetMarkerSize(2.4);
  hcorr_tt_D3->Draw("COLZ,TEXT");
  c_tt_D3->SaveAs("BkgFit_tt_D3.png");

  TCanvas* c_tt_D4 = new TCanvas("c_tt_D4","tt D4",900,900);
  c_tt_D4->Divide(2,2);
  c_tt_D4->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  ttframe_D4->GetYaxis()->SetTitleOffset(2.0);
  ttframe_D4->Draw();
  //chiSquare_tt_D4->Draw("same");
  text_N7_tt_D4->Draw("same");
  c_tt_D4->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_ttframe_D4->GetYaxis()->SetTitleOffset(1.4);
  clone_ttframe_D4->SetMinimum(0.001);
  clone_ttframe_D4->SetMaximum(100000);
  clone_ttframe_D4->remove("background_tt_D4_paramBox");
  clone_ttframe_D4->Draw();
  c_tt_D4->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  framePull_tt_D4->SetYTitle("(Fit-Data)/errData");
  framePull_tt_D4->GetYaxis()->SetTitleOffset(1.4);
  framePull_tt_D4->Draw();
  c_tt_D4->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_tt_D4->SetStats(0);
  hcorr_tt_D4->SetMarkerSize(2.4);
  hcorr_tt_D4->Draw("COLZ,TEXT");
  c_tt_D4->SaveAs("BkgFit_tt_D4.png");





  // TCanvas* c_other_D1 = new TCanvas("c_other_D1","other D1",900,900);
  // c_other_D1->Divide(2,2);
  // c_other_D1->cd(1);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // otherframe_D1->GetYaxis()->SetTitleOffset(2.0);
  // otherframe_D1->Draw();
  // //chiSquare_other_D1->Draw("same");
  // text_N7_other_D1->Draw("same");
  // c_other_D1->cd(2);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(0);
  // gPad->SetLogy();
  // clone_otherframe_D1->GetYaxis()->SetTitleOffset(1.4);
  // clone_otherframe_D1->SetMinimum(0.1);
  // clone_otherframe_D1->SetMaximum(100000);
  // clone_otherframe_D1->remove("bkgMC_other_D1_paramBox");
  // clone_otherframe_D1->Draw();
  // c_other_D1->cd(3);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // framePull_other_D1->SetYTitle("(Fit-Data)/errData");
  // framePull_other_D1->GetYaxis()->SetTitleOffset(1.4);
  // framePull_other_D1->Draw();
  // c_other_D1->cd(4);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(1);
  // hcorr_other_D1->SetStats(0);
  // hcorr_other_D1->SetMarkerSize(2.4);
  // hcorr_other_D1->Draw("COLZ,TEXT");
  // c_other_D1->SaveAs("BkgFit_other_D1.png");

  // TCanvas* c_other_D2 = new TCanvas("c_other_D2","other D2",900,900);
  // c_other_D2->Divide(2,2);
  // c_other_D2->cd(1);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // otherframe_D2->GetYaxis()->SetTitleOffset(2.0);
  // otherframe_D2->Draw();
  // //chiSquare_other_D2->Draw("same");
  // text_N7_other_D2->Draw("same");
  // c_other_D2->cd(2);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(0);
  // gPad->SetLogy();
  // clone_otherframe_D2->GetYaxis()->SetTitleOffset(1.4);
  // clone_otherframe_D2->SetMinimum(0.1);
  // clone_otherframe_D2->SetMaximum(100000);
  // clone_otherframe_D2->remove("bkgMC_other_D2_paramBox");
  // clone_otherframe_D2->Draw();
  // c_other_D2->cd(3);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // framePull_other_D2->SetYTitle("(Fit-Data)/errData");
  // framePull_other_D2->GetYaxis()->SetTitleOffset(1.4);
  // framePull_other_D2->Draw();
  // c_other_D2->cd(4);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(1);
  // hcorr_other_D2->SetStats(0);
  // hcorr_other_D2->SetMarkerSize(2.4);
  // hcorr_other_D2->Draw("COLZ,TEXT");
  // c_other_D2->SaveAs("BkgFit_other_D2.png");

  // TCanvas* c_other_D3 = new TCanvas("c_other_D3","other D3",900,900);
  // c_other_D3->Divide(2,2);
  // c_other_D3->cd(1);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // otherframe_D3->GetYaxis()->SetTitleOffset(2.0);
  // otherframe_D3->Draw();
  // //chiSquare_other_D3->Draw("same");
  // text_N7_other_D3->Draw("same");
  // c_other_D3->cd(2);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(0);
  // gPad->SetLogy();
  // clone_otherframe_D3->GetYaxis()->SetTitleOffset(1.4);
  // clone_otherframe_D3->SetMinimum(0.1);
  // clone_otherframe_D3->SetMaximum(100000);
  // clone_otherframe_D3->remove("bkgMC_other_D3_paramBox");
  // clone_otherframe_D3->Draw();
  // c_other_D3->cd(3);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // framePull_other_D3->SetYTitle("(Fit-Data)/errData");
  // framePull_other_D3->GetYaxis()->SetTitleOffset(1.4);
  // framePull_other_D3->Draw();
  // c_other_D3->cd(4);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(1);
  // hcorr_other_D3->SetStats(0);
  // hcorr_other_D3->SetMarkerSize(2.4);
  // hcorr_other_D3->Draw("COLZ,TEXT");
  // c_other_D3->SaveAs("BkgFit_other_D3.png");

  // TCanvas* c_other_D4 = new TCanvas("c_other_D4","other D4",900,900);
  // c_other_D4->Divide(2,2);
  // c_other_D4->cd(1);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // otherframe_D4->GetYaxis()->SetTitleOffset(2.0);
  // otherframe_D4->Draw();
  // //chiSquare_other_D4->Draw("same");
  // text_N7_other_D4->Draw("same");
  // c_other_D4->cd(2);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(0);
  // gPad->SetLogy();
  // clone_otherframe_D4->GetYaxis()->SetTitleOffset(1.4);
  // clone_otherframe_D4->SetMinimum(0.1);
  // clone_otherframe_D4->SetMaximum(100000);
  // clone_otherframe_D4->remove("bkgMC_other_D4_paramBox");
  // clone_otherframe_D4->Draw();
  // c_other_D4->cd(3);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.05);
  // framePull_other_D4->SetYTitle("(Fit-Data)/errData");
  // framePull_other_D4->GetYaxis()->SetTitleOffset(1.4);
  // framePull_other_D4->Draw();
  // c_other_D4->cd(4);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gStyle->SetOptTitle(1);
  // hcorr_other_D4->SetStats(0);
  // hcorr_other_D4->SetMarkerSize(2.4);
  // hcorr_other_D4->Draw("COLZ,TEXT");
  // c_other_D4->SaveAs("BkgFit_other_D4.png");


  // ===================================================================================

  // Extract and print the correlation matrices

  const TMatrixDSym& cor_tt_D1 = r_tt_D1->correlationMatrix();
  const TMatrixDSym& cor_tt_D2 = r_tt_D2->correlationMatrix();
  const TMatrixDSym& cor_tt_D3 = r_tt_D3->correlationMatrix();
  const TMatrixDSym& cor_tt_D4 = r_tt_D4->correlationMatrix();

  // const TMatrixDSym& cor_other_D1 = r_other_D1->correlationMatrix();
  // const TMatrixDSym& cor_other_D2 = r_other_D2->correlationMatrix();
  // const TMatrixDSym& cor_other_D3 = r_other_D3->correlationMatrix();
  // const TMatrixDSym& cor_other_D4 = r_other_D4->correlationMatrix();
  
  cout << "======================= Correlation Matrices of Fit Results =====================" << endl;
  cout << "tt D1:" << endl;
  cor_tt_D1.Print();
  cout << "tt D2:" << endl;
  cor_tt_D2.Print();
  cout << "tt D3:" << endl;
  cor_tt_D3.Print();
  cout << "tt D4:" << endl;
  cor_tt_D4.Print();
  // cout << "other D1:" << endl;
  // cor_other_D1.Print();
  // cout << "other D2:" << endl;
  // cor_other_D2.Print();
  // cout << "other D3:" << endl;
  // cor_other_D3.Print();
  // cout << "other D4:" << endl;
  // cor_other_D4.Print();

  cout << "======================= Status =====================" << endl;
  cout << "tt D1:  ";;
  cout << r_tt_D1->status() << endl;
  cout << "tt D2:  ";;
  cout << r_tt_D2->status() << endl;
  cout << "tt D3:  ";;
  cout << r_tt_D3->status() << endl;
  cout << "tt D4:  ";;
  cout << r_tt_D4->status() << endl;

  cout << "======================= Minimum Negative Log Likelihood =====================" << endl;
  cout << "tt D1:  ";;
  cout << r_tt_D1->minNll() << endl;
  cout << "tt D2:  ";;
  cout << r_tt_D2->minNll() << endl;
  cout << "tt D3:  ";;
  cout << r_tt_D3->minNll() << endl;
  cout << "tt D4:  ";;
  cout << r_tt_D4->minNll() << endl;

  cout << "======================= Details =====================" << endl;
  cout << "------------------------ D1 ---------------------------" << endl;
  r_tt_D1->Print("V");
  cout << "------------------------ D2 ---------------------------" << endl;
  r_tt_D2->Print("V");
  cout << "------------------------ D3 ---------------------------" << endl;
  r_tt_D3->Print("V");
  cout << "------------------------ D4 ---------------------------" << endl;
  r_tt_D4->Print("V");

  // ===================================================================================

  // Extract and print the post-fit parameters

  RooRealVar* par_p0_tt_D1 = (RooRealVar*) r_tt_D1->floatParsFinal().find("p0_tt_D1");
  RooRealVar* par_p0_tt_D2 = (RooRealVar*) r_tt_D2->floatParsFinal().find("p0_tt_D2");
  RooRealVar* par_p0_tt_D3 = (RooRealVar*) r_tt_D3->floatParsFinal().find("p0_tt_D3");
  RooRealVar* par_p0_tt_D4 = (RooRealVar*) r_tt_D4->floatParsFinal().find("p0_tt_D4");

  RooRealVar* par_p1_tt_D1 = (RooRealVar*) r_tt_D1->floatParsFinal().find("p1_tt_D1");
  RooRealVar* par_p1_tt_D2 = (RooRealVar*) r_tt_D2->floatParsFinal().find("p1_tt_D2");
  RooRealVar* par_p1_tt_D3 = (RooRealVar*) r_tt_D3->floatParsFinal().find("p1_tt_D3");
  RooRealVar* par_p1_tt_D4 = (RooRealVar*) r_tt_D4->floatParsFinal().find("p1_tt_D4");

  RooRealVar* par_p2_tt_D1 = (RooRealVar*) r_tt_D1->floatParsFinal().find("p2_tt_D1");
  RooRealVar* par_p2_tt_D2 = (RooRealVar*) r_tt_D2->floatParsFinal().find("p2_tt_D2");
  RooRealVar* par_p2_tt_D3 = (RooRealVar*) r_tt_D3->floatParsFinal().find("p2_tt_D3");
  RooRealVar* par_p2_tt_D4 = (RooRealVar*) r_tt_D4->floatParsFinal().find("p2_tt_D4");


  // RooRealVar* par_p0_other_D1 = (RooRealVar*) r_other_D1->floatParsFinal().find("p0_other_D1");
  // RooRealVar* par_p0_other_D2 = (RooRealVar*) r_other_D2->floatParsFinal().find("p0_other_D2");
  // RooRealVar* par_p0_other_D3 = (RooRealVar*) r_other_D3->floatParsFinal().find("p0_other_D3");
  // RooRealVar* par_p0_other_D4 = (RooRealVar*) r_other_D4->floatParsFinal().find("p0_other_D4");

  // RooRealVar* par_p1_other_D1 = (RooRealVar*) r_other_D1->floatParsFinal().find("p1_other_D1");
  // RooRealVar* par_p1_other_D2 = (RooRealVar*) r_other_D2->floatParsFinal().find("p1_other_D2");
  // RooRealVar* par_p1_other_D3 = (RooRealVar*) r_other_D3->floatParsFinal().find("p1_other_D3");
  // RooRealVar* par_p1_other_D4 = (RooRealVar*) r_other_D4->floatParsFinal().find("p1_other_D4");

  // RooRealVar* par_p2_other_D1 = (RooRealVar*) r_other_D1->floatParsFinal().find("p2_other_D1");
  // RooRealVar* par_p2_other_D2 = (RooRealVar*) r_other_D2->floatParsFinal().find("p2_other_D2");
  // RooRealVar* par_p2_other_D3 = (RooRealVar*) r_other_D3->floatParsFinal().find("p2_other_D3");
  // RooRealVar* par_p2_other_D4 = (RooRealVar*) r_other_D4->floatParsFinal().find("p2_other_D4");


  cout << "======================= Parameters of Fit Results =====================" << endl;
  cout << "Sample   " << "            p0                 p1                  p2" << endl;

  cout << "tt D1   " << "      " << TString::Format("%0.3f",par_p0_tt_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_tt_D1->getError()) << "    " << TString::Format("%0.3f",par_p1_tt_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_tt_D1->getError()) << "    " << TString::Format("%0.3f",par_p2_tt_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_tt_D1->getError()) << "    " << endl;

  cout << "tt D2   " << "      " << TString::Format("%0.3f",par_p0_tt_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_tt_D2->getError()) << "    " << TString::Format("%0.3f",par_p1_tt_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_tt_D2->getError()) << "    " << TString::Format("%0.3f",par_p2_tt_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_tt_D2->getError()) << "    " << endl;

  cout << "tt D3   " << "      " << TString::Format("%0.3f",par_p0_tt_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_tt_D3->getError()) << "    " << TString::Format("%0.3f",par_p1_tt_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_tt_D3->getError()) << "    " << TString::Format("%0.3f",par_p2_tt_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_tt_D3->getError()) << "    " << endl;

  cout << "tt D4   " << "      " << TString::Format("%0.3f",par_p0_tt_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_tt_D4->getError()) << "    " << TString::Format("%0.3f",par_p1_tt_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_tt_D4->getError()) << "    " << TString::Format("%0.3f",par_p2_tt_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_tt_D4->getError()) << "    " << endl;


  // cout << "other D1" << "      " << TString::Format("%0.3f",par_p0_other_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_other_D1->getError()) << "    " << TString::Format("%0.3f",par_p1_other_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_other_D1->getError()) << "    " << TString::Format("%0.3f",par_p2_other_D1->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_other_D1->getError()) << "    " << endl;

  // cout << "other D2" << "      " << TString::Format("%0.3f",par_p0_other_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_other_D2->getError()) << "    " << TString::Format("%0.3f",par_p1_other_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_other_D2->getError()) << "    " << TString::Format("%0.3f",par_p2_other_D2->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_other_D2->getError()) << "    " << endl;

  // cout << "other D3" << "      " << TString::Format("%0.3f",par_p0_other_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_other_D3->getError()) << "    " << TString::Format("%0.3f",par_p1_other_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_other_D3->getError()) << "    " << TString::Format("%0.3f",par_p2_other_D3->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_other_D3->getError()) << "    " << endl;

  // cout << "other D4" << "      " << TString::Format("%0.3f",par_p0_other_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p0_other_D4->getError()) << "    " << TString::Format("%0.3f",par_p1_other_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p1_other_D4->getError()) << "    " << TString::Format("%0.3f",par_p2_other_D4->getVal()) << " +/- " << TString::Format("%0.3f",par_p2_other_D4->getError()) << "    " << endl;





  // The residual is defined as curve-histogram +/- error_on_histogram
  // The pull is defined as (curve-histogram +/- error_on_histogram) / error_on_histogram
  //   so it is expected that for the case of symmetric errors the error bars on the pull are +/- 1. 

  // I think the pull has the most useful information out of these two

  //new TCanvas("Test","Test",600,600);
  //RooHist* hpull = ttframe_D1->pullHist();
  //RooPlot* frameTest = wspace->var("nj")->frame(Title("Pull Dist"));
  //frameTest->addPlotable(hpull,"P");
  //frameTest->Draw();

  //new TCanvas("Test2","Test2",600,600);
  //RooHist* hresid = ttframe_D1->residHist();
  //RooPlot* frameTest2 = wspace->var("nj")->frame(Title("Resid Dist"));
  //frameTest2->addPlotable(hresid,"P");
  //frameTest2->Draw();




}

