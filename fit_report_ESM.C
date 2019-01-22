// Macro for summarizing and plotting results of full fit
//  but without simultaneous fit to tt or other background MC
// Aron Soha
// April 23, 2018
// September 19, 2018 : Updated to use new DeepESM
// December 21, 2018 : Updated to MVA region naming (e.g. D1 instead of sigD1)

// call with:
// root -l
// .x fit_report_ESM.C("filename")
// or
// .x fit_report_ESM.C  to use default

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH2.h"

using namespace RooFit;

void fit_report_ESM(const string filename = "fitDiagnostics.root", const string outPath = "./", const bool plotErrors = true) {

  TFile* fullfile = TFile::Open(filename.c_str());

  RooWorkspace *w = new RooWorkspace("w","w");
  w->factory("CMS_th1x[0,8]");
  w->var("CMS_th1x")->setBins(8);

  // ===================================================================================
  RooPlot* D1_CMS_th1x_prefit = (RooPlot*)fullfile->Get("D1_CMS_th1x_prefit");
  RooPlot* D1_CMS_th1x_fit_s  = (RooPlot*)fullfile->Get("D1_CMS_th1x_fit_s");
  RooPlot* D1_CMS_th1x_fit_b  = (RooPlot*)fullfile->Get("D1_CMS_th1x_fit_b");

  RooPlot* D2_CMS_th1x_prefit = (RooPlot*)fullfile->Get("D2_CMS_th1x_prefit");
  RooPlot* D2_CMS_th1x_fit_s  = (RooPlot*)fullfile->Get("D2_CMS_th1x_fit_s");
  RooPlot* D2_CMS_th1x_fit_b  = (RooPlot*)fullfile->Get("D2_CMS_th1x_fit_b");

  RooPlot* D3_CMS_th1x_prefit = (RooPlot*)fullfile->Get("D3_CMS_th1x_prefit");
  RooPlot* D3_CMS_th1x_fit_s  = (RooPlot*)fullfile->Get("D3_CMS_th1x_fit_s");
  RooPlot* D3_CMS_th1x_fit_b  = (RooPlot*)fullfile->Get("D3_CMS_th1x_fit_b");

  RooPlot* D4_CMS_th1x_prefit = (RooPlot*)fullfile->Get("D4_CMS_th1x_prefit");
  RooPlot* D4_CMS_th1x_fit_s  = (RooPlot*)fullfile->Get("D4_CMS_th1x_fit_s");
  RooPlot* D4_CMS_th1x_fit_b  = (RooPlot*)fullfile->Get("D4_CMS_th1x_fit_b");

  if (!plotErrors) {
    // Remove all the error bands
    D1_CMS_th1x_prefit->remove(D1_CMS_th1x_prefit->nameOf(1));
    D1_CMS_th1x_fit_s->remove(D1_CMS_th1x_fit_s->nameOf(1));
    D1_CMS_th1x_fit_b->remove(D1_CMS_th1x_fit_b->nameOf(1));
  
    D2_CMS_th1x_prefit->remove(D2_CMS_th1x_prefit->nameOf(1));
    D2_CMS_th1x_fit_s->remove(D2_CMS_th1x_fit_s->nameOf(1));
    D2_CMS_th1x_fit_b->remove(D2_CMS_th1x_fit_b->nameOf(1));
  
    D3_CMS_th1x_prefit->remove(D3_CMS_th1x_prefit->nameOf(1));
    D3_CMS_th1x_fit_s->remove(D3_CMS_th1x_fit_s->nameOf(1));
    D3_CMS_th1x_fit_b->remove(D3_CMS_th1x_fit_b->nameOf(1));
  
    D4_CMS_th1x_prefit->remove(D4_CMS_th1x_prefit->nameOf(1));
    D4_CMS_th1x_fit_s->remove(D4_CMS_th1x_fit_s->nameOf(1));
    D4_CMS_th1x_fit_b->remove(D4_CMS_th1x_fit_b->nameOf(1));
  }

  // Clone the fit result frames, so that I can plot a log version without the fit params
  RooPlot* clone_D1_prefit_frame = (RooPlot*)D1_CMS_th1x_prefit->Clone("D1_prefit_frame");
  RooPlot* clone_D1_fit_frame = (RooPlot*)D1_CMS_th1x_fit_s->Clone("D1_fit_frame");
  RooPlot* clone_D1_fitb_frame = (RooPlot*)D1_CMS_th1x_fit_b->Clone("D1_fitb_frame");

  RooPlot* clone_D2_prefit_frame = (RooPlot*)D2_CMS_th1x_prefit->Clone("D2_prefit_frame");
  RooPlot* clone_D2_fit_frame = (RooPlot*)D2_CMS_th1x_fit_s->Clone("D2_fit_frame");
  RooPlot* clone_D2_fitb_frame = (RooPlot*)D2_CMS_th1x_fit_b->Clone("D2_fitb_frame");

  RooPlot* clone_D3_prefit_frame = (RooPlot*)D3_CMS_th1x_prefit->Clone("D3_prefit_frame");
  RooPlot* clone_D3_fit_frame = (RooPlot*)D3_CMS_th1x_fit_s->Clone("D3_fit_frame");
  RooPlot* clone_D3_fitb_frame = (RooPlot*)D3_CMS_th1x_fit_b->Clone("D3_fitb_frame");

  RooPlot* clone_D4_prefit_frame = (RooPlot*)D4_CMS_th1x_prefit->Clone("D4_prefit_frame");
  RooPlot* clone_D4_fit_frame = (RooPlot*)D4_CMS_th1x_fit_s->Clone("D4_fit_frame");
  RooPlot* clone_D4_fitb_frame = (RooPlot*)D4_CMS_th1x_fit_b->Clone("D4_fitb_frame");

  // ===================================================================================

  // Get the N6 for each fit

  //TPaveText *text_N6_tt_D1 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  //text_N6_tt_D1->AddText(TString::Format("N6 = %.1f",full_ttMC_th1_D1->GetBinContent(1)));
  //text_N6_tt_D1->SetFillStyle(0);
  //text_N6_tt_D1->SetBorderSize(1);
  //text_N6_tt_D1->SetLineColor(0);

  // ===================================================================================

  // Get the chi^2 for each fit

  TPaveText *chiSquare_D1_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D1_prefit->AddText(TString::Format("Chi Square = %0.2f",D1_CMS_th1x_prefit->chiSquare()));
  chiSquare_D1_prefit->SetFillStyle(0);
  chiSquare_D1_prefit->SetBorderSize(1);
  chiSquare_D1_prefit->SetLineColor(0);

  TPaveText *chiSquare_D1_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D1_fit->AddText(TString::Format("Chi Square = %0.2f",D1_CMS_th1x_fit_s->chiSquare()));
  chiSquare_D1_fit->SetFillStyle(0);
  chiSquare_D1_fit->SetBorderSize(1);
  chiSquare_D1_fit->SetLineColor(0);

  TPaveText *chiSquare_D1_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D1_fitb->AddText(TString::Format("Chi Square = %0.2f",D1_CMS_th1x_fit_b->chiSquare()));
  chiSquare_D1_fitb->SetFillStyle(0);
  chiSquare_D1_fitb->SetBorderSize(1);
  chiSquare_D1_fitb->SetLineColor(0);


  TPaveText *chiSquare_D2_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D2_prefit->AddText(TString::Format("Chi Square = %0.2f",D2_CMS_th1x_prefit->chiSquare()));
  chiSquare_D2_prefit->SetFillStyle(0);
  chiSquare_D2_prefit->SetBorderSize(1);
  chiSquare_D2_prefit->SetLineColor(0);

  TPaveText *chiSquare_D2_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D2_fit->AddText(TString::Format("Chi Square = %0.2f",D2_CMS_th1x_fit_s->chiSquare()));
  chiSquare_D2_fit->SetFillStyle(0);
  chiSquare_D2_fit->SetBorderSize(1);
  chiSquare_D2_fit->SetLineColor(0);

  TPaveText *chiSquare_D2_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D2_fitb->AddText(TString::Format("Chi Square = %0.2f",D2_CMS_th1x_fit_b->chiSquare()));
  chiSquare_D2_fitb->SetFillStyle(0);
  chiSquare_D2_fitb->SetBorderSize(1);
  chiSquare_D2_fitb->SetLineColor(0);


  TPaveText *chiSquare_D3_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D3_prefit->AddText(TString::Format("Chi Square = %0.2f",D3_CMS_th1x_prefit->chiSquare()));
  chiSquare_D3_prefit->SetFillStyle(0);
  chiSquare_D3_prefit->SetBorderSize(1);
  chiSquare_D3_prefit->SetLineColor(0);

  TPaveText *chiSquare_D3_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D3_fit->AddText(TString::Format("Chi Square = %0.2f",D3_CMS_th1x_fit_s->chiSquare()));
  chiSquare_D3_fit->SetFillStyle(0);
  chiSquare_D3_fit->SetBorderSize(1);
  chiSquare_D3_fit->SetLineColor(0);

  TPaveText *chiSquare_D3_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D3_fitb->AddText(TString::Format("Chi Square = %0.2f",D3_CMS_th1x_fit_b->chiSquare()));
  chiSquare_D3_fitb->SetFillStyle(0);
  chiSquare_D3_fitb->SetBorderSize(1);
  chiSquare_D3_fitb->SetLineColor(0);


  TPaveText *chiSquare_D4_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D4_prefit->AddText(TString::Format("Chi Square = %0.2f",D4_CMS_th1x_prefit->chiSquare()));
  chiSquare_D4_prefit->SetFillStyle(0);
  chiSquare_D4_prefit->SetBorderSize(1);
  chiSquare_D4_prefit->SetLineColor(0);

  TPaveText *chiSquare_D4_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D4_fit->AddText(TString::Format("Chi Square = %0.2f",D4_CMS_th1x_fit_s->chiSquare()));
  chiSquare_D4_fit->SetFillStyle(0);
  chiSquare_D4_fit->SetBorderSize(1);
  chiSquare_D4_fit->SetLineColor(0);

  TPaveText *chiSquare_D4_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_D4_fitb->AddText(TString::Format("Chi Square = %0.2f",D4_CMS_th1x_fit_b->chiSquare()));
  chiSquare_D4_fitb->SetFillStyle(0);
  chiSquare_D4_fitb->SetBorderSize(1);
  chiSquare_D4_fitb->SetLineColor(0);

  // ===================================================================================

  // Make pull plots of each fit
  // The pull is defined as (curve-histogram +/- error_on_histogram) / error_on_histogram
  //   so it is expected that for the case of symmetric errors the error bars on the pull are +/- 1. 

  RooHist* hpull_D1_prefit = D1_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_D1_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D1_prefit->addPlotable(hpull_D1_prefit,"P");

  RooHist* hpull_D1_fit = D1_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_D1_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D1_fit->addPlotable(hpull_D1_fit,"P");

  RooHist* hpull_D1_fitb = D1_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_D1_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D1_fitb->addPlotable(hpull_D1_fitb,"P");


  RooHist* hpull_D2_prefit = D2_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_D2_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D2_prefit->addPlotable(hpull_D2_prefit,"P");

  RooHist* hpull_D2_fit = D2_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_D2_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D2_fit->addPlotable(hpull_D2_fit,"P");

  RooHist* hpull_D2_fitb = D2_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_D2_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D2_fitb->addPlotable(hpull_D2_fitb,"P");


  RooHist* hpull_D3_prefit = D3_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_D3_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D3_prefit->addPlotable(hpull_D3_prefit,"P");

  RooHist* hpull_D3_fit = D3_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_D3_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D3_fit->addPlotable(hpull_D3_fit,"P");

  RooHist* hpull_D3_fitb = D3_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_D3_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D3_fitb->addPlotable(hpull_D3_fitb,"P");


  RooHist* hpull_D4_prefit = D4_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_D4_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D4_prefit->addPlotable(hpull_D4_prefit,"P");

  RooHist* hpull_D4_fit = D4_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_D4_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D4_fit->addPlotable(hpull_D4_fit,"P");

  RooHist* hpull_D4_fitb = D4_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_D4_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_D4_fitb->addPlotable(hpull_D4_fitb,"P");


  // ===================================================================================

  RooFitResult* fit_s = 0;
  fullfile->GetObject("fit_s",fit_s);
  RooFitResult* fit_b = 0;
  fullfile->GetObject("fit_b",fit_b);

  // Make 2D hists of the correlation matrices
  TH2* hcorr_fit = fit_s->correlationHist();
  TH2* hcorr_fitb = fit_b->correlationHist();


  // ===================================================================================

  // Create the summary canvases

  gStyle->SetPaintTextFormat("2.3f");

  // ------------ D1 -----------

  TCanvas* c_D1_prefit = new TCanvas("c_D1_prefit","Signal D1 prefit",900,900);
  c_D1_prefit->Divide(2,2);
  c_D1_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D1_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  D1_CMS_th1x_prefit->Draw();
  chiSquare_D1_prefit->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_D1_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D1_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D1_prefit_frame->SetMinimum(0.001);
  clone_D1_prefit_frame->SetMaximum(100000);
  //clone_D1_prefit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_D1_prefit_frame->Draw();
  c_D1_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D1_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_D1_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D1_prefit->Draw();
  //c_D1_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D1_prefit->SetStats(0);
  //hcorr_D1_prefit->SetMarkerSize(2.4);
  //hcorr_D1_prefit->Draw("COLZ,TEXT");
  c_D1_prefit->SaveAs((outPath+"fit_report_D1_prefit.png").c_str());

  TCanvas* c_D1_fit = new TCanvas("c_D1_fit","Signal D1 fit",900,900);
  c_D1_fit->Divide(2,2);
  c_D1_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D1_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  D1_CMS_th1x_fit_s->Draw();
  chiSquare_D1_fit->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_D1_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D1_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D1_fit_frame->SetMinimum(0.001);
  clone_D1_fit_frame->SetMaximum(100000);
  //clone_D1_fit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_D1_fit_frame->Draw();
  c_D1_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D1_fit->SetYTitle("(Fit-Data)/errData");
  framePull_D1_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D1_fit->Draw();
  //c_D1_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D1_fit->SetStats(0);
  //hcorr_D1_fit->SetMarkerSize(2.4);
  //hcorr_D1_fit->Draw("COLZ,TEXT");
  c_D1_fit->SaveAs((outPath+"fit_report_D1_fit.png").c_str());

  TCanvas* c_D1_fitb = new TCanvas("c_D1_fitb","D1 fit bkg only",900,900);
  c_D1_fitb->Divide(2,2);
  c_D1_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D1_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  D1_CMS_th1x_fit_b->Draw();
  chiSquare_D1_fitb->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_D1_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D1_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D1_fitb_frame->SetMinimum(0.001);
  clone_D1_fitb_frame->SetMaximum(100000);
  //clone_D1_fit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_D1_fitb_frame->Draw();
  c_D1_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D1_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_D1_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_D1_fitb->Draw();
  //c_D1_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D1_fitb->SetStats(0);
  //hcorr_D1_fitb->SetMarkerSize(2.4);
  //hcorr_D1_fitb->Draw("COLZ,TEXT");
  c_D1_fitb->SaveAs((outPath+"fit_report_D1_fitb.png").c_str());


  // ------------ D2 -----------

  TCanvas* c_D2_prefit = new TCanvas("c_D2_prefit","Signal D2 prefit",900,900);
  c_D2_prefit->Divide(2,2);
  c_D2_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D2_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  D2_CMS_th1x_prefit->Draw();
  chiSquare_D2_prefit->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_D2_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D2_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D2_prefit_frame->SetMinimum(0.001);
  clone_D2_prefit_frame->SetMaximum(100000);
  //clone_D2_prefit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_D2_prefit_frame->Draw();
  c_D2_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D2_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_D2_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D2_prefit->Draw();
  //c_D2_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D2_prefit->SetStats(0);
  //hcorr_D2_prefit->SetMarkerSize(2.4);
  //hcorr_D2_prefit->Draw("COLZ,TEXT");
  c_D2_prefit->SaveAs((outPath+"fit_report_D2_prefit.png").c_str());

  TCanvas* c_D2_fit = new TCanvas("c_D2_fit","Signal D2 fit",900,900);
  c_D2_fit->Divide(2,2);
  c_D2_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D2_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  D2_CMS_th1x_fit_s->Draw();
  chiSquare_D2_fit->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_D2_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D2_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D2_fit_frame->SetMinimum(0.001);
  clone_D2_fit_frame->SetMaximum(100000);
  //clone_D2_fit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_D2_fit_frame->Draw();
  c_D2_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D2_fit->SetYTitle("(Fit-Data)/errData");
  framePull_D2_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D2_fit->Draw();
  //c_D2_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D2_fit->SetStats(0);
  //hcorr_D2_fit->SetMarkerSize(2.4);
  //hcorr_D2_fit->Draw("COLZ,TEXT");
  c_D2_fit->SaveAs((outPath+"fit_report_D2_fit.png").c_str());

  TCanvas* c_D2_fitb = new TCanvas("c_D2_fitb","D2 fit bkg only",900,900);
  c_D2_fitb->Divide(2,2);
  c_D2_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D2_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  D2_CMS_th1x_fit_b->Draw();
  chiSquare_D2_fitb->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_D2_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D2_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D2_fitb_frame->SetMinimum(0.001);
  clone_D2_fitb_frame->SetMaximum(100000);
  //clone_D2_fit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_D2_fitb_frame->Draw();
  c_D2_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D2_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_D2_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_D2_fitb->Draw();
  //c_D2_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D2_fitb->SetStats(0);
  //hcorr_D2_fitb->SetMarkerSize(2.4);
  //hcorr_D2_fitb->Draw("COLZ,TEXT");
  c_D2_fitb->SaveAs((outPath+"fit_report_D2_fitb.png").c_str());


  // ------------ D3 -----------

  TCanvas* c_D3_prefit = new TCanvas("c_D3_prefit","Signal D3 prefit",900,900);
  c_D3_prefit->Divide(2,2);
  c_D3_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D3_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  D3_CMS_th1x_prefit->Draw();
  chiSquare_D3_prefit->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_D3_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D3_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D3_prefit_frame->SetMinimum(0.001);
  clone_D3_prefit_frame->SetMaximum(100000);
  //clone_D3_prefit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_D3_prefit_frame->Draw();
  c_D3_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D3_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_D3_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D3_prefit->Draw();
  //c_D3_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D3_prefit->SetStats(0);
  //hcorr_D3_prefit->SetMarkerSize(2.4);
  //hcorr_D3_prefit->Draw("COLZ,TEXT");
  c_D3_prefit->SaveAs((outPath+"fit_report_D3_prefit.png").c_str());

  TCanvas* c_D3_fit = new TCanvas("c_D3_fit","Signal D3 fit",900,900);
  c_D3_fit->Divide(2,2);
  c_D3_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D3_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  D3_CMS_th1x_fit_s->Draw();
  chiSquare_D3_fit->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_D3_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D3_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D3_fit_frame->SetMinimum(0.001);
  clone_D3_fit_frame->SetMaximum(100000);
  //clone_D3_fit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_D3_fit_frame->Draw();
  c_D3_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D3_fit->SetYTitle("(Fit-Data)/errData");
  framePull_D3_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D3_fit->Draw();
  //c_D3_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D3_fit->SetStats(0);
  //hcorr_D3_fit->SetMarkerSize(2.4);
  //hcorr_D3_fit->Draw("COLZ,TEXT");
  c_D3_fit->SaveAs((outPath+"fit_report_D3_fit.png").c_str());

  TCanvas* c_D3_fitb = new TCanvas("c_D3_fitb","D3 fit bkg only",900,900);
  c_D3_fitb->Divide(2,2);
  c_D3_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D3_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  D3_CMS_th1x_fit_b->Draw();
  chiSquare_D3_fitb->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_D3_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D3_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D3_fitb_frame->SetMinimum(0.001);
  clone_D3_fitb_frame->SetMaximum(100000);
  //clone_D3_fit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_D3_fitb_frame->Draw();
  c_D3_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D3_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_D3_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_D3_fitb->Draw();
  //c_D3_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D3_fitb->SetStats(0);
  //hcorr_D3_fitb->SetMarkerSize(2.4);
  //hcorr_D3_fitb->Draw("COLZ,TEXT");
  c_D3_fitb->SaveAs((outPath+"fit_report_D3_fitb.png").c_str());


  // ------------ D4 -----------

  TCanvas* c_D4_prefit = new TCanvas("c_D4_prefit","Signal D4 prefit",900,900);
  c_D4_prefit->Divide(2,2);
  c_D4_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D4_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  D4_CMS_th1x_prefit->Draw();
  chiSquare_D4_prefit->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_D4_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D4_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D4_prefit_frame->SetMinimum(0.001);
  clone_D4_prefit_frame->SetMaximum(100000);
  //clone_D4_prefit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_D4_prefit_frame->Draw();
  c_D4_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D4_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_D4_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D4_prefit->Draw();
  //c_D4_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D4_prefit->SetStats(0);
  //hcorr_D4_prefit->SetMarkerSize(2.4);
  //hcorr_D4_prefit->Draw("COLZ,TEXT");
  c_D4_prefit->SaveAs((outPath+"fit_report_D4_prefit.png").c_str());

  TCanvas* c_D4_fit = new TCanvas("c_D4_fit","Signal D4 fit",900,900);
  c_D4_fit->Divide(2,2);
  c_D4_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D4_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  D4_CMS_th1x_fit_s->Draw();
  chiSquare_D4_fit->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_D4_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D4_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D4_fit_frame->SetMinimum(0.001);
  clone_D4_fit_frame->SetMaximum(100000);
  //clone_D4_fit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_D4_fit_frame->Draw();
  c_D4_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D4_fit->SetYTitle("(Fit-Data)/errData");
  framePull_D4_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_D4_fit->Draw();
  //c_D4_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D4_fit->SetStats(0);
  //hcorr_D4_fit->SetMarkerSize(2.4);
  //hcorr_D4_fit->Draw("COLZ,TEXT");
  c_D4_fit->SaveAs((outPath+"fit_report_D4_fit.png").c_str());

  TCanvas* c_D4_fitb = new TCanvas("c_D4_fitb","D4 fit bkg only",900,900);
  c_D4_fitb->Divide(2,2);
  c_D4_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  D4_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  D4_CMS_th1x_fit_b->Draw();
  chiSquare_D4_fitb->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_D4_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_D4_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_D4_fitb_frame->SetMinimum(0.001);
  clone_D4_fitb_frame->SetMaximum(100000);
  //clone_D4_fit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_D4_fitb_frame->Draw();
  c_D4_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_D4_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_D4_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_D4_fitb->Draw();
  //c_D4_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_D4_fitb->SetStats(0);
  //hcorr_D4_fitb->SetMarkerSize(2.4);
  //hcorr_D4_fitb->Draw("COLZ,TEXT");
  c_D4_fitb->SaveAs((outPath+"fit_report_D4_fitb.png").c_str());

  TCanvas* c_fit_s_corr = new TCanvas("c_fit_s_corr","signal region correlation matrix",900,900);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_fit->SetStats(0);
  hcorr_fit->SetMarkerSize(1.0);
  //hcorr_fit->Draw("COLZ,TEXT");
  hcorr_fit->Draw("COLZ");
  c_fit_s_corr->SaveAs((outPath+"fit_report_correlation.png").c_str());

  TCanvas* c_fit_b_corr = new TCanvas("c_fit_b_corr","signal region correlation matrix bkg only",900,900);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_fitb->SetStats(0);
  hcorr_fitb->SetMarkerSize(1.0);
  //hcorr_fitb->Draw("COLZ,TEXT");
  hcorr_fitb->Draw("COLZ");
  c_fit_b_corr->SaveAs((outPath+"fit_report_correlation_b.png").c_str());
}
