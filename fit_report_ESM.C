// Macro for summarizing and plotting results of full fit
//  but without simultaneous fit to tt or other background MC
// Aron Soha
// April 23, 2018
// September 19, 2018 : Updated to use new DeepESM

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH2.h"

using namespace RooFit;

void fit_report_ESM() {

 TFile* fullfile = TFile::Open("fitDiagnostics.root");

  RooWorkspace *w = new RooWorkspace("w","w");
  w->factory("CMS_th1x[0,8]");
  w->var("CMS_th1x")->setBins(8);

  // ===================================================================================

  RooPlot* sigD1_CMS_th1x_prefit = 0;
  fullfile->GetObject("sigD1_CMS_th1x_prefit",sigD1_CMS_th1x_prefit);
  RooPlot* sigD1_CMS_th1x_fit_s = 0;
  fullfile->GetObject("sigD1_CMS_th1x_fit_s",sigD1_CMS_th1x_fit_s);
  RooPlot* sigD1_CMS_th1x_fit_b = 0;
  fullfile->GetObject("sigD1_CMS_th1x_fit_b",sigD1_CMS_th1x_fit_b);

  RooPlot* sigD2_CMS_th1x_prefit = 0;
  fullfile->GetObject("sigD2_CMS_th1x_prefit",sigD2_CMS_th1x_prefit);
  RooPlot* sigD2_CMS_th1x_fit_s = 0;
  fullfile->GetObject("sigD2_CMS_th1x_fit_s",sigD2_CMS_th1x_fit_s);
  RooPlot* sigD2_CMS_th1x_fit_b = 0;
  fullfile->GetObject("sigD2_CMS_th1x_fit_b",sigD2_CMS_th1x_fit_b);

  RooPlot* sigD3_CMS_th1x_prefit = 0;
  fullfile->GetObject("sigD3_CMS_th1x_prefit",sigD3_CMS_th1x_prefit);
  RooPlot* sigD3_CMS_th1x_fit_s = 0;
  fullfile->GetObject("sigD3_CMS_th1x_fit_s",sigD3_CMS_th1x_fit_s);
  RooPlot* sigD3_CMS_th1x_fit_b = 0;
  fullfile->GetObject("sigD3_CMS_th1x_fit_b",sigD3_CMS_th1x_fit_b);

  RooPlot* sigD4_CMS_th1x_prefit = 0;
  fullfile->GetObject("sigD4_CMS_th1x_prefit",sigD4_CMS_th1x_prefit);
  RooPlot* sigD4_CMS_th1x_fit_s = 0;
  fullfile->GetObject("sigD4_CMS_th1x_fit_s",sigD4_CMS_th1x_fit_s);
  RooPlot* sigD4_CMS_th1x_fit_b = 0;
  fullfile->GetObject("sigD4_CMS_th1x_fit_b",sigD4_CMS_th1x_fit_b);

  // Clone the fit result frames, so that I can plot a log version without the fit params
  RooPlot* clone_sigD1_prefit_frame = (RooPlot*)sigD1_CMS_th1x_prefit->Clone("sigD1_prefit_frame");
  RooPlot* clone_sigD1_fit_frame = (RooPlot*)sigD1_CMS_th1x_fit_s->Clone("sigD1_fit_frame");
  RooPlot* clone_sigD1_fitb_frame = (RooPlot*)sigD1_CMS_th1x_fit_b->Clone("sigD1_fitb_frame");

  RooPlot* clone_sigD2_prefit_frame = (RooPlot*)sigD2_CMS_th1x_prefit->Clone("sigD2_prefit_frame");
  RooPlot* clone_sigD2_fit_frame = (RooPlot*)sigD2_CMS_th1x_fit_s->Clone("sigD2_fit_frame");
  RooPlot* clone_sigD2_fitb_frame = (RooPlot*)sigD2_CMS_th1x_fit_b->Clone("sigD2_fitb_frame");

  RooPlot* clone_sigD3_prefit_frame = (RooPlot*)sigD3_CMS_th1x_prefit->Clone("sigD3_prefit_frame");
  RooPlot* clone_sigD3_fit_frame = (RooPlot*)sigD3_CMS_th1x_fit_s->Clone("sigD3_fit_frame");
  RooPlot* clone_sigD3_fitb_frame = (RooPlot*)sigD3_CMS_th1x_fit_b->Clone("sigD3_fitb_frame");

  RooPlot* clone_sigD4_prefit_frame = (RooPlot*)sigD4_CMS_th1x_prefit->Clone("sigD4_prefit_frame");
  RooPlot* clone_sigD4_fit_frame = (RooPlot*)sigD4_CMS_th1x_fit_s->Clone("sigD4_fit_frame");
  RooPlot* clone_sigD4_fitb_frame = (RooPlot*)sigD4_CMS_th1x_fit_b->Clone("sigD4_fitb_frame");

  // ===================================================================================

  // Get the N6 for each fit

  //TPaveText *text_N6_tt_D1 = new TPaveText(0.352682,0.785388,0.626392,0.875761,"blNDC");
  //text_N6_tt_D1->AddText(TString::Format("N6 = %.1f",full_ttMC_th1_D1->GetBinContent(1)));
  //text_N6_tt_D1->SetFillStyle(0);
  //text_N6_tt_D1->SetBorderSize(1);
  //text_N6_tt_D1->SetLineColor(0);

  // ===================================================================================

  // Get the chi^2 for each fit

  TPaveText *chiSquare_sigD1_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD1_prefit->AddText(TString::Format("Chi Square = %0.2f",sigD1_CMS_th1x_prefit->chiSquare()));
  chiSquare_sigD1_prefit->SetFillStyle(0);
  chiSquare_sigD1_prefit->SetBorderSize(1);
  chiSquare_sigD1_prefit->SetLineColor(0);

  TPaveText *chiSquare_sigD1_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD1_fit->AddText(TString::Format("Chi Square = %0.2f",sigD1_CMS_th1x_fit_s->chiSquare()));
  chiSquare_sigD1_fit->SetFillStyle(0);
  chiSquare_sigD1_fit->SetBorderSize(1);
  chiSquare_sigD1_fit->SetLineColor(0);

  TPaveText *chiSquare_sigD1_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD1_fitb->AddText(TString::Format("Chi Square = %0.2f",sigD1_CMS_th1x_fit_b->chiSquare()));
  chiSquare_sigD1_fitb->SetFillStyle(0);
  chiSquare_sigD1_fitb->SetBorderSize(1);
  chiSquare_sigD1_fitb->SetLineColor(0);


  TPaveText *chiSquare_sigD2_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD2_prefit->AddText(TString::Format("Chi Square = %0.2f",sigD2_CMS_th1x_prefit->chiSquare()));
  chiSquare_sigD2_prefit->SetFillStyle(0);
  chiSquare_sigD2_prefit->SetBorderSize(1);
  chiSquare_sigD2_prefit->SetLineColor(0);

  TPaveText *chiSquare_sigD2_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD2_fit->AddText(TString::Format("Chi Square = %0.2f",sigD2_CMS_th1x_fit_s->chiSquare()));
  chiSquare_sigD2_fit->SetFillStyle(0);
  chiSquare_sigD2_fit->SetBorderSize(1);
  chiSquare_sigD2_fit->SetLineColor(0);

  TPaveText *chiSquare_sigD2_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD2_fitb->AddText(TString::Format("Chi Square = %0.2f",sigD2_CMS_th1x_fit_b->chiSquare()));
  chiSquare_sigD2_fitb->SetFillStyle(0);
  chiSquare_sigD2_fitb->SetBorderSize(1);
  chiSquare_sigD2_fitb->SetLineColor(0);


  TPaveText *chiSquare_sigD3_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD3_prefit->AddText(TString::Format("Chi Square = %0.2f",sigD3_CMS_th1x_prefit->chiSquare()));
  chiSquare_sigD3_prefit->SetFillStyle(0);
  chiSquare_sigD3_prefit->SetBorderSize(1);
  chiSquare_sigD3_prefit->SetLineColor(0);

  TPaveText *chiSquare_sigD3_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD3_fit->AddText(TString::Format("Chi Square = %0.2f",sigD3_CMS_th1x_fit_s->chiSquare()));
  chiSquare_sigD3_fit->SetFillStyle(0);
  chiSquare_sigD3_fit->SetBorderSize(1);
  chiSquare_sigD3_fit->SetLineColor(0);

  TPaveText *chiSquare_sigD3_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD3_fitb->AddText(TString::Format("Chi Square = %0.2f",sigD3_CMS_th1x_fit_b->chiSquare()));
  chiSquare_sigD3_fitb->SetFillStyle(0);
  chiSquare_sigD3_fitb->SetBorderSize(1);
  chiSquare_sigD3_fitb->SetLineColor(0);


  TPaveText *chiSquare_sigD4_prefit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD4_prefit->AddText(TString::Format("Chi Square = %0.2f",sigD4_CMS_th1x_prefit->chiSquare()));
  chiSquare_sigD4_prefit->SetFillStyle(0);
  chiSquare_sigD4_prefit->SetBorderSize(1);
  chiSquare_sigD4_prefit->SetLineColor(0);

  TPaveText *chiSquare_sigD4_fit = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD4_fit->AddText(TString::Format("Chi Square = %0.2f",sigD4_CMS_th1x_fit_s->chiSquare()));
  chiSquare_sigD4_fit->SetFillStyle(0);
  chiSquare_sigD4_fit->SetBorderSize(1);
  chiSquare_sigD4_fit->SetLineColor(0);

  TPaveText *chiSquare_sigD4_fitb = new TPaveText(0.351476,0.492295,0.701791,0.582667,"blNDC");
  chiSquare_sigD4_fitb->AddText(TString::Format("Chi Square = %0.2f",sigD4_CMS_th1x_fit_b->chiSquare()));
  chiSquare_sigD4_fitb->SetFillStyle(0);
  chiSquare_sigD4_fitb->SetBorderSize(1);
  chiSquare_sigD4_fitb->SetLineColor(0);


  // ===================================================================================

  // Make pull plots of each fit
  // The pull is defined as (curve-histogram +/- error_on_histogram) / error_on_histogram
  //   so it is expected that for the case of symmetric errors the error bars on the pull are +/- 1. 

  RooHist* hpull_sigD1_prefit = sigD1_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_sigD1_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD1_prefit->addPlotable(hpull_sigD1_prefit,"P");

  RooHist* hpull_sigD1_fit = sigD1_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_sigD1_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD1_fit->addPlotable(hpull_sigD1_fit,"P");

  RooHist* hpull_sigD1_fitb = sigD1_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_sigD1_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD1_fitb->addPlotable(hpull_sigD1_fitb,"P");


  RooHist* hpull_sigD2_prefit = sigD2_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_sigD2_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD2_prefit->addPlotable(hpull_sigD2_prefit,"P");

  RooHist* hpull_sigD2_fit = sigD2_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_sigD2_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD2_fit->addPlotable(hpull_sigD2_fit,"P");

  RooHist* hpull_sigD2_fitb = sigD2_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_sigD2_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD2_fitb->addPlotable(hpull_sigD2_fitb,"P");


  RooHist* hpull_sigD3_prefit = sigD3_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_sigD3_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD3_prefit->addPlotable(hpull_sigD3_prefit,"P");

  RooHist* hpull_sigD3_fit = sigD3_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_sigD3_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD3_fit->addPlotable(hpull_sigD3_fit,"P");

  RooHist* hpull_sigD3_fitb = sigD3_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_sigD3_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD3_fitb->addPlotable(hpull_sigD3_fitb,"P");


  RooHist* hpull_sigD4_prefit = sigD4_CMS_th1x_prefit->pullHist();
  RooPlot* framePull_sigD4_prefit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD4_prefit->addPlotable(hpull_sigD4_prefit,"P");

  RooHist* hpull_sigD4_fit = sigD4_CMS_th1x_fit_s->pullHist();
  RooPlot* framePull_sigD4_fit = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD4_fit->addPlotable(hpull_sigD4_fit,"P");

  RooHist* hpull_sigD4_fitb = sigD4_CMS_th1x_fit_b->pullHist();
  RooPlot* framePull_sigD4_fitb = w->var("CMS_th1x")->frame(Title("Pull"));
  framePull_sigD4_fitb->addPlotable(hpull_sigD4_fitb,"P");


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

  TCanvas* c_sigD1_prefit = new TCanvas("c_sigD1_prefit","Signal D1 prefit",900,900);
  c_sigD1_prefit->Divide(2,2);
  c_sigD1_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD1_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  sigD1_CMS_th1x_prefit->Draw();
  chiSquare_sigD1_prefit->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_sigD1_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD1_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD1_prefit_frame->SetMinimum(0.001);
  clone_sigD1_prefit_frame->SetMaximum(100000);
  //clone_sigD1_prefit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_sigD1_prefit_frame->Draw();
  c_sigD1_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD1_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD1_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD1_prefit->Draw();
  //c_sigD1_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD1_prefit->SetStats(0);
  //hcorr_sigD1_prefit->SetMarkerSize(2.4);
  //hcorr_sigD1_prefit->Draw("COLZ,TEXT");
  c_sigD1_prefit->SaveAs("fit_report_sigD1_prefit.png");

  TCanvas* c_sigD1_fit = new TCanvas("c_sigD1_fit","Signal D1 fit",900,900);
  c_sigD1_fit->Divide(2,2);
  c_sigD1_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD1_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  sigD1_CMS_th1x_fit_s->Draw();
  chiSquare_sigD1_fit->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_sigD1_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD1_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD1_fit_frame->SetMinimum(0.001);
  clone_sigD1_fit_frame->SetMaximum(100000);
  //clone_sigD1_fit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_sigD1_fit_frame->Draw();
  c_sigD1_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD1_fit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD1_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD1_fit->Draw();
  //c_sigD1_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD1_fit->SetStats(0);
  //hcorr_sigD1_fit->SetMarkerSize(2.4);
  //hcorr_sigD1_fit->Draw("COLZ,TEXT");
  c_sigD1_fit->SaveAs("fit_report_sigD1_fit.png");

  TCanvas* c_sigD1_fitb = new TCanvas("c_sigD1_fitb","D1 fit bkg only",900,900);
  c_sigD1_fitb->Divide(2,2);
  c_sigD1_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD1_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  sigD1_CMS_th1x_fit_b->Draw();
  chiSquare_sigD1_fitb->Draw("same");  
  //text_N6_tt_D1->Draw("same");
  c_sigD1_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD1_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD1_fitb_frame->SetMinimum(0.001);
  clone_sigD1_fitb_frame->SetMaximum(100000);
  //clone_sigD1_fit_frame->remove("bkgMC_tt_D1_paramBox");
  clone_sigD1_fitb_frame->Draw();
  c_sigD1_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD1_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_sigD1_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD1_fitb->Draw();
  //c_sigD1_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD1_fitb->SetStats(0);
  //hcorr_sigD1_fitb->SetMarkerSize(2.4);
  //hcorr_sigD1_fitb->Draw("COLZ,TEXT");
  c_sigD1_fitb->SaveAs("fit_report_sigD1_fitb.png");


  // ------------ D2 -----------

  TCanvas* c_sigD2_prefit = new TCanvas("c_sigD2_prefit","Signal D2 prefit",900,900);
  c_sigD2_prefit->Divide(2,2);
  c_sigD2_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD2_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  sigD2_CMS_th1x_prefit->Draw();
  chiSquare_sigD2_prefit->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_sigD2_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD2_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD2_prefit_frame->SetMinimum(0.001);
  clone_sigD2_prefit_frame->SetMaximum(100000);
  //clone_sigD2_prefit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_sigD2_prefit_frame->Draw();
  c_sigD2_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD2_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD2_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD2_prefit->Draw();
  //c_sigD2_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD2_prefit->SetStats(0);
  //hcorr_sigD2_prefit->SetMarkerSize(2.4);
  //hcorr_sigD2_prefit->Draw("COLZ,TEXT");
  c_sigD2_prefit->SaveAs("fit_report_sigD2_prefit.png");

  TCanvas* c_sigD2_fit = new TCanvas("c_sigD2_fit","Signal D2 fit",900,900);
  c_sigD2_fit->Divide(2,2);
  c_sigD2_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD2_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  sigD2_CMS_th1x_fit_s->Draw();
  chiSquare_sigD2_fit->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_sigD2_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD2_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD2_fit_frame->SetMinimum(0.001);
  clone_sigD2_fit_frame->SetMaximum(100000);
  //clone_sigD2_fit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_sigD2_fit_frame->Draw();
  c_sigD2_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD2_fit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD2_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD2_fit->Draw();
  //c_sigD2_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD2_fit->SetStats(0);
  //hcorr_sigD2_fit->SetMarkerSize(2.4);
  //hcorr_sigD2_fit->Draw("COLZ,TEXT");
  c_sigD2_fit->SaveAs("fit_report_sigD2_fit.png");

  TCanvas* c_sigD2_fitb = new TCanvas("c_sigD2_fitb","D2 fit bkg only",900,900);
  c_sigD2_fitb->Divide(2,2);
  c_sigD2_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD2_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  sigD2_CMS_th1x_fit_b->Draw();
  chiSquare_sigD2_fitb->Draw("same");  
  //text_N6_tt_D2->Draw("same");
  c_sigD2_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD2_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD2_fitb_frame->SetMinimum(0.001);
  clone_sigD2_fitb_frame->SetMaximum(100000);
  //clone_sigD2_fit_frame->remove("bkgMC_tt_D2_paramBox");
  clone_sigD2_fitb_frame->Draw();
  c_sigD2_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD2_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_sigD2_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD2_fitb->Draw();
  //c_sigD2_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD2_fitb->SetStats(0);
  //hcorr_sigD2_fitb->SetMarkerSize(2.4);
  //hcorr_sigD2_fitb->Draw("COLZ,TEXT");
  c_sigD2_fitb->SaveAs("fit_report_sigD2_fitb.png");


  // ------------ D3 -----------

  TCanvas* c_sigD3_prefit = new TCanvas("c_sigD3_prefit","Signal D3 prefit",900,900);
  c_sigD3_prefit->Divide(2,2);
  c_sigD3_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD3_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  sigD3_CMS_th1x_prefit->Draw();
  chiSquare_sigD3_prefit->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_sigD3_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD3_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD3_prefit_frame->SetMinimum(0.001);
  clone_sigD3_prefit_frame->SetMaximum(100000);
  //clone_sigD3_prefit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_sigD3_prefit_frame->Draw();
  c_sigD3_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD3_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD3_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD3_prefit->Draw();
  //c_sigD3_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD3_prefit->SetStats(0);
  //hcorr_sigD3_prefit->SetMarkerSize(2.4);
  //hcorr_sigD3_prefit->Draw("COLZ,TEXT");
  c_sigD3_prefit->SaveAs("fit_report_sigD3_prefit.png");

  TCanvas* c_sigD3_fit = new TCanvas("c_sigD3_fit","Signal D3 fit",900,900);
  c_sigD3_fit->Divide(2,2);
  c_sigD3_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD3_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  sigD3_CMS_th1x_fit_s->Draw();
  chiSquare_sigD3_fit->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_sigD3_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD3_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD3_fit_frame->SetMinimum(0.001);
  clone_sigD3_fit_frame->SetMaximum(100000);
  //clone_sigD3_fit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_sigD3_fit_frame->Draw();
  c_sigD3_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD3_fit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD3_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD3_fit->Draw();
  //c_sigD3_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD3_fit->SetStats(0);
  //hcorr_sigD3_fit->SetMarkerSize(2.4);
  //hcorr_sigD3_fit->Draw("COLZ,TEXT");
  c_sigD3_fit->SaveAs("fit_report_sigD3_fit.png");

  TCanvas* c_sigD3_fitb = new TCanvas("c_sigD3_fitb","D3 fit bkg only",900,900);
  c_sigD3_fitb->Divide(2,2);
  c_sigD3_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD3_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  sigD3_CMS_th1x_fit_b->Draw();
  chiSquare_sigD3_fitb->Draw("same");  
  //text_N6_tt_D3->Draw("same");
  c_sigD3_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD3_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD3_fitb_frame->SetMinimum(0.001);
  clone_sigD3_fitb_frame->SetMaximum(100000);
  //clone_sigD3_fit_frame->remove("bkgMC_tt_D3_paramBox");
  clone_sigD3_fitb_frame->Draw();
  c_sigD3_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD3_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_sigD3_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD3_fitb->Draw();
  //c_sigD3_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD3_fitb->SetStats(0);
  //hcorr_sigD3_fitb->SetMarkerSize(2.4);
  //hcorr_sigD3_fitb->Draw("COLZ,TEXT");
  c_sigD3_fitb->SaveAs("fit_report_sigD3_fitb.png");


  // ------------ D4 -----------

  TCanvas* c_sigD4_prefit = new TCanvas("c_sigD4_prefit","Signal D4 prefit",900,900);
  c_sigD4_prefit->Divide(2,2);
  c_sigD4_prefit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD4_CMS_th1x_prefit->GetYaxis()->SetTitleOffset(2.0);
  sigD4_CMS_th1x_prefit->Draw();
  chiSquare_sigD4_prefit->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_sigD4_prefit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD4_prefit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD4_prefit_frame->SetMinimum(0.001);
  clone_sigD4_prefit_frame->SetMaximum(100000);
  //clone_sigD4_prefit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_sigD4_prefit_frame->Draw();
  c_sigD4_prefit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD4_prefit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD4_prefit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD4_prefit->Draw();
  //c_sigD4_prefit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD4_prefit->SetStats(0);
  //hcorr_sigD4_prefit->SetMarkerSize(2.4);
  //hcorr_sigD4_prefit->Draw("COLZ,TEXT");
  c_sigD4_prefit->SaveAs("fit_report_sigD4_prefit.png");

  TCanvas* c_sigD4_fit = new TCanvas("c_sigD4_fit","Signal D4 fit",900,900);
  c_sigD4_fit->Divide(2,2);
  c_sigD4_fit->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD4_CMS_th1x_fit_s->GetYaxis()->SetTitleOffset(2.0);
  sigD4_CMS_th1x_fit_s->Draw();
  chiSquare_sigD4_fit->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_sigD4_fit->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD4_fit_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD4_fit_frame->SetMinimum(0.001);
  clone_sigD4_fit_frame->SetMaximum(100000);
  //clone_sigD4_fit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_sigD4_fit_frame->Draw();
  c_sigD4_fit->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD4_fit->SetYTitle("(Fit-Data)/errData");
  framePull_sigD4_fit->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD4_fit->Draw();
  //c_sigD4_fit->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD4_fit->SetStats(0);
  //hcorr_sigD4_fit->SetMarkerSize(2.4);
  //hcorr_sigD4_fit->Draw("COLZ,TEXT");
  c_sigD4_fit->SaveAs("fit_report_sigD4_fit.png");

  TCanvas* c_sigD4_fitb = new TCanvas("c_sigD4_fitb","D4 fit bkg only",900,900);
  c_sigD4_fitb->Divide(2,2);
  c_sigD4_fitb->cd(1);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  sigD4_CMS_th1x_fit_b->GetYaxis()->SetTitleOffset(2.0);
  sigD4_CMS_th1x_fit_b->Draw();
  chiSquare_sigD4_fitb->Draw("same");  
  //text_N6_tt_D4->Draw("same");
  c_sigD4_fitb->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(0);
  gPad->SetLogy();
  clone_sigD4_fitb_frame->GetYaxis()->SetTitleOffset(1.4);
  clone_sigD4_fitb_frame->SetMinimum(0.001);
  clone_sigD4_fitb_frame->SetMaximum(100000);
  //clone_sigD4_fit_frame->remove("bkgMC_tt_D4_paramBox");
  clone_sigD4_fitb_frame->Draw();
  c_sigD4_fitb->cd(3);
  gPad->SetLeftMargin(0.25);
  gPad->SetRightMargin(0.05);
  framePull_sigD4_fitb->SetYTitle("(Fit-Data)/errData");
  framePull_sigD4_fitb->GetYaxis()->SetTitleOffset(1.4);
  framePull_sigD4_fitb->Draw();
  //c_sigD4_fitb->cd(4);
  //gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.15);
  //gStyle->SetOptTitle(1);
  //hcorr_sigD4_fitb->SetStats(0);
  //hcorr_sigD4_fitb->SetMarkerSize(2.4);
  //hcorr_sigD4_fitb->Draw("COLZ,TEXT");
  c_sigD4_fitb->SaveAs("fit_report_sigD4_fitb.png");



  TCanvas* c_fit_s_corr = new TCanvas("c_fit_s_corr","signal region correlation matrix",900,900);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_fit->SetStats(0);
  hcorr_fit->SetMarkerSize(1.0);
  hcorr_fit->Draw("COLZ,TEXT");
  c_fit_s_corr->SaveAs("fit_report_correlation.png");


  TCanvas* c_fit_b_corr = new TCanvas("c_fit_b_corr","signal region correlation matrix bkg only",900,900);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gStyle->SetOptTitle(1);
  hcorr_fitb->SetStats(0);
  hcorr_fitb->SetMarkerSize(1.0);
  hcorr_fitb->Draw("COLZ,TEXT");
  c_fit_b_corr->SaveAs("fit_report_correlation_b.png");

}
