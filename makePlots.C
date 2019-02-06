// Aron Soha
// October 2017
//
// Compile this macro using:
//  root -l
//  .L makePlots.C
// and then call using
//  makePlots();
//
// or compile and call in one step:
//  root -l
//  .x makePlots.C+

#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TPaveText.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLine.h"
#include "TLegend.h"
//#include "/uscms/home/soha/scripts/CMS_lumi.h"

void makePlots(const string today = "Jan17_2019", const string filedir = "fit_results_v5_Jan17_2019", const string year = "2017", const string model = "RPV", const string fitType = "AsymptoticLimits") {

  // =============================================================
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.03);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
  
  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

  // =============================================================


  // Settings for CMS_lumi macro
  //gROOT->LoadMacro("/uscms/home/soha/scripts/CMS_lumi.C");
  bool writeExtraText = false;
  TString extraText = "Preliminary";
  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  int iPos = 11;

 // **** Set the following each time before running ****

  //string channel = "RPV";
  //string channel = "SYY";
  //string channel = "SHH";

  //string channel = "RPVL";
  //string channel = "SYYL";
  //string channel = "SHHL";

  //string channel = "RPV2L";
  //string channel = "SYY2L";
  //string channel = "SHH2L";

  bool STOP_PAIRS = true;
  bool SYST_HALF = false;
  bool SYST_DOUBLE = false;
  
  string date = "Jan17_19";
  string ssave = "---------------------------------------------------------------";
  string ssave_base = "./sigBrLim_"+model+"_"+year+"_";

  bool DRAW_OBS = true;
  bool DRAW_LOGOS = false;

  // ****************************************************

  int W = 800;
  int H = 600;
  int W_ref = 800; 
  int H_ref = 600; 
  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  gStyle->SetCanvasDefH      (H);
  gStyle->SetCanvasDefW      (W);
  gStyle->SetTitleOffset( 1, "y");

  // *****
  // Extract limit results from set of root files produced by Higgs Combine tool
  //std::vector<double> xpoints = {300,350,400,450,500,550,600,650,700,750,800,850,900};  // mass points
  std::vector<double> xpoints = {350,450,550,650,750,850};  // mass points
  const int npoints = xpoints.size();

  // Arrays for storing results
  // The following are the values of r from the fitter, where r is
  //  the number of signal events / number of expected signal events
  std::vector<double> limits_obs(npoints,0);
  std::vector<double> limits_obsErr(npoints,0);
  std::vector<double> limits_m2s(npoints,0);
  std::vector<double> limits_m1s(npoints,0);
  std::vector<double> limits_mean(npoints,0);
  std::vector<double> limits_p1s(npoints,0);
  std::vector<double> limits_p2s(npoints,0);

  // Loop over mass points
  for (int i=0; i<npoints; i++) 
  {
    const string& mass = std::to_string(int(xpoints[i]));
    const string& fitter_file = filedir+"/"+model+"_"+mass+"_"+year+"/higgsCombine"+year+"."+fitType+".mH"+mass+".MODEL"+model+".root";
    
    std::cout << fitter_file << std::endl;
    // Load the root file and read the tree and leaves
    TFile *f = new TFile(fitter_file.c_str());
    TTreeReader reader("limit");
    TTreeReaderValue<double> lim(reader,"limit");
    TTreeReaderValue<double> lim_err(reader,"limitErr");

    reader.Next();
    limits_m2s[i] = *lim;
    reader.Next();
    limits_m1s[i] = *lim;
    reader.Next();
    limits_mean[i] = *lim;
    reader.Next();
    limits_p1s[i] = *lim;
    reader.Next();
    limits_p2s[i] = *lim;
    reader.Next();
    limits_obs[i] = *lim;
    limits_obsErr[i] = *lim_err;

    f->Close();
  }

  // Define the cross section times branchingn ratio (called sigBr here):
  std::vector<double> sigBr;
  std::vector<double> sigBr1SPpercent;
  std::vector<double> sigBr1SMpercent;
  
  // cross sections and uncertainties from
  //  https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVstopsbottom
  //const std::vector<double>& stop_pair_Br =           {8.51615, 3.78661,  1.83537, 0.948333, 0.51848, 0.296128, 0.174599, 0.107045, 0.0670476, 0.0431418, 0.0283338, 0.0189612, 0.0128895};
  //const std::vector<double>& stop_pair_Br1SPpercent = {13.9223, 13.6877, 13.6985, 13.4559,  13.3797, 13.2687,  13.2074,  12.9232,  13.3429,   13.7455,   14.171,    14.702,    15.2026};
  //const std::vector<double>& stop_pair_Br1SMpercent = {13.9223, 13.6877, 13.6985, 13.4559,  13.3797, 13.2687,  13.2074,  12.9232,  13.3429,   13.7455,   14.171,    14.702,    15.2026};
  const std::vector<double>& stop_pair_Br =           {3.78661,  0.948333, 0.296128, 0.107045, 0.0431418, 0.0189612};
  const std::vector<double>& stop_pair_Br1SPpercent = {13.6877, 13.4559,  13.2687,  12.9232,  13.7455,   14.702};
  const std::vector<double>& stop_pair_Br1SMpercent = {13.6877, 13.4559,  13.2687,  12.9232,  13.7455,   14.702};

  // For stop pair production
  if (STOP_PAIRS) 
  {
      sigBr = stop_pair_Br;
      sigBr1SPpercent = stop_pair_Br1SPpercent;
      sigBr1SMpercent = stop_pair_Br1SMpercent;
  }

  std::vector<double> sigBr1SP(npoints,0);
  std::vector<double> sigBr1SM(npoints,0);
  for (int i=0; i<npoints; i++) {
    sigBr1SP[i] = sigBr[i]*sigBr1SPpercent[i]/100.0;
    sigBr1SM[i] = sigBr[i]*sigBr1SMpercent[i]/100.0;
  }

  bool projectingRLimitLogY = true;
  double projectingXmin = 300, projectingXmax = 900;
  double projectingRLimitYmin = 0.01, projectingRLimitYmax = 100;
  std::string projectingRLimitXYtitles = ";m_{#tilde{t}} [GeV]; 95% CL limit on #sigma#timesBr [pb]";

  if (STOP_PAIRS) ssave = ssave_base+today+"_CLs";
  if (SYST_HALF) ssave = ssave_base+"half_"+today+"_CLs";
  if (SYST_DOUBLE) ssave = ssave_base+"double_"+today+"_CLs";

  std::vector<double> limits_exp(npoints,0);
  for(int n=0; n<npoints; n++)
  {
    limits_m2s[n]=limits_m2s[n]*sigBr[n];
    limits_m1s[n]=limits_m1s[n]*sigBr[n];
    limits_mean[n]=limits_mean[n]*sigBr[n];
    limits_p1s[n]=limits_p1s[n]*sigBr[n];
    limits_p2s[n]=limits_p2s[n]*sigBr[n];

    limits_exp[n]=limits_mean[n];
    limits_obs[n]=limits_obs[n]*sigBr[n];
  }

  TPaveText* pt = nullptr;
  if (DRAW_LOGOS) 
  {
    pt = new TPaveText(0.46, 0.75, 0.6, 0.95, "ndc");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
  } 
  else 
  {
    pt = new TPaveText(0.46, 0.75, 0.6, 0.95, "ndc");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
  }


  if (model=="RPV")
    pt->AddText("pp #rightarrow #tilde{t}#tilde{t} #rightarrow (t #tilde{#chi}^{0}_{1})(t #tilde{#chi}^{0}_{1}) #rightarrow (t jjj)(t jjj), Lepton");
  else if (model=="SYY")
    pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SYY coupling, Lepton");
  else if (model=="SHH")
    pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SHH coupling, Lepton");

  // if (channel=="RPV")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t} #rightarrow (t #tilde{#chi}^{0}_{1})(t #tilde{#chi}^{0}_{1}) #rightarrow (t jjj)(t jjj)");
  // else if (channel=="RPVL")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t} #rightarrow (t #tilde{#chi}^{0}_{1})(t #tilde{#chi}^{0}_{1}) #rightarrow (t jjj)(t jjj), Lepton");
  // else if (channel=="RPV2L")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t} #rightarrow (t #tilde{#chi}^{0}_{1})(t #tilde{#chi}^{0}_{1}) #rightarrow (t jjj)(t jjj), 2 Leptons");
  // else if (channel=="SYY")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SYY coupling");
  // else if (channel=="SYYL")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SYY coupling, Lepton");
  // else if (channel=="SHH")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SHH coupling");
  // else if (channel=="SHHL")
  //   pt->AddText("pp #rightarrow #tilde{t}#tilde{t}, SHH coupling, Lepton");

  cout << "npoints = " << npoints << endl;
  for (int n=0; n<npoints; n++)
  {
    cout << "limitx_m2s = " << limits_m2s[n] << endl;
    cout << "limitx_m1s = " << limits_m1s[n] << endl;
    cout << "limitx_exp = " << limits_exp[n] << endl;
    cout << "limitx_p1s = " << limits_p1s[n] << endl;
    cout << "limitx_p2s = " << limits_p2s[n] << endl;
    cout << "limitx_obs = " << limits_obs[n] << endl;
    cout << "xpoints = " << xpoints[n] << endl;
    cout << endl;
  }

  // PlotWithBelts lb(limits_m1s, limits_p1s, limits_m2s, limits_p2s,limits_exp, limits_obs, npoints, xpoints, ssave+"_1", pt, projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
  //lb.plot();

  TCanvas *cCanvas = new TCanvas(ssave.c_str(),"Canvas");
  TString stmp = "hframe"; stmp += ssave;
  TH1F *hframe= new TH1F(stmp, projectingRLimitXYtitles.c_str(), 1000, projectingXmin, projectingXmax);
  hframe->SetMinimum(projectingRLimitYmin);
  hframe->SetMaximum(projectingRLimitYmax);
  hframe->SetStats(0);
  hframe->SetFillStyle(1);
  hframe->Draw(" ");
        
  cCanvas->SetLogy(projectingRLimitLogY);


  TGraph *grMean = new TGraph(npoints, xpoints.data(), limits_exp.data()); 
  TGraph *grYellow = new TGraph(2*npoints);
  for(int n=0; n<npoints; n++){
    grYellow->SetPoint(n, xpoints[n], limits_p2s[n]);
    grYellow->SetPoint(npoints+n, xpoints[npoints-n-1], limits_m2s[npoints-n-1]);
  }
  grYellow->SetFillColor(kYellow);
  grYellow->SetLineColor(kYellow);
  grYellow->Draw("f");

  TGraph *grGreen = new TGraph(2*npoints);
  for(int n=0; n<npoints; n++)
  {
    grGreen->SetPoint(n, xpoints[n], limits_p1s[n]);
    grGreen->SetPoint(npoints+n, xpoints[npoints-n-1], limits_m1s[npoints-n-1]);
  }
  grGreen->SetFillColor(kGreen);
  grGreen->SetLineColor(kGreen);
  grGreen->Draw("f");

  TLine *lineOne = new TLine(projectingXmin,1, projectingXmax, 1);
  lineOne->SetLineWidth(2);
  lineOne->SetLineStyle(1);
  lineOne->SetLineColor(kBlack);
  lineOne->Draw("same");

  grMean->SetLineWidth(2);
  grMean->SetLineStyle(2);
  grMean->SetLineColor(kBlue);
  grMean->SetMarkerSize(0);
  grMean->Draw("lp");

  TGraph* grObs = nullptr;
  if (DRAW_OBS) 
  {
      grObs=new TGraph(npoints, xpoints.data(), limits_obs.data());
      grObs->SetMarkerStyle(20);
      grObs->SetMarkerColor(kBlack);
      grObs->SetLineWidth(2);
      grObs->SetLineWidth(1);
      grObs->SetLineColor(kBlack);
      grObs->Draw("lp");
  }

  pt->Draw();


  TPaveText *br1 = new TPaveText(0.207, 0.258, 0.357, 0.258+0.066, "ndc");
  br1->SetBorderSize(0);
  br1->SetFillStyle(0);
  br1->SetTextAlign(12);
  br1->SetTextFont(42);
  br1->SetTextSize(0.035);
  if (model=="RPV")
    br1->AddText("B(#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}) = 1.0");
  else if (model=="SYY")
    br1->AddText("B(#tilde{t} #rightarrow t#tilde{S}g) = 1.0");
  else if (model=="SHH")
    br1->AddText("B(#tilde{t} #rightarrow t#tilde{S}) = 1.0");
  // if (channel=="RPV" || channel=="RPVL" || channel=="RPV2L")
  //   br1->AddText("B(#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}) = 1.0");
  // else if (channel=="SYY" || channel=="SYYL")
  //   br1->AddText("B(#tilde{t} #rightarrow t#tilde{S}g) = 1.0");
  // else if (channel=="SHH" || channel=="SHHL")
  //   br1->AddText("B(#tilde{t} #rightarrow t#tilde{S}) = 1.0");
  br1->Draw("same");

  TPaveText *br2 = new TPaveText(0.207, 0.204, 0.357, 0.204+0.066, "ndc");
  br2->SetBorderSize(0);
  br2->SetFillStyle(0);
  br2->SetTextAlign(12);
  br2->SetTextFont(42);
  br2->SetTextSize(0.035);
  if (model=="RPV")
    br2->AddText("B(#tilde{#chi}^{0}_{1} #rightarrow jjj) = 1.0");
  else if (model=="SYY")
    br2->AddText("B(#tilde{S} #rightarrow S#tilde{G}) = 1.0");
  else if (model=="SHH")
    br2->AddText("B(#tilde{S} #rightarrow S#tilde{G}) = 1.0");
  // if (channel=="RPV" || channel=="RPVL" || channel=="RPV2L")
  //   br2->AddText("B(#tilde{#chi}^{0}_{1} #rightarrow jjj) = 1.0");
  // else if (channel=="SYY" || channel=="SYYL")
  //   br2->AddText("B(#tilde{S} #rightarrow S#tilde{G}) = 1.0");
  // else if (channel=="SHH" || channel=="SHHL")
  //   br2->AddText("B(#tilde{S} #rightarrow S#tilde{G}) = 1.0");
  br2->Draw("same");

  TPaveText *br3 = new TPaveText(0.207, 0.150, 0.357, 0.150+0.066, "ndc");
  br3->SetBorderSize(0);
  br3->SetFillStyle(0);
  br3->SetTextAlign(12);
  br3->SetTextFont(42);
  br3->SetTextSize(0.035);
  if (model=="SYY")
    br3->AddText("B(S #rightarrow gg) = 1.0");
  if (model=="SHH")
    br3->AddText("B(S #rightarrow bb) = 1.0");
  // if (channel=="SYY" || channel=="SYYL")
  //   br3->AddText("B(S #rightarrow gg) = 1.0");
  // if (channel=="SHH" || channel=="SHHL")
  //   br3->AddText("B(S #rightarrow bb) = 1.0");
  br3->Draw("same");

  if (DRAW_LOGOS) 
  {
    cCanvas->SetLeftMargin( L/W );
    cCanvas->SetRightMargin( R/W );
    cCanvas->SetTopMargin( T/H );
    cCanvas->SetBottomMargin( B/H );
  }

  TGraphAsymmErrors *grTheoryErr = new TGraphAsymmErrors(npoints,xpoints.data(),sigBr.data(),nullptr,nullptr,sigBr1SM.data(),sigBr1SP.data());
  grTheoryErr->SetLineColor(2);
  grTheoryErr->SetLineWidth(2);
  //grTheoryErr->SetMarkerColor(2);
  //grTheoryErr->SetMarkerStyle(20);
  grTheoryErr->SetFillColor(42);
  TGraph *grTheory = new TGraph(npoints,xpoints.data(),sigBr.data());
  grTheory->SetLineColor(2);
  grTheory->SetLineWidth(2);

  TLegend *legend = new TLegend(0.45, 0.62, 0.849, 0.814);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextAlign(12);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(grGreen,"Expected #pm 1 #sigma", "f");
  legend->AddEntry(grYellow,"Expected #pm 2 #sigma", "f");
  if(DRAW_OBS) legend->AddEntry(grObs,"Observed limit", "lp");

  if (model=="RPV")
    legend->AddEntry(grTheoryErr,"RPV stop (with xsec uncertainty)", "lf");
  else if (model=="SYY")
    legend->AddEntry(grTheoryErr,"SYY (with xsec uncertainty)", "lf");
  else if (model=="SHH")
    legend->AddEntry(grTheoryErr,"SHH (with xsec uncertainty)", "lf");
  // if (channel=="RPV" || channel=="RPVL" || channel=="RPV2L")
  //   legend->AddEntry(grTheoryErr,"RPV stop (with xsec uncertainty)", "lf");
  // else if (channel=="SYY" || channel=="SYYL")
  //   legend->AddEntry(grTheoryErr,"SYY (with xsec uncertainty)", "lf");
  // else if (channel=="SHH" || channel=="SHHL")
  //   legend->AddEntry(grTheoryErr,"SHH (with xsec uncertainty)", "lf");
  legend->Draw();

  grTheoryErr->Draw("3,same");
  grTheory->Draw("l,same");

  // redraw mean, so that it appears over the signal lines
  grMean->Draw("lp");

  // redraw obs, so that it appears over the expected lines
  if (DRAW_OBS) grObs->Draw("lp");

  lineOne->Delete();
  //lb.getLine()->SetLineColor(kRed);
  //lb.getLine()->Draw();

  cCanvas->cd();
  //gr = new TGraph(npoints,xpoints,limits_obs);
  //gr->SetLineWidth(2);
  //gr->SetMarkerStyle(1);
  //gr->Draw("ALP,same");

  if (DRAW_LOGOS) 
  {
    // ALS:    CMS_lumi(cCanvas,iPeriod,iPos);
    cCanvas->Update();
    cCanvas->RedrawAxis();
    cCanvas->GetFrame()->Draw();
  }

  TPaveText *cmstext=NULL;
  if (!DRAW_LOGOS) 
  {
    TPaveText* cmstext = new TPaveText(0.308789, 0.958188, 0.806516, 0.996516, "ndc");
    cmstext->SetBorderSize(0);
    cmstext->SetFillStyle(0);
    cmstext->SetTextAlign(12);
    cmstext->SetTextFont(42);
    cmstext->SetTextSize(0.035);

    if (year=="2016") 
    {
      if (STOP_PAIRS)  cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 fb^{-1}");
      if (SYST_HALF)   cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 fb^{-1}, HALF SYST");
      if (SYST_DOUBLE) cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 fb^{-1}, DOUBLE SYST");
    }
    else if (year=="2017") 
    {
      if (STOP_PAIRS) cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=41.5 fb^{-1}");
      if (SYST_HALF)  cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=41.5 fb^{-1}, HALF SYST");
      if (SYST_DOUBLE)cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=41.5 fb^{-1}, DOUBLE SYST");
    }
    else if (year=="Combo")
    {
      if (STOP_PAIRS) cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 + 41.5 fb^{-1}");
      if (SYST_HALF)  cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 + 41.5 fb^{-1}, HALF SYST");
      if (SYST_DOUBLE)cmstext->AddText("CMS Preliminary, #sqrt{s}=13 TeV, L_{Int}=35.9 + 41.5 fb^{-1}, DOUBLE SYST");
    }
    cmstext->Draw("same");
  }

  std::string seps = filedir+"/"+ssave+".eps";
  std::string sgif = filedir+"/"+ssave+".gif";
  std::string sroot = filedir+"/"+ssave+".root";
  std::string spdf = filedir+"/"+ssave+".pdf";
  cCanvas->Print(sroot.c_str());
  cCanvas->Print(seps.c_str());
  cCanvas->Print(sgif.c_str());
  cCanvas->Print(spdf.c_str());
}
