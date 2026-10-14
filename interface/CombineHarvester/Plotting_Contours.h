#ifndef ICHiggsTauTau_Core_Plotting_Countours_h
#define ICHiggsTauTau_Core_Plotting_Countours_h

// Mostly copied from here:
// https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/plotting/contours2D.cxx
// with minor modifications
#include "Plotting.h"

#include <iostream>
#include <math.h>
#include "TGraph.h"
#include "TTree.h"
#include "TCut.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

// ---------------------------------------------------------------------------
// Function declarations
// ---------------------------------------------------------------------------
TGraph *bestFit(TTree *t, TString x, TString y, TCut cut);

TH2 *treeToHist2D(TTree *t, TString x, TString y, TString name, TCut cut,
                  double xmin, double xmax, double ymin, double ymax, int xbins,
                  int ybins);

TList *contourFromTH2(TH2 *h2in, double threshold, int minPoints = 20);

TH2D *frameTH2D(TH2D *in, double threshold);

void styleMultiGraph(TList *tmg, int lineColor, int lineWidth, int lineStyle);

void styleMultiGraphMarker(TList *tmg, int markerColor, int markerSize,
                           int markerStyle);

void drawContours(TPad *pad, TString file, TString x_title, TString y_title,
                  TLegend *leg);
// ---------------------------------------------------------------------------
// Function implementations
// ---------------------------------------------------------------------------
/** Make a 2D contour plot from the output of MultiDimFit
 *
 * Inputs:
 *  - gFile should point to the TFile containing the 'limit' TTree
 *  - xvar should be the variable to use on the X axis, with xbins bins in the [xmin, xmax] range
 *  - yvar should be the variable to use on the Y axis, with ybins bins in the [ymin, ymax] range
 *  - (smx, smy) are the coordinates at which to put a diamond representing the SM expectation
 *  - if fOut is not null, then the output objects will be saved to fOut:
 *     - the 2D histogram will be saved as a TH2 with name name+"_h2d"
 *     - the 68% CL contour will be saved as a TList of TGraphs with name name+"_c68"
 *     - the 95% CL contour will be saved as a TList of TGraphs with name name+"_c95"
 *     - the 99.7% CL contour will be saved as a TList of TGraphs with name name+"_c997"
 *     - the best fit point will be saved as a TGraph with name name+"_best"
 *
 * Notes:
 *     - it's up to you to make sure that the binning you use for this plot matches the one used
 *       when running MultiDimFit (but you can just plot a subset of the points; e.g. if you had
 *       100x100 points in [-1,1]x[-1,1] you can make a 50x50 plot for [0,1]x[0,1])
 *     - the 99.7 contour is not plotted by default
 *     - the SM marker is not saved
*/
void contour2D(TString filein, TString fileout, TString xvar, int xbins, float xmin, float xmax, TString yvar,
               int ybins, float ymin, float ymax, TString name = "contour2D") {
  TFile fin(filein);
  fin.cd();
  TTree *tree = (TTree *)gFile->Get("limit");
  TH2 *hist2d = treeToHist2D(tree, xvar, yvar, "h2d", "", xmin, xmax, ymin,
                             ymax, xbins, ybins);
  hist2d->SetContour(200);
  hist2d->GetZaxis()->SetRangeUser(0, 21);
  TGraph *fit = bestFit(tree, xvar, yvar, "");
  TList *c68 = contourFromTH2(hist2d, 2.30);
  TList *c95 = contourFromTH2(hist2d, 5.99);
  TList *c997 = contourFromTH2(hist2d, 11.83);
  TFile fout(fileout, "RECREATE");
  fout.cd();
  hist2d->SetName(name + "_h2d");
  fout.WriteTObject(hist2d, 0);
  fit->SetName(name + "_best");
  fout.WriteTObject(fit, 0);
  c68->SetName(name + "_c68");
  fout.WriteTObject(c68, 0, "SingleKey");
  c95->SetName(name + "_c95");
  fout.WriteTObject(c95, 0, "SingleKey");
  c997->SetName(name + "_c997");
  fout.WriteTObject(c997, 0, "SingleKey");
  c68->At(1)->Write("second");
  fout.Close();
}



TGraph *bestFit(TTree *t, TString x, TString y, TCut cut) {
  int nfind = t->Draw(y + ":" + x, cut + "deltaNLL == 0");
  if (nfind == 0) {
    TGraph *gr0 = new TGraph(1);
    gr0->SetPoint(0, -999, -999);
    gr0->SetMarkerStyle(34);
    gr0->SetMarkerSize(2.0);
    return gr0;
  } else {
    TGraph *gr0 = (TGraph *)gROOT->FindObject("Graph")->Clone();
    gr0->SetMarkerStyle(34);
    gr0->SetMarkerSize(2.0);
    if (gr0->GetN() > 1) gr0->Set(1);
    return gr0;
  }
}

TH2 *treeToHist2D(TTree *t, TString x, TString y, TString name, TCut cut,
                  double xmin, double xmax, double ymin, double ymax, int xbins,
                  int ybins) {
  t->Draw(Form("2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", y.Data(),
               x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax),
          cut + "deltaNLL != 0", "PROF");
  TH2 *prof = (TH2 *)gROOT->FindObject(name + "_prof");
  TH2D *h2d = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
  for (int ix = 1; ix <= xbins; ++ix) {
    for (int iy = 1; iy <= ybins; ++iy) {
      double z = prof->GetBinContent(ix, iy);
      if (z == TMath::Infinity()) z = 999;
      if (z != z)
        z = (name.Contains("bayes") ? 0 : 999);  // protect agains NANs
      h2d->SetBinContent(ix, iy, z);
    }
  }
  h2d->GetXaxis()->SetTitle(x);
  h2d->GetYaxis()->SetTitle(y);
  h2d->SetDirectory(0);
  return h2d;
}

TList *contourFromTH2(TH2 *h2in, double threshold, int minPoints) {
  std::cout << "Getting contour at threshold " << threshold << " from "
            << h2in->GetName() << std::endl;
  // http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
  Double_t contours[1];
  contours[0] = threshold;
  if (h2in->GetNbinsX() * h2in->GetNbinsY() > 10000) minPoints = 50;
  if (h2in->GetNbinsX() * h2in->GetNbinsY() <= 100) minPoints = 10;

  TH2D *h2 = frameTH2D((TH2D *)h2in, threshold);

  h2->SetContour(1, contours);

  // Draw contours as filled regions, and Save points
  h2->Draw("CONT Z LIST");
  gPad->Update();  // Needed to force the plotting and retrieve the contours in
                   // TGraphs

  // Get Contours
  TObjArray *conts =
      (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  TList *contLevel = NULL;

  if (conts == NULL || conts->GetSize() == 0) {
    printf("*** No Contours Were Extracted!\n");
    return 0;
  }

  TList *ret = new TList();
  for (int i = 0; i < conts->GetSize(); i++) {
    contLevel = (TList *)conts->At(i);
    // printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
    for (int j = 0, n = contLevel->GetSize(); j < n; ++j) {
      TGraph *gr1 = (TGraph *)contLevel->At(j);
      // printf("\t Graph %d has %d points\n", j, gr1->GetN());
      if (gr1->GetN() > minPoints) ret->Add(gr1->Clone());
      // break;
    }
  }
  return ret;
}

TH2D *frameTH2D(TH2D *in, double threshold) {
  // NEW LOGIC:
  //   - pretend that the center of the last bin is on the border if the frame
  //   - add one tiny frame with huge values
  double frameValue = 1000;
  if (TString(in->GetName()).Contains("bayes")) frameValue = -1000;

  Double_t xw = in->GetXaxis()->GetBinWidth(1);
  Double_t yw = in->GetYaxis()->GetBinWidth(1);

  Int_t nx = in->GetNbinsX();
  Int_t ny = in->GetNbinsY();

  Double_t x0 = in->GetXaxis()->GetXmin();
  Double_t x1 = in->GetXaxis()->GetXmax();

  Double_t y0 = in->GetYaxis()->GetXmin();
  Double_t y1 = in->GetYaxis()->GetXmax();
  Double_t xbins[999], ybins[999];
  double eps = 0.1;
  double mult = 5.;

  xbins[0] = x0 - eps * xw - xw * mult;
  xbins[1] = x0 + eps * xw - xw * mult;
  for (int ix = 2; ix <= nx; ++ix) xbins[ix] = x0 + (ix - 1) * xw;
  xbins[nx + 1] = x1 - eps * xw + 0.5 * xw * mult;
  xbins[nx + 2] = x1 + eps * xw + xw * mult;

  ybins[0] = y0 - eps * yw - yw * mult;
  ybins[1] = y0 + eps * yw - yw * mult;
  for (int iy = 2; iy <= ny; ++iy) ybins[iy] = y0 + (iy - 1) * yw;
  ybins[ny + 1] = y1 - eps * yw + yw * mult;
  ybins[ny + 2] = y1 + eps * yw + yw * mult;

  TH2D *framed =
      new TH2D(Form("%s framed", in->GetName()),
               Form("%s framed", in->GetTitle()), nx + 2, xbins, ny + 2, ybins);

  // Copy over the contents
  for (int ix = 1; ix <= nx; ix++) {
    for (int iy = 1; iy <= ny; iy++) {
      framed->SetBinContent(1 + ix, 1 + iy, in->GetBinContent(ix, iy));
    }
  }
  // Frame with huge values
  nx = framed->GetNbinsX();
  ny = framed->GetNbinsY();
  for (int ix = 1; ix <= nx; ix++) {
    framed->SetBinContent(ix, 1, frameValue);
    framed->SetBinContent(ix, ny, frameValue);
  }
  for (int iy = 2; iy <= ny - 1; iy++) {
    framed->SetBinContent(1, iy, frameValue);
    framed->SetBinContent(nx, iy, frameValue);
  }

  return framed;
}

// void styleMultiGraph(TList *tmg, int lineColor, int lineWidth, int lineStyle) {
//   for (int i = 0; i < tmg->GetSize(); ++i) {
//     TGraph *g = (TGraph *)tmg->At(i);
//     g->SetLineColor(lineColor);
//     g->SetLineWidth(lineWidth);
//     g->SetLineStyle(lineStyle);
//   }
// }

// void styleMultiGraphMarker(TList *tmg, int markerColor, int markerSize,
//                            int markerStyle) {
//   for (int i = 0; i < tmg->GetSize(); ++i) {
//     TGraph *g = (TGraph *)tmg->At(i);
//     g->SetMarkerColor(markerColor);
//     g->SetMarkerSize(markerSize);
//     g->SetMarkerStyle(markerStyle);
//   }
// }

void drawContours(TPad *pad, TString file, TString x_title, TString y_title,
                  TLegend *leg) {
  TFile infile(file);
  TH2D *h2d = (TH2D *)gDirectory->Get("contour2D_h2d");

  h2d->GetXaxis()->SetTitle(x_title);
  h2d->GetYaxis()->SetTitle(y_title);

  pad->cd();
  h2d->Draw("AXIS");

  TGraph *bestfit = (TGraph *)gDirectory->Get("contour2D_best");
  std::vector<TGraph *> c68;
  std::vector<TGraph *> c95;

  TList * l68 = (TList *)gDirectory->Get("contour2D_c68");
  for (int i = 0; i < l68->GetSize(); ++i) {
    c68.push_back((TGraph *)l68->At(i));
  }

  TList * l95 = (TList *)gDirectory->Get("contour2D_c95");
  for (int i = 0; i < l95->GetSize(); ++i) {
    c95.push_back((TGraph *)l95->At(i));
  }


  // TGraph *c68 = (TGraph *)((TList *)gDirectory->Get("contour2D_c68"))->At(0);
  // TGraph *c95 = (TGraph *)((TList *)gDirectory->Get("contour2D_c95"))->At(0);

  for (unsigned i = 0; i < c95.size(); ++i) {
    c95[i]->SetLineStyle(1);
    c95[i]->SetLineColor(kBlack);
    // c95[i]->SetLineWidth(3);
    c95[i]->SetFillColor(kBlue - 10);
    c95[i]->Draw("F SAME");
  }

  for (unsigned i = 0; i < c68.size(); ++i) {
    c68[i]->SetLineStyle(1);
    c68[i]->SetLineColor(kBlack);
    // c68[i]->SetLineWidth(3);
    c68[i]->SetFillColor(kBlue - 8);
    c68[i]->Draw("F SAME");
  }

  for (unsigned i = 0; i < c95.size(); ++i) c95[i]->Draw("L SAME");
  for (unsigned i = 0; i < c68.size(); ++i) c68[i]->Draw("L SAME");

  bestfit->SetMarkerStyle(34);
  bestfit->SetMarkerSize(3.0);
  bestfit->SetMarkerColor(kBlack);
  // bestfit->Draw("P SAME");
  if (leg) {
    leg->AddEntry(c95[0], "95% CL", "F");
    leg->AddEntry(c68[0], "68% CL", "F");
    }
  // leg->AddEntry(bestfit, "Best fit", "P");

  // leg->Draw("SAME");

  // TLatex *title_latex = new TLatex();
  // title_latex->SetNDC();
  // title_latex->SetTextSize(0.035);
  // title_latex->SetTextFont(62);
  // title_latex->SetTextAlign(31);
  // double height = 0.94;
  // title_latex->DrawLatex(0.95, height,
  //                        "19.7 fb^{-1} (8 TeV) + 4.9 fb^{-1} (7 TeV)");
  // title_latex->SetTextAlign(11);
  // title_latex->DrawLatex(0.17, height, "H#rightarrow#tau#tau");
  // title_latex->SetTextSize(0.08);
  // title_latex->DrawLatex(0.22, 0.84, "#mu#tau_{h}");
  // canv.SaveAs(output);
}

#endif
