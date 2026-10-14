#ifndef Core_Plotting_h
#define Core_Plotting_h

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "TStyle.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TList.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TColor.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "THStack.h"



/**
 * @name Pad and axis histogram creation
 */
/**@{*/

/**
 * Just creates a single pad filling the entire canvas
 *
 * Exists just to provide a similar output interface to TwoPadSplit()
 *
 * @return A vector containing a single pointer to the newly created TPad
 */
std::vector<TPad*> OnePad();

/**
 * Create two pads, split horizontally, on the current canvas split
 *
 * Note that both pads will actually cover the full area of the canvas. We
 * create the appearance of two separated pads by adjusting the top and bottom
 * margins, and then making the background of each pad transparent.
 *
 * @param split_point Fraction of the canvas that should be occuplied by the
 * lower pad.
 * @param gap_low The top margin of the lower pad (with respect to the split
 * point), given as a fraction of the input canvas height
 * @param gap_high The bottom margin of the upper pad (with respect to the
 * split point), given as a fraction of the input canvas height
 *
 * @return A vector of pointers to the new pads: [0] = upper, [1] = lower. The
 * current pad (`gPad`) will be set to the upper one.
 */
std::vector<TPad*> TwoPadSplit(double split_point, double gap_low,
                               double gap_high);

/**
 * Create an empty TH1 from another TH1 for drawing the axes
 *
 * A copy of the input TH1 is made, and all bin contents and errors are
 * removed.
 *
 * @param src TH1 input
 * @param xmin adjust the minimum of the drawn range of the x-axis
 * @param xmax adjust the maximum of the drawn range of the x-axis (only if
 * `xmax > xmin`)
 *
 * @return An empty TH1 with axes built from the input TH1
 */
TH1* CreateAxisHist(TH1* src, double xmin = 0, double xmax = -1);

/**
 * Create multiple axis TH1s from another TH1
 *
 * Just calls CreateAxisHist(TH1*, double, double) multiple times.
 *
 * @param n Number of histograms to make
 * @param src TH1 input
 * @param xmin adjust the minimum of the drawn range of the x-axis
 * @param xmax adjust the maximum of the drawn range of the x-axis (only if
 * `xmax > xmin`)
 */
std::vector<TH1*> CreateAxisHists(unsigned n, TH1* src, double xmin = 0,
                                  double xmax = -1);

/**
 * Create an empty TH1 from a TGraph for drawing the axes
 *
 * This function just draws the input TGraph on a temporary TCanvas and
 * extracts the internal TH1 that ROOT creates to do this. In this process the
 * x- and y-axis ranges are computing automatically.
 *
 * @param src TGraph input
 * @param at_limits By default the x-axis range will go beyond the min and max
 * x-values found in the TGraph. Set this option to true to fix the range to
 * these min/max values.
 *
 * @return An empty TH1 with axes built from the input TGraph
 */
TH1* CreateAxisHist(TGraph* src, bool at_limits);

/**
 * Create multiple axis TH1s from a TGraph
 *
 * Just calls CreateAxisHist(TGraph*, bool) multiple times.
 *
 * @param n Number of histograms to make
 * @param src TGraph input
 * @param at_limits By default the x-axis range will go beyond the min and max
 * x-values found in the TGraph. Set this option to true to fix the range to
 * these min/max values.
 */
std::vector<TH1*> CreateAxisHists(unsigned n, TGraph* src,  bool at_limits);


/**
 * Finds the TH1 used to draw the axes on a given TPad
 *
 * Should be able to find a TH1 regardless of whether a TH1, THStack or TGraph
 * was used to define the axes.
 *
 * @param pad The TPad to search
 *
 * @return A pointer to the axis TH1 if it exists, a NULL pointer otherwise
 */
TH1* GetAxisHist(TPad *pad);
/**@}*/


/**
 * @name Object extraction and manipulation
 */
/**@{*/

/**
 * Create a new histogram by dividing one by the other
 *
 * The errors on the two inputs are optionally propagated. Note that these are
 * just combined in quadrature, thus assuming the inputs are uncorrelated.
 * This does not hold for efficiency calculations.
 *
 * @param num Numerator histogram
 * @param den Denominator histogram
 * @param num_err Propagate errors on the numerator
 * @param den_err Propagate errors on the denominator
 *
 * @return Pointer to the newly created TH1
 */
TH1* MakeRatioHist(TH1* num, TH1* den, bool num_err, bool den_err);


/**
 * Create a TGraph from entries in a TTree
 *
 * Note that the graph entries will be filled in whichever order they are
 * found in the TTree. It may be desirable to call `TGraph::Sort` on the
 * resulting object.
 *
 * @param tree Input TTree
 * @param xvar Branch or expression for the x-values
 * @param yvar Branch or expression for the y-values
 * @param selection Optional cut string to apply to each entry
 */
TGraph TGraphFromTree(TTree* tree, TString const& xvar, TString const& yvar,
                      TString const& selection = "");

/**
 * Create a TGraph2D from entries in a TTree
 *
 * Note that the graph entries will be filled in whichever order they are
 * found in the TTree. It may be desirable to call `TGraph::Sort` on the
 * resulting object.
 *
 * @param tree Input TTree
 * @param xvar Branch or expression for the x-values
 * @param yvar Branch or expression for the y-values
 * @param zvar Branch or expression for the y-values
 * @param selection Optional cut string to apply to each entry
 */
TGraph2D TGraph2DFromTree(TTree* tree, TString const& xvar, TString const& yvar,
                          TString const& zvar, TString const& selection = "");

/**
 * Shift all the graph y-values upwards such that there are no negative values
 * and the minimum point is at zero.
 *
 * Useful for fixing graphs of log-likelihood scans extracted from `combine`
 * where the true minimum was not found correctly in the initial fit.
 */
void ReZeroTGraph(TGraph *gr);
/**@}*/



/**
 * @name Plot decoration (colours, styles, text and legends)
 */
/**@{*/

/**
 * Create a transparent version of a colour
 *
 * @param color ROOT colour code
 * @param alpha Transparency (0 = fully transparent, 1 = opaque)
 *
 * @return The code for the new TColor
 */
int CreateTransparentColor(int color, float alpha);

/**
 * Set a few style options for a two-pad setup used to show a data-MC
 * comparison and ratio plot.
 *
 * Use this version of the function if the axis histograms haven't been drawn
 * yet. This function will:
 *
 *    - remove the x-axis title and labels on the main plot
 *    - fix the axis tick-lengths to be consistent between both plots
 *    - enable x- and y-axis ticks on the ratio plot
 *    - Set the range and title options for the ratio plot y-axis
 *
 * @param pads Vector of TPad pointers where [0] = upper pad, [1] = lower pad
 * @param upper Pointer to the upper pad axis TH1
 * @param lower Pointer to the lower pad axis TH1
 * @param y_title Title for the y-axis of the ratio plot
 * @param y_centered Centre (vertically) the y-axis title on the ratio plot.
 * @param y_min Set the y-axis minimum in the ratio plot
 * @param y_max Set the y-axis maximum in the ratio plot
 */
void SetupTwoPadSplitAsRatio(std::vector<TPad*> const& pads, TH1* upper,
                             TH1* lower, TString y_title, bool y_centered,
                             float y_min, float y_max);

/**
 * Set a few style options for a two-pad setup used to show a data-MC
 * comparison and ratio plot.
 *
 * This function should be called after the axis histograms have been created
 * in both pads (e.g. after some TH1, TGraph etc. has been drawn). It will:
 *
 *  - remove the x-axis title and labels on the main plot
 *  - fix the axis tick-lengths to be consistent between both plots
 *  - enable x- and y-axis ticks on the ratio plot
 *  - Set the range and title options for the ratio plot y-axis
 *
 * @param pads Vector of TPad pointers where [0] = upper pad, [1] = lower pad
 * @param y_title Title for the y-axis of the ratio plot
 * @param y_centered Centre (vertically) the y-axis title on the ratio plot.
 * @param y_min Set the y-axis minimum in the ratio plot
 * @param y_max Set the y-axis maximum in the ratio plot
 */
void SetupTwoPadSplitAsRatio(std::vector<TPad*> const& pads, TString y_title,
                             bool y_centered, float y_min, float y_max);

/**
 * Sets standard x- and y-axis titles with given units
 *
 * Will create an x-axis title of the form `<variable> (<units>)` and a y-axis
 * title of the form `Events / <bin_width> <units>`.
 *
 * @param xaxis Target x-axis
 * @param yaxis Target y-axis
 * @param var String of the variable name
 * @param units String of the units, e.g. `GeV`. Leave empty if the variable
 * is dimensionless
 */
void StandardAxes(TAxis* xaxis, TAxis* yaxis, TString var, TString units);

/**
 * Sets standard x- and y-axis titles when histograms are drawn divided by bin
 * width
 *
 * Will create an x-axis title of the form `<variable> (<units>)` and a y-axis
 * title of the form `dN/d<var> (1/<units>)`.
 *
 * @param xaxis Target x-axis
 * @param yaxis Target y-axis
 * @param var String of the variable name
 * @param units String of the units, e.g. `GeV`
 */
void UnitAxes(TAxis* xaxis, TAxis* yaxis, TString var, TString units);

/**
 * Create a legend with fixed height, width and positioning
 *
 * Legend will be created in one three fixed positions:
 *
 *    1. Top-left corner
 *    2. Top centered
 *    3. Top-right corner
 *    4. Bottom-left corner
 *    5. Bottom centered
 *    6. Bottom-right corner
 *
 * @param width Legend width as a fraction of the pad
 * @param height Legend height as a fraction of the pad
 * @param pos Position index (see above)
 * @param offset Position offset from edges of frame
 *
 * @return Pointer to newly created TLegend
 */
TLegend* PositionedLegend(double width, double height, int pos, double offset);

/**
 * Use an existing TLine to draw a new horizontal line across the current
 * frame
 *
 * Will find the TH1 used to draw the axes on the pad, get the min/max
 * x-values and then draw a TLine at the specified y-value
 *
 * @param pad TPad where the line should be drawn
 * @param line Existing TLine already customised with the desired colour,
 * width and style
 * @param yval The y-value (given in the axis units) where the line should be
 * drawn
 *
 * @return A pointer to the newly created TLine
 */
TLine* DrawHorizontalLine(TPad *pad, TLine* line, double yval);

/**
 * Draw text in the top-margin region of a TPad
 *
 * Specified text will be latex-rendered in one of three positions:
 *
 *    1. Left aligned with left side of frame
 *    2. Centre aligned along top of frame
 *    3. Right aligned with right side of fram
 *
 * Note also that the text size will be set automatically to fill the top
 * margin areas.
 *
 * @param pad Input Tpad
 * @param text Text to draw
 * @param align Alignment option (see above)
 */
void DrawTitle(TPad* pad, TString text, int align);

/**
 * Draw the CMS logo and subtitle in the new style
 *
 * Adapted from the example [here](https://ghm.web.cern.ch/ghm/plots/).
 *
 * Position setting:
 *
 *     iPos=0  : out of frame
 *     iPos=11 : top-left, left-aligned
 *     iPos=33 : top-right, right-aligned
 *     iPos=22 : center, centered
 *
 * mode generally :
 *
 *     iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
 *
 * @param pad The target TPad
 * @param cmsText Usually just "CMS"
 * @param extraText The subtitle, e.g. "Preliminary"
 * @param iPosX Sets the position and alignment of the text (see above)
 * @param relPosX Sets the x-position of the text relative to the frame
 * @param relPosY Sets the y-position of the text relative to the frame
 * @param relExtraDY Sets the position of the *extraText* relative to the
 * *cmsText*
 */
void DrawCMSLogo(TPad* pad, TString cmsText, TString extraText, int iPosX,
                 float relPosX, float relPosY, float relExtraDY);

/**
 * Call DrawCMSLogo with some sensible defaults for the positioning parameters
 *
 * @param pad The target TPad
 * @param cmsText Usually just "CMS"
 * @param extraText The subtitle, e.g. "Preliminary"
 * @param iPosX Sets the position and alignment of the text
 */
void DrawCMSLogo(TPad* pad, TString cmsText, TString extraText, int iPosX);
/**@}*/



/**
 * @name Axis adjustment and overlap checks
 */
/**@{*/


/**
 * Find the maximum value of all drawn objects in a given x-axis range
 *
 * Currently only objects derived from TH1 and TGraph are examined.
 *
 * @param pad TPad to examine
 * @param x_min Minimum x-axis value to start the scan
 * @param x_max Maximum x-axis alue to conclude the scan
 *
 * @return Value on y-axis where the maximum is found
 */
double GetPadYMax(TPad* pad, double x_min, double x_max);

/**
 * Call GetPadYMax(TPad*, double, double) with the full range of the
 * given TPad
 */
double GetPadYMax(TPad* pad);

/**
 * Adjusts the y-axis maximum on the pad such that the specified y-value is
 * positioned a fixed fraction from this new maximum.
 *
 * Should work for both linear and log-scale y-axes. This function is most
 * useful in conjunction with GetPadYMax() to supply the y-axis value to fix,
 * thus ensuring that a fixed fraction at the top of the frame is free of any
 * histogram or graph.
 *
 * @param pad TPad to adjust
 * @param fix_y The position on the y-axis that will be fixed
 * @param fraction The fraction of the y-axis that *fix_y* should be
 * positioned (measured from the top of the frame)
 */
void FixTopRange(TPad *pad, double fix_y, double fraction);


/**
 * Modify the pad y-axis range to ensure there is at least a given gap between
 * a particular TBox and the highest drawn element in the x-axis range of this
 * box.
 *
 * Most useful after drawing a legend somewhere near the top of the frame,
 * which may have overlapped with some of the histograms or graphs that had
 * already been drawn. This function will check for the overlap, then adjust
 * the y-axis range to ensure there is a gap. Note that if a sufficient (or
 * larger) gap already exists, the y-axis range will not be adujusted.
 *
 * @param pad TPad to modify
 * @param box Pointer to the reference TBox (e.g. could be a TLegend)
 * @param frac Desired fraction of the frame height between the bottom of the
 * TBox and the top of the drawn objects.
 */
void FixBoxPadding(TPad *pad, TBox *box, double frac);

/**
 * Just re-draws the axes on the current TPad
 */
void FixOverlay();
/**@}*/



//----------------------------------------------------------------------------
// Pad and axis histogram creation
//----------------------------------------------------------------------------
std::vector<TPad*> OnePad() {
  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->Draw();
  pad->cd();
  std::vector<TPad*> result;
  result.push_back(pad);
  return result;
}

std::vector<TPad*> TwoPadSplit(double split_point, double gap_low,
                               double gap_high) {
  TPad* upper = new TPad("upper", "upper", 0., 0., 1., 1.);
  upper->SetBottomMargin(split_point + gap_high);
  upper->SetFillStyle(4000);
  upper->Draw();
  TPad* lower = new TPad("lower", "lower", 0., 0., 1., 1.);
  lower->SetTopMargin(1 - split_point + gap_low);
  lower->SetFillStyle(4000);
  lower->Draw();
  upper->cd();
  std::vector<TPad*> result;
  result.push_back(upper);
  result.push_back(lower);
  return result;
}

TH1* CreateAxisHist(TH1* src, double xmin, double xmax) {
  TH1 *res = reinterpret_cast<TH1F*>(src->Clone());
  res->Reset();
  if (xmax > xmin) {
    res->SetAxisRange(xmin, xmax);
  }
  return res;
}

std::vector<TH1*> CreateAxisHists(unsigned n, TH1* src, double xmin,
                                  double xmax) {
  std::vector<TH1*> result;
  for (unsigned i = 0; i < n; ++i) {
    result.push_back(CreateAxisHist(src, xmin, xmax));
  }
  return result;
}

TH1* CreateAxisHist(TGraph* src, bool at_limits) {
  TVirtualPad *backup = gPad;
  TCanvas tmp;
  tmp.cd();
  src->Draw("AP");
  TH1 * result = (TH1*)src->GetHistogram()->Clone();
  if (at_limits) {
    double min = 0;
    double max = 0;
    double x = 0;
    double y = 0;
    src->GetPoint(0, x, y);
    min = x;
    max = x;
    for (int i = 1; i < src->GetN(); ++i) {
      src->GetPoint(i, x, y);
      if (x < min) min = x;
      if (x > max) max = x;
    }
    result->GetXaxis()->SetLimits(min, max);
  }
  gPad = backup;
  return result;
}

std::vector<TH1*> CreateAxisHists(unsigned n, TGraph* src,  bool at_limits) {
  std::vector<TH1*> result;
  for (unsigned i = 0; i < n; ++i) {
    result.push_back(CreateAxisHist(src, at_limits));
  }
  return result;
}

TH1* GetAxisHist(TPad *pad) {
  TList* pad_obs = pad->GetListOfPrimitives();
  if (!pad_obs) return NULL;
  TIter next(pad_obs);
  TObject* obj;
  TH1 *hobj = 0;
  while ((obj = next())) {
    if (obj->InheritsFrom(TH1::Class())) {
      hobj = reinterpret_cast<TH1*>(obj);
      break;
    }
    if (obj->InheritsFrom(TMultiGraph::Class())) {
      TMultiGraph *mg = reinterpret_cast<TMultiGraph*>(obj);
      if (mg) return mg->GetHistogram();
    }
    if (obj->InheritsFrom(TGraph::Class())) {
      TGraph *g = reinterpret_cast<TGraph*>(obj);
      if (g) return g->GetHistogram();
    }
    if (obj->InheritsFrom(THStack::Class())) {
      THStack *hs = reinterpret_cast<THStack*>(obj);
      if (hs) return hs->GetHistogram();
    }
  }
  return hobj;
}


//----------------------------------------------------------------------------
// Object extraction and manipulation
//----------------------------------------------------------------------------
TH1* MakeRatioHist(TH1* num, TH1* den, bool num_err, bool den_err) {
  TH1* result = reinterpret_cast<TH1*>(num->Clone());
  if (!num_err) {
    for (int i = 1; i <= result->GetNbinsX(); ++i) result->SetBinError(i, 0.);
  }
  TH1* den_fix = reinterpret_cast<TH1*>(den->Clone());
  if (!den_err) {
    for (int i = 1; i <= den_fix->GetNbinsX(); ++i) den_fix->SetBinError(i, 0.);
  }
  result->Divide(den_fix);
  delete den_fix;
  return result;
}

TGraph TGraphFromTree(TTree * tree, TString const& xvar, TString const& yvar,
                     TString const& selection) {
  tree->Draw(xvar+":"+yvar, selection, "goff");
  TGraph gr(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2());
  return gr;
}

TGraph2D TGraph2DFromTree(TTree* tree, TString const& xvar, TString const& yvar,
                          TString const& zvar, TString const& selection) {
  tree->Draw(xvar+":"+yvar+":"+zvar, selection, "goff");
  TGraph2D gr(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(), tree->GetV3());
  return gr;
}

void ReZeroTGraph(TGraph *gr) {
  unsigned n = gr->GetN();
  double min = 0.;
  double x = 0;
  double y = 0;
  for (unsigned i = 0; i < n; ++i) {
    gr->GetPoint(i, x, y);
    if (y < min) min = y;
  }
  for (unsigned i = 0; i < n; ++i) {
    gr->GetPoint(i, x, y);
    y = y - min;
    gr->SetPoint(i, x, y);
  }
}

//----------------------------------------------------------------------------
// Plot decoration (colours, styles, text and legends)
//----------------------------------------------------------------------------
int CreateTransparentColor(int color, float alpha) {
  TColor* adapt = gROOT->GetColor(color);
  int new_idx = gROOT->GetListOfColors()->GetSize() + 1;
  TColor* trans = new TColor(new_idx, adapt->GetRed(), adapt->GetGreen(),
                             adapt->GetBlue(), "", alpha);
  trans->SetName(Form("userColor%i", new_idx));
  return new_idx;
}

void SetupTwoPadSplitAsRatio(std::vector<TPad*> const& pads,
                      TString y_title, bool y_centered, float y_min,
                      float y_max) {
  TH1 *upper = GetAxisHist(pads[0]);
  TH1 *lower = GetAxisHist(pads[1]);
  SetupTwoPadSplitAsRatio(pads, upper, lower, y_title, y_centered, y_min,
                          y_max);
}

void SetupTwoPadSplitAsRatio(std::vector<TPad*> const& pads, TH1* upper,
                             TH1* lower, TString y_title, bool y_centered,
                             float y_min, float y_max) {
  upper->GetXaxis()->SetTitle("");
  upper->GetXaxis()->SetLabelSize(0);
  double upper_h = 1. - pads[0]->GetTopMargin() - pads[0]->GetBottomMargin();
  double lower_h = 1. - pads[1]->GetTopMargin() - pads[1]->GetBottomMargin();
  lower->GetYaxis()->SetTickLength(gStyle->GetTickLength() *
                                   (upper_h / lower_h));
  pads[1]->SetTickx(1);
  pads[1]->SetTicky(1);
  lower->GetYaxis()->SetTitle(y_title);
  lower->GetYaxis()->CenterTitle(y_centered);
  if (y_max > y_min) {
    lower->SetMinimum(y_min);
    lower->SetMaximum(y_max);
  }
}

void StandardAxes(TAxis* xaxis, TAxis* yaxis, TString var, TString units) {
  double width = xaxis->GetBinWidth(1);
  TString w_label = TString::Format("%.1f", width);
  if (units == "") {
    xaxis->SetTitle(var);
    yaxis->SetTitle("Events / " + w_label);
  } else {
    xaxis->SetTitle(var + " (" + units + ")");
    yaxis->SetTitle("Events / " + w_label + " " + units);
  }
}


void UnitAxes(TAxis* xaxis, TAxis* yaxis, TString var, TString units) {
  xaxis->SetTitle(var + " (" + units + ")");
  yaxis->SetTitle("dN/d"+ var + " (1/" + units + ")");
}

TLegend* PositionedLegend(double width, double height, int pos, double offset) {
  double o = offset;
  double w = width;
  double h = height;
  double l = gPad->GetLeftMargin();
  double t = gPad->GetTopMargin();
  double b = gPad->GetBottomMargin();
  double r = gPad->GetRightMargin();
  TLegend* leg = 0;
  if (pos == 1) {
    leg = new TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, "", "NBNDC");
  }
  if (pos == 2) {
    double c = l + 0.5 * (1 - l - r);
    leg = new TLegend(c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, "",
                      "NBNDC");
  }
  if (pos == 3) {
    leg = new TLegend(1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, "",
                      "NBNDC");
  }
  if (pos == 4) {
    leg = new TLegend(l + o, b + o, l + o + w, b + o + h, "", "NBNDC");
  }
  if (pos == 5) {
    double c = l + 0.5 * (1 - l - r);
    leg = new TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, "",
                      "NBNDC");
  }
  if (pos == 6) {
    leg = new TLegend(1 - r - o - w, b + o, 1 - r - o, b + o + h, "",
                      "NBNDC");
  }
  return leg;
}

TLine* DrawHorizontalLine(TPad *pad, TLine* line, double yval) {
  TH1 *axis = GetAxisHist(pad);
  std::cout << "pad: " << pad << "\n";
  double xmin = axis->GetXaxis()->GetXmin();
  double xmax = axis->GetXaxis()->GetXmax();
  return line->DrawLine(xmin, yval, xmax, yval);
}

void DrawTitle(TPad* pad, TString text, int align) {
  TVirtualPad *pad_backup = gPad;
  pad->cd();
  float t = pad->GetTopMargin();
  float l = pad->GetLeftMargin();
  float r = pad->GetRightMargin();

  float pad_ratio = (static_cast<float>(pad->GetWh()) * pad->GetAbsHNDC()) /
                    (static_cast<float>(pad->GetWw()) * pad->GetAbsWNDC());
  if (pad_ratio < 1.) pad_ratio = 1.;

  float textSize = 0.6;
  float textOffset = 0.2;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  latex.SetTextSize(textSize * t * pad_ratio);
  double y_off = 1 - t + textOffset * t;
  if (align == 1) latex.SetTextAlign(11);
  if (align == 1) latex.DrawLatex(l, y_off, text);
  if (align == 2) latex.SetTextAlign(21);
  if (align == 2) latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text);
  if (align == 3) latex.SetTextAlign(31);
  if (align == 3) latex.DrawLatex(1 - r, y_off, text);
  pad_backup->cd();
}

void DrawCMSLogo(TPad* pad, TString cmsText, TString extraText, int iPosX,
                 float relPosX, float relPosY, float relExtraDY) {
  TVirtualPad *pad_backup = gPad;
  pad->cd();
  float cmsTextFont = 61;  // default is helvetic-bold

  bool writeExtraText = extraText.Length() > 0;
  float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextOffset = 0.2;
  float cmsTextSize = 0.8;
  // float cmsTextOffset    = 0.1;  // only used in outOfFrame version

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize = 0.76;

  TString lumi_13TeV = "20.1 fb^{-1}";
  TString lumi_8TeV = "19.7 fb^{-1}";
  TString lumi_7TeV = "5.1 fb^{-1}";

  bool outOfFrame = false;
  if (iPosX / 10 == 0) {
    outOfFrame = true;
  }
  int alignY_ = 3;
  int alignX_ = 2;
  if (iPosX / 10 == 0) alignX_ = 1;
  if (iPosX == 0) alignX_ = 1;
  if (iPosX == 0) alignY_ = 1;
  if (iPosX / 10 == 1) alignX_ = 1;
  if (iPosX / 10 == 2) alignX_ = 2;
  if (iPosX / 10 == 3) alignX_ = 3;
  if (iPosX == 0) relPosX = 0.14;
  int align_ = 10 * alignX_ + alignY_;

  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  float extraTextSize = extraOverCmsTextSize * cmsTextSize;
  float pad_ratio = (static_cast<float>(pad->GetWh()) * pad->GetAbsHNDC()) /
                    (static_cast<float>(pad->GetWw()) * pad->GetAbsWNDC());
  if (pad_ratio < 1.) pad_ratio = 1.;


  if (outOfFrame) {
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11);
    latex.SetTextSize(cmsTextSize * t * pad_ratio);
    latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText);
  }

  float posX_ = 0;
  if (iPosX % 10 <= 1) {
    posX_ = l + relPosX * (1 - l - r);
  } else if (iPosX % 10 == 2) {
    posX_ = l + 0.5 * (1 - l - r);
  } else if (iPosX % 10 == 3) {
    posX_ = 1 - r - relPosX * (1 - l - r);
  }
  float posY_ = 1 - t - relPosY * (1 - t - b);
  if (!outOfFrame) {
    latex.SetTextFont(cmsTextFont);
    latex.SetTextSize(cmsTextSize * t * pad_ratio);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, cmsText);
    if (writeExtraText) {
      latex.SetTextFont(extraTextFont);
      latex.SetTextAlign(align_);
      latex.SetTextSize(extraTextSize * t * pad_ratio);
      latex.DrawLatex(posX_, posY_ - relExtraDY * cmsTextSize * t, extraText);
    }
  } else if (writeExtraText) {
    if (iPosX == 0) {
      posX_ = l + relPosX * (1 - l - r);
      posY_ = 1 - t + lumiTextOffset * t;
    }
    latex.SetTextFont(extraTextFont);
    latex.SetTextSize(extraTextSize * t * pad_ratio);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, extraText);
  }
  pad_backup->cd();
}

void DrawCMSLogo(TPad* pad, TString cmsText, TString extraText, int iPosX) {
  DrawCMSLogo(pad, cmsText, extraText, iPosX, 0.045, 0.035, 1.2);
}

//----------------------------------------------------------------------------
// Axis adjustment and overlap checks
//----------------------------------------------------------------------------
double GetPadYMax(TPad* pad, double x_min, double x_max) {
  TList *pad_obs = pad->GetListOfPrimitives();
  if (!pad_obs) return 0.;
  TIter next(pad_obs);
  double h_max = -99999;
  TObject *obj;
  while ((obj = next())) {
    if (obj->InheritsFrom(TH1::Class())) {
      TH1 *hobj = reinterpret_cast<TH1*>(obj);
      for (int j = 1; j < hobj->GetNbinsX(); ++j) {
        if (hobj->GetBinLowEdge(j) + hobj->GetBinWidth(j) < x_min ||
            hobj->GetBinLowEdge(j) > x_max)
          continue;
        if (hobj->GetBinContent(j) + hobj->GetBinError(j) > h_max) {
          h_max = hobj->GetBinContent(j) + hobj->GetBinError(j);
        }
      }
    }
    if (obj->InheritsFrom(TGraph::Class())) {
      TGraph *gobj = reinterpret_cast<TGraph*>(obj);
      unsigned n = gobj->GetN();
      double x = 0;
      double y = 0;
      for (unsigned k = 0; k < n; ++k) {
        gobj->GetPoint(k, x, y);
        if (x < x_min || x > x_max) continue;
        if (y > h_max) h_max = y;
      }
    }
  }
  return h_max;
}

double GetPadYMax(TPad* pad) {
  TList *pad_obs = pad->GetListOfPrimitives();
  if (!pad_obs) return 0.;
  double xmin = GetAxisHist(pad)->GetXaxis()->GetXmin();
  double xmax = GetAxisHist(pad)->GetXaxis()->GetXmax();
  return GetPadYMax(pad, xmin, xmax);
}

void FixTopRange(TPad *pad, double fix_y, double fraction) {
  TH1* hobj = GetAxisHist(pad);
  double ymin = hobj->GetMinimum();
  hobj->SetMaximum((fix_y - fraction * ymin) / (1. - fraction));
  if (gPad->GetLogy()) {
    if (ymin == 0.) {
      std::cout
          << "Can't adjust log-scale y-axis range if the minimum is zero!\n";
      return;
    }
    double max =
        (std::log10(fix_y) - fraction * std::log10(ymin)) / (1 - fraction);
    max = std::pow(10, max);
    hobj->SetMaximum(max);
  }
}

void FixBoxPadding(TPad *pad, TBox *box, double frac) {
  // Get the bounds of the box - these are in the normalised
  // Pad co-ordinates.
  double p_x1 = box->GetX1();
  double p_x2 = box->GetX2();
  double p_y1 = box->GetY1();
  // double p_y2 = box->GetY2();

  // Convert to normalised co-ordinates in the frame
  double f_x1 = (p_x1 - pad->GetLeftMargin()) /
                (1. - pad->GetLeftMargin() - pad->GetRightMargin());
  double f_x2 = (p_x2 - pad->GetLeftMargin()) /
                (1. - pad->GetLeftMargin() - pad->GetRightMargin());
  double f_y1 = (p_y1 - pad->GetBottomMargin()) /
                (1. - pad->GetTopMargin() - pad->GetBottomMargin());

  // Extract histogram providing the frame and axes
  TH1 *hobj = GetAxisHist(pad);

  double xmin = hobj->GetBinLowEdge(hobj->GetXaxis()->GetFirst());
  double xmax = hobj->GetBinLowEdge(hobj->GetXaxis()->GetLast()+1);
  double ymin = hobj->GetMinimum();
  double ymax = hobj->GetMaximum();

  // Convert box bounds to x-axis values
  double a_x1 = xmin + (xmax - xmin) * f_x1;
  double a_x2 = xmin + (xmax - xmin) * f_x2;

  // Get the histogram maximum in this range, given as y-axis value
  double a_max_h = GetPadYMax(pad, a_x1, a_x2);

  // Convert this to a normalised frame value
  double f_max_h = (a_max_h - ymin) / (ymax - ymin);
  if (gPad->GetLogy()) {
    f_max_h = (log10(a_max_h) - log10(ymin)) / (log10(ymax) - log10(ymin));
  }

  if (f_y1 - f_max_h < frac) {
    double f_target = 1. - (f_y1 - frac);
    FixTopRange(pad, a_max_h, f_target);
  }
}

void FixOverlay() { gPad->GetFrame()->Draw(); gPad->RedrawAxis(); }

#endif
