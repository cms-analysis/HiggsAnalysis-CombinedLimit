#ifndef Core_Plotting_Style_h
#define Core_Plotting_Style_h

#include "TStyle.h"
#include "TROOT.h"

/**
 * Sets the semi-official CMS plotting global style
 *
 * This is copied verbatim from the latest guidlines
 * [here](https://ghm.web.cern.ch/ghm/plots/). You should call this function
 * at the start of your macro/program. But note, the function ModTDRStyle()
 * might be more useful - this will call SetTDRStyle() first then adjust a few
 * more settings to give a better plotting style.
 */
void SetTDRStyle();

/**
 * Sets an improved plotting style, using the CMS default as a base
 *
 * Main changes are to be more consistent with margin sizes when stretching
 * the TCanvas horizontally or vertically, and to fix some style options
 * omitted from SetTDRStyle(). The canvas width and height, and the pad
 * margins are passed as parameters.
 */
void ModTDRStyle(int width, int height, double t, double b, double l, double r);

/**
 * Apply the modified plotting style using some sensible defaults for the
 * canvas size and pad margins.
 */
void ModTDRStyle();

/**
 * Modify the global style to enable/disable the drawing of axis gridlines
 */
void TdrGrid(bool grid_on);

//----------------------------------------------------------------------------
// Implementation
//----------------------------------------------------------------------------
void SetTDRStyle() {
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600);  // Height of canvas
  gStyle->SetCanvasDefW(600);  // Width of canvas
  gStyle->SetCanvasDefX(0);    // POsition on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  // For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  // gStyle->SetErrorMarker(20);
  // gStyle->SetErrorX(0.);

  gStyle->SetMarkerStyle(20);

  // For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  // For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);  // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.02);

  // For the Global title:
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // gStyle->SetTitleH(0); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // Another way to set the size?
  // gStyle->SetTitleXSize(Float_t size = 0.02);
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  // Postscript options:
  gStyle->SetPaperSize(20., 20.);
  // gStyle->SetLineScalePS(Float_t scale = 3);
  // gStyle->SetLineStyleString(Int_t i, const char* text);
  // gStyle->SetHeaderPS(const char* header);
  // gStyle->SetTitlePS(const char* pstitle);

  // gStyle->SetBarOffset(Float_t baroff = 0.5);
  // gStyle->SetBarWidth(Float_t barwidth = 0.5);
  // gStyle->SetPaintTextFormat(const char* format = "g");
  // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // gStyle->SetTimeOffset(Double_t toffset);
  // gStyle->SetHistMinimumZero(kTRUE);

  gStyle->SetHatchesLineWidth(5);
  gStyle->SetHatchesSpacing(0.05);
}

void ModTDRStyle(int width, int height, double t, double b, double l,
                 double r) {
  SetTDRStyle();

  // Set the default canvas width and height in pixels
  gStyle->SetCanvasDefW(width);
  gStyle->SetCanvasDefH(height);

  // Set the default margins.
  // These are given as fractions of the pad height for `Top` and `Bottom` and
  // the pad width for `Left` and `Right`. But we want to specify all of these
  // as fractions of the shortest length.
  double def_w = static_cast<double>(gStyle->GetCanvasDefW());
  double def_h = static_cast<double>(gStyle->GetCanvasDefH());

  double scale_h = (def_h > def_w) ? def_w / def_h : 1.;
  double scale_w = (def_w > def_h) ? def_h / def_w : 1.;

  double def_min = def_h < def_w ? def_h : def_w;

  gStyle->SetPadTopMargin(t * scale_h);      // default 0.05
  gStyle->SetPadBottomMargin(b * scale_h);   // default 0.13
  gStyle->SetPadLeftMargin(l * scale_w);     // default 0.16
  gStyle->SetPadRightMargin(r * scale_w);    // default 0.02
  // But note the new CMS style sets these:
  // 0.08, 0.12, 0.12, 0.04

  // Set number of axis tick divisions
  gStyle->SetNdivisions(506, "XYZ");            // default 510

  // Some marker properties not set in the default tdr style
  gStyle->SetMarkerColor(kBlack);
  gStyle->SetMarkerSize(1.0);


  gStyle->SetLabelOffset(0.007, "YZ");
  // This is an adhoc adjustment to scale the x-axis label
  // offset when we strect plot vertically
  // Will also need to increase if first x-axis label has more than one digit
  gStyle->SetLabelOffset(0.005 * (3. - 2. / scale_h), "X");

  // In this next part we do a slightly involved calculation to set the axis
  // title offsets, depending on the values of the TPad dimensions and margins.
  // This is to try and ensure that regardless of how these pad values are set,
  // the axis titles will be located towards the edges of the canvas and not get
  // pushed off the edge - which can often happen if a fixed value is used.
  double title_size = 0.05;
  double title_px = title_size * def_min;
  double label_size = 0.04;
  gStyle->SetTitleSize(title_size, "XYZ");
  gStyle->SetLabelSize(label_size, "XYZ");

  gStyle->SetTitleXOffset(0.5 * scale_h *
                          (1.2 * (def_h * b * scale_h - 0.6 * title_px)) /
                          title_px);
  gStyle->SetTitleYOffset(0.5 * scale_w *
                          (1.2 * (def_w * l * scale_w - 0.6 * title_px)) /
                          title_px);

  // Only draw ticks where we have an axis
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetTickLength(0.02, "XYZ");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendFillColor(0);
  gStyle->SetFillColor(0);

  gROOT->ForceStyle();
}

void ModTDRStyle() { ModTDRStyle(600, 600, 0.06, 0.12, 0.16, 0.04); }

void TdrGrid(bool grid_on) {
  gStyle->SetPadGridX(grid_on);
  gStyle->SetPadGridY(grid_on);
}

#endif
