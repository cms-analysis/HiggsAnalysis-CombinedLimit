#include "CombineHarvester/CombineTools/interface/SOverBTools.h"
#include "TH1.h"
#include "TAxis.h"

namespace ch {

SOverBInfo::SOverBInfo(TH1F const* sig, TH1F const* bkg, unsigned steps,
                       double frac) {
  double xmin = sig->GetXaxis()->GetXmin();
  double xmax = sig->GetXaxis()->GetXmax();
  double step_size = (xmax - xmin) / static_cast<double>(steps);
  double sig_tot = sig->Integral();
  double lower_limit = 0;
  double upper_limit = 0;
  double ofrac = (1. - frac) / 2.;
  for (unsigned j = 0; j < steps; ++j) {
    double integral = ch::IntegrateFloatRange(
        sig, xmin, xmin + (step_size * static_cast<double>(j)));
    if (integral / sig_tot > ofrac) {
      lower_limit = xmin + (step_size * static_cast<double>(j));
      break;
    }
  }
  for (unsigned j = 0; j < steps; ++j) {
    double integral = ch::IntegrateFloatRange(
        sig, xmax - (step_size * static_cast<double>(j)), xmax);
    if (integral / sig_tot > ofrac) {
      upper_limit = xmax - (step_size * static_cast<double>(j));
      break;
    }
  }
  x_lo = lower_limit;
  x_hi = upper_limit;
  s = ch::IntegrateFloatRange(sig, lower_limit, upper_limit);
  b = ch::IntegrateFloatRange(bkg, lower_limit, upper_limit);
}

double IntegrateFloatRange(TH1F const* hist, double xmin, double xmax) {
  TAxis const* axis = hist->GetXaxis();
  int bmin = axis->FindFixBin(xmin);
  int bmax = axis->FindFixBin(xmax);
  double integral = hist->Integral(bmin, bmax);
  integral -= hist->GetBinContent(bmin) * (xmin - axis->GetBinLowEdge(bmin)) /
              axis->GetBinWidth(bmin);
  integral -= hist->GetBinContent(bmax) * (axis->GetBinUpEdge(bmax) - xmax) /
              axis->GetBinWidth(bmax);
  return integral;
}
}
