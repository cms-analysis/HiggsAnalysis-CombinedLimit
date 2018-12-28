#include <HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_1D_fast.h>

HZZ4L_RooSpinZeroPdf_1D_fast::HZZ4L_RooSpinZeroPdf_1D_fast() :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(),
  fai1("fai1", "fai1", this)
{}

HZZ4L_RooSpinZeroPdf_1D_fast::HZZ4L_RooSpinZeroPdf_1D_fast(
  const char *name, const char *title,
  RooAbsReal& in_fai1,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
) :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(name, title, inObsList, inCoefList),
  fai1("fai1", "fai1", this, in_fai1)
{}

HZZ4L_RooSpinZeroPdf_1D_fast::HZZ4L_RooSpinZeroPdf_1D_fast(
  const HZZ4L_RooSpinZeroPdf_1D_fast& other, const char* name
) :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(other, name),
  fai1("fai1", this, other.fai1)
{}

double HZZ4L_RooSpinZeroPdf_1D_fast::a1Val() const {
  return sqrt(1-abs(fai1));
}
double HZZ4L_RooSpinZeroPdf_1D_fast::ai1Val() const {
  return copysign(sqrt(abs(fai1)), fai1);
}

bool HZZ4L_RooSpinZeroPdf_1D_fast::isAnomalousCouplingValueValid() const {
  return abs(fai1) <= 1;
}

ClassImp(HZZ4L_RooSpinZeroPdf_1D_fast);
