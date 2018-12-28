#include <HiggsAnalysis/CombinedLimit/interface/VVHZZ4L_RooSpinZeroPdf_1D_fast.h>

VVHZZ4L_RooSpinZeroPdf_1D_fast::VVHZZ4L_RooSpinZeroPdf_1D_fast() :
  VVHZZ4L_RooSpinZeroPdf_1D_fast_base(),
  fai1("fai1", "fai1", this)
{}

VVHZZ4L_RooSpinZeroPdf_1D_fast::VVHZZ4L_RooSpinZeroPdf_1D_fast(
  const char *name, const char *title,
  RooAbsReal& in_fai1,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
) :
  VVHZZ4L_RooSpinZeroPdf_1D_fast_base(name, title, inObsList, inCoefList),
  fai1("fai1", "fai1", this, in_fai1)
{}

VVHZZ4L_RooSpinZeroPdf_1D_fast::VVHZZ4L_RooSpinZeroPdf_1D_fast(
  const VVHZZ4L_RooSpinZeroPdf_1D_fast& other, const char* name
) :
  VVHZZ4L_RooSpinZeroPdf_1D_fast_base(other, name),
  fai1("fai1", this, other.fai1)
{}

double VVHZZ4L_RooSpinZeroPdf_1D_fast::a1Val() const {
  return sqrt(1-abs(fai1));
}
double VVHZZ4L_RooSpinZeroPdf_1D_fast::ai1Val() const {
  return copysign(sqrt(abs(fai1)), fai1);
}

bool VVHZZ4L_RooSpinZeroPdf_1D_fast::isAnomalousCouplingValueValid() const {
  return abs(fai1) <= 1;
}

ClassImp(VVHZZ4L_RooSpinZeroPdf_1D_fast);
