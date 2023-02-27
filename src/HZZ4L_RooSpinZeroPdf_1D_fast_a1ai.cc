#include <HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_1D_fast_a1ai.h>

HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::HZZ4L_RooSpinZeroPdf_1D_fast_a1ai() :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(),
  a1("a1", "a1", this),
  ai1("ai1", "ai1", this)
{}

HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(
  const char *name, const char *title,
  RooAbsReal& in_a1,
  RooAbsReal& in_ai1,
  const RooArgList& inObsList,
  const RooArgList& inCoefList
) :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(name, title, inObsList, inCoefList),
  a1("a1", "a1", this, in_a1),
  ai1("ai1", "ai1", this, in_ai1)
{}

HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::HZZ4L_RooSpinZeroPdf_1D_fast_a1ai(
  const HZZ4L_RooSpinZeroPdf_1D_fast_a1ai& other, const char* name
) :
  HZZ4L_RooSpinZeroPdf_1D_fast_base(other, name),
  a1("a1", this, other.a1),
  ai1("ai1", this, other.ai1)
{}

double HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::a1Val() const {
  return a1;
}
double HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::ai1Val() const {
  return ai1;
}

bool HZZ4L_RooSpinZeroPdf_1D_fast_a1ai::isAnomalousCouplingValueValid() const {
  return true;
}

ClassImp(HZZ4L_RooSpinZeroPdf_1D_fast_a1ai);
