#include "../interface/CMSHggFormula.h"
#include "RooConstVar.h"

CMSHggFormulaA1::CMSHggFormulaA1(const char* name, const char* title, RooAbsReal& p0, RooAbsReal& p1,
                               RooAbsReal& p2, RooAbsReal& p3, RooArgList const& terms,
                               std::vector<double> const& coeffs)
    : RooAbsReal(name, title),
      p0_("p0", "p0", this, p0),
      p1_("p1", "p1", this, p1),
      p2_("p2", "p2", this, p2),
      p3_("p3", "p3", this, p3),
      terms_("terms", "terms", this),
      coeffs_(coeffs) {
  for (RooAbsArg* a : terms) {
    RooAbsReal* rar = dynamic_cast<RooAbsReal*>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Component ") + a->GetName() +
                                  " of CMSHggFormulaA1 is a " + a->ClassName());
    }
    terms_.add(*rar);
    vterms_.push_back(rar);
  }
  if (terms_.getSize() != int(coeffs.size())) {
    throw std::invalid_argument("Terms and coeffs not of equal size in CMSHggFormulaA1");
  }
}

CMSHggFormulaA1::CMSHggFormulaA1(const CMSHggFormulaA1& other, const char* name)
    : RooAbsReal(other, name),
      p0_("p0", this, other.p0_),
      p1_("p1", this, other.p1_),
      p2_("p2", this, other.p2_),
      p3_("p3", this, other.p3_),
      terms_("terms", this, other.terms_),
      coeffs_(other.coeffs_),
      vterms_(other.vterms_) {}

Double_t CMSHggFormulaA1::evaluate() const {
  if (vterms_.empty()) {
    vterms_.resize(terms_.getSize());
    for (int i = 0; i < terms_.getSize(); ++i) {
      vterms_[i] = dynamic_cast<RooAbsReal*>(terms_.at(i));
    }
  }
  // (@0+@1)*(1.+@2+@3+@4*@5+@6*@7+@8*@9+@10*@11+@12*@13+@14*@15+@16*@17+@18*@19+@20*@21+@22*@23)
  double ret = 1. + p2_ + p3_;
  for (unsigned i = 0; i < coeffs_.size(); ++i) {
    ret += coeffs_[i] * vterms_[i]->getVal();
  }
  return (p0_ + p1_) * ret;
}

CMSHggFormulaA2::CMSHggFormulaA2(const char* name, const char* title, RooAbsReal& p0, double const& p1,
                               RooAbsReal& p2, RooAbsReal& p3, RooArgList const& terms,
                               std::vector<double> const& coeffs)
    : RooAbsReal(name, title),
      p0_("p0", "p0", this, p0),
      p1_(p1),
      p2_("p2", "p2", this, p2),
      p3_("p3", "p3", this, p3),
      terms_("terms", "terms", this),
      coeffs_(coeffs) {
  for (RooAbsArg *a : terms) {
    RooAbsReal* rar = dynamic_cast<RooAbsReal*>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Component ") + a->GetName() +
                                  " of CMSHggFormulaA2 is a " + a->ClassName());
    }
    terms_.add(*rar);
    vterms_.push_back(rar);
  }
  if (terms_.getSize() != int(coeffs.size())) {
    throw std::invalid_argument("Terms and coeffs not of equal size in CMSHggFormulaA2");
  }
}

CMSHggFormulaA2::CMSHggFormulaA2(const CMSHggFormulaA2& other, const char* name)
    : RooAbsReal(other, name),
      p0_("p0", this, other.p0_),
      p1_(other.p1_),
      p2_("p2", this, other.p2_),
      p3_("p3", this, other.p3_),
      terms_("terms", this, other.terms_),
      coeffs_(other.coeffs_),
      vterms_(other.vterms_) {}

Double_t CMSHggFormulaA2::evaluate() const {
  if (vterms_.empty()) {
    vterms_.resize(terms_.getSize());
    for (int i = 0; i < terms_.getSize(); ++i) {
      vterms_[i] = dynamic_cast<RooAbsReal*>(terms_.at(i));
    }
  }
  // (@0+@1)*(1.+@2+@3+@4*@5+@6*@7+@8*@9+@10*@11+@12*@13+@14*@15+@16*@17+@18*@19+@20*@21+@22*@23)
  double ret = 1. + p2_ + p3_;
  for (unsigned i = 0; i < coeffs_.size(); ++i) {
    ret += coeffs_[i] * vterms_[i]->getVal();
  }
  return (p0_ + p1_) * ret;
}

CMSHggFormulaB1::CMSHggFormulaB1(const char* name, const char* title, RooAbsReal& p0,
                               RooArgList const& terms, std::vector<double> const& coeffs)
    : RooAbsReal(name, title),
      p0_("p0", "p0", this, p0),
      terms_("terms", "terms", this),
      coeffs_(coeffs) {
  for (RooAbsArg *a : terms) {
    RooAbsReal* rar = dynamic_cast<RooAbsReal*>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Component ") + a->GetName() +
                                  " of CMSHggFormulaB1 is a " + a->ClassName());
    }
    terms_.add(*rar);
    vterms_.push_back(rar);
  }
  if (terms_.getSize() != int(coeffs.size())) {
    throw std::invalid_argument("Terms and coeffs not of equal size in CMSHggFormulaB1");
  }
}

CMSHggFormulaB1::CMSHggFormulaB1(const CMSHggFormulaB1& other, const char* name)
    : RooAbsReal(other, name),
      p0_("p0", this, other.p0_),
      terms_("terms", this, other.terms_),
      coeffs_(other.coeffs_),
      vterms_(other.vterms_) {}

Double_t CMSHggFormulaB1::evaluate() const {
  if (vterms_.empty()) {
    vterms_.resize(terms_.getSize());
    for (int i = 0; i < terms_.getSize(); ++i) {
      vterms_[i] = dynamic_cast<RooAbsReal*>(terms_.at(i));
    }
  }
  // @0*TMath::Max(1.e-2,(1.+@1*@2+@3*@4+@5*@6+@7*@8+@9*@10+@11*@12+@13*@14+@15*@16+@17*@18+@19*@20+@21*@22+@23*@24+@25*@26+@27*@28+@29*@30+@31*@32+@33*@34+@35*@36+@37*@38+@39*@40+@41*@42+@43*@44))
  double ret = 1.;
  for (unsigned i = 0; i < coeffs_.size(); ++i) {
    ret += coeffs_[i] * vterms_[i]->getVal();
  }
  return p0_ * TMath::Max(1.e-2, ret);
}

CMSHggFormulaB2::CMSHggFormulaB2(const char* name, const char* title, double const& p0,
                               RooArgList const& terms, std::vector<double> const& coeffs)
    : RooAbsReal(name, title),
      p0_(p0),
      terms_("terms", "terms", this),
      coeffs_(coeffs) {
  for (RooAbsArg* a : terms) {
    RooAbsReal* rar = dynamic_cast<RooAbsReal*>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Component ") + a->GetName() +
                                  " of CMSHggFormulaB2 is a " + a->ClassName());
    }
    terms_.add(*rar);
    vterms_.push_back(rar);
  }
  if (terms_.getSize() != int(coeffs.size())) {
    throw std::invalid_argument("Terms and coeffs not of equal size in CMSHggFormulaB2");
  }
}

CMSHggFormulaB2::CMSHggFormulaB2(const CMSHggFormulaB2& other, const char* name)
    : RooAbsReal(other, name),
      p0_(other.p0_),
      terms_("terms", this, other.terms_),
      coeffs_(other.coeffs_),
      vterms_(other.vterms_) {}

Double_t CMSHggFormulaB2::evaluate() const {
  if (vterms_.empty()) {
    vterms_.resize(terms_.getSize());
    for (int i = 0; i < terms_.getSize(); ++i) {
      vterms_[i] = dynamic_cast<RooAbsReal*>(terms_.at(i));
    }
  }
  // @0*TMath::Max(1.e-2,(1.+@1*@2+@3*@4+@5*@6+@7*@8+@9*@10+@11*@12+@13*@14+@15*@16+@17*@18+@19*@20+@21*@22+@23*@24+@25*@26+@27*@28+@29*@30+@31*@32+@33*@34+@35*@36+@37*@38+@39*@40+@41*@42+@43*@44))
  double ret = 1.;
  for (unsigned i = 0; i < coeffs_.size(); ++i) {
    ret += coeffs_[i] * vterms_[i]->getVal();
  }
  return p0_ * TMath::Max(1.e-2, ret);
}

CMSHggFormulaC1::CMSHggFormulaC1(const char* name, const char* title, RooArgList const& terms,
                               std::vector<double> const& coeffs)
    : RooAbsReal(name, title),
      terms_("terms", "terms", this),
      coeffs_(coeffs) {
  for (RooAbsArg* a : terms) {
    RooAbsReal* rar = dynamic_cast<RooAbsReal*>(a);
    if (!rar) {
      throw std::invalid_argument(std::string("Component ") + a->GetName() +
                                  " of CMSHggFormulaC1 is a " + a->ClassName());
    }
    terms_.add(*rar);
    vterms_.push_back(rar);
  }
  if (terms_.getSize() != int(coeffs.size())) {
    throw std::invalid_argument("Terms and coeffs not of equal size in CMSHggFormulaC1");
  }
}

CMSHggFormulaC1::CMSHggFormulaC1(const CMSHggFormulaC1& other, const char* name)
    : RooAbsReal(other, name),
      terms_("terms", this, other.terms_),
      coeffs_(other.coeffs_),
      vterms_(other.vterms_) {}

Double_t CMSHggFormulaC1::evaluate() const {
  if (vterms_.empty()) {
    vterms_.resize(terms_.getSize());
    for (int i = 0; i < terms_.getSize(); ++i) {
      vterms_[i] = dynamic_cast<RooAbsReal*>(terms_.at(i));
    }
  }
  // (1.+@0*@1+@2*@3+@4*@5+@6*@7+@8*@9+@10*@11+@12*@13) @0[C],@1[V],@2[C],@3[V],@4[C],@5[V],@6[C],@7[V],@8[C],@9[V],@10[C],@11[V],@12[C],@13[V]
  double ret = 1.;
  for (unsigned i = 0; i < coeffs_.size(); ++i) {
    ret += coeffs_[i] * vterms_[i]->getVal();
  }
  return ret;
}

CMSHggFormulaD1::CMSHggFormulaD1(const char* name, const char* title, RooAbsReal& p0, RooAbsReal& p1)
    : RooAbsReal(name, title),
      p0_("p0", "p0", this, p0),
      p1_("p1", "p1", this, p1) { }

CMSHggFormulaD1::CMSHggFormulaD1(const CMSHggFormulaD1& other, const char* name)
    : RooAbsReal(other, name),
      p0_("p0", this, other.p0_),
      p1_("p1", this, other.p1_) {}

Double_t CMSHggFormulaD1::evaluate() const {
  // TMath::Min(@0+@1,1.0) @0[V],@1[V]
  return TMath::Min(p0_+p1_,1.0);
}

CMSHggFormulaD2::CMSHggFormulaD2(const char* name, const char* title, RooAbsReal& p0, double const& p1)
    : RooAbsReal(name, title),
      p0_("p0", "p0", this, p0),
      p1_(p1) { }

CMSHggFormulaD2::CMSHggFormulaD2(const CMSHggFormulaD2& other, const char* name)
    : RooAbsReal(other, name),
      p0_("p0", this, other.p0_),
      p1_(other.p1_) {}

Double_t CMSHggFormulaD2::evaluate() const {
  // TMath::Min(@0+@1,1.0) @0[V],@1[V]
  return TMath::Min(p0_+p1_,1.0);
}
