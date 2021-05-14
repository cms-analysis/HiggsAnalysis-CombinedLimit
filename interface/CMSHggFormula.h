#ifndef HiggsAnalysis_CombinedLimit_CMSHggFormula_h
#define HiggsAnalysis_CombinedLimit_CMSHggFormula_h
#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include <vector>

class CMSHggFormulaA1 : public RooAbsReal {
 public:
  CMSHggFormulaA1() {}
  CMSHggFormulaA1(const char* name, const char* title, RooAbsReal & p0, RooAbsReal & p1,
                 RooAbsReal & p2, RooAbsReal & p3, RooArgList const& terms,
                 std::vector<double> const& coeffs);
  CMSHggFormulaA1(const CMSHggFormulaA1& other, const char* name = 0);
  virtual ~CMSHggFormulaA1() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaA1(*this, newname); }

 protected:
  RooRealProxy p0_;
  RooRealProxy p1_;
  RooRealProxy p2_;
  RooRealProxy p3_;
  RooListProxy terms_;
  std::vector<double> coeffs_;
  mutable std::vector<RooAbsReal*> vterms_; //! not to be serialized
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaA1,1)

};

class CMSHggFormulaA2 : public RooAbsReal {
 public:
  CMSHggFormulaA2() {}
  CMSHggFormulaA2(const char* name, const char* title, RooAbsReal & p0, double const& p1,
                 RooAbsReal & p2, RooAbsReal & p3, RooArgList const& terms,
                 std::vector<double> const& coeffs);
  CMSHggFormulaA2(const CMSHggFormulaA2& other, const char* name = 0);
  virtual ~CMSHggFormulaA2() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaA2(*this, newname); }

 protected:
  RooRealProxy p0_;
  double p1_;
  RooRealProxy p2_;
  RooRealProxy p3_;
  RooListProxy terms_;
  std::vector<double> coeffs_;
  mutable std::vector<RooAbsReal*> vterms_; //! not to be serialized
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaA2,1)

};

class CMSHggFormulaB1 : public RooAbsReal {
 public:
  CMSHggFormulaB1() {}
  CMSHggFormulaB1(const char* name, const char* title, RooAbsReal & p0, RooArgList const& terms,
                 std::vector<double> const& coeffs);
  CMSHggFormulaB1(const CMSHggFormulaB1& other, const char* name = 0);
  virtual ~CMSHggFormulaB1() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaB1(*this, newname); }

 protected:
  RooRealProxy p0_;
  RooListProxy terms_;
  std::vector<double> coeffs_;
  mutable std::vector<RooAbsReal*> vterms_; //! not to be serialized
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaB1,1)

};

class CMSHggFormulaB2 : public RooAbsReal {
 public:
  CMSHggFormulaB2() {}
  CMSHggFormulaB2(const char* name, const char* title, double const& p0, RooArgList const& terms,
                 std::vector<double> const& coeffs);
  CMSHggFormulaB2(const CMSHggFormulaB2& other, const char* name = 0);
  virtual ~CMSHggFormulaB2() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaB2(*this, newname); }

 protected:
  double p0_;
  RooListProxy terms_;
  std::vector<double> coeffs_;
  mutable std::vector<RooAbsReal*> vterms_; //! not to be serialized
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaB2,1)

};

class CMSHggFormulaC1 : public RooAbsReal {
 public:
  CMSHggFormulaC1() {}
  CMSHggFormulaC1(const char* name, const char* title, RooArgList const& terms,
                 std::vector<double> const& coeffs);
  CMSHggFormulaC1(const CMSHggFormulaC1& other, const char* name = 0);
  virtual ~CMSHggFormulaC1() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaC1(*this, newname); }

 protected:
  RooListProxy terms_;
  std::vector<double> coeffs_;
  mutable std::vector<RooAbsReal*> vterms_; //! not to be serialized
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaC1,1)

};

class CMSHggFormulaD1 : public RooAbsReal {
 public:
  CMSHggFormulaD1() {}
  CMSHggFormulaD1(const char* name, const char* title, RooAbsReal & p0, RooAbsReal & p1);
  CMSHggFormulaD1(const CMSHggFormulaD1& other, const char* name = 0);
  virtual ~CMSHggFormulaD1() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaD1(*this, newname); }

 protected:
  RooRealProxy p0_;
  RooRealProxy p1_;
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaD1,1)

};

class CMSHggFormulaD2 : public RooAbsReal {
 public:
  CMSHggFormulaD2() {}
  CMSHggFormulaD2(const char* name, const char* title, RooAbsReal & p0, double const& p1);
  CMSHggFormulaD2(const CMSHggFormulaD2& other, const char* name = 0);
  virtual ~CMSHggFormulaD2() {}
  virtual TObject* clone(const char* newname) const { return new CMSHggFormulaD2(*this, newname); }

 protected:
  RooRealProxy p0_;
  double p1_;
  virtual Double_t evaluate() const;

 private:
  ClassDef(CMSHggFormulaD2,1)

};

#endif
