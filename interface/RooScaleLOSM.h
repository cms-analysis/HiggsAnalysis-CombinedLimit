#ifndef HiggsAnalysis_CombinedLimit_RooScaleLOSM_h
#define HiggsAnalysis_CombinedLimit_RooScaleLOSM_h

#include <cmath>
#include <complex>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

typedef std::complex<double> complexD;

class RooScaleLOSM: public RooAbsReal
{
  public:
    RooScaleLOSM(){};
    RooScaleLOSM(const char *name, const char *title, RooAbsReal &mH);
    ~RooScaleLOSM() override{};

    TObject* clone(const char *newname) const override = 0;

  protected:
    complexD f(double tau) const;
    inline complexD AmpSpinOneHalf(double tau) const;
    inline complexD AmpSpinOne(double tau) const;
    Double_t evaluate() const override = 0;

    RooRealProxy mH_;
    static const double mt_, mW_;
    complexD At_, Ab_, AW_;

    double C_SM_;

  private:

    ClassDefOverride(RooScaleLOSM,1)
};

class RooScaleHGamGamLOSM: public RooScaleLOSM
{
  public:
    RooScaleHGamGamLOSM(){};
    RooScaleHGamGamLOSM(const char *name, const char *title,
    		RooAbsReal &mH, RooAbsReal &ct, RooAbsReal &cW, RooAbsReal &mb, RooAbsReal &cb);
    ~RooScaleHGamGamLOSM() override{};

    TObject* clone(const char *newname) const override;

  protected:
      Double_t evaluate() const override;
      RooRealProxy ct_, cW_, mb_, cb_;

  private:

    ClassDefOverride(RooScaleHGamGamLOSM,1)
};

class RooScaleHGluGluLOSM: public RooScaleLOSM
{
  public:
    RooScaleHGluGluLOSM(){};
    RooScaleHGluGluLOSM(const char *name, const char *title,
    		RooAbsReal &mH, RooAbsReal &ct, RooAbsReal &mb, RooAbsReal &cb);
    ~RooScaleHGluGluLOSM() override{};

    TObject* clone(const char *newname) const override;

  protected:
      Double_t evaluate() const override;
      RooRealProxy ct_, mb_, cb_;

  private:

    ClassDefOverride(RooScaleHGluGluLOSM,1)
};

class RooScaleHGamGamLOSMPlusX: public RooScaleHGamGamLOSM
{
  public:
    RooScaleHGamGamLOSMPlusX(){};
    RooScaleHGamGamLOSMPlusX(const char *name, const char *title,
    		RooAbsReal &mH, RooAbsReal &ct, RooAbsReal &cW, RooAbsReal &mb, RooAbsReal &cb, RooAbsReal &X);
    ~RooScaleHGamGamLOSMPlusX() override{};

    TObject* clone(const char *newname) const override;

  protected:
      Double_t evaluate() const override;
      RooRealProxy X_;

  private:

    ClassDefOverride(RooScaleHGamGamLOSMPlusX,1)
};

class RooScaleHGluGluLOSMPlusX: public RooScaleHGluGluLOSM
{
  public:
    RooScaleHGluGluLOSMPlusX(){};
    RooScaleHGluGluLOSMPlusX(const char *name, const char *title,
    		RooAbsReal &mH, RooAbsReal &ct, RooAbsReal &mb, RooAbsReal &cb, RooAbsReal &X);
    ~RooScaleHGluGluLOSMPlusX() override{};

    TObject* clone(const char *newname) const override;

  protected:
      Double_t evaluate() const override;
      RooRealProxy X_;

  private:

    ClassDefOverride(RooScaleHGluGluLOSMPlusX,1)
};

#endif
