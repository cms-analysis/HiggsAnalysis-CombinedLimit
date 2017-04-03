#ifndef HiggsAnalysis_CombinedLimit_FastTplFunc
#define HiggsAnalysis_CombinedLimit_FastTplFunc

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "FastTemplate.h"
#include <iostream>
#include <cstring>
#include <cassert>

template <typename T> class FastTemplateFunc_t : public RooAbsReal{
protected:
  RooListProxy obsList;

public:
  FastTemplateFunc_t() : RooAbsReal(), obsList("obsList", "obsList", this){}
  FastTemplateFunc_t(const char *name, const char *title, RooArgList& inObsList) : RooAbsReal(name, title), obsList("obsList", "obsList", this){ setProxyList(obsList, inObsList); }
  FastTemplateFunc_t(const FastTemplateFunc_t& other, const char* name=0) : RooAbsReal(other, name), obsList("obsList", this, other.obsList){}
  virtual inline ~FastTemplateFunc_t(){}

  void setProxyList(RooListProxy& proxyList, RooArgList& varList){
    TIterator* varIter = varList.createIterator();
    RooAbsArg* var;
    while ((var = (RooAbsArg*)varIter->Next())) {
      if (!dynamic_cast<RooAbsReal*>(var)) {
        assert(0);
      }
      proxyList.add(*var);
    }
    delete varIter;
  }

  virtual TObject* clone(const char* newname) const = 0;
  virtual Double_t evaluate() const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const = 0;

private:
  ClassDef(FastTemplateFunc_t, 1)

};
typedef FastTemplateFunc_t<Float_t> FastTemplateFunc_f;
typedef FastTemplateFunc_t<Double_t> FastTemplateFunc_d;

template <typename T, typename U=Double_t> class FastHistoFunc_t : public FastTemplateFunc_t<T>{
protected:
  FastHisto_t<T,U> tpl;
  T fullIntegral;

public:
  FastHistoFunc_t() : FastTemplateFunc_t<T>(), fullIntegral(0){}
  FastHistoFunc_t(const char *name, const char *title, RooArgList& inObsList, FastHisto_t<T,U>& inTpl) : FastTemplateFunc_t<T>(name, title, inObsList), tpl(inTpl), fullIntegral(tpl.IntegralWidth()){
    if ((this->obsList).getSize()!=1){
      std::cerr << "FastHistoFunc_t::FastHistoFunc_t(" << this->GetName() << "): obsList.size()!=1!" << std::endl;
      assert(0);
    }
  }
  FastHistoFunc_t(const FastHistoFunc_t& other, const char* name=0) : FastTemplateFunc_t<T>(other, name), tpl(other.tpl), fullIntegral(other.fullIntegral){}
  ~FastHistoFunc_t(){}
  TObject* clone(const char* newname) const { return new FastHistoFunc_t(*this, newname); }

  FastHisto_t<T,U> getHistogram() const{ return tpl; }
  const Int_t getFullIntegralCode() const{ return 2; }

  Double_t evaluate() const{
    T x = (T)(dynamic_cast<RooAbsReal*>((this->obsList).at(0))->getVal());
    Double_t value=tpl.GetAt(x);
    return value;
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const{
    Int_t code=1;
    const Int_t code_prime[1]={ 2 };
    for (int ic=0; ic<(this->obsList).getSize(); ic++){
      if (ic>=1){ code=0; break; }
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(ic));
      if (var!=0){
        if (this->matchArgs(allVars, analVars, RooArgSet(*var))) code*=code_prime[ic];
      }
    }
    if (code==1) code=0;
    return code;
  }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const{
    if (code==0) return evaluate();
    U xmin = (U)(dynamic_cast<RooRealVar*>((this->obsList).at(0))->getMin(rangeName));
    U xmax = (U)(dynamic_cast<RooRealVar*>((this->obsList).at(0))->getMax(rangeName));
    if (
      std::abs((Float_t)(xmin/tpl.GetXmin()-U(1)))<1e-5
      &&
      std::abs((Float_t)(xmax/tpl.GetXmax()-U(1)))<1e-5
      ) return (Double_t)fullIntegral;
    else{
      int binmin = tpl.FindBin(xmin);
      int binmax = tpl.FindBin(xmax);
      return (Double_t)tpl.IntegralWidth(binmin, binmax);
    }
  }

private:
  ClassDef(FastHistoFunc_t, 1)

};
typedef FastHistoFunc_t<Float_t> FastHistoFunc_f;
typedef FastHistoFunc_t<Double_t> FastHistoFunc_d;

template <typename T, typename U=Double_t> class FastHisto2DFunc_t : public FastTemplateFunc_t<T>{
protected:
  FastHisto2D_t<T,U> tpl;
  T fullIntegral;

public:
  FastHisto2DFunc_t() : FastTemplateFunc_t<T>(), fullIntegral(0){}
  FastHisto2DFunc_t(const char *name, const char *title, RooArgList& inObsList, FastHisto2D_t<T,U>& inTpl) : FastTemplateFunc_t<T>(name, title, inObsList), tpl(inTpl), fullIntegral(tpl.IntegralWidth()){
    if ((this->obsList).getSize()!=2){
      std::cerr << "FastHisto2DFunc_t::FastHisto2DFunc_t(" << this->GetName() << "): obsList.size()!=2!" << std::endl;
      assert(0);
    }
  }
  FastHisto2DFunc_t(const FastHisto2DFunc_t& other, const char* name=0) : FastTemplateFunc_t<T>(other, name), tpl(other.tpl), fullIntegral(other.fullIntegral){}
  ~FastHisto2DFunc_t(){}
  TObject* clone(const char* newname) const { return new FastHisto2DFunc_t(*this, newname); }

  FastHisto2D_t<T,U> getHistogram() const{ return tpl; }
  const Int_t getFullIntegralCode() const{ return /*2*3*/6; }

  Double_t evaluate() const{
    U x = (U)(dynamic_cast<RooAbsReal*>((this->obsList).at(0))->getVal());
    U y = (U)(dynamic_cast<RooAbsReal*>((this->obsList).at(1))->getVal());
    Double_t value=tpl.GetAt(x, y);
    return value;
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const{
    Int_t code=1;
    const Int_t code_prime[2]={ 2, 3 };
    for (int ic=0; ic<(this->obsList).getSize(); ic++){
      if (ic>=2){ code=0; break; }
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(ic));
      if (var!=0){
        if (this->matchArgs(allVars, analVars, RooArgSet(*var))) code*=code_prime[ic];
      }
    }
    if (code==1) code=0;
    return code;
  }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const{
    if (code==0) return evaluate();
    U xmin, xmax, ymin, ymax;
    if (code%2==0){
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(0));
      xmin = (U)(var->getMin(rangeName));
      xmax = (U)(var->getMax(rangeName));
    }
    else{
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(0));
      xmin = (U)(var->getVal());
      xmax = xmin;
    }
    if (code%3==0){
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(1));
      ymin = (U)(var->getMin(rangeName));
      ymax = (U)(var->getMax(rangeName));
    }
    else{
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(1));
      ymin = (U)(var->getVal());
      ymax = ymin;
    }
    if (
      code%6==0
      &&
      fabs((Float_t)(xmin/tpl.GetXmin()-U(1)))<1e-5
      &&
      fabs((Float_t)(xmax/tpl.GetXmax()-U(1)))<1e-5
      &&
      fabs((Float_t)(ymin/tpl.GetYmin()-U(1)))<1e-5
      &&
      fabs((Float_t)(ymax/tpl.GetYmax()-U(1)))<1e-5
      ) return (Double_t)fullIntegral;
    else{
      int xbinmin = tpl.FindBinX(xmin);
      int xbinmax = tpl.FindBinX(xmax);
      int ybinmin = tpl.FindBinY(ymin);
      int ybinmax = tpl.FindBinY(ymax);
      T result = tpl.IntegralWidth(xbinmin, xbinmax, ybinmin, ybinmax);
      if (code%2!=0) result /= tpl.GetBinWidthX(xbinmin);
      if (code%3!=0) result /= tpl.GetBinWidthY(ybinmin);
      return (Double_t)result;
    }
  }

private:
  ClassDef(FastHisto2DFunc_t, 1)

};
typedef FastHisto2DFunc_t<Float_t> FastHisto2DFunc_f;
typedef FastHisto2DFunc_t<Double_t> FastHisto2DFunc_d;

template <typename T, typename U=Double_t> class FastHisto3DFunc_t : public FastTemplateFunc_t<T>{
protected:
  FastHisto3D_t<T,U> tpl;
  T fullIntegral;

public:
  FastHisto3DFunc_t() : FastTemplateFunc_t<T>(), fullIntegral(0){}
  FastHisto3DFunc_t(const char *name, const char *title, RooArgList& inObsList, FastHisto3D_t<T,U>& inTpl) : FastTemplateFunc_t<T>(name, title, inObsList), tpl(inTpl), fullIntegral(tpl.IntegralWidth()){
    if ((this->obsList).getSize()!=3){
      std::cerr << "FastHisto3DFunc_t::FastHisto3DFunc_t(" << this->GetName() << "): obsList.size()!=3!" << std::endl;
      assert(0);
    }
  }
  FastHisto3DFunc_t(const FastHisto3DFunc_t& other, const char* name=0) : FastTemplateFunc_t<T>(other, name), tpl(other.tpl), fullIntegral(other.fullIntegral){}
  ~FastHisto3DFunc_t(){}
  TObject* clone(const char* newname) const { return new FastHisto3DFunc_t(*this, newname); }

  FastHisto3D_t<T,U> getHistogram() const{ return tpl; }
  const Int_t getFullIntegralCode() const{ return /*2*3*5*/30; }

  Double_t evaluate() const{
    U x = (U)(dynamic_cast<RooAbsReal*>((this->obsList).at(0))->getVal());
    U y = (U)(dynamic_cast<RooAbsReal*>((this->obsList).at(1))->getVal());
    U z = (U)(dynamic_cast<RooAbsReal*>((this->obsList).at(2))->getVal());
    Double_t value=tpl.GetAt(x, y, z);
    return value;
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const{
    Int_t code=1;
    const Int_t code_prime[3]={ 2, 3, 5 };
    for (int ic=0; ic<(this->obsList).getSize(); ic++){
      if (ic>=3){ code=0; break; }
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(ic));
      if (var!=0){
        if (this->matchArgs(allVars, analVars, RooArgSet(*var))) code*=code_prime[ic];
      }
    }
    if (code==1) code=0;
    return code;
  }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const{
    if (code==0) return evaluate();
    U xmin, xmax, ymin, ymax, zmin, zmax;
    if (code%2==0){
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(0));
      xmin = (U)(var->getMin(rangeName));
      xmax = (U)(var->getMax(rangeName));
    }
    else{
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(0));
      xmin = (U)(var->getVal());
      xmax = xmin;
    }
    if (code%3==0){
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(1));
      ymin = (U)(var->getMin(rangeName));
      ymax = (U)(var->getMax(rangeName));
    }
    else{
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(1));
      ymin = (U)(var->getVal());
      ymax = ymin;
    }
    if (code%5==0){
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(2));
      zmin = (U)(var->getMin(rangeName));
      zmax = (U)(var->getMax(rangeName));
    }
    else{
      RooRealVar* var = dynamic_cast<RooRealVar*>((this->obsList).at(2));
      zmin = (U)(var->getVal());
      zmax = zmin;
    }
    if (
      code%30==0
      &&
      fabs((Float_t)(xmin/tpl.GetXmin()-U(1)))<1e-5
      &&
      fabs((Float_t)(xmax/tpl.GetXmax()-U(1)))<1e-5
      &&
      fabs((Float_t)(ymin/tpl.GetYmin()-U(1)))<1e-5
      &&
      fabs((Float_t)(ymax/tpl.GetYmax()-U(1)))<1e-5
      &&
      fabs((Float_t)(zmin/tpl.GetZmin()-U(1)))<1e-5
      &&
      fabs((Float_t)(zmax/tpl.GetZmax()-U(1)))<1e-5
      ) return (Double_t)fullIntegral;
    else{
      int xbinmin = tpl.FindBinX(xmin);
      int xbinmax = tpl.FindBinX(xmax);
      int ybinmin = tpl.FindBinY(ymin);
      int ybinmax = tpl.FindBinY(ymax);
      int zbinmin = tpl.FindBinZ(zmin);
      int zbinmax = tpl.FindBinZ(zmax);
      T result = tpl.IntegralWidth(xbinmin, xbinmax, ybinmin, ybinmax, zbinmin, zbinmax);
      if (code%2!=0) result /= tpl.GetBinWidthX(xbinmin);
      if (code%3!=0) result /= tpl.GetBinWidthY(ybinmin);
      if (code%5!=0) result /= tpl.GetBinWidthZ(zbinmin);
      return (Double_t)result;
    }
  }

private:
  ClassDef(FastHisto3DFunc_t, 1)

};
typedef FastHisto3DFunc_t<Float_t> FastHisto3DFunc_f;
typedef FastHisto3DFunc_t<Double_t> FastHisto3DFunc_d;


#endif
