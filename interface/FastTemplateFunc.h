#ifndef HiggsAnalysis_CombinedLimit_FastTplFunc
#define HiggsAnalysis_CombinedLimit_FastTplFunc

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "FastTemplate.h"
#include <cstring>

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
template <typename T> class FastHistoFunc_t : public FastTemplateFunc_t<T>{
protected:
  FastHisto_t<T> tpl;
  T fullIntegral;

public:
  FastHistoFunc_t() : FastTemplateFunc_t<T>(){}
  FastHistoFunc_t(const char *name, const char *title, RooArgList& inObsList, FastHisto_t<T>& inTpl) : FastTemplateFunc_t<T>(name, title, inObsList), tpl(inTpl), fullIntegral(tpl.IntegralWidth()){}
  FastHistoFunc_t(const FastHistoFunc_t& other, const char* name=0) : FastTemplateFunc_t<T>(other, name), tpl(other.tpl), fullIntegral(other.fullIntegral){}
  ~FastHistoFunc_t(){}
  TObject* clone(const char* newname) const { return new FastHistoFunc_t(*this, newname); }

  Double_t evaluate() const{
    T x = (T)((this->obsList).at(0)->getVal());
    Double_t value=tpl.GetAt(x);
    return value;
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const{
    Int_t code=1;
    const Int_t code_prime[1]={ 2 };
    for (int ic=0; ic<(this->obsList).getSize(); ic++){
      if (ic>=1){ code=0; break; }
      if (dynamic_cast<RooRealVar*>((this->obsList).at(ic))!=0){
        RooRealVar* var = (this->obsList).at(ic);
        if (this->matchArgs(allVars, analVars, RooArgSet(*var))) code*=code_prime[ic];
      }
    }
    if (code==1) code=0;
    return code;
  }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const{
    if (code==0) return evaluate();
    T xmin = (T)((this->obsList).at(0)->getMin(rangeName));
    T xmax = (T)((this->obsList).at(0)->getMax(rangeName));
    if (
      fabs((Float_t)(xmin/tpl.GetXmin()-1.))<1e-5
      &&
      fabs((Float_t)(xmax/tpl.GetXmax()-1.))<1e-5
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

typedef FastTemplateFunc_t<Float_t> FastTemplateFunc_f;
typedef FastTemplateFunc_t<Double_t> FastTemplateFunc_d;

#endif
