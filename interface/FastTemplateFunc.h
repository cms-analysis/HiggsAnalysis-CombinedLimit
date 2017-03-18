#ifndef HiggsAnalysis_CombinedLimit_FastTplFunc
#define HiggsAnalysis_CombinedLimit_FastTplFunc

#include "RooAbsReal.h"
#include "FastTemplate.h"

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

public:
  FastHistoFunc_t() : FastTemplateFunc_t<T>(){}
  FastHistoFunc_t(const char *name, const char *title, RooArgList& inObsList, FastHisto_t<T>& inTpl) : FastTemplateFunc_t<T>(name, title, inObsList), tpl(inTpl){}
  FastHistoFunc_t(const FastHistoFunc_t& other, const char* name=0) : FastTemplateFunc_t<T>(other, name), tpl(other.tpl){}
  ~FastHistoFunc_t(){}
  TObject* clone(const char* newname) const { return new FastHistoFunc_t(*this, newname); }

  Double_t evaluate() const{
    Double_t value=0;
    return value;
  }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const{
    Int_t code=1;
    return code;
  }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const{
    Double_t value=0;
    return value;
  }

private:
  ClassDef(FastHistoFunc_t, 1)

};

typedef FastTemplateFunc_t<Float_t> FastTemplateFunc_f;
typedef FastTemplateFunc_t<Double_t> FastTemplateFunc_d;

#endif
