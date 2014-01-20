#ifndef ROO_MULTIPDF
#define ROO_MULTIPDF

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooAbsCategory.h"
#include "RooCategory.h"
#include "RooCategoryProxy.h"
#include "RooArgProxy.h"
#include "RooAbsProxy.h"
#include "RooFormulaVar.h"
#include "RooLinkedList.h"
#include "RooConstVar.h"


#include "TIterator.h"
#include "RooListProxy.h"

#include <iostream>
#include <vector>

class RooAbsArg;
class RooAbsPdf;
class RooAbsReal;
class RooRealProxy;
class RooArgList;

using namespace std;

class RooMultiPdf : public RooAbsPdf {
public:
  enum PenatlyScheme{PVAL, AIC};
  RooMultiPdf() {} ;
  RooMultiPdf(const char *name, const char *title, RooCategory &, const RooArgList& _c);
  RooMultiPdf(const RooMultiPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooMultiPdf(*this,newname); }
  inline virtual ~RooMultiPdf() { }

/*
  RooAbsReal* createNLL(RooAbsData& data, const RooCmdArg& arg1=RooCmdArg::none(),  const RooCmdArg& arg2=RooCmdArg::none(),  
                                const RooCmdArg& arg3=RooCmdArg::none(),
const RooCmdArg& arg4=RooCmdArg::none(), const RooCmdArg&
arg5=RooCmdArg::none(),  
                                 const RooCmdArg& arg6=RooCmdArg::none(),
const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg&
arg8=RooCmdArg::none());

  RooAbsReal* createNLL(RooAbsData &data,const RooLinkedList&);
*/

//}
/*
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
*/
  bool checkIndexDirty() const;
  double getCorrection() const;
  RooAbsPdf *getCurrentPdf() const;
  int getNumPdfs(){return nPdfs;};
  void setCorrectionFactor(PenatlyScheme penal);
  Double_t getValV(const RooArgSet *nset = 0) const;
protected:
  RooListProxy c;
  RooListProxy corr;
  RooCategoryProxy x;
  //RooFormulaVar *cval;
 // RooRealProxy nllcorr;
//  RooAbsCatgeory *fIndex_r;

  int fIndex; // sigh, there should be a better way than this
  int nPdfs;
  mutable Int_t _oldIndex;

  Double_t evaluate() const;
  Double_t getLogVal(const RooArgSet *set = 0) const;
  //std::string createCorrectionString();	// should only do this once really
  double cFactor;

private:
  ClassDef(RooMultiPdf,1) // Multi PDF
};
#endif
