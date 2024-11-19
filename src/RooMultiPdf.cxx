//
// BEGIN_HTML
// Multiple function p.d.f
// END_HTML
//

#include "../interface/RooMultiPdf.h"

#include "RooConstVar.h"

ClassImp(RooMultiPdf)


//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const char *name, const char *title, RooCategory& _x, const RooArgList& _c) : 
  RooAbsPdf(name, title),
  c("_pdfs","The list of pdfs",this),
  corr("_corrs","The list of correction factors",this),
  x("_index","the pdf index",this,_x) 
{
  int count=0;

  c.add(_c);
  for (RooAbsArg *pdf : c) {
	// This is done by the user BUT is there a way to do it at construction?
	_x.defineType(("_pdf" + std::to_string(count)).c_str(), count);
  std::unique_ptr<RooArgSet> variables(pdf->getVariables());
  std::unique_ptr<RooAbsCollection> nonConstVariables(variables->selectByAttrib("Constant", false));
	// Isn't there a better wat to hold on to these values?
	RooConstVar *tmp = new RooConstVar((std::string{"const"} + pdf->GetName()).c_str(),"",nonConstVariables->size());
	corr.addOwned(*tmp);
	count++;
  }
 _oldIndex=fIndex;
}


//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const RooMultiPdf& other, const char* name) :
 RooAbsPdf(other, name),c("_pdfs",this,other.c),corr("_corrs",this,other.corr),x("_index",this,other.x)
{
 fIndex=other.fIndex;
 _oldIndex=fIndex;
 cFactor=other.cFactor; // correction to 2*NLL by default is -> 2*0.5 per param
}

//_____________________________________________________________________________
Double_t RooMultiPdf::evaluate() const{
  double val = getCurrentPdf()->getVal(c.nset());
  _oldIndex=x;
  return val;
}

//_____________________________________________________________________________
Double_t  RooMultiPdf::getLogVal(const RooArgSet* nset) const {
  double logval = getCurrentPdf()->getLogVal(nset);
  _oldIndex=x;
  return logval;	
}
