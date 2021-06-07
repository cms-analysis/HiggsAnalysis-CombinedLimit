//
// BEGIN_HTML
// Multiple function p.d.f
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include <stdexcept>

ClassImp(RooMultiPdf)


//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const char *name, const char *title, RooCategory& _x, const RooArgList& _c) : 
  RooAbsPdf(name, title),  //Why is this here? just to use the names to be used? 
  c("_pdfs","The list of pdfs",this),
  x("_index","the pdf index",this,_x) 
{
  TIterator *pdfIter=_c.createIterator(); 
  int count=0;

  RooAbsPdf *fPdf;
  while ( (fPdf = (RooAbsPdf*) pdfIter->Next()) ){
	c.add(*fPdf);
	// This is done by the user BUT is there a way to do it at construction?
	_x.defineType(Form("_pdf%d",count),count);//(fPdf->getParameters())->getSize());
  std::unique_ptr<RooArgSet> variables(fPdf->getVariables());
  std::unique_ptr<RooAbsCollection> nonConstVariables(variables->selectByAttrib("Constant", false));
	// Isn't there a better wat to hold on to these values?
	RooConstVar *tmp = new RooConstVar(Form("const%s",fPdf->GetName()),"",nonConstVariables->getSize());
	corr.add(*tmp);
	count++;
  }
  nPdfs=c.getSize();
  cFactor=0.5; // correction to 2*NLL by default is -> 2*0.5 per param
 _oldIndex=fIndex;
 
}


//_____________________________________________________________________________
RooMultiPdf::RooMultiPdf(const RooMultiPdf& other, const char* name) :
 RooAbsPdf(other, name),c("_pdfs",this,RooListProxy()),x("_index",this,other.x)
{

 fIndex=other.fIndex;
 nPdfs=other.nPdfs;

 TIterator *pdfIter=(other.c).createIterator();

 RooAbsPdf *fPdf;
 while ( (fPdf = (RooAbsPdf*) pdfIter->Next()) ){
	c.add(*fPdf);
  std::unique_ptr<RooArgSet> variables(fPdf->getVariables());
  std::unique_ptr<RooAbsCollection> nonConstVariables(variables->selectByAttrib("Constant", false));

	RooConstVar *tmp = new RooConstVar(Form("const%s",fPdf->GetName())
		,"",nonConstVariables->getSize());
	corr.add(*tmp);
 }

 _oldIndex=fIndex;
  cFactor=other.cFactor; // correction to 2*NLL by default is -> 2*0.5 per param
}

bool RooMultiPdf::checkIndexDirty() const {
  return _oldIndex!=x;  
}
//_____________________________________________________________________________
void RooMultiPdf::setCorrectionFactor(PenatlyScheme penal){
  if ( penal==AIC ){
  	cFactor=1.0;
  } else if ( penal==PVAL ){
	cFactor=0.5;
  }
}
//_____________________________________________________________________________
void RooMultiPdf::setCorrectionFactor(double penal){
  cFactor=penal;
}
//_____________________________________________________________________________
double RooMultiPdf::getCorrection() const {

  double val = ((RooAbsReal*)corr.at(x))->getVal(); 
  return cFactor*val;  //PVAL correction
}
//_____________________________________________________________________________
RooAbsPdf* RooMultiPdf::getCurrentPdf() const {

  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  return cPdf; 
}
RooAbsPdf* RooMultiPdf::getPdf(int index) const {
    
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(index));
  return cPdf;
}

int RooMultiPdf::getCurrentIndex() const {
    Int_t index = x;
    return index;
}

//_____________________________________________________________________________
Double_t RooMultiPdf::getValV(const RooArgSet* nset) const {
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  double val = cPdf->getVal(nset);
  _oldIndex=x;
  return val;
}

//_____________________________________________________________________________
Double_t RooMultiPdf::evaluate() const{
  // This is dangerous since if the underlying pdf is a RooAddPdf the meaning of the 
  // coefficients depends on the normalization set, and we don't really know
  // how this information is propagated.
  // So, we just forward the getVal which is anyway the contract for RooMultiPdf.
  throw std::invalid_argument("RooMultiPdf::evaluate() called\n");
}

//_____________________________________________________________________________
Double_t  RooMultiPdf::getLogVal(const RooArgSet* nset) const {
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  double logval = cPdf->getLogVal(nset);
  _oldIndex=x;
  return logval;	
}

