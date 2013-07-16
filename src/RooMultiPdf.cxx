//
// BEGIN_HTML
// Multiple function p.d.f
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "../interface/RooMultiPdf.h"
#include "RooRealVar.h"

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
	// Isn't there a better wat to hold on to these values?
	RooConstVar *tmp = new RooConstVar(Form("const%s",fPdf->GetName()),"",fPdf->getVariables()->getSize());
	corr.add(*tmp);
	count++;
  }
  nPdfs=c.getSize();

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
	RooConstVar *tmp = new RooConstVar(Form("const%s",fPdf->GetName())
		,"",fPdf->getVariables()->getSize());
	corr.add(*tmp);
 }

 _oldIndex=fIndex;
}

bool RooMultiPdf::checkIndexDirty() const {
  
  return _oldIndex!=x;
  
}
//_____________________________________________________________________________
double RooMultiPdf::getCorrection() const {

  double val = ((RooAbsReal*)corr.at(x))->getVal(); 
  return 0.5*val;  //PVAL correction
}
//_____________________________________________________________________________
//_____________________________________________________________________________
Double_t RooMultiPdf::evaluate() const{
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  double val = cPdf->getVal();
  _oldIndex=x;
  return val;
}


