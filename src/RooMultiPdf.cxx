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
#include "RooAddPdf.h"

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
	RooConstVar *tmp = new RooConstVar(Form("const%s",fPdf->GetName())
		,"",fPdf->getVariables()->getSize());
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
double RooMultiPdf::getCorrection() const {

  double val = ((RooAbsReal*)corr.at(x))->getVal(); 
  return cFactor*val;  //PVAL correction
}
//_____________________________________________________________________________
RooAbsPdf* RooMultiPdf::getCurrentPdf() const {

  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  return cPdf; 
}
//_____________________________________________________________________________
Double_t RooMultiPdf::evaluate() const{

  double val=0;
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  if (cPdf->IsA()->InheritsFrom(RooAddPdf::Class()))
  {
	const RooAddPdf *aPdf = dynamic_cast<const RooAddPdf*>(cPdf);
	val = aPdf->evaluate();
  } else {
   	val = cPdf->getVal();
  }
//   val = cPdf->getVal();
  _oldIndex=x;
  return val;

}

//_____________________________________________________________________________
Double_t  RooMultiPdf::getLogVal(const RooArgSet* nset) const{

  double logval=0;
  RooAbsPdf *cPdf = ((RooAbsPdf*)c.at(x)); 
  if (cPdf->IsA()->InheritsFrom(RooAddPdf::Class()))
  {
	const RooAddPdf *aPdf = dynamic_cast<const RooAddPdf*>(cPdf);
	logval = aPdf->getLogVal(nset);
  } else {
	logval = cPdf->getLogVal(nset);
  }
  _oldIndex=x;
  return logval;//cPdf->getLogVal(nset);	
}

