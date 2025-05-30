//
// BEGIN_HTML
// Multiple function p.d.f
// END_HTML
//

#include <ROOT/RConfig.hxx> // for ROOT_VERSION

#if ROOT_VERSION_CODE < ROOT_VERSION(6,37,00)

#include "../interface/RooMultiPdfCombine.h"

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

#endif // ROOT_VERSION_CODE < ROOT_VERSION(6,37,00)

#include "../interface/CachingMultiPdf.h"
#include "vectorized.h"
#include "../interface/utils.h"

// Uncomment do do regression testing wrt uncached multipdf
//#define CachingMultiPdf_VALIDATE

cacheutils::CachingMultiPdf::CachingMultiPdf(const RooMultiPdf &pdf, const RooArgSet &obs) :
    pdf_(&pdf)
{
    //std::cout << "Making a CachingMultiPdf for " << pdf.GetName() << " with " <<  pdf_->getNumPdfs() << " pdfs." << std::endl;
    for (int i = 0, n = pdf_->getNumPdfs(); i < n; ++i) {
        cachingPdfs_.emplace_back(makeCachingPdf(pdf_->getPdf(i), &obs));
        //cachingPdfs_.push_back(new CachingPdf(pdf_->getPdf(i), &obs));
        //std::cout << "      MultiPdfAdding " <<  pdf.GetName() << "[" << i << "]: " << pdf_->getPdf(i)->ClassName() << " " << pdf_->getPdf(i)->GetName() << " using " << typeid(cachingPdfs_.back()).name() << std::endl;
    }
#ifdef CachingMultiPdf_VALIDATE
    cachingPdfs_.push_back(new CachingPdf((RooAbsPdf*)&pdf,&obs));
#endif

}

cacheutils::CachingMultiPdf::~CachingMultiPdf()
{
}

const std::vector<Double_t> & cacheutils::CachingMultiPdf::eval(const RooAbsData &data)
{
    const std::vector<Double_t> & ret = cachingPdfs_[pdf_->getCurrentIndex()]->eval(data);
#ifdef CachingMultiPdf_VALIDATE
    cachingPdfs_.back().setDataDirty();
    const std::vector<Double_t> & chk = cachingPdfs_.back()->eval(data);
    double normratio = 0.0;
    for (unsigned int i = 0, n = ret.size(); i < n; ++i) {
        if (chk[i] > 1e-6) normratio = max(normratio, abs(ret[i]/chk[i]-1));
    }  
    printf("%-70s: normalization match %+13.10f, selected pdf %s\n", pdf_->GetName(), normratio, pdf_->getCurrentPdf()->ClassName()); fflush(stdout);
    if (normratio > 1e-5 && typeid(*pdf_->getCurrentPdf()) == typeid(RooAddPdf) ) {
        RooAddPdf *add = (RooAddPdf*)pdf_->getCurrentPdf();
        for (unsigned int i = 0, n = add->pdfList().getSize(); i < n; ++i) {
            printf("%-70s: RooAdd[%d] is a %s\n", pdf_->GetName(), i, add->pdfList().at(i)->ClassName());
        }
    }
    if (normratio > 1e-5) {
        double sum = 0.0, sumc = 0;
        std::cout << "Before normalization " << std::endl;
        for (unsigned int i = 0, n = ret.size(); i < n; ++i) {
            printf("%s[%6d]   %18.10f  %18.10f  %+13.10f %+13.10f\n", pdf_->GetName(), i, ret[i], chk[i], ret[i]-chk[i], ret[i]/chk[i]-1);
            sum += ret[i]; sumc += chk[i];
        }
        printf("%s[%6s]   %18.10f  %18.10f  %+13.10f %+13.10f\n", pdf_->GetName(), "TOTAL", sum, sumc, sum-sumc, sum/sumc-1);
    }
    return chk;
#endif
    return ret;
}

void cacheutils::CachingMultiPdf::setDataDirty()
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setDataDirty();
    }
}

void cacheutils::CachingMultiPdf::setIncludeZeroWeights(bool includeZeroWeights) 
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setIncludeZeroWeights(includeZeroWeights);
    }
}



cacheutils::CachingAddPdf::CachingAddPdf(const RooAddPdf &pdf, const RooArgSet &obs) :
    pdf_(&pdf)
{
    const RooArgList & pdfs   = pdf.pdfList();
    const RooArgList & coeffs = pdf.coefList();
    //std::cout << "Making a CachingAddPdf for " << pdf.GetName() << " with " << pdfs.getSize() << " pdfs, " << coeffs.getSize() << " coeffs." << std::endl;
    for (int i = 0, n = pdfs.getSize(); i < n; ++i) {
        RooAbsPdf *pdfi = (RooAbsPdf*) pdfs.at(i);
        cachingPdfs_.emplace_back(makeCachingPdf(pdfi, &obs));
        //std::cout << "      AddPdfAdding " << pdf.GetName() << "[" << i << "]: " << pdfi->ClassName() << " " << pdfi->GetName() << " using " << typeid(cachingPdfs_.back()).name() << std::endl;
    }
    for (int i = 0, n = coeffs.getSize(); i < n; ++i) {
        coeffs_.push_back((RooAbsReal*)coeffs.at(i));
    }
}

cacheutils::CachingAddPdf::~CachingAddPdf()
{
}

const std::vector<Double_t> & cacheutils::CachingAddPdf::eval(const RooAbsData &data)
{
    double coefSum = 0.0;
    std::vector<double> coefVals(cachingPdfs_.size());
    for (int i = 0, n = coeffs_.size(); i < n; ++i) {
        coefVals[i] =  coeffs_[i]->getVal();
        coefSum     += coefVals[i];
    }
    if (coeffs_.size() != coefVals.size()) {
        coefVals.back() = 1.0 - coefSum; 
    } else if (coefSum != 1.0) {
        for (int i = 0, n = coeffs_.size(); i < n; ++i) {
            coefVals[i] /= coefSum;
        }
    }
    const std::vector<Double_t> & one = cachingPdfs_.front()->eval(data);
    unsigned int size = one.size();
    work_.resize(size);
    std::fill(work_.begin(), work_.end(), 0.0);
    for (int i = 0, n =  coefVals.size(); i < n; ++i) {
        vectorized::mul_add(size, coefVals[i], &(i ? cachingPdfs_[i]->eval(data) : one)[0], &work_[0]);
    }
    //std::cout << "Evaluated " << pdf_->GetName() << " of type CachingAddPdf" << std::endl;
    return work_;
}

void cacheutils::CachingAddPdf::setDataDirty()
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setDataDirty();
    }
}

void cacheutils::CachingAddPdf::setIncludeZeroWeights(bool includeZeroWeights) 
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setIncludeZeroWeights(includeZeroWeights);
    }
}

cacheutils::CachingProduct::CachingProduct(const RooProduct &pdf, const RooArgSet &obs) :
    pdf_(&pdf)
{
#if ROOT_VERSION_CODE < ROOT_VERSION(6,28,0)
    const RooArgList & pdfs = utils::factors(pdf);
#else
    const RooArgList & pdfs = pdf.realComponents();
#endif
    //std::cout << "Making a CachingProduct for " << pdf.GetName() << " with " << pdfs.getSize() << " pdfs, " << coeffs.getSize() << " coeffs." << std::endl;
    for (int i = 0, n = pdfs.getSize(); i < n; ++i) {
        RooAbsReal *pdfi = (RooAbsReal*) pdfs.at(i);
        cachingPdfs_.emplace_back(makeCachingPdf(pdfi, &obs));
    }
}

cacheutils::CachingProduct::~CachingProduct()
{
}

const std::vector<Double_t> & cacheutils::CachingProduct::eval(const RooAbsData &data)
{
    const std::vector<Double_t> & one = cachingPdfs_.front()->eval(data);
    unsigned int size = one.size();
    work_.resize(size);
    std::copy(one.begin(), one.end(), work_.begin());
    for (int i = 1, n =  cachingPdfs_.size(); i < n; ++i) {
        vectorized::mul_inplace(size, &cachingPdfs_[i]->eval(data)[0], &work_[0]);
    }
    return work_;
}

void cacheutils::CachingProduct::setDataDirty()
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setDataDirty();
    }
}

void cacheutils::CachingProduct::setIncludeZeroWeights(bool includeZeroWeights) 
{
    for (auto &pdf : cachingPdfs_) {
        pdf->setIncludeZeroWeights(includeZeroWeights);
    }
}

