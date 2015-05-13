#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include <stdexcept>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooProduct.h>
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"
#include <HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h>
#include <HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h>
#include <HiggsAnalysis/CombinedLimit/interface/VectorizedGaussian.h>
#include <HiggsAnalysis/CombinedLimit/interface/VectorizedCB.h>
#include <HiggsAnalysis/CombinedLimit/interface/VectorizedSimplePdfs.h>
#include <HiggsAnalysis/CombinedLimit/interface/CachingMultiPdf.h>
#include "vectorized.h"

namespace cacheutils {
    typedef OptimizedCachingPdfT<FastVerticalInterpHistPdf,FastVerticalInterpHistPdfV> CachingHistPdf;
    typedef OptimizedCachingPdfT<FastVerticalInterpHistPdf2,FastVerticalInterpHistPdf2V> CachingHistPdf2;
    typedef OptimizedCachingPdfT<RooGaussian,VectorizedGaussian> CachingGaussPdf;
    typedef OptimizedCachingPdfT<RooCBShape,VectorizedCBShape> CachingCBPdf;
    typedef OptimizedCachingPdfT<RooExponential,VectorizedExponential> CachingExpoPdf;
    typedef OptimizedCachingPdfT<RooPower,VectorizedPower> CachingPowerPdf;

    class ReminderSum : public RooAbsReal {
        public:
            ReminderSum() {}
            ReminderSum(const char *name, const char *title, const RooArgList& sumSet) ;
            ReminderSum(const ReminderSum &other, const char* name = 0) :
                RooAbsReal(other, name),
                list_("deps",this,other.list_),
                terms_(other.terms_) {} 
            ~ReminderSum() {}
            virtual TObject* clone(const char* newname) const { return new ReminderSum(*this,newname); }
        private:
            RooListProxy list_;
            std::vector<const RooAbsReal *> terms_;
            Double_t evaluate() const;
    };
}

//---- Uncomment this to get a '.' printed every some evals
//#define TRACE_NLL_EVALS

//---- Uncomment this to get a total of the evals done
//#define TRACE_NLL_EVAL_COUNT

//---- Uncomment this and run with --perfCounters to get cache statistics
//#define DEBUG_CACHE

//---- Uncomment to dump PDF values inside CachingAddNLL
//#define LOG_ADDPDFS

//---- Uncomment to enable Kahan's summation (if enabled at runtime with --X-rtd = ...
// http://en.wikipedia.org/wiki/Kahan_summation_algorithm
//#define ADDNLL_KAHAN_SUM
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"

//std::map<std::string,double> cacheutils::CachingAddNLL::offsets_;
bool cacheutils::CachingSimNLL::noDeepLEE_ = false;
bool cacheutils::CachingSimNLL::hasError_  = false;
bool cacheutils::CachingSimNLL::optimizeContraints_  = true;

//#define DEBUG_TRACE_POINTS
#ifdef DEBUG_TRACE_POINTS
namespace { 
    template<unsigned int>
    void tracePoint(const RooAbsCollection &point) {
        static const RooAbsCollection *lastPoint = 0;
        static std::vector<double> values;
        if (&point != lastPoint) {
            std::cout << "Arrived in a completely new point. " << std::endl;
            values.resize(point.getSize());
            RooLinkedListIter iter = point.iterator();
            for (RooAbsArg *a  = (RooAbsArg*)(iter.Next()); a != 0; a  = (RooAbsArg*)(iter.Next())) {
                RooRealVar *rrv = dynamic_cast<RooRealVar *>(a); if (!rrv) continue;
                values.push_back(rrv->getVal());
            }
            lastPoint = &point;
        } else {
            std::cout << "Moved: ";
            RooLinkedListIter iter = point.iterator();
            int i = 0;
            for (RooAbsArg *a  = (RooAbsArg*)(iter.Next()); a != 0; a  = (RooAbsArg*)(iter.Next())) {
                RooRealVar *rrv = dynamic_cast<RooRealVar *>(a); if (!rrv) continue;
                if (values[i] != rrv->getVal()) std::cout << a->GetName() << " " << values[i] << " => " << rrv->getVal() << "    "; 
                values[i++] = rrv->getVal();
            }
            std::cout << std::endl;
        }
    }
}
#define TRACE_POINT2(x,i)  ::tracePoint<i>(x);
#define TRACE_POINT(x)  ::tracePoint<0>(x);
#define TRACE_NLL(x)    std::cout << x << std::endl;
#else
#define TRACE_POINT2(x,i)
#define TRACE_POINT(x) 
#define TRACE_NLL(x) 
#endif

#ifdef TRACE_NLL_EVAL_COUNT
namespace { unsigned long CachingSimNLLEvalCount = 0; }
#endif

cacheutils::ArgSetChecker::ArgSetChecker(const RooAbsCollection &set) 
{
    std::auto_ptr<TIterator> iter(set.createIterator());
    for (RooAbsArg *a  = dynamic_cast<RooAbsArg *>(iter->Next()); 
                    a != 0; 
                    a  = dynamic_cast<RooAbsArg *>(iter->Next())) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv) { // && !rrv->isConstant()) { 
            vars_.push_back(rrv);
            vals_.push_back(rrv->getVal());
            continue;
        }
        RooCategory *cat =  dynamic_cast<RooCategory *>(a);
        if (cat) {
            cats_.push_back(cat);
            states_.push_back(cat->getIndex()); 
        }
    }
}

bool 
cacheutils::ArgSetChecker::changed(bool updateIfChanged) 
{
    bool changed = false;
    std::vector<RooRealVar *>::const_iterator it = vars_.begin(), ed = vars_.end();
    std::vector<double>::iterator itv = vals_.begin();
    for ( ; it != ed; ++it, ++itv) {
        double val = (*it)->getVal();
        if (val != *itv) { 
            //std::cerr << "var::CachingPdf " << (*it)->GetName() << " changed: " << *itv << " -> " << val << std::endl;
            changed = true; 
            if (updateIfChanged) { *itv = val; }
            else break;
        }
    }
    std::vector<RooCategory *>::const_iterator itc = cats_.begin(), edc = cats_.end();
    std::vector<int>::iterator itvc = states_.begin();
    for ( ; itc != edc; ++itc, ++itvc) {
        int val = (*itc)->getIndex();
        if (val != *itvc) { 
            //std::cerr << "cat::CachingPdf " << (*itc)->GetName() << " changed: " << *itvc << " -> " << val << std::endl;
            changed = true; 
            if (updateIfChanged) { *itvc = val; }
            else break;
        }
    }

    return changed;
}

cacheutils::ValuesCache::ValuesCache(const RooAbsCollection &params, int size) :
    size_(1),
    maxSize_(size)
{
    assert(size <= MaxItems_);
    items[0] = new Item(params);
}
cacheutils::ValuesCache::ValuesCache(const RooAbsReal &pdf, const RooArgSet &obs, int size) :
    size_(1),
    maxSize_(size)
{
    assert(size <= MaxItems_);
    std::auto_ptr<RooArgSet> params(pdf.getParameters(obs));
    //std::cout << "Parameters for pdf " << pdf.GetName() << " (" << pdf.ClassName() << "):"; params->Print("");
    items[0] = new Item(*params);
}


cacheutils::ValuesCache::~ValuesCache() 
{
    for (int i = 0; i < size_; ++i) delete items[i];
}

void cacheutils::ValuesCache::clear() 
{
    for (int i = 0; i < size_; ++i) items[i]->good = false;
}

std::pair<std::vector<Double_t> *, bool> cacheutils::ValuesCache::get() 
{
    int found = -1; bool good = false;
    for (int i = 0; i < size_; ++i) {
        if (items[i]->good) {
            // valid entry, check if fresh
            if (!items[i]->checker.changed()) {
#ifdef DEBUG_CACHE
                PerfCounter::add(i == 0 ? "ValuesCache::get hit first" : "ValuesCache::get hit other");
#endif
                // fresh: done! 
                found = i; 
                good = true; 
                break;
            } 
        } else if (found == -1) {
            // invalid entry, can be replaced
            found = i;
#ifdef DEBUG_CACHE
            PerfCounter::add("ValuesCache::get hit invalid");
#endif
        }
    } 
    if (found == -1) {
        // all entries are valid but old 
#ifdef DEBUG_CACHE
        PerfCounter::add("ValuesCache::get miss");
#endif
        if (size_ < maxSize_) {
            // if I can, make a new entry
            items[size_] = new Item(items[0]->checker); // create a new item, copying the ArgSetChecker from the first one
            found = size_; 
            size_++;
        } else {
            // otherwise, pick the last one
            found = size_-1;
        }
    }
    // make sure new entry is the first one
    if (found != 0) {
        // remember what found is pointing to
        Item *f = items[found];
        // shift the other items down one place
        while (found > 0) { items[found] = items[found-1]; --found; } 
        // and put found on top
        items[found] = f;
    }
    if (!good) items[found]->checker.changed(true); // store new values in cache sentry
    items[found]->good = true;                      // mark this as valid entry
    return std::pair<std::vector<Double_t> *, bool>(&items[found]->values, good);
}

cacheutils::CachingPdf::CachingPdf(RooAbsReal *pdf, const RooArgSet *obs) :
    obs_(obs),
    pdfOriginal_(pdf),
    pdfPieces_(),
    pdf_(utils::fullCloneFunc(pdf, pdfPieces_)),
    lastData_(0),
    cache_(*pdf_,*obs_),
    includeZeroWeights_(false)
{
}

cacheutils::CachingPdf::CachingPdf(const CachingPdf &other) :
    obs_(other.obs_),
    pdfOriginal_(other.pdfOriginal_),
    pdfPieces_(),
    pdf_(utils::fullCloneFunc(pdfOriginal_, pdfPieces_)),
    lastData_(0),
    cache_(*pdf_,*obs_),
    includeZeroWeights_(other.includeZeroWeights_)
{
}

cacheutils::CachingPdf::~CachingPdf() 
{
}

const std::vector<Double_t> & 
cacheutils::CachingPdf::eval(const RooAbsData &data) 
{
#ifdef DEBUG_CACHE
    PerfCounter::add("CachingPdf::eval called");
#endif
    bool newdata = (lastData_ != &data);
    if (newdata) newData_(data);
    std::pair<std::vector<Double_t> *, bool> hit = cache_.get();
    if (!hit.second) { 
        realFill_(data, *hit.first);
    } 
    return *hit.first;
}

void
cacheutils::CachingPdf::newData_(const RooAbsData &data) 
{
    lastData_ = &data;
    pdf_->optimizeCacheMode(*data.get());
    pdf_->attachDataSet(data);
    const_cast<RooAbsData*>(lastData_)->setDirtyProp(false);
    cache_.clear();
    nonZeroW_.resize(data.numEntries());
    nonZeroWEntries_ = 0;
    for (unsigned int i = 0, n = nonZeroW_.size(); i < n; ++i) {
        data.get(i);
        if (data.weight() > 0 || includeZeroWeights_) {
            nonZeroWEntries_++;
            nonZeroW_[i] = 1; 
        } else {
            nonZeroW_[i] = 0;
        }
    }
}


void
cacheutils::CachingPdf::realFill_(const RooAbsData &data, std::vector<Double_t> &vals) 
{
#ifdef DEBUG_CACHE
    PerfCounter::add("CachingPdf::realFill_ called");
#endif
    //std::cout << "CachingPdf::realFill_ called for " << pdf_->GetName() << " (" << pdf_->ClassName() << ")\n";
    //utils::printPdf((RooAbsPdf*)pdf_);
    
    int n = data.numEntries();
    vals.resize(nonZeroWEntries_); // should be a no-op if size is already >= n.
    std::vector<Double_t>::iterator itv = vals.begin();
    for (int i = 0; i < n; ++i) {
        if (!nonZeroW_[i]) continue;
        data.get(i);
        *itv = pdf_->getVal(obs_); ++itv;
        //std::cout << " at i = " << i << " pdf = " << *itv << std::endl;
        TRACE_NLL("PDF value for " << pdf_->GetName() << " is " << *itv << " at this point.") 
        TRACE_POINT2(*obs_,1)
    }
}


template <typename PdfT, typename VPdfT>
void
cacheutils::OptimizedCachingPdfT<PdfT,VPdfT>::newData_(const RooAbsData &data) 
{
    lastData_ = &data;
    pdf_->optimizeCacheMode(*data.get());
    pdf_->attachDataSet(data);
    const_cast<RooAbsData*>(lastData_)->setDirtyProp(false);
    cache_.clear();
    delete vpdf_;
    vpdf_ = new VPdfT(static_cast<const PdfT &>(*pdf_), data, includeZeroWeights_);
}

template <typename PdfT, typename VPdfT>
void
cacheutils::OptimizedCachingPdfT<PdfT,VPdfT>::realFill_(const RooAbsData &data, std::vector<Double_t> &vals) 
{
    vpdf_->fill(vals);
}


cacheutils::ReminderSum::ReminderSum(const char *name, const char *title, const RooArgList& sumSet) :
    RooAbsReal(name,title),
    list_("deps","",this)
{
    RooLinkedListIter iter(sumSet.iterator());
    for (RooAbsReal *rar = (RooAbsReal *) iter.Next(); rar != 0; rar = (RooAbsReal *) iter.Next()) {
        list_.add(*rar);
        terms_.push_back(rar);
    }
}
Double_t cacheutils::ReminderSum::evaluate() const {
    Double_t ret = 1.0;
    for (auto x : terms_) ret -= x->getVal();
    return ret;
}

cacheutils::CachingAddNLL::CachingAddNLL(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data, bool includeZeroWeights) :
    RooAbsReal(name, title),
    pdf_(pdf),
    params_("params","parameters",this),
    includeZeroWeights_(includeZeroWeights),
    zeroPoint_(0)
{
    if (pdf == 0) throw std::invalid_argument(std::string("Pdf passed to ")+name+" is null");
    setData(*data);
    setup_();
}

cacheutils::CachingAddNLL::CachingAddNLL(const CachingAddNLL &other, const char *name) :
    RooAbsReal(name ? name : (TString("nll_")+other.pdf_->GetName()).Data(), ""),
    pdf_(other.pdf_),
    params_("params","parameters",this),
    includeZeroWeights_(other.includeZeroWeights_),
    zeroPoint_(0)
{
    setData(*other.data_);
    setup_();
}

cacheutils::CachingAddNLL::~CachingAddNLL() 
{
    for (int i = 0, n = integrals_.size(); i < n; ++i) delete integrals_[i];
    integrals_.clear();
}

cacheutils::CachingAddNLL *
cacheutils::CachingAddNLL::clone(const char *name) const 
{
    return new cacheutils::CachingAddNLL(*this, name);
}

void cacheutils::CachingAddNLL::addPdfs_(RooAddPdf *addpdf, bool recursive, const RooArgList & basecoeffs)
{
    int npdf = addpdf->pdfList().getSize();
    //std::cout << "Unpacking RooAddPdf " << addpdf->GetName() << " with " << npdf << " components:" << std::endl;
    RooAbsReal *lastcoeff = 0;
    if (npdf == addpdf->coefList().getSize()) {
        lastcoeff =  dynamic_cast<RooAbsReal*>(addpdf->coefList().at(npdf-1));
    } else {
        prods_.push_back(new ReminderSum((std::string("reminder_of_")+addpdf->GetName()).c_str(),"", addpdf->coefList()));
        lastcoeff = & prods_.back(); 
    }
    //std::cout << "   Last coefficient is a " << lastcoeff->ClassName() << " aka " << typeid(*lastcoeff).name() << ": "; lastcoeff->Print("");
    for (int i = 0; i < npdf; ++i) {
        RooAbsReal * coeff = (i < npdf-1 ? dynamic_cast<RooAbsReal*>(addpdf->coefList().at(i)) : lastcoeff);
        RooAbsPdf  * pdfi  = dynamic_cast<RooAbsPdf *>(addpdf->pdfList().at(i));
        if (recursive && typeid(*pdfi) == typeid(RooAddPdf)) {
            RooAddPdf *apdfi = static_cast<RooAddPdf*>(pdfi);
            RooArgList list(*coeff);
            if (basecoeffs.getSize()) list.add(basecoeffs);
            //std::cout << "    Invoking recursive unpack on " << addpdf->GetName() << "[" << i << "]: RooAddPdf " << apdfi->GetName() << std::endl;
            addPdfs_(apdfi, recursive, list);
            continue;
        } 
        if (basecoeffs.getSize() == 0) {
            coeffs_.push_back(coeff);
        } else {
            RooArgList list(*coeff);
            list.add(basecoeffs);
            prods_.push_back(new RooProduct("","",list));
            coeffs_.push_back(&prods_.back());
            //std::cout << "Coefficient of " << pdfi->GetName() << std::endl; prods_.back().Print("");
        }
        const RooArgSet *obs = data_->get();
        pdfs_.push_back(makeCachingPdf(pdfi,obs));
        pdfs_.back().setIncludeZeroWeights(includeZeroWeights_);
        //std::cout << "    Adding " << addpdf->GetName() << "[" << i << "]: " << pdfi->ClassName() << " " << pdfi->GetName() << " using " << typeid(pdfs_.back()).name() << std::endl;
    }
}

cacheutils::CachingPdfBase *
cacheutils::makeCachingPdf(RooAbsReal *pdf, const RooArgSet *obs) {
    static bool histNll  = runtimedef::get("ADDNLL_HISTNLL");
    static bool gaussNll  = runtimedef::get("ADDNLL_GAUSSNLL");
    static bool multiNll  = runtimedef::get("ADDNLL_MULTINLL");
    static bool cbNll  = runtimedef::get("ADDNLL_CBNLL");

    if (histNll && typeid(*pdf) == typeid(FastVerticalInterpHistPdf)) {
        return new CachingHistPdf(pdf, obs);
    } else if (histNll && typeid(*pdf) == typeid(FastVerticalInterpHistPdf2)) {
        return new CachingHistPdf2(pdf, obs);
    } else if (gaussNll && typeid(*pdf) == typeid(RooGaussian)) {
        return new CachingGaussPdf(pdf, obs);
    } else if (cbNll && typeid(*pdf) == typeid(RooCBShape)) {
        return new CachingCBPdf(pdf, obs);
    } else if (gaussNll && typeid(*pdf) == typeid(RooExponential)) {
	std::auto_ptr<RooArgSet> params(pdf->getParameters(obs));
	if(params->getSize()!=1) {return new CachingPdf(pdf,obs);}
        return new CachingExpoPdf(pdf, obs);
    } else if (gaussNll && typeid(*pdf) == typeid(RooPower)) {
        return new CachingPowerPdf(pdf, obs);
    } else if (multiNll && typeid(*pdf) == typeid(RooMultiPdf)) {
        return new CachingMultiPdf(static_cast<RooMultiPdf&>(*pdf), *obs);
    } else if (multiNll && typeid(*pdf) == typeid(RooAddPdf)) {
        return new CachingAddPdf(static_cast<RooAddPdf&>(*pdf), *obs);
    } else {
        return new CachingPdf(pdf, obs);
    }

}

void
cacheutils::CachingAddNLL::setup_() 
{
    fastExit_ = !runtimedef::get("NO_ADDNLL_FASTEXIT");
    for (int i = 0, n = integrals_.size(); i < n; ++i) delete integrals_[i];
    integrals_.clear(); pdfs_.clear(); coeffs_.clear(); prods_.clear();
    RooAddPdf *addpdf = 0;
    RooRealSumPdf *sumpdf = 0;
    if ((addpdf = dynamic_cast<RooAddPdf *>(pdf_)) != 0) {
        isRooRealSum_ = false; basicIntegrals_ = 0;
        addPdfs_(addpdf, runtimedef::get("ADDNLL_RECURSIVE"), RooArgList());
    } else if ((sumpdf = dynamic_cast<RooRealSumPdf *>(pdf_)) != 0) {
        const RooArgSet *obs = data_->get();
        isRooRealSum_ = true;  basicIntegrals_ = canBasicIntegrals_;
        int npdf = sumpdf->coefList().getSize();
        coeffs_.reserve(npdf);
        pdfs_.reserve(npdf);
        integrals_.reserve(npdf);
        for (int i = 0; i < npdf; ++i) {
            RooAbsReal * coeff = dynamic_cast<RooAbsReal*>(sumpdf->coefList().at(i));
            RooAbsReal * funci = dynamic_cast<RooAbsReal*>(sumpdf->funcList().at(i));
            static int tryfactor = runtimedef::get("ADDNLL_ROOREALSUM_FACTOR"); 
            RooProduct *prodi = 0;
            if (tryfactor && ((prodi = dynamic_cast<RooProduct *>(funci)) != 0)) {
                RooArgList newcoeffs(*coeff), newfuncs; 
                utils::factorizeFunc(*obs, *funci, newfuncs, newcoeffs);
                if (newcoeffs.getSize() > 1) {
                    prods_.push_back(new RooProduct("","",newcoeffs));
                    coeff = &prods_.back();
                }
                if (newfuncs.getSize() > 1) {
                    prods_.push_back(new RooProduct("","",newfuncs));
                    funci = &prods_.back();
                } else {
                    funci =  dynamic_cast<RooAbsReal *>(newfuncs.first());
                }
            }
            coeffs_.push_back(coeff);
            pdfs_.push_back(new CachingPdf(funci, obs));
            pdfs_.back().setIncludeZeroWeights(includeZeroWeights_);
            integrals_.push_back(funci->createIntegral(*obs));
        }
    } else {
        std::string errmsg = "ERROR: CachingAddNLL: Pdf ";
        errmsg += pdf_->GetName();
        errmsg += " is neither a RooAddPdf nor a RooRealSumPdf, but a ";
        errmsg += pdf_->ClassName();
        throw std::invalid_argument(errmsg);
    }

    std::auto_ptr<RooArgSet> params(pdf_->getParameters(*data_));
    std::auto_ptr<TIterator> iter(params->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        //if (rrv != 0 && !rrv->isConstant()) params_.add(*rrv);
        if (rrv != 0) params_.add(*rrv);
    }

    multiPdfs_.clear();
    for (auto itp = pdfs_.begin(), edp = pdfs_.end(); itp != edp; ++itp) {
        bool isMultiPdf = itp->pdf()->IsA()->InheritsFrom(RooMultiPdf::Class());
        if (isMultiPdf) {
            const RooMultiPdf *mpdf = dynamic_cast<const RooMultiPdf*>((*itp).pdf());
            multiPdfs_.push_back(std::make_pair(mpdf, &*itp));
        }
    }
 
}

void
cacheutils::CachingAddNLL::setIncludeZeroWeights(bool includeZeroWeights) 
{
    includeZeroWeights_ = includeZeroWeights;
    for (CachingPdfBase &pdf : pdfs_) {
        pdf.setIncludeZeroWeights(includeZeroWeights_);
    }
}

Double_t 
cacheutils::CachingAddNLL::evaluate() const 
{
#ifdef DEBUG_CACHE
    PerfCounter::add("CachingAddNLL::evaluate called");
#endif

    // For multi pdf's need to reset the cache if index changed before evaluations
    // unless they're being properly treated in the CachingPdf
    static bool multiNll  = runtimedef::get("ADDNLL_MULTINLL");
    if (!multiNll && !multiPdfs_.empty()) {
        for (std::vector<std::pair<const RooMultiPdf*,CachingPdfBase*> >::iterator itp = multiPdfs_.begin(), edp = multiPdfs_.end(); itp != edp; ++itp) {
		bool hasChangedPdf = itp->first->checkIndexDirty();
		if (hasChangedPdf) itp->second->setDataDirty();
        }
    }

    std::fill( partialSum_.begin(), partialSum_.end(), 0.0 );

    std::vector<RooAbsReal*>::iterator  itc = coeffs_.begin(), edc = coeffs_.end();
    boost::ptr_vector<CachingPdfBase>::iterator   itp = pdfs_.begin();//,   edp = pdfs_.end();
    std::vector<Double_t>::const_iterator itw, bgw = weights_.begin();//,    edw = weights_.end();
    std::vector<Double_t>::iterator       its, bgs = partialSum_.begin(), eds = partialSum_.end();
    double sumCoeff = 0;
    bool allBasicIntegralsOk = (basicIntegrals_ == 1);
    //std::cout << "Performing evaluation of " << GetName() << std::endl;
    for ( ; itc != edc; ++itp, ++itc ) {
        // get coefficient
        Double_t coeff = (*itc)->getVal();
        if (isRooRealSum_ && basicIntegrals_ < 2) {
            sumCoeff += coeff * integrals_[itc - coeffs_.begin()]->getVal();
            //std::cout << "  coefficient = " << coeff << ", integral = " << integrals_[itc - coeffs_.begin()]->getVal() << std::endl;
        } else {
            sumCoeff += coeff;
        }
        // get vals
        const std::vector<Double_t> &pdfvals = itp->eval(*data_);
        if (basicIntegrals_) {
            double integral = (binWidths_.size() > 1) ? 
                                    vectorized::dot_product(pdfvals.size(), &pdfvals[0], &binWidths_[0]) :
                                    binWidths_.front() * std::accumulate(pdfvals.begin(), pdfvals.end(), 0.0);
            if (basicIntegrals_ == 1) {
                double refintegral = integrals_[itc - coeffs_.begin()]->getVal();
                if (refintegral > 0) {
                    if (std::abs((integral - refintegral)/refintegral) > 1e-5) {
                        printf("integrals don't match: %+10.6f  %+10.6f  %10.7f %s\n", refintegral, integral, refintegral ? std::abs((integral - refintegral)/refintegral) : 0,  itp->pdf()->GetName());
                        allBasicIntegralsOk = false;
                        basicIntegrals_ = 0; // don't waste time on this anymore
                    }
                }
            } else {
                sumCoeff += coeff*(integral - 1.0); // I had added up coeff before.
            }
        }
#ifdef LOG_ADDPDFS
        printf("%s coefficient %s (%s) = %20.15f\n", itp->pdf()->GetName(), (*itc)->GetName(), (*itc)->ClassName(), coeff);
        //(*itc)->Print("");
        for (unsigned int i = 0, n = pdfvals.size(); i < n; ++i) {
            if (i%84==0) printf("%-80s[%3d] = %20.15f\n", itp->pdf()->GetName(), i, pdfvals[i]);
        }
#endif
        // update running sum
        //    std::vector<Double_t>::const_iterator itv = pdfvals.begin();
        //    for (its = bgs; its != eds; ++its, ++itv) {
        //         *its += coeff * (*itv); // sum (n_i * pdf_i)
        //    }
        // vectorize to make it faster
        vectorized::mul_add(pdfvals.size(), coeff, &pdfvals[0], &partialSum_[0]);
    }
    // if all basic integrals evaluated ok, use them
    if (allBasicIntegralsOk) basicIntegrals_ = 2;
    // then get the final nll
    double ret = 0;
    for (its = bgs; its != eds ; ++its) {
        if (!isnormal(*its) || *its <= 0) {
            if ((weights_[its-bgs] == 0) && (*its == 0)) {
                // this special case we don't care, as zero is fine and it will be multiplied by zero afterwards,
                // we just need to protect it for the logarithm
                *its = 1.0; // arbitrary number, to avoid log(0)
                continue;
            } else if (weights_[its-bgs] == 0) {
                // this is a special case we should in principle care, even if it does not alter the likelihood
                // since it's multiplied by zero. However, normally RooFit ignores errors in zero-weight bins,
                // so we comply to his policy (but we issue a warning, and we protect the logarithm)
                static int nwarn = 0;
                if (++nwarn < 100) {
                    std::cout << "WARNING: underflow to " << *its << " in " << pdf_->GetName() << " for zero-entry bin " << its-bgs << std::endl;
                }
                *its = 1.0; // arbitrary number, to avoid bad logs
                continue;
            }
            std::cout << "WARNING: underflow to " << *its << " in " << pdf_->GetName() << " for bin " << its-bgs << ", weight " << weights_[its-bgs] << std::endl; 
            if (!CachingSimNLL::noDeepLEE_) logEvalError("Number of events is negative or error"); else CachingSimNLL::hasError_ = true;
            if (fastExit_) { return 9e9; }
            else *its = 1;
        }
    }
    #ifndef ADDNLL_KAHAN_SUM
    // Do the reduction 
    //      for ( its = bgs, itw = bgw ; its != eds ; ++its, ++itw ) {
    //         ret += (*itw) * log( ((*its) / sumCoeff) );
    //      }
    ret += vectorized::nll_reduce(partialSum_.size(), &partialSum_[0], &weights_[0], sumCoeff, &workingArea_[0]);
    #else
    double compensation = 0;
    static bool do_kahan = runtimedef::get("ADDNLL_KAHAN_SUM");
    for ( its = bgs, itw = bgw ; its != eds ; ++its, ++itw ) {
        double thispiece = (*itw) * log( ((*its) / sumCoeff) );
        if (do_kahan) {
            double kahan_y = thispiece  - compensation;
            double kahan_t = ret + kahan_y;
            double kahan_d = (kahan_t - ret);
            compensation = kahan_d - kahan_y;
            ret  = kahan_t;
        } else {
            ret += thispiece;
        }
        ret += thispiece;
    }
    #endif
    // then flip sign
    ret = -ret;
    // std::cout << "AddNLL for " << pdf_->GetName() << ": " << ret << std::endl;
    // and add extended term: expected - observed*log(expected);
    static bool expEventsNoNorm = runtimedef::get("ADDNLL_ROOREALSUM_NONORM");
    double expectedEvents = (isRooRealSum_ && !expEventsNoNorm ? pdf_->getNorm(data_->get()) : sumCoeff);
    if (expectedEvents <= 0) {
        if (!CachingSimNLL::noDeepLEE_) logEvalError("Expected number of events is negative"); else CachingSimNLL::hasError_ = true;
        expectedEvents = 1e-6;
    }
    //ret += expectedEvents - UInt_t(sumWeights_) * log(expectedEvents); // no, doesn't work with Asimov dataset
    ret += expectedEvents - sumWeights_ * log(expectedEvents);
    ret += zeroPoint_;

    // multipdfs want to add a correction factor to the NLL
    if (!multiPdfs_.empty()) {
        double correctionFactor = 0;
        for (std::vector<std::pair<const RooMultiPdf*,CachingPdfBase*> >::iterator itp = multiPdfs_.begin(), edp = multiPdfs_.end(); itp != edp; ++itp) {
            correctionFactor += itp->first->getCorrection();
        }
        // Add correction 
        ret+=correctionFactor;
    }

    TRACE_NLL("AddNLL for " << pdf_->GetName() << ": " << ret)
    return ret;
}

void 
cacheutils::CachingAddNLL::setData(const RooAbsData &data) 
{
    //std::cout << "Setting data for pdf " << pdf_->GetName() << std::endl;
    //utils::printRAD(&data);
    data_ = &data;
    setValueDirty();
    sumWeights_ = 0.0;
    weights_.clear(); weights_.reserve(data.numEntries());
    #ifdef ADDNLL_KAHAN_SUM
    double compensation = 0;
    #endif
    for (int i = 0, n = data.numEntries(); i < n; ++i) {
        data.get(i);
        double w = data.weight();
        if (w || includeZeroWeights_) weights_.push_back(w); 
        #ifdef ADDNLL_KAHAN_SUM
        static bool do_kahan = runtimedef::get("ADDNLL_KAHAN_SUM");
        if (do_kahan) {
            double kahan_y = w - compensation;
            double kahan_t = sumWeights_ + kahan_y;
            double kahan_d = (kahan_t - sumWeights_);
            compensation = kahan_d - kahan_y;
            sumWeights_  = kahan_t;
        } else {
            sumWeights_ += w;
        }
        #else
        sumWeights_ += w;
        #endif
    }
    partialSum_.resize(weights_.size());
    workingArea_.resize(weights_.size());
    for (auto itp = pdfs_.begin(), edp = pdfs_.end(); itp != edp; ++itp) {
        itp->setDataDirty();
    }
    binWidths_.clear(); canBasicIntegrals_ = 0;
    if (dynamic_cast<RooRealSumPdf *>(pdf_) != 0 && runtimedef::get("ADDNLL_ROOREALSUM_BASICINT") > 0) {
        const RooArgSet *obs = data_->get();
        RooRealVar *xvar = dynamic_cast<RooRealVar *>(obs->first());
        if (obs->getSize() == 1 && xvar != 0 && xvar->numBins() == data_->numEntries()) {
            const RooAbsBinning &bins = xvar->getBinning();
            binWidths_.resize(xvar->numBins());
            bool all_equal = true;
            canBasicIntegrals_ = runtimedef::get("ADDNLL_ROOREALSUM_BASICINT");
            for (unsigned int ibin = 0, nbin = binWidths_.size(); ibin < nbin; ++ibin) {
                double bc = bins.binCenter(ibin), dc = data_->get(ibin)->getRealValue(xvar->GetName());
                //printf("bin %3d: center %+8.5f ( data %+8.5f , diff %+8.5f ), width %8.5f, data weight %10.5f, channel %s\n", ibin, bc, dc, abs(dc-bc)/bins.binWidth(ibin), bins.binWidth(ibin), data_->weight(), pdf_->GetName());
                binWidths_[ibin] = bins.binWidth(ibin);
                if (std::abs(bc-dc) > 1e-5*binWidths_[ibin]) {
                    printf("channel %s, for observable %s, bin %d mismatch: binning %+8.5f ( data %+8.5f , diff %+7.8f of width %8.5f\n",
                        pdf_->GetName(), xvar->GetName(), ibin, bc, dc, std::abs(bc-dc)/binWidths_[ibin], binWidths_[ibin]);
                    canBasicIntegrals_ = 0; 
                    break;
                }
                if ((ibin > 0) && (binWidths_[ibin] != binWidths_[ibin-1])) all_equal = false;
            }
            if (all_equal) binWidths_.resize(1);
        } else {
            printf("channel %s (binned likelihood? %d), can't do binned intergals. nobs %d, obs %s, nbins %d, ndata %d\n", pdf_->GetName(), pdf_->getAttribute("BinnedLikelihood"), obs->getSize(), (xvar ? xvar->GetName() : "<nil>"), (xvar ? xvar->numBins() : -999), data_->numEntries());
        }
    }
}

RooArgSet* 
cacheutils::CachingAddNLL::getObservables(const RooArgSet* depList, Bool_t valueOnly) const 
{
    return new RooArgSet();
}

RooArgSet* 
cacheutils::CachingAddNLL::getParameters(const RooArgSet* depList, Bool_t stripDisconnected) const 
{
    return new RooArgSet(params_); 
}


cacheutils::CachingSimNLL::CachingSimNLL(RooSimultaneous *pdf, RooAbsData *data, const RooArgSet *nuis) :
    pdfOriginal_(pdf),
    dataOriginal_(data),
    nuis_(nuis),
    params_("params","parameters",this)
{
    setup_();
}

cacheutils::CachingSimNLL::CachingSimNLL(const CachingSimNLL &other, const char *name) :
    pdfOriginal_(other.pdfOriginal_),
    dataOriginal_(other.dataOriginal_),
    nuis_(other.nuis_),
    params_("params","parameters",this)
{
    setup_();
}

cacheutils::CachingSimNLL *
cacheutils::CachingSimNLL::clone(const char *name) const 
{
    return new cacheutils::CachingSimNLL(*this, name);
}

cacheutils::CachingSimNLL::~CachingSimNLL()
{
    std::vector<bool>::const_iterator ito = constrainPdfsFastOwned_.begin();
    for (std::vector<SimpleGaussianConstraint*>::iterator it = constrainPdfsFast_.begin(), ed = constrainPdfsFast_.end(); it != ed; ++it, ++ito) {
        if (*ito) { delete *it; }
    }
    #ifdef TRACE_NLL_EVAL_COUNT
        std::cout << "CachingSimNLLEvalCount: " << ::CachingSimNLLEvalCount << std::endl;
    #endif
}

void
cacheutils::CachingSimNLL::setup_() 
{
    // Allow runtime-flag to switch off logEvalErrors
    noDeepLEE_ = runtimedef::get("SIMNLL_NO_LEE");

    //RooAbsPdf *pdfclone = runtimedef::get("SIMNLL_CLONE") ? pdfOriginal_  : utils::fullClonePdf(pdfOriginal_, piecesForCloning_);
    RooAbsPdf *pdfclone = pdfOriginal_; // never clone

    //---- Instead of getting the parameters here, we get them from the individual constraint terms and single pdfs ----
    //---- This seems to save memory.
    //std::auto_ptr<RooArgSet> params(pdfclone->getParameters(*dataOriginal_));
    //params_.add(*params);

    RooArgList constraints;
    factorizedPdf_.reset(dynamic_cast<RooSimultaneous *>(utils::factorizePdf(*dataOriginal_->get(), *pdfclone, constraints)));
    
    RooSimultaneous *simpdf = factorizedPdf_.get();
    constrainPdfs_.clear(); 
    if (constraints.getSize()) {
        int FastConstraints = optimizeContraints_ && runtimedef::get("SIMNLL_FASTGAUSS");
        for (int i = 0, n = constraints.getSize(); i < n; ++i) {
            RooAbsPdf *pdfi = dynamic_cast<RooAbsPdf*>(constraints.at(i));
            if (optimizeContraints_ && typeid(*pdfi) == typeid(SimpleGaussianConstraint)) {
                constrainPdfsFast_.push_back(static_cast<SimpleGaussianConstraint *>(pdfi));
                constrainPdfsFastOwned_.push_back(false);
                constrainZeroPointsFast_.push_back(0);
            } else if (FastConstraints && typeid(*pdfi) == typeid(RooGaussian)) {
                constrainPdfsFast_.push_back(new SimpleGaussianConstraint(dynamic_cast<const RooGaussian&>(*pdfi)));
                constrainPdfsFastOwned_.push_back(true);
                constrainZeroPointsFast_.push_back(0);
            } else {
                constrainPdfs_.push_back(pdfi);
                constrainZeroPoints_.push_back(0);
            }
            //std::cout << "Constraint pdf: " << constraints.at(i)->GetName() << std::endl;
            std::auto_ptr<RooArgSet> params(pdfi->getParameters(*dataOriginal_));
            params_.add(*params, false);
        }
    } else {
        std::cerr << "PDF didn't factorize!" << std::endl;
        std::cout << "Parameters: " << std::endl;
        std::auto_ptr<RooArgSet> params(pdfclone->getParameters(*dataOriginal_));
        params->Print("V");
        std::cout << "Obs: " << std::endl;
        dataOriginal_->get()->Print("V");
        factorizedPdf_.release();
        simpdf = dynamic_cast<RooSimultaneous *>(pdfclone);
    }

    
    std::auto_ptr<RooAbsCategoryLValue> catClone((RooAbsCategoryLValue*) simpdf->indexCat().Clone());
    pdfs_.resize(catClone->numBins(NULL), 0);
    //dataSets_.reset(dataOriginal_->split(pdfOriginal_->indexCat(), true));
    datasets_.resize(pdfs_.size(), 0);
    splitWithWeights(*dataOriginal_, simpdf->indexCat(), true);
    //std::cout << "Pdf " << simpdf->GetName() <<" is a SimPdf over category " << catClone->GetName() << ", with " << pdfs_.size() << " bins" << std::endl;
    for (int ib = 0, nb = pdfs_.size(); ib < nb; ++ib) {
        catClone->setBin(ib);
        RooAbsPdf *pdf = simpdf->getPdf(catClone->getLabel());
        if (pdf != 0) {
            RooAbsData *data = (RooAbsData *) datasets_[ib]; //dataSets_->FindObject(catClone->getLabel());
            //RooAbsData *data = (RooAbsData *) dataSets_->FindObject(catClone->getLabel());
            //std::cout << "   bin " << ib << " (label " << catClone->getLabel() << ") has pdf " << pdf->GetName() << " of type " << pdf->ClassName() << " and " << (data ? data->numEntries() : -1) << " dataset entries" << std::endl;
            if (data == 0) { throw std::logic_error("Error: no data"); }
            bool includeZeroWeights = (runtimedef::get("ADDNLL_ROOREALSUM_BASICINT") && runtimedef::get("ADDNLL_ROOREALSUM_KEEPZEROS") && (dynamic_cast<RooRealSumPdf*>(pdf)!=0));
            pdfs_[ib] = new CachingAddNLL(catClone->getLabel(), "", pdf, data, includeZeroWeights);
            params_.add(pdfs_[ib]->params(), /*silent=*/true); 
        } else { 
            pdfs_[ib] = 0; 
            //std::cout << "   bin " << ib << " (label " << catClone->getLabel() << ") has no pdf" << std::endl;
        }
    }   

    setValueDirty();
}

Double_t 
cacheutils::CachingSimNLL::evaluate() const 
{
    TRACE_POINT(params_)
#ifdef TRACE_NLL_EVAL_COUNT
    ::CachingSimNLLEvalCount++;
#endif
#ifdef DEBUG_CACHE
    PerfCounter::add("CachingSimNLL::evaluate called");
#endif
    double ret = 0;
    for (std::vector<CachingAddNLL*>::const_iterator it = pdfs_.begin(), ed = pdfs_.end(); it != ed; ++it) {
        if (*it != 0) {
            double nllval = (*it)->getVal();
            // what sanity check could I put here?
            ret += nllval;
        }
    }
    if (!constrainPdfs_.empty() || !constrainPdfsFast_.empty()) {
        /// ============= GENERIC CONSTRAINTS  =========
        std::vector<double>::const_iterator itz = constrainZeroPoints_.begin();
        for (std::vector<RooAbsPdf *>::const_iterator it = constrainPdfs_.begin(), ed = constrainPdfs_.end(); it != ed; ++it, ++itz) { 
            double pdfval = (*it)->getVal(nuis_);
            if (!isnormal(pdfval) || pdfval <= 0) {
                if (!noDeepLEE_) logEvalError((std::string("Constraint pdf ")+(*it)->GetName()+" evaluated to zero, negative or error").c_str());
                pdfval = 1e-9;
            }
            ret -= (log(pdfval) + *itz);
        }
        /// ============= FAST GAUSSIAN CONSTRAINTS  =========
        itz = constrainZeroPointsFast_.begin();
        for (std::vector<SimpleGaussianConstraint*>::const_iterator it = constrainPdfsFast_.begin(), ed = constrainPdfsFast_.end(); it != ed; ++it, ++itz) { 
            double logpdfval = (*it)->getLogValFast();
            //std::cout << "pdf " << (*it)->GetName() << " = " << logpdfval << std::endl;
            ret -= (logpdfval + *itz);
        }
    }
#ifdef TRACE_NLL_EVALS
    static unsigned long _trace_ = 0; _trace_++;
    if (_trace_ % 10 == 0)  { putchar('.'); fflush(stdout); }
    //if (_trace_ % 250 == 0) { printf("               NLL % 10.4f after %10lu evals.\n", ret, _trace_); fflush(stdout); }
#endif
    TRACE_NLL("SimNLL for " << GetName() << ": " << ret)
    return ret;
}

void 
cacheutils::CachingSimNLL::setData(const RooAbsData &data) 
{
    dataOriginal_ = &data;
    //std::cout << "combined data has " << data.numEntries() << " dataset entries (sumw " << data.sumEntries() << ", weighted " << data.isWeighted() << ")" << std::endl;
    //utils::printRAD(&data);
    //dataSets_.reset(dataOriginal_->split(pdfOriginal_->indexCat(), true));
    splitWithWeights(*dataOriginal_, pdfOriginal_->indexCat(), true);
    for (int ib = 0, nb = pdfs_.size(); ib < nb; ++ib) {
        CachingAddNLL *canll = pdfs_[ib];
        if (canll == 0) continue;
        RooAbsData *data = datasets_[ib];
        //RooAbsData *data = (RooAbsData *) dataSets_->FindObject(canll->GetName());
        if (data == 0) { throw std::logic_error("Error: no data"); }
        //std::cout << "   bin " << ib << " (label " << canll->GetName() << ") has pdf " << canll->pdf()->GetName() << " of type " << canll->pdf()->ClassName() <<
        //             " and " << (data ? data->numEntries() : -1) << " dataset entries (sumw " << data->sumEntries() << ", weighted " << data->isWeighted() << ")" << std::endl;
        canll->setData(*data);
    }
}

void cacheutils::CachingSimNLL::splitWithWeights(const RooAbsData &data, const RooAbsCategory& splitCat, Bool_t createEmptyDataSets) {
    RooCategory *cat = dynamic_cast<RooCategory *>(data.get()->find(splitCat.GetName()));
    if (cat == 0) throw std::logic_error("Error: no category");
    std::auto_ptr<RooAbsCategoryLValue> catClone((RooAbsCategoryLValue*) splitCat.Clone());
    int nb = cat->numBins((const char *)0), ne = data.numEntries();
    RooArgSet obs(*data.get()); obs.remove(*cat, true, true);
    RooRealVar weight("_weight_","",1);
    RooArgSet obsplus(obs); obsplus.add(weight);
    if (nb != int(datasets_.size())) throw std::logic_error("Number of categories changed"); // this can happen due to bugs in RooDataSet
    std::vector<int> includeZeroWeights(nb,0); bool includeZeroWeightsAny = false;
    if (runtimedef::get("ADDNLL_ROOREALSUM_BASICINT") && runtimedef::get("ADDNLL_ROOREALSUM_KEEPZEROS")) {
        for (int ib = 0; ib < nb; ++ib) {
            catClone->setBin(ib);
            RooAbsPdf *pdf = factorizedPdf_->getPdf(catClone->getLabel());
            if (dynamic_cast<RooRealSumPdf*>(pdf)!=0) {
                includeZeroWeights[ib] = 1;
                includeZeroWeightsAny = true;
            }
        }
    }
    for (int ib = 0; ib < nb; ++ib) {
        if (datasets_[ib] == 0) {
            catClone->setBin(ib);
            RooAbsPdf *pdf = pdfOriginal_->getPdf(catClone->getLabel());
            if (pdf) {
                std::auto_ptr<RooArgSet> myobs(pdf->getObservables(obs));
                myobs->add(weight);
                //std::cout << "Observables for bin " << ib << ":"; myobs->Print("");
                datasets_[ib] = new RooDataSet("", "", *myobs, "_weight_");
            } else {
                datasets_[ib] = new RooDataSet("", "", obsplus, "_weight_");
            }
        } else {
            datasets_[ib]->reset();
        }
    }
    //utils::printRDH((RooAbsData*)&data);
    for (int i = 0; i < ne; ++i) {
        data.get(i); if (data.weight() == 0 && !includeZeroWeightsAny) continue;
        int ib = cat->getBin();
        //std::cout << "Event " << i << " of weight " << data.weight() << " is in bin " << ib << " label " << cat->getLabel() << std::endl;
        if (data.weight() > 0 || includeZeroWeights[ib]) datasets_[ib]->add(obs, data.weight());
    }
}

void cacheutils::CachingSimNLL::setZeroPoint() {
    for (std::vector<CachingAddNLL*>::const_iterator it = pdfs_.begin(), ed = pdfs_.end(); it != ed; ++it) {
        if (*it != 0) (*it)->setZeroPoint();
    }
    std::vector<double>::iterator itz = constrainZeroPoints_.begin();
    for (std::vector<RooAbsPdf *>::const_iterator it = constrainPdfs_.begin(), ed = constrainPdfs_.end(); it != ed; ++it, ++itz) {
        double pdfval = (*it)->getVal(nuis_);
        if (isnormal(pdfval) || pdfval > 0) *itz = -log(pdfval);
    }
    itz = constrainZeroPointsFast_.begin();
    for (std::vector<SimpleGaussianConstraint*>::const_iterator it = constrainPdfsFast_.begin(), ed = constrainPdfsFast_.end(); it != ed; ++it, ++itz) {
        double logpdfval = (*it)->getLogValFast();
        *itz = -logpdfval;
    }
    setValueDirty();
}

void cacheutils::CachingSimNLL::clearZeroPoint() {
    for (std::vector<CachingAddNLL*>::const_iterator it = pdfs_.begin(), ed = pdfs_.end(); it != ed; ++it) {
        if (*it != 0) (*it)->clearZeroPoint();
    }
    std::fill(constrainZeroPoints_.begin(), constrainZeroPoints_.end(), 0.0);
    std::fill(constrainZeroPointsFast_.begin(), constrainZeroPointsFast_.end(), 0.0);
    setValueDirty();
}

RooArgSet* 
cacheutils::CachingSimNLL::getObservables(const RooArgSet* depList, Bool_t valueOnly) const 
{
    return new RooArgSet();
}

RooArgSet* 
cacheutils::CachingSimNLL::getParameters(const RooArgSet* depList, Bool_t stripDisconnected) const 
{
    return new RooArgSet(params_); 
}
