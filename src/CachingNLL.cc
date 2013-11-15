#include "../interface/CachingNLL.h"
#include "../interface/utils.h"
#include <stdexcept>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooProduct.h>
#include "../interface/ProfilingTools.h"
#include <../interface/RooMultiPdf.h>
#include <../interface/VerticalInterpHistPdf.h>
#include <../interface/VectorizedGaussian.h>
#include "vectorized.h"

namespace cacheutils {
    typedef OptimizedCachingPdfT<FastVerticalInterpHistPdf,FastVerticalInterpHistPdfV> CachingHistPdf;
    typedef OptimizedCachingPdfT<FastVerticalInterpHistPdf2,FastVerticalInterpHistPdf2V> CachingHistPdf2;
    typedef OptimizedCachingPdfT<RooGaussian,VectorizedGaussian> CachingGaussPdf;

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

//---- Uncomment this and run with --perfCounters to get cache statistics
//#define DEBUG_CACHE

//---- Uncomment to enable Kahan's summation (if enabled at runtime with --X-rtd = ...
// http://en.wikipedia.org/wiki/Kahan_summation_algorithm
//#define ADDNLL_KAHAN_SUM
#include "../interface/ProfilingTools.h"

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
        }
    }
}

bool 
cacheutils::ArgSetChecker::changed(bool updateIfChanged) 
{
    std::vector<RooRealVar *>::const_iterator it = vars_.begin(), ed = vars_.end();
    std::vector<double>::iterator itv = vals_.begin();
    bool changed = false;
    for ( ; it != ed; ++it, ++itv) {
        double val = (*it)->getVal();
        if (val != *itv) { 
            //std::cerr << "var::CachingPdfable " << (*it)->GetName() << " changed: " << *itv << " -> " << val << std::endl;
            changed = true; 
            if (updateIfChanged) { *itv = val; }
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
    cache_(*pdf_,*obs_)
{
}

cacheutils::CachingPdf::CachingPdf(const CachingPdf &other) :
    obs_(other.obs_),
    pdfOriginal_(other.pdfOriginal_),
    pdfPieces_(),
    pdf_(utils::fullCloneFunc(pdfOriginal_, pdfPieces_)),
    lastData_(0),
    cache_(*pdf_,*obs_)
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
        if (data.weight() > 0) {
            nonZeroWEntries_++;
            nonZeroW_[i] = (data.weight() > 0);
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
    vpdf_ = new VPdfT(static_cast<const PdfT &>(*pdf_), data);
}

template <typename PdfT, typename VPdfT>
void
cacheutils::OptimizedCachingPdfT<PdfT,VPdfT>::realFill_(const RooAbsData &data, std::vector<Double_t> &vals) 
{
    vpdf_->fill(vals);
}


cacheutils::ReminderSum::ReminderSum(const char *name, const char *title, const RooArgList& sumSet) :
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

cacheutils::CachingAddNLL::CachingAddNLL(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data) :
    RooAbsReal(name, title),
    pdf_(pdf),
    params_("params","parameters",this),
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
    bool histNll  = runtimedef::get("ADDNLL_HISTNLL");
    bool gaussNll  = runtimedef::get("ADDNLL_GAUSSNLL");
    int npdf = addpdf->pdfList().getSize();
    //std::cout << "Unpacking RooAddPdf " << addpdf->GetName() << " with " << npdf << " components:" << std::endl;
    RooAbsReal *lastcoeff = 0;
    if (npdf == addpdf->coefList().getSize()) {
        lastcoeff =  dynamic_cast<RooAbsReal*>(addpdf->coefList().at(npdf-1));
    } else {
        prods_.push_back(new ReminderSum("","", addpdf->coefList()));
        lastcoeff = & prods_.back(); 
    }
    for (int i = 0; i < npdf; ++i) {
        RooAbsReal * coeff = (i < npdf-1 ? dynamic_cast<RooAbsReal*>(addpdf->coefList().at(i)) : lastcoeff);
        RooAbsPdf  * pdfi  = dynamic_cast<RooAbsPdf *>(addpdf->pdfList().at(i));
        if (recursive && typeid(*pdfi) == typeid(RooAddPdf)) {
            RooAddPdf *apdfi = static_cast<RooAddPdf*>(pdfi);
            RooArgList list(*coeff);
            if (basecoeffs.getSize()) list.add(basecoeffs);
            //std::cout << "    Invoking recursive unpack on " << i << ": RooAddPdf " << apdfi->GetName() << std::endl;
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
        }
        const RooArgSet *obs = data_->get();
        //std::cout << "    Adding " << i << ": " << pdfi->ClassName() << " " << pdfi->GetName() << std::endl;
        if (histNll && typeid(*pdfi) == typeid(FastVerticalInterpHistPdf)) {
            pdfs_.push_back(new CachingHistPdf(pdfi, obs));
        } else if (histNll && typeid(*pdfi) == typeid(FastVerticalInterpHistPdf2)) {
            pdfs_.push_back(new CachingHistPdf2(pdfi, obs));
        } else if (gaussNll && typeid(*pdfi) == typeid(RooGaussian)) {
            pdfs_.push_back(new CachingGaussPdf(pdfi, obs));
        } else {
            pdfs_.push_back(new CachingPdf(pdfi, obs));
        }
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
        isRooRealSum_ = false;
        addPdfs_(addpdf, runtimedef::get("ADDNLL_RECURSIVE"), RooArgList());
    } else if ((sumpdf = dynamic_cast<RooRealSumPdf *>(pdf_)) != 0) {
        const RooArgSet *obs = data_->get();
        isRooRealSum_ = true;
        int npdf = sumpdf->coefList().getSize();
        coeffs_.reserve(npdf);
        pdfs_.reserve(npdf);
        integrals_.reserve(npdf);
        for (int i = 0; i < npdf; ++i) {
            RooAbsReal * coeff = dynamic_cast<RooAbsReal*>(sumpdf->coefList().at(i));
            RooAbsReal * funci = dynamic_cast<RooAbsReal*>(sumpdf->funcList().at(i));
            coeffs_.push_back(coeff);
            pdfs_.push_back(new CachingPdf(funci, obs));
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

    // For multi pdf's need to reset the cache if index changed before evaluations
    multiPdfs_.clear();
    for (auto itp = pdfs_.begin(), edp = pdfs_.end(); itp != edp; ++itp) {
	bool isMultiPdf = itp->pdf()->IsA()->InheritsFrom(RooMultiPdf::Class());
	if (isMultiPdf) {
            const RooMultiPdf *mpdf = dynamic_cast<const RooMultiPdf*>((*itp).pdf());
            multiPdfs_.push_back(std::make_pair(mpdf, &*itp));
	}
    }
 
}

Double_t 
cacheutils::CachingAddNLL::evaluate() const 
{
#ifdef DEBUG_CACHE
    PerfCounter::add("CachingAddNLL::evaluate called");
#endif

    // For multi pdf's need to reset the cache if index changed before evaluations
    if (!multiPdfs_.empty()) {
        for (std::vector<std::pair<const RooMultiPdf*,CachingPdf*> >::iterator itp = multiPdfs_.begin(), edp = multiPdfs_.end(); itp != edp; ++itp) {
		bool hasChangedPdf = itp->first->checkIndexDirty();
		if (hasChangedPdf) itp->second->setDataDirty();
        }
    }

    std::fill( partialSum_.begin(), partialSum_.end(), 0.0 );

    std::vector<RooAbsReal*>::iterator  itc = coeffs_.begin(), edc = coeffs_.end();
    boost::ptr_vector<CachingPdf>::iterator   itp = pdfs_.begin();//,   edp = pdfs_.end();
    std::vector<Double_t>::const_iterator itw, bgw = weights_.begin();//,    edw = weights_.end();
    std::vector<Double_t>::iterator       its, bgs = partialSum_.begin(), eds = partialSum_.end();
    double sumCoeff = 0;
    //std::cout << "Performing evaluation of " << GetName() << std::endl;
    for ( ; itc != edc; ++itp, ++itc ) {
        // get coefficient
        Double_t coeff = (*itc)->getVal();
        if (isRooRealSum_) {
            sumCoeff += coeff * integrals_[itc - coeffs_.begin()]->getVal();
            //std::cout << "  coefficient = " << coeff << ", integral = " << integrals_[itc - coeffs_.begin()]->getVal() << std::endl;
        } else {
            sumCoeff += coeff;
        }
        // get vals
        const std::vector<Double_t> &pdfvals = itp->eval(*data_);
        // update running sum
        //    std::vector<Double_t>::const_iterator itv = pdfvals.begin();
        //    for (its = bgs; its != eds; ++its, ++itv) {
        //         *its += coeff * (*itv); // sum (n_i * pdf_i)
        //    }
        // vectorize to make it faster
        vectorized::mul_add(pdfvals.size(), coeff, &pdfvals[0], &partialSum_[0]);
    }
    // then get the final nll
    double ret = 0;
    for (its = bgs; its != eds ; ++its) {
        if (!isnormal(*its) || *its <= 0) {
            std::cerr << "WARNING: underflow to " << *its << " in " << GetName() << std::endl; 
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
    double expectedEvents = (isRooRealSum_ ? pdf_->getNorm(data_->get()) : sumCoeff);
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
        for (std::vector<std::pair<const RooMultiPdf*,CachingPdf*> >::iterator itp = multiPdfs_.begin(), edp = multiPdfs_.end(); itp != edp; ++itp) {
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
        if (w) weights_.push_back(w); 
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
            pdfs_[ib] = new CachingAddNLL(catClone->getLabel(), "", pdf, data);
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
        data.get(i); if (data.weight() == 0) continue;
        int ib = cat->getBin();
        //std::cout << "Event " << i << " of weight " << data.weight() << " is in bin " << ib << " label " << cat->getLabel() << std::endl;
        if (data.weight() > 0) datasets_[ib]->add(obs, data.weight());
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
