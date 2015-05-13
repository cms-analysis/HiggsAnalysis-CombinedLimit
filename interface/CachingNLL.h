#ifndef HiggsAnalysis_CombinedLimit_CachingNLL_h
#define HiggsAnalysis_CombinedLimit_CachingNLL_h

#include <memory>
#include <map>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooRealSumPdf.h>
#include <RooProdPdf.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooSetProxy.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooGaussian.h>
#include <RooProduct.h>
#include "HiggsAnalysis/CombinedLimit/interface/SimpleGaussianConstraint.h"
#include <boost/ptr_container/ptr_vector.hpp>

class RooMultiPdf;

// Part zero: ArgSet checker
namespace cacheutils {
    class ArgSetChecker {
        public:
            ArgSetChecker() {}
            ArgSetChecker(const RooAbsCollection &set) ;
            bool changed(bool updateIfChanged=false) ;
        private:
            std::vector<RooRealVar *> vars_;
            std::vector<double> vals_;
            std::vector<RooCategory *> cats_;
            std::vector<int> states_;
    };

// Part zero point five: Cache of pdf values for different parameters
    class ValuesCache {
        public:
            ValuesCache(const RooAbsReal &pdf, const RooArgSet &obs, int size=MaxItems_);
            ValuesCache(const RooAbsCollection &params, int size=MaxItems_);
            ~ValuesCache();
            // search for the item corresponding to the current values of the parameters.
            // if available, return (&values, true)
            // if not available, return (&room, false)
            // and it will be up to the caller code to fill the room the new item
            std::pair<std::vector<Double_t> *, bool> get(); 
            void clear();
        private:
            struct Item {
                Item(const RooAbsCollection &set)   : checker(set),   good(false) {}
                Item(const ArgSetChecker    &check) : checker(check), good(false) {}
                std::vector<Double_t> values;
                ArgSetChecker         checker;
                bool                  good;
            };
            int size_, maxSize_;
            enum { MaxItems_ = 3 };
            Item *items[MaxItems_];
    };
// Part one: cache all values of a pdf
class CachingPdfBase {
    public:
        CachingPdfBase() {}
        virtual ~CachingPdfBase() {}
        virtual const std::vector<Double_t> & eval(const RooAbsData &data) = 0;
        virtual const RooAbsReal *pdf() const = 0;
        virtual void  setDataDirty() = 0;
        virtual void  setIncludeZeroWeights(bool includeZeroWeights) = 0;
};
class CachingPdf : public CachingPdfBase {
    public:
        CachingPdf(RooAbsReal *pdf, const RooArgSet *obs) ;
        CachingPdf(const CachingPdf &other) ;
        virtual ~CachingPdf() ;
        virtual const std::vector<Double_t> & eval(const RooAbsData &data) ;
        const RooAbsReal *pdf() const { return pdf_; }
        virtual void  setDataDirty() { lastData_ = 0; }
        virtual void  setIncludeZeroWeights(bool includeZeroWeights) { includeZeroWeights_ = includeZeroWeights;  setDataDirty(); }
    protected:
        const RooArgSet *obs_;
        RooAbsReal *pdfOriginal_;
        RooArgSet  pdfPieces_;
        RooAbsReal *pdf_;
        const RooAbsData *lastData_;
        ValuesCache cache_;
        std::vector<uint8_t> nonZeroW_;
        unsigned int         nonZeroWEntries_;
        bool                 includeZeroWeights_;
        virtual void newData_(const RooAbsData &data) ;
        virtual void realFill_(const RooAbsData &data, std::vector<Double_t> &values) ;
};

template <typename PdfT, typename VPdfT> 
class OptimizedCachingPdfT : public CachingPdf {
    public:
        OptimizedCachingPdfT(RooAbsReal *pdf, const RooArgSet *obs) :
            CachingPdf(pdf,obs), vpdf_(0) {}
        OptimizedCachingPdfT(const OptimizedCachingPdfT &other) : 
            CachingPdf(other), vpdf_(0) {}
        virtual ~OptimizedCachingPdfT() { delete vpdf_; }
    protected:
        virtual void realFill_(const RooAbsData &data, std::vector<Double_t> &values) ;
        virtual void newData_(const RooAbsData &data) ;
        VPdfT *vpdf_;
};

CachingPdfBase * makeCachingPdf(RooAbsReal *pdf, const RooArgSet *obs) ;

class CachingAddNLL : public RooAbsReal {
    public:
        CachingAddNLL(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data, bool includeZeroWeights = false) ;
        CachingAddNLL(const CachingAddNLL &other, const char *name = 0) ;
        virtual ~CachingAddNLL() ;
        virtual CachingAddNLL *clone(const char *name = 0) const ;
        virtual Double_t evaluate() const ;
        virtual Bool_t isDerived() const { return kTRUE; }
        virtual Double_t defaultErrorLevel() const { return 0.5; }
        void setData(const RooAbsData &data) ;
        virtual RooArgSet* getObservables(const RooArgSet* depList, Bool_t valueOnly = kTRUE) const ;
        virtual RooArgSet* getParameters(const RooArgSet* depList, Bool_t stripDisconnected = kTRUE) const ;
        double  sumWeights() const { return sumWeights_; }
        const RooAbsPdf *pdf() const { return pdf_; }
        void setZeroPoint() { zeroPoint_ = -this->getVal(); setValueDirty(); }
        void clearZeroPoint() { zeroPoint_ = 0.0; setValueDirty();  }
        /// note: setIncludeZeroWeights(true) won't have effect unless you also re-call setData
        virtual void  setIncludeZeroWeights(bool includeZeroWeights) ;
        RooSetProxy & params() { return params_; }
    private:
        void setup_();
        void addPdfs_(RooAddPdf *addpdf, bool recursive, const RooArgList & basecoeffs) ;
        RooAbsPdf *pdf_;
        RooSetProxy params_;
        const RooAbsData *data_;
        std::vector<Double_t>  weights_, binWidths_;
        double               sumWeights_;
        bool includeZeroWeights_;
        mutable std::vector<RooAbsReal*> coeffs_;
        mutable boost::ptr_vector<CachingPdfBase>  pdfs_;
        mutable boost::ptr_vector<RooAbsReal>  prods_;
        mutable std::vector<RooAbsReal*> integrals_;
        mutable std::vector<std::pair<const RooMultiPdf*,CachingPdfBase*> > multiPdfs_;
        mutable std::vector<Double_t> partialSum_;
        mutable std::vector<Double_t> workingArea_;
        mutable bool isRooRealSum_, fastExit_;
        mutable int canBasicIntegrals_, basicIntegrals_;
        double zeroPoint_;
};

class CachingSimNLL  : public RooAbsReal {
    public:
        CachingSimNLL(RooSimultaneous *pdf, RooAbsData *data, const RooArgSet *nuis=0) ;
        CachingSimNLL(const CachingSimNLL &other, const char *name = 0) ;
        ~CachingSimNLL() ;
        virtual CachingSimNLL *clone(const char *name = 0) const ;
        virtual Double_t evaluate() const ;
        virtual Bool_t isDerived() const { return kTRUE; }
        virtual Double_t defaultErrorLevel() const { return 0.5; }
        void setData(const RooAbsData &data) ;
        virtual RooArgSet* getObservables(const RooArgSet* depList, Bool_t valueOnly = kTRUE) const ;
        virtual RooArgSet* getParameters(const RooArgSet* depList, Bool_t stripDisconnected = kTRUE) const ;
        void splitWithWeights(const RooAbsData &data, const RooAbsCategory& splitCat, Bool_t createEmptyDataSets) ;
        static void setNoDeepLogEvalError(bool noDeep) { noDeepLEE_ = noDeep; }
        void setZeroPoint() ; 
        void clearZeroPoint() ;
        static void forceUnoptimizedConstraints() { optimizeContraints_ = false; }
        friend class CachingAddNLL;
    private:
        void setup_();
        RooSimultaneous   *pdfOriginal_;
        const RooAbsData  *dataOriginal_;
        const RooArgSet   *nuis_;
        RooSetProxy        params_;
        RooArgSet piecesForCloning_;
        std::auto_ptr<RooSimultaneous>  factorizedPdf_;
        std::vector<RooAbsPdf *>        constrainPdfs_;
        std::vector<SimpleGaussianConstraint *>  constrainPdfsFast_;
        std::vector<bool>                        constrainPdfsFastOwned_;
        std::vector<CachingAddNLL*>     pdfs_;
        std::auto_ptr<TList>            dataSets_;
        std::vector<RooDataSet *>       datasets_;
        static bool noDeepLEE_;
        static bool hasError_;
        static bool optimizeContraints_;
        std::vector<double> constrainZeroPoints_;
        std::vector<double> constrainZeroPointsFast_;
};

}
#endif
