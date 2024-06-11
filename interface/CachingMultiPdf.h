#ifndef CachingMultiPdf_h
#define CachingMultiPdf_h

#include "RooMultiPdf.h"
#include "CachingNLL.h"
#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <vector>

namespace cacheutils {
    class CachingMultiPdf : public CachingPdfBase {
        public:
            CachingMultiPdf(const RooMultiPdf &pdf, const RooArgSet &obs) ;
            ~CachingMultiPdf() override ;
            const std::vector<Double_t> & eval(const RooAbsData &data) override ;
            const RooAbsReal *pdf() const override { return pdf_; }
            void  setDataDirty() override ;
            void  setIncludeZeroWeights(bool includeZeroWeights) override ;
        protected:
            const RooMultiPdf * pdf_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfs_;
    };

    class CachingAddPdf : public CachingPdfBase {
        public:
            CachingAddPdf(const RooAddPdf &pdf, const RooArgSet &obs) ;
            ~CachingAddPdf() override ;
            const std::vector<Double_t> & eval(const RooAbsData &data) override ;
            const RooAbsReal *pdf() const override { return pdf_; }
            void  setDataDirty() override ;
            void  setIncludeZeroWeights(bool includeZeroWeights) override ;
        protected:
            const RooAddPdf * pdf_;
            std::vector<const RooAbsReal *> coeffs_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfs_;
            std::vector<Double_t> work_;
    };

    class CachingProduct : public CachingPdfBase {
        public:
            CachingProduct(const RooProduct &pdf, const RooArgSet &obs) ;
            ~CachingProduct() override ;
            const std::vector<Double_t> & eval(const RooAbsData &data) override ;
            const RooAbsReal *pdf() const override { return pdf_; }
            void  setDataDirty() override ;
            void  setIncludeZeroWeights(bool includeZeroWeights) override ;
        protected:
            const RooProduct * pdf_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfs_;
            std::vector<Double_t> work_;
    };



} // namespace

#endif
