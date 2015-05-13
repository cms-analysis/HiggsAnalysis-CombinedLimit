#ifndef CachingMultiPdf_h
#define CachingMultiPdf_h

#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <vector>

namespace cacheutils {
    class CachingMultiPdf : public CachingPdfBase {
        public:
            CachingMultiPdf(const RooMultiPdf &pdf, const RooArgSet &obs) ;
            ~CachingMultiPdf() ;
            virtual const std::vector<Double_t> & eval(const RooAbsData &data) ;
            const RooAbsReal *pdf() const { return pdf_; }
            virtual void  setDataDirty() ;
            virtual void  setIncludeZeroWeights(bool includeZeroWeights) ;
        protected:
            const RooMultiPdf * pdf_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfs_;
    };

    class CachingAddPdf : public CachingPdfBase {
        public:
            CachingAddPdf(const RooAddPdf &pdf, const RooArgSet &obs) ;
            ~CachingAddPdf() ;
            virtual const std::vector<Double_t> & eval(const RooAbsData &data) ;
            const RooAbsReal *pdf() const { return pdf_; }
            virtual void  setDataDirty() ;
            virtual void  setIncludeZeroWeights(bool includeZeroWeights) ;
        protected:
            const RooAddPdf * pdf_;
            std::vector<const RooAbsReal *> coeffs_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfs_;
            std::vector<Double_t> work_;
    };

} // namespace

#endif
