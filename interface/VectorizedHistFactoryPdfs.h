#ifndef VectorizedHistFactoryPdfs_h
#define VectorizedHistFactoryPdfs_h

#include <RooHistFunc.h>
#include <RooStats/HistFactory/ParamHistFunc.h>
#include <RooStats/HistFactory/PiecewiseInterpolation.h>
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"

namespace cacheutils {
    class VectorizedHistFunc : public CachingPdfBase {
        public:
            VectorizedHistFunc(const RooHistFunc &pdf, bool includeZeroWeights=false) ;
            virtual ~VectorizedHistFunc() {}
            virtual const std::vector<Double_t> & eval(const RooAbsData &data) ;
            virtual const RooAbsReal *pdf() const { return pdf_; };
            virtual void  setDataDirty() { data_ = 0; }
            virtual void  setIncludeZeroWeights(bool includeZeroWeights) { includeZeroWeights_ = includeZeroWeights; }
        private:
            const RooHistFunc * pdf_;
            const RooAbsData * data_;
            bool includeZeroWeights_;
            std::vector<Double_t> yvals_;
            void fill() ;
    };

    class VectorizedParamHistFunc {
        public:
            VectorizedParamHistFunc(const ParamHistFunc &pdf, const RooAbsData &data, bool includeZeroWeights=false) ;
            void fill(std::vector<Double_t> &out) const ;
        private:
            std::vector<const RooRealVar *> yvars_;
    };

    class CachingPiecewiseInterpolation : public CachingPdfBase {
        public:
            CachingPiecewiseInterpolation(const PiecewiseInterpolation &pdf, const RooArgSet &obs) ;
            ~CachingPiecewiseInterpolation() ;
            virtual const std::vector<Double_t> & eval(const RooAbsData &data) ;
            const RooAbsReal *pdf() const { return pdf_; }
            virtual void  setDataDirty() ;
            virtual void  setIncludeZeroWeights(bool includeZeroWeights) ;
        protected:
            const PiecewiseInterpolation * pdf_;
            std::vector<const RooAbsReal *> coeffs_;
            std::vector<int>                codes_;
            bool positiveDefinite_;
            std::unique_ptr<CachingPdfBase>    cachingPdfNominal_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfsHi_;
            boost::ptr_vector<CachingPdfBase>  cachingPdfsLow_;
            std::vector<Double_t> work_;
    };
}


#endif
