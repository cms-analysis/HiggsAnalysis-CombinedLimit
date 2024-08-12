#ifndef VectorizedHistFactoryPdfs_h
#define VectorizedHistFactoryPdfs_h

#include <RooHistFunc.h>
#include <RooStats/HistFactory/ParamHistFunc.h>
#include <RooStats/HistFactory/PiecewiseInterpolation.h>
#include "CachingNLL.h"

namespace cacheutils {
    class VectorizedHistFunc : public CachingPdfBase {
        public:
            VectorizedHistFunc(const RooHistFunc &pdf, bool includeZeroWeights=false) ;
            ~VectorizedHistFunc() override {}
            const std::vector<Double_t> & eval(const RooAbsData &data) override ;
            const RooAbsReal *pdf() const override { return pdf_; };
            void  setDataDirty() override { data_ = 0; }
            void  setIncludeZeroWeights(bool includeZeroWeights) override { includeZeroWeights_ = includeZeroWeights; }
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
            std::vector<RooAbsReal *> yvars_;
    };

    class CachingPiecewiseInterpolation : public CachingPdfBase {
        public:
            CachingPiecewiseInterpolation(const PiecewiseInterpolation &pdf, const RooArgSet &obs) ;
            ~CachingPiecewiseInterpolation() override ;
            const std::vector<Double_t> & eval(const RooAbsData &data) override ;
            const RooAbsReal *pdf() const override { return pdf_; }
            void  setDataDirty() override ;
            void  setIncludeZeroWeights(bool includeZeroWeights) override ;
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
