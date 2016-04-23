#ifndef VectorizedExponential_h
#define VectorizedExponential_h

#include <RooExponential.h>
#include <RooAbsData.h>
#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include <vector>

class VectorizedExponential {
    public:
        VectorizedExponential(const RooExponential &pdf, const RooAbsData &data, bool includeZeroWeights=false) ;
        void fill(std::vector<Double_t> &out) const ;
    private:
        const RooRealVar * x_;
        const RooAbsReal * lambda_;
        std::vector<Double_t> xvals_;
        mutable std::vector<Double_t> work_;
};

class VectorizedPower {
    public:
        VectorizedPower(const RooPower &pdf, const RooAbsData &data, bool includeZeroWeights=false) ;
        void fill(std::vector<Double_t> &out) const ;
    private:
        const RooRealVar * x_;
        const RooAbsReal * exponent_;
        std::vector<Double_t> xvals_;
        mutable std::vector<Double_t> work_;
};

#endif
