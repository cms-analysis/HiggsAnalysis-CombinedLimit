#ifndef VectorizedGaussian_h
#define VectorizedGaussian_h

#include <RooGaussian.h>
#include <RooAbsData.h>
#include <vector>

class VectorizedGaussian {
    class Worker : public RooGaussian {
        public:
            Worker(const RooGaussian &g) : RooGaussian(g, "") {}
            const RooAbsReal & xvar()    const { return x.arg(); }
            const RooAbsReal & meanvar() const { return mean.arg(); }
            const RooAbsReal & sigvar()  const { return sigma.arg(); }
    };
    public:
        VectorizedGaussian(const RooGaussian &gaus, const RooAbsData &data, bool includeZeroWeights=false) ;
        void fill(std::vector<Double_t> &out) const ;
    private:
        const RooRealVar * x_;
        const RooAbsReal * mean_;
        const RooAbsReal * sigma_;
        std::vector<Double_t> xvals_;
        mutable std::vector<Double_t> work_, work2_;
};

#endif
