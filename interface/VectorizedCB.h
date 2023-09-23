#ifndef VectorizedCBShape_h
#define VectorizedCBShape_h

#include <RooCBShape.h>
#include <RooAbsData.h>
#include <vector>
#include <cmath>

class VectorizedCBShape {
    class Worker : public RooCBShape {
        public:
            Worker(const RooCBShape &g) : RooCBShape(g, "") {}
            const RooAbsReal & alphavar() const { return alpha.arg(); }
            const RooAbsReal & mvar()     const { return m.arg(); }
            const RooAbsReal & m0var()    const { return m0.arg(); }
            const RooAbsReal & nvar()     const { return n.arg(); }
            const RooAbsReal & sigmavar() const { return sigma.arg(); }
    };
    public:
        VectorizedCBShape(const RooCBShape &gaus, const RooAbsData &data, bool includeZeroWeights=false) ;
        void fill(std::vector<Double_t> &out) const ;
        double getIntegral() const ;

    private:
        const RooRealVar * x_;
        const RooAbsReal * m0_, * n_, *alpha_, *sigma_;
        std::vector<Double_t> xvals_;
        mutable std::vector<Double_t> work1_, work2_;
        std::vector<int> xindex_;
        std::vector<Double_t> xsorted_;
        enum WorkingMode { 
            Plain=1, FastPlain=2, Sorted=11, FastSorted=12, PreSort=21, FastPreSort=22
        } mode_;
        inline bool hasPreSort() const { return mode_ == PreSort || mode_ == FastPreSort; }
        inline bool hasSort() const { return mode_ == Sorted || mode_ == FastSorted || mode_ == PreSort || mode_ == FastPreSort; }
        inline bool hasFast() const { return mode_ == FastPlain || mode_ == FastSorted || mode_ == FastPreSort; }
        void cbGauss(double* __restrict__ t, unsigned int n, double norm, double* __restrict__ out,  double* __restrict__ work2) const ;
        void cbCB(double* __restrict__ t, unsigned int n, double norm, double* __restrict__ out,  double* __restrict__ work2) const ;
};

#endif
