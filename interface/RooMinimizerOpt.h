#ifndef HiggsAnalysis_CombinedLimit_RooMinimizerOpt
#define HiggsAnalysis_CombinedLimit_RooMinimizerOpt

#if defined(ROO_MINIMIZER) || defined(ROO_MINIMIZER_FCN)
   #error "You cannot include RooMinimizer.h or RooMinimizerFcn.h before RooMinimizerOpt.h"
#else
   #define private protected
   #include <RooMinimizer.h>
   #undef protected
#endif

class RooMinimizerOpt : public RooMinimizer {
    public:
        RooMinimizerOpt(RooAbsReal& function) ;
        Double_t edm();
        Int_t minimize(const char* type, const char* alg=0) ;
        Int_t improve() ;
        Int_t migrad() ;
        Int_t hesse() ;
        Int_t minos() ;
        Int_t minos(const RooArgSet& minosParamList) ;

};

class RooMinimizerFcnOpt : public RooMinimizerFcn {
    public: 
        RooMinimizerFcnOpt(RooAbsReal *funct, RooMinimizer *context,  bool verbose = false);
        RooMinimizerFcnOpt(const RooMinimizerFcnOpt &other) ;
        virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const;
        Bool_t Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, Bool_t optConst, Bool_t verbose);
        void initStdVects() const ;
    protected:
        virtual double DoEval(const double * x) const;
        mutable std::vector<RooRealVar *> _vars;
        mutable std::vector<double>       _vals;
        mutable std::vector<bool  >       _hasOptimzedBounds;
        struct OptBound {
            double softMin, softMax, hardMin, hardMax;
            double transform(double x) const {
                if (x < softMin) {
                    double dx = softMin-x; // > 0
                    return hardMin + ( (softMin-hardMin) + dx ) * std::exp ( -2*dx/(softMin-hardMin) );
                } else if (x > softMax) { 
                    double dx = x-softMax;
                    return hardMax - ( (hardMax-softMax) + dx ) * std::exp ( -2*dx/(hardMax-softMax) );
                } else {
                    return x;
                }
            }
        };
        mutable std::vector<OptBound> _optimzedBounds;
};

#endif
