#ifndef SimplePoissonConstraint_h
#define SimplePoissonConstraint_h

#include <cmath>
#include <RooPoisson.h>

class SimplePoissonConstraint : public RooPoisson {
    public:
        SimplePoissonConstraint() {} ;
        SimplePoissonConstraint(const char *name, const char *title,
                RooAbsReal& _x, RooAbsReal& _mean, bool noRounding = false ):
            RooPoisson(name,title,_x,_mean, noRounding) { init(); }
        SimplePoissonConstraint(const SimplePoissonConstraint& other, const char* name=0) :
            RooPoisson(other, name) { init(); }
        SimplePoissonConstraint(const RooPoisson &g) : RooPoisson(g, "") { init(); }

        virtual TObject* clone(const char* newname) const { return new SimplePoissonConstraint(*this,newname); }
        inline virtual ~SimplePoissonConstraint() { }

        const RooAbsReal & getMean() const { return mean.arg(); }

        double getLogValFast() const { 
            if (_valueDirty) {
                Double_t expected = mean;
                Double_t observed = x;
                if (std::abs(observed)<1e-10) {
                    _value = (std::abs(expected)<1e-10) ? 0 : -1*expected;
                } else {
                    if(observed<1000000) {
                        _value = - ( - observed * log(expected) + expected + logGamma_ );
                    } else {
                        //if many observed events, use Gauss approximation
                        Double_t sigma_square = expected;
                        Double_t diff = observed - expected;
                        _value = log(sigma_square)/2 - (diff*diff)/(2*sigma_square);
                    }
                }
                _valueDirty = false;
            }
            return _value;
        }

        static RooPoisson * make(RooPoisson &c) ;
    private:
        double logGamma_;
        void init() ;

        ClassDef(SimplePoissonConstraint,1) // Poisson PDF with fast log
};

#endif
