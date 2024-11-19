#ifndef SimpleGaussianConstraint_h
#define SimpleGaussianConstraint_h

#include <RooGaussian.h>

class SimpleGaussianConstraint : public RooGaussian {
    public:
        SimpleGaussianConstraint() {} ;
        SimpleGaussianConstraint(const char *name, const char *title,
                RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma):
            RooGaussian(name,title,_x,_mean,_sigma) { init(); }
        SimpleGaussianConstraint(const SimpleGaussianConstraint& other, const char* name=0) :
            RooGaussian(other, name) { init(); }
        SimpleGaussianConstraint(const RooGaussian &g) : RooGaussian(g, "") { init(); }

        TObject* clone(const char* newname) const override { return new SimpleGaussianConstraint(*this,newname); }

#if ROOT_VERSION_CODE < ROOT_VERSION(6,26,0)
        // function was upstreamed to RooGaussian in ROOT 6.26
        const RooAbsReal & getX() const { return x.arg(); }
#endif

        double getLogValFast() const { 
            if (_valueDirty) {
                Double_t arg = x - mean;  
                //Double_t sig = sigma ;
                //return -0.5*arg*arg/(sig*sig);
                _value = scale_*arg*arg;
                _valueDirty = false;
            }
            return _value;
        }

        // RooFit should make no attempt to normalize this constraint, as the
        // "getLogValFast()" function that combined CachingNLL is calling also
        // doesn't do any normalization.
        bool selfNormalized() const override { return true; }

        static RooGaussian * make(RooGaussian &c) ;

    private:
        double scale_;
        void init() ;

        ClassDefOverride(SimpleGaussianConstraint,1) // Gaussian PDF with fast log
};

#endif
