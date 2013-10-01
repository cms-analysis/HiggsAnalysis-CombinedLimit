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

        virtual TObject* clone(const char* newname) const { return new SimpleGaussianConstraint(*this,newname); }
        inline virtual ~SimpleGaussianConstraint() { }

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

    private:
        double scale_;
        void init() ;

        ClassDef(SimpleGaussianConstraint,1) // Gaussian PDF with fast log
};

#endif
