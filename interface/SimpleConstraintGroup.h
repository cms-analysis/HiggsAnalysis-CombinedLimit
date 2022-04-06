#ifndef SimpleGaussianConstraintGroup_h
#define SimpleGaussianConstraintGroup_h

#include "HiggsAnalysis/CombinedLimit/interface/SimpleGaussianConstraint.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimplePoissonConstraint.h"
#include "RooSetProxy.h"

class SimpleConstraintGroup : public RooAbsReal {
    public:
        SimpleConstraintGroup() ;
        SimpleConstraintGroup(const char *name, const char *title) ;
        SimpleConstraintGroup(const SimpleConstraintGroup & other, const char *newname = 0) ;
        void add(const SimpleGaussianConstraint * gaus) ;
        void add(const SimplePoissonConstraint * pois) ;
        void setZeroPoint() ;
        void clearZeroPoint() ;
        void updateZeroPoint() { clearZeroPoint(); setZeroPoint(); }
 
        TObject * clone(const char *newname) const { return new SimpleConstraintGroup(*this, newname); }

        unsigned int size() const { return _gaus.size() + _pois.size(); }
    protected:
        Double_t evaluate() const;

    private:
        RooSetProxy _deps;
        std::vector<const SimpleGaussianConstraint *> _gaus;
        std::vector<const SimplePoissonConstraint *> _pois;
        std::vector<double> _gaus0, _pois0;

        ClassDef(SimpleConstraintGroup,1) // group of constraints
};

#endif
