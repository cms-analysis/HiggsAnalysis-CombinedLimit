#ifndef HiggsAnalysis_CombinedLimit_SimpleTaylorExpansion1D_h
#define HiggsAnalysis_CombinedLimit_SimpleTaylorExpansion1D_h

#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooRealProxy.h>


//_________________________________________________
/*
BEGIN_HTML
<p>
SimpleTaylorExpansion1D is a class encoding a fixed simple taylor expansion of a function around a given point and with a fixed step size</p>
END_HTML
*/
//
class SimpleTaylorExpansion1D : public RooAbsReal {

   public:
      enum { MaxOrder = 4 };

      SimpleTaylorExpansion1D() {}
      SimpleTaylorExpansion1D(const char *name, const char *title, RooAbsReal & func, RooRealVar & xvar, double dx, int order=2) ;
      SimpleTaylorExpansion1D(const SimpleTaylorExpansion1D &other, const char *newname=0) ;

      ~SimpleTaylorExpansion1D() ;

      TObject * clone(const char *newname) const ;

    protected:
        Double_t evaluate() const;

    private:

        RooRealProxy x_;
        double x0_, ci_[MaxOrder+1];


  ClassDef(SimpleTaylorExpansion1D,1) // Simple Taylor Expansion in 1D 
};

#endif
