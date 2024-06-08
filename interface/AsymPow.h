#ifndef HiggsAnalysis_CombinedLimit_AsymPow_h
#define HiggsAnalysis_CombinedLimit_AsymPow_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>

#include "CombineCodegenImpl.h"

//_________________________________________________
/*
BEGIN_HTML
<p>
AsymPow is helper class for implementing asymmetric log-normal errors. 
It has two parameters <i>kappa<sub>Low</sub></i>, <i>kappa<sub>High</sub></i> and one variable (<i>theta</i>).
<ul>
<li>for <i>theta &gt; 0</i>, it evaluates to <b>pow</b>(<i>kappa<sub>High</sub></i>, <i>theta</i>). </li>
<li>for <i>theta &lt; 0</i>, it evaluates to <b>pow</b>(<i>kappa<sub>Low</sub></i>, &minus;<i>theta</i>). </li>
</ul>
</p>
END_HTML
*/
//
class AsymPow : public RooAbsReal {

   public:
      AsymPow() {}
      AsymPow(const char *name, const char *title, RooAbsReal &kappaLow, RooAbsReal &kappaHigh, RooAbsReal &theta) ;
      AsymPow(const AsymPow &other, const char *newname=0) ;

      TObject * clone(const char *newname) const override { return new AsymPow(*this, newname); }

      COMBINE_DECLARE_TRANSLATE;

      RooAbsReal const &kappaLow() const { return kappaLow_.arg(); }
      RooAbsReal const &kappaHigh() const { return kappaHigh_.arg(); }
      RooAbsReal const &theta() const { return theta_.arg(); }

    protected:
        Double_t evaluate() const override;

    private:
        RooRealProxy kappaLow_, kappaHigh_;
        RooRealProxy theta_;

  ClassDefOverride(AsymPow,1) // Asymmetric power	
};

COMBINE_DECLARE_CODEGEN_IMPL(AsymPow);

#endif
