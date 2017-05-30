//---------------------------------------------------------------------------
#ifndef ROOPARAMETRICSHAPEBINPDF
#define ROOPARAMETRICSHAPEBINPDF
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include <TF1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/WrappedTF1.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooParametricShapeBinPdf : public RooAbsPdf
{
public:
   RooParametricShapeBinPdf() {} ;
   RooParametricShapeBinPdf(const char *name, const char *title,  const char *formula, 
		  RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape );
   RooParametricShapeBinPdf(const RooParametricShapeBinPdf& other,
      const char* name = 0);
   void setTH1Binning(const TH1& _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooParametricShapeBinPdf(*this,newname); }
   inline virtual ~RooParametricShapeBinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy x;        // dependent variable
   RooListProxy pars;
   TF1 * myfunc;
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Int_t nPars;
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooParametricShapeBinPdf,1) // RooParametricShapeBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

