//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooDijetBinPdf_h
#define HiggsAnalysis_CombinedLimit_RooDijetBinPdf_h
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
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooDijetBinPdf : public RooAbsPdf
{
public:
   RooDijetBinPdf() {} ;
   RooDijetBinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		  RooAbsReal& _p2, RooAbsReal& _p3,
		  RooAbsReal& _sqrts, RooAbsReal& _meff, RooAbsReal& _seff);
   RooDijetBinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		  RooAbsReal& _p2, RooAbsReal& _p3,
		  RooAbsReal& _sqrts);
   RooDijetBinPdf(const RooDijetBinPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooDijetBinPdf(*this,newname); }
   inline virtual ~RooDijetBinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy sqrts;        // sqrts
   RooRealProxy meff;        // meff
   RooRealProxy seff;        // seff
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooDijetBinPdf,1) // RazorDijetBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class DijetFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = pow(1-x/p[0],p[1])/pow(x/p[0],p[2]+p[3]*log(x/p[0]));
     double eff = 1.;
     if (p[4]>0 && p[5]>0) eff = 0.5 * (1.0 + TMath::Erf((x - p[4])/p[5])) ;
     return pdf*eff;
   }
   
   double DoEval(double x) const
   {
     double pdf = pow(1-x/pars[0],pars[1])/pow(x/pars[0],pars[2]+pars[3]*log(x/pars[0]));
     double eff = 1.;
     if (pars[4]>0 && pars[5]>0) eff = 0.5 * (1.0 + TMath::Erf((x - pars[4])/pars[5]));
     return pdf*eff;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new DijetFunction();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 6;
   }
};
