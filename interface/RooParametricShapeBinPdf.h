//---------------------------------------------------------------------------
#ifndef ROOPARAMETRICSHAPEBINPDF
#define ROOPARAMETRICSHAPEBINPDF
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>

//---------------------------------------------------------------------------
class RooParametricShapeBinPdf : public RooAbsPdf
{
public:
   RooParametricShapeBinPdf() {} ;
   RooParametricShapeBinPdf(const char *name, const char *title,  RooAbsReal& _pdf, 
		  RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape );
   RooParametricShapeBinPdf(const RooParametricShapeBinPdf& other,
      const char* name = 0);
   void setTH1Binning(const TH1& _Hnominal);
   RooAbsPdf* getPdf() const;
   RooAbsReal* getIntegral(int index) const;
   virtual TObject* clone(const char* newname) const { return new RooParametricShapeBinPdf(*this,newname); }
   inline virtual ~RooParametricShapeBinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy x;        // dependent variable
   RooListProxy pars;
   RooRealProxy mypdf;   
   RooListProxy myintegrals;
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min

   Double_t evaluate() const;
private:
   virtual RooPlot* plotOn(RooPlot* frame, 
              const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
              const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
              const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
              const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none(),
              const RooCmdArg& arg9=RooCmdArg::none(), const RooCmdArg& arg10=RooCmdArg::none()
                 ) const {
     return mypdf.arg().plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10) ;
     }
     virtual RooPlot* plotOn(RooPlot* frame, RooLinkedList& cmdList) const {
       return mypdf.arg().plotOn(frame,cmdList) ;
     }
     
   ClassDef(RooParametricShapeBinPdf,2) // RooParametricShapeBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

