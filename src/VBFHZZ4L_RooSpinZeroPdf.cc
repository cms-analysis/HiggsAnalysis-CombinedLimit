/*
Unlike for HZZ4l_RooSpinZeroPdf, T1...T5 have to be normalized
so that a1^4*T1 + a1^3*ai*T2 + a1^2*ai^2+T3 + a1*ai^3*T4 + ai^4*T5 = total pdf.
This means that for example for L1, T5's normalization should be tiny, because then it gets
multiplied by a very large g1_prime2^4.
*/

#include "Riostream.h"
#include <HiggsAnalysis/CombinedLimit/interface/VBFHZZ4L_RooSpinZeroPdf.h>
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;

ClassImp(VBFHZZ4L_RooSpinZeroPdf)

  VBFHZZ4L_RooSpinZeroPdf::VBFHZZ4L_RooSpinZeroPdf(const char *name, const char *title,
                                                   RooAbsReal& _kd,
                                                   RooAbsReal& _kdint,
                                                   RooAbsReal& _ksmd,
                                                   RooAbsReal& _a1,
                                                   RooAbsReal& _ai,
                                                   const RooArgList& inCoefList):
   RooAbsPdf(name,title),
   kd("kd","kd",this,_kd),
   kdint("kdint","kdint",this,_kdint),
   ksmd("ksmd","ksmd",this,_ksmd),
   a1("a1","a1",this,_a1),
   ai("ai","ai",this,_ai),
  _coefList("coefList","List of funcficients",this)
 {
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :VBFHZZ4L_RooSpinZeroPdf(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;

  Integral_T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))-> analyticalIntegral(1000);
  Integral_T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))-> analyticalIntegral(1000);
  Integral_T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))-> analyticalIntegral(1000);
  Integral_T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(3))-> analyticalIntegral(1000);
  Integral_T5 = dynamic_cast<const RooHistFunc*>(_coefList.at(4))-> analyticalIntegral(1000);
// _coefIter = _coefList.createIterator() ;
 }


 VBFHZZ4L_RooSpinZeroPdf::VBFHZZ4L_RooSpinZeroPdf(const VBFHZZ4L_RooSpinZeroPdf& other, const char* name) :
   RooAbsPdf(other,name),
   kd("kd",this,other.kd),
   kdint("kdint",this,other.kdint),
   ksmd("ksmd",this,other.ksmd),
   a1("a1",this,other.a1),
   ai("ai",this,other.ai),
  _coefList("coefList",this,other._coefList)
 {
         Integral_T1 = other.Integral_T1;
         Integral_T2 = other.Integral_T2;
         Integral_T3 = other.Integral_T3;
         Integral_T4 = other.Integral_T4;
         Integral_T5 = other.Integral_T5;
 }


 Double_t VBFHZZ4L_RooSpinZeroPdf::evaluate() const
 {
   double value = 0.;


   Double_t T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->getVal();
   Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();
   Double_t T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(3))->getVal();
   Double_t T5 = dynamic_cast<const RooHistFunc*>(_coefList.at(4))->getVal();

   value = (
              a1*a1*a1*a1               * T1
            + a1*a1*a1              *ai * T2
            + a1*a1              *ai*ai * T3
            + a1              *ai*ai*ai * T4
            +               ai*ai*ai*ai * T5
           );

   if ( value <= 0.) return 1.0e-200;

   return value ;

 }

Int_t VBFHZZ4L_RooSpinZeroPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg(), *kdint.absArg(), *ksmd.absArg()))) return 4 ;

  return 0 ;

}

Double_t VBFHZZ4L_RooSpinZeroPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

     case 4:
       {
         double Int_T1  = Integral_T1;
         double Int_T2  = Integral_T2;
         double Int_T3  = Integral_T3;
         double Int_T4  = Integral_T4;
         double Int_T5  = Integral_T5;

         double integral = (
                            a1*a1*a1*a1               * Int_T1
                          + a1*a1*a1              *ai * Int_T2
                          + a1*a1              *ai*ai * Int_T3
                          + a1              *ai*ai*ai * Int_T4
                          +               ai*ai*ai*ai * Int_T5
                         );

         if (integral <= 0.) return 1.0e-200;
         return integral;
       }

     }

   assert(0) ;
   return 0 ;
}

