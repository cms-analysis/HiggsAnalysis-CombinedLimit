#include "Riostream.h" 
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_phase.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;

ClassImp(HZZ4L_RooSpinZeroPdf_phase) 

  HZZ4L_RooSpinZeroPdf_phase::HZZ4L_RooSpinZeroPdf_phase(const char *name, const char *title, 
					     RooAbsReal& _kd,
					     RooAbsReal& _kdint,
					     RooAbsReal& _ksmd,
					     RooAbsReal& _fai,
					     RooAbsReal& _phi,
					     const RooArgList& inCoefList): 
   RooAbsPdf(name,title), 
   kd("kd","kd",this,_kd),
   kdint("kdint","kdint",this,_kdint),
   ksmd("ksmd","ksmd",this,_ksmd),
   fai("fai","fai",this,_fai),
   phi("phi","phi",this,_phi),
  _coefList("coefList","List of funcficients",this) 
  
 { 
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :HZZ4L_RooSpinZeroPdf_phase(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;
  
// _coefIter = _coefList.createIterator() ;
 } 


 HZZ4L_RooSpinZeroPdf_phase::HZZ4L_RooSpinZeroPdf_phase(const HZZ4L_RooSpinZeroPdf_phase& other, const char* name) :  
   RooAbsPdf(other,name), 
   kd("kd",this,other.kd),
   kdint("kdint",this,other.kdint),
   ksmd("ksmd",this,other.ksmd),
   fai("fai",this,other.fai),
   phi("phi",this,other.phi),
  _coefList("coefList",this,other._coefList)

 { 
 // _coefIter = _coefList.createIterator() ;
 } 


 Double_t HZZ4L_RooSpinZeroPdf_phase::evaluate() const 
 { 
   double value = 0.;

	
   Double_t T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->getVal();
   Double_t T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();
   Double_t T5 = dynamic_cast<const RooHistFunc*>(_coefList.at(3))->getVal();

   
   value = (1.-fabs(fai)) * T1 + fabs(fai) * T2 + sqrt((1.-fabs(fai))*fabs(fai)) * (cos(phi)*T4 +sin(phi)*T5); 
   
   if ( value <= 0.) return 1.0e-200;
   
   return value ; 
   
 } 

Int_t HZZ4L_RooSpinZeroPdf_phase::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg(), *kdint.absArg(), *ksmd.absArg()))) return 4 ;
  //if (matchArgs(allVars,analVars,kd)) return 1 ;
  //if (matchArgs(allVars,analVars,kdint)) return 2 ;
  //if (matchArgs(allVars,analVars,ksmd)) return 3 ;

  return 0 ;

}

Double_t HZZ4L_RooSpinZeroPdf_phase::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

       // integrate out kd, depend on kdint
//     case 1: 
//       {
//
//	 int biny = histo0->GetYaxis()->FindBin(kdint);
//	 
//	 double Int_T1 = histo0->Integral(1, nbinsx, biny, biny, 1, nbinsz);
//	 double Int_T2 = histo1->Integral(1, nbinsx, biny, biny, 1, nbinsz);
//	 double Int_T4 = histo2->Integral(1, nbinsx, biny, biny, 1, nbinsz);
//	 // something related to phase factor, this is by guess
//
//	 double mysgn = 1.;
//	 if(fai < 0.) 
//	   {
//	     mysgn = -1.;
//	   }
//
//	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2  + mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; ;
//
//	 integral = integral * dx; 
//	  
//	 return integral; 
//	 
//       }
//       
//       // integrate out  kdint, depend on kd
//     case 2: 
//       {
//	 
//	 int binx = histo0->GetXaxis()->FindBin(kd);
//
//	 double Int_T1 = histo0->Integral(binx, binx, 1, nbinsy, 1, nbinsz);
//	 double Int_T2 = histo1->Integral(binx, binx, 1, nbinsy, 1, nbinsz);
//	 double Int_T4 = histo2->Integral(binx, binx, 1, nbinsy, 1, nbinsz);
//
//	 double mysgn = 1.;
//	 if(fai < 0.) 
//	   {
//	     mysgn = -1.;
//	   }
//
//	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2 +  mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; 
//	 
//	 // something related to phase factor, this is by guess
//	 integral = integral * dy;
//
//	 return integral;
//       }
//       
//     case 3: 
//       {
//	 
//	 int binz = histo0->GetZaxis()->FindBin(ksmd);
//
//	 double Int_T1 = histo0->Integral(1, nbinsx, 1, nbinsy, binz, binz);
//	 double Int_T2 = histo1->Integral(1, nbinsx, 1, nbinsy, binz, binz);
//	 double Int_T4 = histo2->Integral(1, nbinsx, 1, nbinsy, binz, binz);
//
//	 double mysgn = 1.;
//	 if(fai < 0.) 
//	   {
//	     mysgn = -1.;
//	   }
//
//	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2 +  mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; 
//	 
//	 // something related to phase factor, this is by guess
//	 integral = integral * dz;
//
//	 return integral;
//       }
     case 4: 
       {
 double Int_T1  = dynamic_cast<const RooHistFunc*>(_coefList.at(0))-> analyticalIntegral(1000);
 double Int_T2  = dynamic_cast<const RooHistFunc*>(_coefList.at(1))-> analyticalIntegral(1000);
 double Int_T4  = dynamic_cast<const RooHistFunc*>(_coefList.at(2))-> analyticalIntegral(1000);
 double Int_T5  = dynamic_cast<const RooHistFunc*>(_coefList.at(3))-> analyticalIntegral(1000);



	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2 + sqrt((1.-fabs(fai))*fabs(fai)) *( cos(phi)* Int_T4 + sin(phi) * Int_T5) ;

	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

