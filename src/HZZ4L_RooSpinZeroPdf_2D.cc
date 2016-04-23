#include "Riostream.h" 
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_2D.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"

using namespace TMath;

ClassImp(HZZ4L_RooSpinZeroPdf_2D) 

  HZZ4L_RooSpinZeroPdf_2D::HZZ4L_RooSpinZeroPdf_2D(const char *name, const char *title, 
					     RooAbsReal& _kd,
					     RooAbsReal& _kdint,
					     RooAbsReal& _ksmd,
					     RooAbsReal& _fai1,
					     RooAbsReal& _fai2,
					     RooAbsReal& _phi1,
					     RooAbsReal& _phi2,
					     const RooArgList& inCoefList): 
   RooAbsPdf(name,title), 
   kd("kd","kd",this,_kd),
   kdint("kdint","kdint",this,_kdint),
   ksmd("ksmd","ksmd",this,_ksmd),
   fai1("fai1","fai1",this,_fai1),
   fai2("fai2","fai2",this,_fai2),
   phi1("phi1","phi1",this,_phi1),
   phi2("phi2","phi2",this,_phi2),
  _coefList("coefList","List of funcficients",this) 
  
 { 
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :HZZ4L_RooSpinZeroPdf_2D(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;
  
// _coefIter = _coefList.createIterator() ;
 } 


 HZZ4L_RooSpinZeroPdf_2D::HZZ4L_RooSpinZeroPdf_2D(const HZZ4L_RooSpinZeroPdf_2D& other, const char* name) :  
   RooAbsPdf(other,name), 
   kd("kd",this,other.kd),
   kdint("kdint",this,other.kdint),
   ksmd("ksmd",this,other.ksmd),
   fai1("fai1",this,other.fai1),
   fai2("fai2",this,other.fai2),
   phi1("phi1",this,other.phi1),
   phi2("phi2",this,other.phi2),
  _coefList("coefList",this,other._coefList)

 { 
 // _coefIter = _coefList.createIterator() ;
 } 


 Double_t HZZ4L_RooSpinZeroPdf_2D::evaluate() const 
 { 
   double value = 0.;

// pure terms	
   Double_t T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->getVal();
   Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();

// interference cos term 
   Double_t T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(3))->getVal();
   Double_t T5 = dynamic_cast<const RooHistFunc*>(_coefList.at(4))->getVal();
   Double_t T6 = dynamic_cast<const RooHistFunc*>(_coefList.at(5))->getVal();

// interference sin term 
   Double_t T7 = dynamic_cast<const RooHistFunc*>(_coefList.at(6))->getVal();
   Double_t T8 = dynamic_cast<const RooHistFunc*>(_coefList.at(7))->getVal();
   Double_t T9 = dynamic_cast<const RooHistFunc*>(_coefList.at(8))->getVal();
   
   value = (1.-fai1 - fai2) * T1 + fai1 * T2 + fai2 * T3 + sqrt((1.-fai1- fai2)*fai1) * (cos(phi1)*T4 +sin(phi1)*T7) + sqrt((1.-fai1- fai2)*fai2) * (cos(phi2)*T5 +sin(phi2)*T8) + sqrt(fai1*fai2) * (cos(phi1-phi2)*T6 +sin(phi1-phi2)*T9); 
   
   if ( value <= 0.) return 1.0e-200;
   
   return value ; 
   
 } 

Int_t HZZ4L_RooSpinZeroPdf_2D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg(), *kdint.absArg(), *ksmd.absArg()))) return 4 ;
  //if (matchArgs(allVars,analVars,kd)) return 1 ;
  //if (matchArgs(allVars,analVars,kdint)) return 2 ;
  //if (matchArgs(allVars,analVars,ksmd)) return 3 ;

  return 0 ;

}

Double_t HZZ4L_RooSpinZeroPdf_2D::analyticalIntegral(Int_t code, const char* rangeName) const
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
 double Int_T3  = dynamic_cast<const RooHistFunc*>(_coefList.at(2))-> analyticalIntegral(1000);

 double Int_T4  = dynamic_cast<const RooHistFunc*>(_coefList.at(3))-> analyticalIntegral(1000);
 double Int_T5  = dynamic_cast<const RooHistFunc*>(_coefList.at(4))-> analyticalIntegral(1000);
 double Int_T6  = dynamic_cast<const RooHistFunc*>(_coefList.at(5))-> analyticalIntegral(1000);

 double Int_T7  = dynamic_cast<const RooHistFunc*>(_coefList.at(6))-> analyticalIntegral(1000);
 double Int_T8  = dynamic_cast<const RooHistFunc*>(_coefList.at(7))-> analyticalIntegral(1000);
 double Int_T9  = dynamic_cast<const RooHistFunc*>(_coefList.at(8))-> analyticalIntegral(1000);


   double integral = (1.-fai1 - fai2) * Int_T1 + fai1 * Int_T2 + fai2 * Int_T3 + sqrt((1.-fai1- fai2)*fai1) * (cos(phi1)*Int_T4 +sin(phi1)*Int_T7) + sqrt((1.-fai1- fai2)*fai2) * (cos(phi2)*Int_T5 +sin(phi2)*Int_T8) + sqrt(fai1*fai2) * (cos(phi1-phi2)*Int_T6 +sin(phi1-phi2)*Int_T9); 
	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

