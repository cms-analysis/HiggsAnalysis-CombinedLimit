#include "Riostream.h" 
#include "../interface/HZZ4L_RooHighmass.h"
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"

using namespace TMath;

ClassImp(HZZ4L_RooHighmass) 

  HZZ4L_RooHighmass::HZZ4L_RooHighmass(const char *name, const char *title, 
					     RooAbsReal& _mass,
					     RooAbsReal& _dbkg,
					     RooAbsReal& _coupl,
					     const RooArgList& inCoefList): 
   RooAbsPdf(name,title), 
   mass("mass","mass",this,_mass),
   dbkg("dbkg","dbkg",this,_dbkg),
   coupl("coupl","coupl",this,_coupl),
 // _normSet("normSet","List of observables",this), 
  _coefList("coefList","List of funcficients",this) 
  
 { 
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :HZZ4L_RooHighmass(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;
  
 } 


 HZZ4L_RooHighmass::HZZ4L_RooHighmass(const HZZ4L_RooHighmass& other, const char* name) :  
   RooAbsPdf(other,name), 
//	 _normSet(other._normSet),
   mass("mass",this,other.mass),
   dbkg("dbkg",this,other.dbkg),
   coupl("coupl",this,other.coupl),
  //_normSet("normSet",this,other._normSet),
  _coefList("coefList",this,other._coefList)

 { 
  _coefIter = _coefList.createIterator() ;
 } 


 Double_t HZZ4L_RooHighmass::evaluate() const 
 { 
   double value = 0.;

   Double_t T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->getVal();
//   RooArgSet *nset = new RooArgSet(mass.arg());
   Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();

   value = T2 + coupl*T1 + sqrt(coupl)*T3;
   
   if ( value <= 0.) { return 1.0e-200;}
   
   return value ; 
   
 } 

Int_t HZZ4L_RooHighmass::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*mass.absArg(),*dbkg.absArg()))) return 4 ;

  return 0 ;

}

Double_t HZZ4L_RooHighmass::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

     case 4: 
       {
 double Int_T1  = dynamic_cast<const RooHistFunc*>(_coefList.at(0))-> analyticalIntegral(1000);
 double Int_T2  = dynamic_cast<const RooHistFunc*>(_coefList.at(1))-> analyticalIntegral(1000);
 Double_t T3_int= dynamic_cast<const RooHistFunc*>(_coefList.at(2))->analyticalIntegral(1000);


	 double mysgn = 1.;
	 if(coupl < 0.) 
	   {
	     mysgn = -1.;
	   }

	 double integral =  coupl*Int_T1 + Int_T2 + mysgn*sqrt(fabs(coupl)) * T3_int; ;
	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

