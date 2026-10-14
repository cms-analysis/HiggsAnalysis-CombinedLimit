#include "../interface/HZZ4LRooPdfs.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"
#include "RooMath.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"
#include "TComplex.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>

using namespace TMath;

namespace RooFit{
	
	void readFile();
	
	const Double_t FracEventsNoBrem_4mu = 0.703618;
	const Double_t FracEventsNoBrem_4e = 0.583196;
	const Double_t FracEventsNoBrem_2e2mu = 0.641297;
	
	Double_t Lgg_7(Double_t mHstar);
	Double_t HiggsWidth(Int_t ID,Double_t mHrequested);
        Double_t HiggsWidthSM4(Int_t ID,Double_t mHrequested);

	Double_t pdf1(double mHstar,double mHreq);
        Double_t pdf1SM4(double mHstar,double mHreq);
        Double_t pdf1Param(double mHstar,double mHreq,double x,double widthSF=1.);
	Double_t rho(double r, TString proc);
	Double_t Sigma(double mHreq, TString proc);
	
	Double_t N(Double_t mH, TString proc);
	Double_t sigma_CB(Double_t mH, TString proc);
	Double_t mean(Double_t mH, TString proc);
	
	Double_t scratchMass;
	Double_t BR[26][217];
        Double_t BR_SM4[18][102];

        void readFile() {
          static constexpr std::array<std::array<Double_t, 217>, 26> sm = {{
#include "../data/hzz4l/HiggsWidthSM.inc"
          }};
          static constexpr std::array<std::array<Double_t, 102>, 18> sm4 = {{
#include "../data/hzz4l/HiggsWidthSM4.inc"
          }};

          std::size_t channel = 0;
          for (auto& destination : BR) {
            // The original initializer left SM channel 12 untouched.
            if (channel != 12) {
              const auto& row = sm.at(channel);
              std::copy_n(row.begin(), row.size(), std::begin(destination));
            }
            ++channel;
          }
          channel = 0;
          for (auto& destination : BR_SM4) {
            const auto& row = sm4.at(channel++);
            std::copy_n(row.begin(), row.size(), std::begin(destination));
          }
        }
	Double_t HiggsWidth(Int_t ID, Double_t mHreq){
		
		Double_t TotalWidth = 0;
		Double_t PartialWidth = 0;
		Double_t Width = 0;
		Int_t i = 0;
		Double_t closestMass = 0;
		Double_t tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
		Double_t deltaY2, tmpLow2, tmpHigh2, slope2, step;
		
		// If ID is unavailable return -1
		if(ID > 25 || ID < 0){return -1;}
		
		// If mH is out of range return -1
		// else find what array number to read
		if( mHreq < 90 || mHreq > 1000){return -1;}
		else{
			
			//Find index and closest higgs mass for which we have numbers 
			if(mHreq <= 110 ){step = 5; i = (Int_t)((mHreq - 90)/step); closestMass = (Int_t)(step*i + 90);}
			if(mHreq > 110 && mHreq <= 140 ){step = 0.5; i = (Int_t)(4 + (mHreq - 110)/step); closestMass = (step*(i-4) + 110);}
			if(mHreq > 140 && mHreq <= 160 ){step = 1; i = (Int_t)(64 + (mHreq - 140)/step); closestMass = (Int_t)(step*(i-64) + 140);}
			if(mHreq > 160 && mHreq <= 290 ){step = 2; i = (Int_t)(84 + (mHreq - 160)/step); closestMass = (Int_t)(step*(i-84) + 160);}
			if(mHreq > 290 && mHreq <= 350 ){step = 5; i = (Int_t)(149 + (mHreq - 290)/step); closestMass = (Int_t)(step*(i-149) + 290);}
			if(mHreq > 350 && mHreq <= 400 ){step = 10; i = (Int_t)(161 + (mHreq - 350)/step); closestMass = (Int_t)(step*(i-161) + 350);}
			if(mHreq > 400){step = 20; i = (Int_t)(166 + (mHreq - 400)/step); closestMass = (Int_t)(step*(i-166) + 400);}
			
			
			tmpLow1 = BR[ID][i]*BR[0][i];
			tmpHigh1 = BR[ID][i+1]*BR[0][i+1];
			
			
			tmpLow2 = BR[0][i];
			tmpHigh2 = BR[0][i+1];
			deltaX = mHreq - closestMass;
			deltaY1 = tmpHigh1 - tmpLow1;
			slope1 = deltaY1/step;
			
			
			deltaY2 = tmpHigh2 - tmpLow2;
			slope2 = deltaY2/step;
			
			// For partial widths                 
			if(deltaX == 0){ PartialWidth = tmpLow1;
			TotalWidth = tmpLow2;}
			else{ PartialWidth = slope1*deltaX + tmpLow1;
			TotalWidth = slope2*deltaX + tmpLow2;}
			
			// For total width
			if( ID == 0 ){ Width = TotalWidth; }
			else{ Width = PartialWidth;}
			
		}
		
		return Width;
		
	}
	




// HiggsWidth takes process ID and higgs mass mH
double HiggsWidthSM4(int ID, double mHrequested){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                   H->4e/4mu = 12               */
  /*                    H->2e2mu = 13               */
  /*             H->4l(e/mu/tau) = 14               */
  /*                       H->4q = 15               */
  /*                     H->2l2q = 16               */
  /*                       H->4f = 17               */
  /**************************************************/



  double TotalWidth = 0;
  double PartialWidth = 0;
  double Width = 0;
  int i = 0;
  double closestMass = 0;
  double tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
  double deltaY2, tmpLow2, tmpHigh2, slope2, step;


  // If ID is unavailable return -1                                           
  if(ID > 17 || ID < 0){return 0;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mHrequested < 100 || mHrequested > 1000){return 0;}
  else{

    //Find index and closest higgs mass for which we have numbers
    if(mHrequested <= 140){step = 10; i = (int)((mHrequested - 100)/step); closestMass = (int)(step*i + 100);}
    if(mHrequested > 140 && mHrequested <= 150 ){step = 5; i = (int)(4 + (mHrequested-140)/step); closestMass = (step*(i-4) + 140);}
    if(mHrequested > 150 && mHrequested <= 190 ){step = 1; i = (int)(6 + (mHrequested-150)/step); closestMass = (int)(step*(i-6) + 150);}
    if(mHrequested > 190 && mHrequested <= 200 ){step = 5; i = (int)(46 + (mHrequested-190)/step); closestMass = (int)(step*(i-46) + 190);}
    if(mHrequested > 200 && mHrequested <= 330 ){step = 10; i = (int)(48 + (mHrequested-200)/step); closestMass = (int)(step*(i-48) + 200);}
    if(mHrequested > 330 && mHrequested <= 335 ){step = 5; i = (int)(61 + (mHrequested-330)/step); closestMass = (int)(step*(i-61) + 330);}
    if(mHrequested > 335 && mHrequested <= 340 ){step = 1; i = (int)(62 + (mHrequested-335)/step); closestMass = (int)(step*(i-62) + 335);}
    if(mHrequested > 340 && mHrequested <= 345 ){step = 5; i = (int)(67 + (mHrequested-340)/step); closestMass = (int)(step*(i-67) + 340);}
    if(mHrequested > 345 && mHrequested <= 355 ){step = 1; i = (int)(68 + (mHrequested-345)/step); closestMass = (int)(step*(i-68) + 345);}
    if(mHrequested > 355 && mHrequested <= 360 ){step = 5; i = (int)(78 + (mHrequested-355)/step); closestMass = (int)(step*(i-78) + 355);}
    if(mHrequested > 360 && mHrequested <= 370 ){step = 10; i = (int)(79 + (mHrequested-360)/step); closestMass = (int)(step*(i-79) + 360);}
    if(mHrequested > 370 && mHrequested <= 380 ){step = 5; i = (int)(80 + (mHrequested-370)/step); closestMass = (int)(step*(i-80) + 370);}
    if(mHrequested > 380 && mHrequested <= 500 ){step = 10; i = (int)(82 + (mHrequested-380)/step); closestMass = (int)(step*(i-82) + 380);}
    if(mHrequested > 500 && mHrequested <= 800 ){step = 50; i = (int)(94 + (mHrequested-500)/step); closestMass = (int)(step*(i-94) + 500);}
    if(mHrequested > 800 && mHrequested <= 1000 ){step = 200; i = (int)(100 + (mHrequested-800)/step); closestMass = (int)(step*(i-100) + 800);}



      tmpLow1 = BR_SM4[ID][i]*BR_SM4[0][i];
      tmpHigh1 = BR_SM4[ID][i+1]*BR_SM4[0][i+1];

      tmpLow2 = BR_SM4[0][i];
      tmpHigh2 = BR_SM4[0][i+1];
      deltaX = mHrequested - closestMass;

      deltaY1 = tmpHigh1 - tmpLow1;
      slope1 = deltaY1/step;


      deltaY2 = tmpHigh2 - tmpLow2;
      slope2 = deltaY2/step;


      // For partial widths                                                                                                                 
      if(deltaX == 0){ PartialWidth = tmpLow1;
	TotalWidth = tmpLow2;}
      else{ PartialWidth = slope1*deltaX + tmpLow1;
	TotalWidth = slope2*deltaX + tmpLow2;}

      // For total width  
      if( ID == 0 ){ Width = TotalWidth; }
      else{ Width = PartialWidth;}

  }

  return Width;

} 



	Double_t Lgg_7( Double_t mHstar ){
		
		mHstar = mHstar/7000;
		Double_t u = log(1000*mHstar);
#if 1 	// evaluate a polynomial the proper way using Horner's method
        // http://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation
                Double_t Z =  ((((((((((-6.4e-09 * u - 7e-09)*u - 7.15e-08)*u -5.5e-09)*u - 7.5e-08)*u - 1.6e-05)*u - 1e-05)* u - 5e-05)* u - 5.5e-04)*u -.33)*u-2)*u;
#else   // unreasonably slow implementation
		Double_t z1 = (-2)*u;
		Double_t z2 = (-.33)*pow(u,2);
		Double_t z3 = (-5.5e-04)*pow(u,3);
		Double_t z4 = (-5e-05)*pow(u,4);
		Double_t z5 = (-1e-05)*pow(u,5);
		Double_t z6 = (-1.6e-05)*pow(u,6);
		Double_t z7 = (-7.5e-08)*pow(u,7);
		Double_t z8 = (-5.5e-09)*pow(u,8);
		Double_t z9 = (-7.15e-08)*pow(u,9);
		Double_t z10 = (-7e-09)*pow(u,10);
		Double_t z11 = (-6.4e-09)*pow(u,11);
		Double_t Z = z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11;
          //printf("Zold %.9f Znew %.9f reldiff %.9f \n", Z, Znew, std::abs((Znew-Z)/(Znew+Z)));
#endif

		Double_t zz = exp(Z);
		
		Double_t lum = 1.2e+10*zz;
		
		return lum;
		
	}
	
	// rho tail function
	// r is ratio = m4l/mH* 
	Double_t rho(Double_t r, TString proc)
	{
		
		Double_t a = 0,b = 0,c = 0;
		if(proc == "4mu"){ a=.608398; b=.06931; c=.10224;}
		if(proc == "4e"){ a=1.901196; b=.074754; c=.091285;}
		if(proc == "2e2mu"){ a=1.243678; b=.072587; c=.096782;}
		
		double T = (a + (b/(1-r))) * exp(-(1-r)/c);
		
		
		return T;
		
	}
	
	
	
	
	Double_t pdf1(Double_t mHstar, Double_t mHreq)
	{
		
		Double_t Gamma_gg = HiggsWidth(7,mHstar);
		Double_t Gamma_ZZ = HiggsWidth(11,mHstar);
		Double_t Gamma_TOT = HiggsWidth(0,mHstar);
		
		Double_t ggLum = Lgg_7(mHstar);
		
		Double_t pdf = Gamma_gg*Gamma_ZZ*ggLum*mHstar/( pow(pow(mHstar,2)-pow(mHreq,2),2 )+pow(mHstar,2)*pow(Gamma_TOT,2) );
		
		return pdf;
		
	}

	Double_t pdf1SM4(Double_t mHstar, Double_t mHreq)
	{
		
		Double_t Gamma_gg = HiggsWidthSM4(7,mHstar);
		Double_t Gamma_ZZ = HiggsWidthSM4(11,mHstar);
		Double_t Gamma_TOT = HiggsWidthSM4(0,mHstar);
		
		Double_t ggLum = Lgg_7(mHstar);
		
		Double_t pdf = Gamma_gg*Gamma_ZZ*ggLum*mHstar/( pow(pow(mHstar,2)-pow(mHreq,2),2 )+pow(mHstar,2)*pow(Gamma_TOT,2) );
		
		return pdf;
		
	}


  Double_t pdf1Param(Double_t mHstar, Double_t mHreq, Double_t x, Double_t widthSF)
	{
		
		Double_t Gamma_gg = HiggsWidth(7,mHstar);
		Double_t Gamma_ZZ = HiggsWidth(11,mHstar);
		Double_t Gamma_TOT = widthSF*HiggsWidth(0,mHstar);
		Double_t Gamma_TOTFixed = HiggsWidth(0,mHreq);
		
		Double_t ggLum = Lgg_7(mHstar);
		
		//Double_t pdf = Gamma_gg*Gamma_ZZ*ggLum/( pow(pow(mHstar,2)-pow(mHreq,2),2 )+pow(mHstar,2)*pow(Gamma_TOT,2) );
		Double_t numerator = Gamma_gg*pow(Gamma_ZZ, x)*pow(Gamma_TOTFixed, (1.-x))*ggLum*mHstar;
		Double_t pdf = numerator/( pow(pow(mHstar,2)-pow(mHreq,2),2 )+pow(mHstar,2)*pow(Gamma_TOT,2) );
		
		return pdf;
		
	}
	
	
	Double_t Sigma(Double_t mH, TString proc)
	{
		
		Double_t a = 0,b = 0;
		
		if(proc == "4mu")
		{
			a = .09942;
			b = .00003913;
		}
		if(proc == "2e2mu")
		{
			a = .1599;
			b = .00003396;
		}
		if(proc == "4e")
		{
			a = .2236;
			b = .00004084;
		}
		
		double sig = mH*sqrt((pow(a,2)/mH) + pow(b*mH,2));
		
		return sig;
		
	}
	
	
	Double_t sigma_CB(Double_t mH, TString proc)
	{
		
		Double_t a = 0;
		Double_t b = 0;
		Double_t c = 0;
		Double_t d = 0;
		if( proc == "4mu"  ){ a = 0.006688; b = 1.106e-05; c = 0; d = 0; }
		if( proc == "4e"   ){ a = 0.01826; b = -3.192e-05; c = 2.836e-09; d = 3.47e-11; }
		if( proc == "2e2mu"){ a = 0.01524; b = -4.259e-05; c = 1.081e-07; d = -8.564e-11; }
		
		Double_t s = mH*(a + b*mH + c*mH*mH + d*mH*mH*mH);
		
		return s;
		
	}
	
	
	Double_t mean(Double_t mH, TString proc)
	{
		
		Double_t a = 0;
		
		if( proc == "4mu"  ){ a = -0.0001532;}
		if( proc == "4e"   ){ a = -0.002602; }
		if( proc == "2e2mu"){ a = -0.00162; }
		
		Double_t m = a*mH;
		
		return m;
		
	}
	
	Double_t N(double mH, TString proc)
	{
		
		Double_t a = 0;
		Double_t b = 0;
		if( proc == "4mu"  ){ a = -0.04191; b = 50.69; }
		if( proc == "4e"   ){ a = 0.05952;  b = 15.38; }
		if( proc == "2e2mu"){ a = 0.01461;  b = 28.42; }
		
		Double_t N = a + (b/mH);
		
		return N;
		
	}
	
	
}




/**********RooqqZZPdf*********/

ClassImp(RooqqZZPdf)

RooqqZZPdf::RooqqZZPdf(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a1,
					   RooAbsReal& _a2,
					   RooAbsReal& _a3,
					   RooAbsReal& _b1,
					   RooAbsReal& _b2,
					   RooAbsReal& _b3,
					   RooAbsReal& _frac) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
b1("b1","b1",this,_b1),
b2("b2","b2",this,_b2),
b3("b3","b3",this,_b3),
frac("frac","frac",this,_frac)
{
}


RooqqZZPdf::RooqqZZPdf(const RooqqZZPdf& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
b1("b1",this,other.b1),
b2("b2",this,other.b2),
b3("b3",this,other.b3),
frac("frac",this,other.frac)
{
}


Double_t RooqqZZPdf::evaluate() const
{
	
	//std::cout << "a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << std::endl;                                                                    
	//std::cout << "b1 = " << b1 << ", b2 = " << b2 << ", b3 = " << b3 << std::endl;                                                                    
	
	double signa = 0.;
	if ((m4l-a1) > 0) { signa = 1.; }
	else if ((m4l-a1) < 0) { signa = -1.; }
	else { signa = 0.; }
	
	double signb = 0.;
	if ((m4l-b1) > 0) { signb = 1.; }
	else if ((m4l-b1) < 0) { signb = -1.; }
	else { signb = 0.; }
	
	double bkglo = (0.5 + 0.5*signa * TMath::Erf(TMath::Abs(m4l-a1)/a2)) * exp(-1.*m4l/a3);
	double bkghi = (0.5 + 0.5*signb * TMath::Erf(TMath::Abs(m4l-b1)/b2)) * exp(-1.*m4l/b3);
	double total = bkglo*frac + (1-frac)*bkghi;
	
	double dynamicKqq = 1. + m4l*0.001074 + m4l*m4l*-7.851e-07;
	
	double totalNLO = total*dynamicKqq;
	return totalNLO ;
}


/************RooggZZPdf***********/

ClassImp(RooggZZPdf)

RooggZZPdf::RooggZZPdf(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a1,
					   RooAbsReal& _a2,
					   RooAbsReal& _a3,
					   RooAbsReal& _b1,
					   RooAbsReal& _b2,
					   RooAbsReal& _b3,
					   RooAbsReal& _frac) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
b1("b1","b1",this,_b1),
b2("b2","b2",this,_b2),
b3("b3","b3",this,_b3),
frac("frac","frac",this,_frac)
{
}


RooggZZPdf::RooggZZPdf(const RooggZZPdf& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
b1("b1",this,other.b1),
b2("b2",this,other.b2),
b3("b3",this,other.b3),
frac("frac",this,other.frac)
{
}

Double_t RooggZZPdf::evaluate() const
{
	
	//std::cout << "a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << std::endl;                                                                    
	//std::cout << "b1 = " << b1 << ", b2 = " << b2 << ", b3 = " << b3 << std::endl;                                                                    
	
	// ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE                                                                                             
	double signa = 0.;
	if ((m4l-a1) > 0) { signa = 1.; }
	else if ((m4l-a1) < 0) { signa = -1.; }
	else { signa = 0.; }
	
	double signb = 0.;
	if ((m4l-b1) > 0) { signb = 1.; }
	else if ((m4l-b1) < 0) { signb = -1.; }
	else { signb = 0.; }
	
	double bkglo = (0.5 + 0.5*signa * TMath::Erf(TMath::Abs(m4l-a1)/a2)) * exp(-1.*m4l/a3);
	double bkghi = (0.5 + 0.5*signb * TMath::Erf(TMath::Abs(m4l-b1)/b2)) * exp(-1.*m4l/b3);
	double total = bkglo*frac + (1-frac)*bkghi;
	
	double m_a = 0.17;
	double m_b = -2.95e-04;
	double m_c = 1.55e-07;
	double m_d = 169.2;
	double m_e = 45.94;	
	double dynamicKgg = (m_a + m_b*m4l + m_c*m4l*m4l)*(0.5 + 0.5*TMath::Erf( (m4l-m_d)/m_e ) );
	double totalNLO = total*dynamicKgg;
	
	return totalNLO ;
}

//// ------- v2 below ------

ClassImp(RooqqZZPdf_v2)

RooqqZZPdf_v2::RooqqZZPdf_v2(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a0,
					   RooAbsReal& _a1,
					   RooAbsReal& _a2,
					   RooAbsReal& _a3,
					   RooAbsReal& _a4,
					   RooAbsReal& _a5,
					   RooAbsReal& _a6,
					   RooAbsReal& _a7,
					   RooAbsReal& _a8,
					   RooAbsReal& _a9,
					   RooAbsReal& _a10,
					   RooAbsReal& _a11,
					   RooAbsReal& _a12,
					   RooAbsReal& _a13
					   ) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a0("a0","a0",this,_a0),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
a4("a4","a4",this,_a4),
a5("a5","a5",this,_a5),
a6("a6","a6",this,_a6),
a7("a7","a7",this,_a7),
a8("a8","a8",this,_a8),
a9("a9","a9",this,_a9),
a10("a10","a10",this,_a10),
a11("a11","a11",this,_a11),
a12("a12","a12",this,_a12),
a13("a13","a13",this,_a13)

{
}


RooqqZZPdf_v2::RooqqZZPdf_v2(const RooqqZZPdf_v2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a0("a0",this,other.a0),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
a4("a4",this,other.a4),
a5("a5",this,other.a5),
a6("a6",this,other.a6),
a7("a7",this,other.a7),
a8("a8",this,other.a8),
a9("a9",this,other.a9),
a10("a10",this,other.a10),
a11("a11",this,other.a11),
a12("a12",this,other.a12),
a13("a13",this,other.a13)

{
}


Double_t RooqqZZPdf_v2::evaluate() const
{
	
	double ZZ = (.5+.5*TMath::Erf((m4l-a0)/a1))*(a3/(1+exp((m4l-a0)/a2)))+
    (.5+.5*TMath::Erf((m4l-a4)/a5))*(a7/(1+exp((m4l-a4)/a6))+a9/(1+exp((m4l-a4)/a8)))
    +(.5+.5*TMath::Erf((m4l-a10)/a11))*(a13/(1+exp((m4l-a10)/a12)) );
	
	return ZZ;
}


#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
std::unique_ptr<RooAbsArg> RooqqZZPdf_v2::compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext & ctx) const
{
    return RooGenericPdf{GetName(),
        "(.5+.5*TMath::Erf((x[0]-x[1])/x[2]))*(x[4]/(1+exp((x[0]-x[1])/x[3])))+(.5+.5*TMath::Erf((x[0]-x[5])/x[6]))*(x[8]/(1+exp((x[0]-x[5])/x[7]))+x[10]/(1+exp((x[0]-x[5])/x[9])))+(.5+.5*TMath::Erf((x[0]-x[11])/x[12]))*(x[14]/(1+exp((x[0]-x[11])/x[13])) )",
        {*m4l, *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8, *a9, *a10, *a11, *a12, *a13}
    }.compileForNormSet(normSet, ctx);
}
#endif


/************RooggZZPdf_v2***********/

ClassImp(RooggZZPdf_v2)

RooggZZPdf_v2::RooggZZPdf_v2(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a0,
					   RooAbsReal& _a1,
					   RooAbsReal& _a2,
					   RooAbsReal& _a3,
					   RooAbsReal& _a4,
					   RooAbsReal& _a5,
					   RooAbsReal& _a6,
					   RooAbsReal& _a7,
					   RooAbsReal& _a8,
					   RooAbsReal& _a9
					   ) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a0("a0","a0",this,_a0),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
a4("a4","a4",this,_a4),
a5("a5","a5",this,_a5),
a6("a6","a6",this,_a6),
a7("a7","a7",this,_a7),
a8("a8","a8",this,_a8),
a9("a9","a9",this,_a9)

{
}


RooggZZPdf_v2::RooggZZPdf_v2(const RooggZZPdf_v2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a0("a0",this,other.a0),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
a4("a4",this,other.a4),
a5("a5",this,other.a5),
a6("a6",this,other.a6),
a7("a7",this,other.a7),
a8("a8",this,other.a8),
a9("a9",this,other.a9)
{
}


Double_t RooggZZPdf_v2::evaluate() const
{
	
	double ZZ = (.5+.5*TMath::Erf((m4l-a0)/a1))*(a3/(1+exp((m4l-a0)/a2)))+
    (.5+.5*TMath::Erf((m4l-a4)/a5))*(a7/(1+exp((m4l-a4)/a6))+a9/(1+exp((m4l-a4)/a8)));
	
	
	return ZZ;
}

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
std::unique_ptr<RooAbsArg> RooggZZPdf_v2::compileForNormSet(RooArgSet const &normSet, RooFit::Detail::CompileContext & ctx) const
{
    return RooGenericPdf{GetName(),
        "(.5+.5*TMath::Erf((x[0]-x[1])/x[2]))*(x[4]/(1+exp((x[0]-x[1])/x[3]))) + (.5+.5*TMath::Erf((x[0]-x[5])/x[6]))*(x[8]/(1+exp((x[0]-x[5])/x[7]))+x[10]/(1+exp((x[0]-x[5])/x[9])))",
        {*m4l, *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8, *a9}
    }.compileForNormSet(normSet, ctx);
}
#endif


//// ------- RooVBFZZPdf ------

ClassImp(RooVBFZZPdf)

RooVBFZZPdf::RooVBFZZPdf(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a0,
					   RooAbsReal& _a1,
					   RooAbsReal& _a2,
					   RooAbsReal& _a3,
					   RooAbsReal& _a4,
					   RooAbsReal& _a5,
					   RooAbsReal& _a6,
					   RooAbsReal& _a7,
					   RooAbsReal& _a8,
					   RooAbsReal& _a9,
					   RooAbsReal& _a10,
					   RooAbsReal& _a11,
					   RooAbsReal& _a12,
			                   RooAbsReal& _a13,
			                   RooAbsReal& _a14,
			                   RooAbsReal& _a15,
			                   RooAbsReal& _a16
					   ) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a0("a0","a0",this,_a0),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
a4("a4","a4",this,_a4),
a5("a5","a5",this,_a5),
a6("a6","a6",this,_a6),
a7("a7","a7",this,_a7),
a8("a8","a8",this,_a8),
a9("a9","a9",this,_a9),
a10("a10","a10",this,_a10),
a11("a11","a11",this,_a11),
a12("a12","a12",this,_a12),
a13("a13","a13",this,_a13),
a14("a14","a14",this,_a14),
a15("a15","a15",this,_a15),
a16("a16","a16",this,_a16)

{
}


RooVBFZZPdf::RooVBFZZPdf(const RooVBFZZPdf& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a0("a0",this,other.a0),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
a4("a4",this,other.a4),
a5("a5",this,other.a5),
a6("a6",this,other.a6),
a7("a7",this,other.a7),
a8("a8",this,other.a8),
a9("a9",this,other.a9),
a10("a10",this,other.a10),
a11("a11",this,other.a11),
a12("a12",this,other.a12),
a13("a13",this,other.a13),
a14("a14",this,other.a14),
a15("a15",this,other.a15),
a16("a16",this,other.a16)

{
}


Double_t RooVBFZZPdf::evaluate() const
{
	
  double ZZ = (.5+.5*TMath::Erf((m4l-a0)/a1))*(a3/(1+exp((m4l-a0)/a2))) +
    (.5+.5*TMath::Erf((m4l-a4)/a5))*(a7/(1+exp((m4l-a4)/a6))+a9/(1+exp((m4l-a4)/a8))) + 
    (.5+.5*TMath::Erf((m4l-a10)/a11))*(a13/(1+exp((m4l-a10)/a12))) +
    a14*exp(-fabs((m4l-a16)/a15));


// // Decent
//   double ZZ = (a3/(1+exp((m4l-a0)/a2))) +
//     (.5+.5*TMath::Erf((m4l-a4)/a5))*(a7/(1+exp((m4l-a4)/a6))+a9/(1+exp((m4l-a4)/a8))) + 
//     (.5+.5*TMath::Erf((m4l-a10)/a11))*(a13/(1+exp((m4l-a10)/a12))) +
//     a14*exp(-fabs((m4l-a16)/a15));




	
	return ZZ;
}



//// ------- RooVBFZZPdf_v2 ------

ClassImp(RooVBFZZPdf_v2)

RooVBFZZPdf_v2::RooVBFZZPdf_v2(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _a4,
					   RooAbsReal& _a5,
					   RooAbsReal& _a6,
					   RooAbsReal& _a7,
					   RooAbsReal& _a8,
					   RooAbsReal& _a9,
					   RooAbsReal& _a10,
					   RooAbsReal& _a11,
					   RooAbsReal& _a12,
			                   RooAbsReal& _a13,
			                   RooAbsReal& _a14,
			                   RooAbsReal& _a15,
			                   RooAbsReal& _a16
					   ) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
a4("a4","a4",this,_a4),
a5("a5","a5",this,_a5),
a6("a6","a6",this,_a6),
a7("a7","a7",this,_a7),
a8("a8","a8",this,_a8),
a9("a9","a9",this,_a9),
a10("a10","a10",this,_a10),
a11("a11","a11",this,_a11),
a12("a12","a12",this,_a12),
a13("a13","a13",this,_a13),
a14("a14","a14",this,_a14),
a15("a15","a15",this,_a15),
a16("a16","a16",this,_a16)

{
}


RooVBFZZPdf_v2::RooVBFZZPdf_v2(const RooVBFZZPdf_v2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
a4("a4",this,other.a4),
a5("a5",this,other.a5),
a6("a6",this,other.a6),
a7("a7",this,other.a7),
a8("a8",this,other.a8),
a9("a9",this,other.a9),
a10("a10",this,other.a10),
a11("a11",this,other.a11),
a12("a12",this,other.a12),
a13("a13",this,other.a13),
a14("a14",this,other.a14),
a15("a15",this,other.a15),
a16("a16",this,other.a16)

{
}


Double_t RooVBFZZPdf_v2::evaluate() const
{
	
  double ZZ = (.5+.5*TMath::Erf((m4l-a4)/a5))*(a7/(1+exp((m4l-a4)/a6))+a9/(1+exp((m4l-a4)/a8))) + 
    + (.5+.5*TMath::Erf((m4l-a10)/a11))*(a13/(1+exp((m4l-a10)/a12)));
    //+ a14*exp(-fabs((m4l-a16)/a15));
	
	return ZZ;
}




ClassImp(RooBetaFunc_v2) 

RooBetaFunc_v2::RooBetaFunc_v2(){}

RooBetaFunc_v2::RooBetaFunc_v2(const char *name, const char *title, 
							   RooAbsReal& _mZstar,	     
							   RooAbsReal& _mZ,
							   RooAbsReal& _m0,
							   RooAbsReal& _mZZ,	     
							   RooAbsReal& _Gamma,
							   RooAbsReal& _Gamma0,
							   RooAbsReal& _a0,  // mZZ distribution vars
							   RooAbsReal& _a1,
							   RooAbsReal& _a2,
							   RooAbsReal& _a3,
							   RooAbsReal& _f,
							   RooAbsReal& _f0
							   ) :
RooAbsPdf(name,title), 
mZstar("mZstar","mZstar",this,_mZstar),
mZ("mZ","mZ",this,_mZ),
m0("m0","m0",this,_m0),
mZZ("mZZ","mZZ",this,_mZZ),
Gamma("Gamma","Gamma",this,_Gamma),
Gamma0("Gamma0","Gamma0",this,_Gamma0),
a0("a0","a0",this,_a0),
a1("a1","a1",this,_a1),
a2("a2","a2",this,_a2),
a3("a3","a3",this,_a3),
f("f","f",this,_f),
f0("f0","f0",this,_f0)
{ 
} 


RooBetaFunc_v2::RooBetaFunc_v2(const RooBetaFunc_v2& other, const char* name) :  
RooAbsPdf(other,name), 
mZstar("mZstar",this,other.mZstar),
mZ("mZ",this,other.mZ),
m0("m0",this,other.m0),
mZZ("mZZ",this,other.mZZ),
Gamma("Gamma",this,other.Gamma),
Gamma0("Gamma0",this,other.Gamma0),
a0("a0",this,other.a0),
a1("a1",this,other.a1),
a2("a2",this,other.a2),
a3("a3",this,other.a3),
f("f",this,other.f),
f0("f0",this,other.f0)
{ 
} 



double RooBetaFunc_v2::evaluate() const 
{ 
	
	double mZDistribution, acceptance;
	
	//================= mZ Distribution =================
	double beta= (1-(mZstar-mZ)*(mZstar-mZ)/(mZZ*mZZ))*(1-(mZstar+mZ)*(mZstar+mZ)/(mZZ*mZZ));
	if(beta<0) return 0.0;
	else mZDistribution = beta;
	
	//================= acceptance (chebychev polynomials) 
	acceptance = a0+a1*mZstar+a2*(2*mZstar*mZstar-1)+a3*(4*mZstar*mZstar*mZstar-3*mZstar);
	
	//================= on-shell Z contamination ========
	
	double onShellZ = exp(-(mZstar-mZ)*(mZstar-mZ)/(2*Gamma*Gamma))/sqrt(2*3.1415*Gamma*Gamma);
	//double errorFunc = TMath::Erf((mZstar - m0)/Gamma);
	double errorFunc = exp(-(mZstar-m0)*(mZstar-m0)/(2*Gamma0*Gamma0))/sqrt(2*3.1415*Gamma0*Gamma0);
	
	//std::cout << mZZDistribution << " : " << mZDistribution << " : " << acceptance << std::endl;
	double final =  mZDistribution*((1-f)*(f0*acceptance+(1-f0)*errorFunc)+f*onShellZ);
	if (final <= 0) final = 1e-12;
	return final;
}

ClassImp(Roo4lMasses2D_Bkg) 

Roo4lMasses2D_Bkg::Roo4lMasses2D_Bkg(){}

Roo4lMasses2D_Bkg::Roo4lMasses2D_Bkg(const char *name, const char *title, 
									 RooAbsReal& _mZstar,	     
									 RooAbsReal& _mZZ,	     
									 RooAbsReal& _channelVal
									 ) :
RooAbsPdf(name,title), 
mZstar("mZstar","mZstar",this,_mZstar),
mZZ("mZZ","mZZ",this,_mZZ),
channelVal("channelVal","channelVal",this,_channelVal)
{ 
} 


Roo4lMasses2D_Bkg::Roo4lMasses2D_Bkg(const Roo4lMasses2D_Bkg& other, const char* name) :  
RooAbsPdf(other,name), 
mZstar("mZstar",this,other.mZstar),
mZZ("mZZ",this,other.mZZ),
channelVal("channelVal",this,other.channelVal)
{ 
} 


double Roo4lMasses2D_Bkg::UnitStep(double arg) const
{
	if(arg<0.0){
		return 0.0;
	}
	else{
		return 1.0;
	}
}


double Roo4lMasses2D_Bkg::evaluate() const 
{ 
	
	double Gamma=0, mZ=0, a0=0, a1=0, a2=0, a3=0, Gamma0=0, m0=0, f=0, f0=0;
	double a0_bkgd=0, a1_bkgd=0, a2_bkgd=0, a3_bkgd=0, a4_bkgd=0, a5_bkgd=0, a6_bkgd=0, a7_bkgd=0;
	double a8_bkgd=0, a9_bkgd=0, a10_bkgd=0, a11_bkgd=0, a12_bkgd=0, a13_bkgd=0;
	
	a2 = 0.; a3 = 0.;
	if (channelVal == 1.){
		mZ = 91.2;
		Gamma0 = 24.3063 + -0.169814*mZZ + 0.00274393*(mZZ - 108.86)*(mZZ - 108.86);
		m0 = -273.264 + 6.81815*mZZ + -0.0520652*TMath::Power(mZZ,2) + 0.000129716*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = -34.7502 + 1.15173*mZZ +-0.0152287*TMath::Power(mZZ,2) + 0.000100538*TMath::Power(mZZ,3)+ -3.317e-07*TMath::Power(mZZ,4)+ 4.37829e-10*TMath::Power(mZZ,5);
		a2 = 0.0710702 + -0.00232026*mZZ + 3.00042e-05*TMath::Power(mZZ,2) + -1.92831e-07*TMath::Power(mZZ,3)+ 6.18009e-10*TMath::Power(mZZ,4)+ -7.92751e-13*TMath::Power(mZZ,5);
		f = 4.81681e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);
		/*
		 Gamma = 1.50320e+01;
		 Gamma0 = 83.4831 + -0.711931*mZZ + 0.0367013*(mZZ - 104.704)*(mZZ - 104.704);
		 m0 = 146.261 + -1.12857*mZZ + 0.0401877*(mZZ - 106.72)*(mZZ - 106.72);
		 a0 = -0.0348702 + 0.000251279*mZZ + 5.63152e-07*(mZZ - 106.548)*(mZZ - 106.548);
		 a1 = 0.0013887 + -9.69053e-06*mZZ + 1.62389e-11*(mZZ - 105.303)*(mZZ - 105.303)*(mZZ - 105.303)*(mZZ - 105.303);
		 f = 4.02164e-01*(TMath::Erf( (mZZ-1.55799e+02)/1.36356e+01 ) + 1.);
		 f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);	
		 */
		
		a0_bkgd = 116.19;
		a1_bkgd = 4.52219;
		a2_bkgd = 105.493;
		a3_bkgd = 0.0823749;
		a4_bkgd = 180.89;
		a5_bkgd = 9.58208;
		a6_bkgd = 30.445;
		a7_bkgd = 0.182481;
		a8_bkgd = 56.1627;
		a9_bkgd = 0.0957905;
		a10_bkgd = 119.406;
		a11_bkgd = -99.674;
		a12_bkgd = 9995.69;
		a13_bkgd = 0.110563;
	}
	else if (channelVal == 2.){
		mZ = 91.2;
		Gamma0 = 64.5791 + -0.465027*mZZ + 0.00680112*(mZZ - 109.083)*(mZZ - 109.083);
		m0 = -4528.09 + 98.4925*mZZ + -0.70193*TMath::Power(mZZ,2) + 0.00164952*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = 0.0549555 + -0.000367875*mZZ + 6.16188e-10*(mZZ - 106.028)*(mZZ - 106.028)*(mZZ - 106.028)*(mZZ - 106.028);
		a2 = -0.000222052 + 1.51606e-06*mZZ + -9.57278e-08*(mZZ - 144.667)*(mZZ - 144.667);
		f = 4.82988e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);	
		
		/*
		 Gamma = 1.50320e+01;
		 Gamma0 = -56.1597 + 0.542356*mZZ + -0.00250288*(mZZ - 107.336)*(mZZ - 107.336);
		 m0 = -2.14585 + 0.188968*mZZ + 0.0210833*(mZZ - 108.32)*(mZZ - 108.32);
		 a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		 a1 = -0.000170976 + 1.83011e-06*mZZ + -1.05582e-11*(mZZ - 105.893)*(mZZ - 105.893)*(mZZ - 105.893)*(mZZ - 105.893);
		 f = 3.53555e-01*(TMath::Erf( (mZZ-1.53281e+02)/1.75311e+00 ) + 1.);
		 f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);
		 */
		a0_bkgd = 120.096;
		a1_bkgd = 15.6766;
		a2_bkgd = 149.141;
		a3_bkgd = 0.0361427;
		a4_bkgd = 181.587;
		a5_bkgd = 10.4505;
		a6_bkgd = 9.17084;
		a7_bkgd = 0.0345372;
		a8_bkgd = 38.655;
		a9_bkgd = 0.114464;
		a10_bkgd = 87.2772;
		a11_bkgd = -6.94956;
		a12_bkgd = 3.35268;
		a13_bkgd = 0.0320278;
		/*
		 a0_bkgd = 114.507;
		 a1_bkgd = 16.8809;
		 a2_bkgd = 116.481;
		 a3_bkgd = 0.0386084;
		 a4_bkgd = 182.841;
		 a5_bkgd = 15.2248;
		 a6_bkgd = 27.79;
		 a7_bkgd = 0.133664;
		 a8_bkgd = 59.2288;
		 a9_bkgd = 0.0400656;
		 a10_bkgd = 93.8226;
		 a11_bkgd = -2.64585;
		 a12_bkgd = 5066.99;
		 a13_bkgd = 0.125328;	
		 */
	}
	else if (channelVal == 3.){
		mZ = 91.2;
		Gamma0 = 36.9986 + -0.250136*mZZ + 0.00474031*(mZZ - 108.732)*(mZZ - 108.732);
		m0 = -3266.9 + 73.5969*mZZ + -0.544448*TMath::Power(mZZ,2) + 0.00133127*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = 0.0206819 + -0.000168968*mZZ + 4.26934e-11*(mZZ - 40.4456)*(mZZ - 40.4456)*(mZZ - 40.4456)*(mZZ - 40.4456);
		a2 = -1.271e-05 + 8.742e-08*mZZ + -6.34276e-08*(mZZ - 144.032)*(mZZ - 144.032);
		f = 4.82988e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);	
		/*
		 Gamma = 1.50320e+01;
		 mZ = 91.2;
		 Gamma0 = -78.9298 + 0.744417*mZZ + -0.0105399*(mZZ - 109.541)*(mZZ - 109.541);
		 m0 = -76.0686 + 0.840304*mZZ + -0.00797722*(mZZ - 105.677)*(mZZ - 105.677);
		 a0 = 0.378445 + -0.00310967*mZZ + 5.78471e-05*(mZZ - 107.739)*(mZZ - 107.739);
		 a1 = 0.00109279 + -8.90024e-06*mZZ + 5.21042e-12*(mZZ - 103.517)*(mZZ - 103.517)*(mZZ - 103.517)*(mZZ - 103.517);
		 f = 4.81681e-01*(TMath::Erf( (mZZ-1.30465e+02)/8.69644e-01 ) + 1.);
		 f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);
		 */
		a0_bkgd = 114.504;
		a1_bkgd = 14.3013;
		a2_bkgd = 101.43;
		a3_bkgd = 0.0455842;
		a4_bkgd = 179.976;
		a5_bkgd = 9.89769;
		a6_bkgd = 28.5377;
		a7_bkgd = 0.0692236;
		a8_bkgd = 51.5593;
		a9_bkgd = 0.0561965;
		a10_bkgd = 101.547;
		a11_bkgd = -1.71367;
		a12_bkgd = 7513.75;
		a13_bkgd = 0.543254;
	}
	else{
		std::cout << "Invalid parameters for this PDF!" << std::endl;
	}
	
	
	//  ZZ shape
	double ZZshape = (.5+.5*TMath::Erf((mZZ-a0_bkgd)/a1_bkgd))*(a3_bkgd/(1+exp((mZZ-a0_bkgd)/a2_bkgd)))+(.5+.5*TMath::Erf((mZZ-a4_bkgd)/a5_bkgd))*(a7_bkgd/(1+exp((mZZ-a4_bkgd)/a6_bkgd))+a9_bkgd/(1+exp((mZZ-a4_bkgd)/a8_bkgd)))+(.5+.5*TMath::Erf((mZZ-a10_bkgd)/a11_bkgd))*(a13_bkgd/(1+exp((mZZ-a10_bkgd)/a12_bkgd)) );
	
	
	double mZDistribution, acceptance;
	
	//================= mZ Distribution =================
	double beta= (1-(mZstar-mZ)*(mZstar-mZ)/(mZZ*mZZ))*(1-(mZstar+mZ)*(mZstar+mZ)/(mZZ*mZZ));
	if(beta<0) return 0.0000001;
	else mZDistribution = beta;
	
	//================= acceptance (chebychev polynomials) 
	acceptance = a0+a1*mZstar+a2*(2*mZstar*mZstar-1)+a3*(4*mZstar*mZstar*mZstar-3*mZstar);
	
	//================= on-shell Z contamination ========
	
	double onShellZ = exp(-(mZstar-mZ)*(mZstar-mZ)/(2*Gamma*Gamma))/sqrt(2*3.1415*Gamma*Gamma);
	//double errorFunc = TMath::Erf((mZstar - m0)/Gamma);
	double errorFunc = exp(-(mZstar-m0)*(mZstar-m0)/(2*Gamma0*Gamma0))/sqrt(2*3.1415*Gamma0*Gamma0);
	
	//std::cout << mZZDistribution << " : " << mZDistribution << " : " << acceptance << std::endl;
	double final =  mZDistribution*((1-f)*(f0*acceptance+(1-f0)*errorFunc)+f*onShellZ);
	if (final <= 0) final = 1e-12;
	
	///*
	double integral;
	double E = 2.71828183;
	
	// ========================================================
	// integral
	///*
	double mZstarFix=12.;
	double integralT1Lo=((1 - f)*f0*((a1*Power(mZstarFix,6))/6. + (2*a2*Power(mZstarFix,7))/7. + (a0 - a2)*mZstarFix*Power(Power(mZ,2) - Power(mZZ,2),2) + (a1*Power(mZstarFix,2)*Power(Power(mZ,2) - Power(mZZ,2),2))/2. - (a1*Power(mZstarFix,4)*(Power(mZ,2) + Power(mZZ,2)))/2. + (Power(mZstarFix,5)*(a0 - a2*(1 + 4*Power(mZ,2) + 4*Power(mZZ,2))))/5. + (2*Power(mZstarFix,3)*(-(a0*(Power(mZ,2) + Power(mZZ,2))) + a2*(Power(mZ,4) + Power(mZZ,2) + Power(mZZ,4) + Power(mZ,2)*(1 - 2*Power(mZZ,2)))))/3.))/Power(mZZ,4);
	mZstarFix=mZZ-mZ;
	double integralT1Hi=((1 - f)*f0*((a1*Power(mZstarFix,6))/6. + (2*a2*Power(mZstarFix,7))/7. + (a0 - a2)*mZstarFix*Power(Power(mZ,2) - Power(mZZ,2),2) + (a1*Power(mZstarFix,2)*Power(Power(mZ,2) - Power(mZZ,2),2))/2. - (a1*Power(mZstarFix,4)*(Power(mZ,2) + Power(mZZ,2)))/2. + (Power(mZstarFix,5)*(a0 - a2*(1 + 4*Power(mZ,2) + 4*Power(mZZ,2))))/5. + (2*Power(mZstarFix,3)*(-(a0*(Power(mZ,2) + Power(mZZ,2))) + a2*(Power(mZ,4) + Power(mZZ,2) + Power(mZZ,4) + Power(mZ,2)*(1 - 2*Power(mZZ,2)))))/3.))/Power(mZZ,4);
	
	//std::cout << "IntegralHI: " << integralT1Hi << ", IntegralLO: " << integralT1Lo << std::endl;
	
	mZstarFix=12.;
	double t2loA = (3*Power(Gamma0,4) + Power(m0,4) + Power(mZ,4) - 2*Power(Gamma0,2)*Power(mZZ,2) + Power(mZZ,4) + Power(m0,2)*(6*Power(Gamma0,2) - 2*Power(mZ,2) - 2*Power(mZZ,2)) - 2*Power(mZ,2)*(Power(Gamma0,2) + Power(mZZ,2)));
	double integralT2Lo=((-1 + f)*(-1 + f0)*(-1.*((Power(Gamma0,2)*(Power(m0,3) + Power(m0,2)*mZstarFix + mZstarFix*(3*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) + m0*(5*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(m0 - mZstarFix,2)/(2.*Power(Gamma0,2)))) - Gamma0*t2loA*sqrt(TMath::Pi()/2.)*TMath::Erf((m0 - mZstarFix)/(sqrt(2)*Gamma0))))/Power(mZZ,4);
	integralT2Lo/=sqrt(2*TMath::Pi()*Gamma0*Gamma0);
	mZstarFix=mZZ-mZ;
	double t2hiA = (3*Power(Gamma0,4) + Power(m0,4) + Power(mZ,4) - 2*Power(Gamma0,2)*Power(mZZ,2) + Power(mZZ,4) + Power(m0,2)*(6*Power(Gamma0,2) - 2*Power(mZ,2) - 2*Power(mZZ,2)) - 2*Power(mZ,2)*(Power(Gamma0,2) + Power(mZZ,2)));
	double integralT2Hi=((-1 + f)*(-1 + f0)*(-1.*((Power(Gamma0,2)*(Power(m0,3) + Power(m0,2)*mZstarFix + mZstarFix*(3*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) + m0*(5*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(m0 - mZstarFix,2)/(2.*Power(Gamma0,2)))) - Gamma0*t2hiA*sqrt(TMath::Pi()/2.)*TMath::Erf((m0 - mZstarFix)/(sqrt(2)*Gamma0))))/Power(mZZ,4);
	integralT2Hi/=sqrt(2*TMath::Pi()*Gamma0*Gamma0);
	
	mZstarFix=12.;
	double t3loA = (3*Power(Gamma,4) + 4*Power(Gamma,2)*Power(mZ,2) - 2*(Power(Gamma,2) + 2*Power(mZ,2))*Power(mZZ,2) + Power(mZZ,4));
	double integralT3Lo=(f*((Power(Gamma,2)*(Power(mZ,3) + Power(mZ,2)*mZstarFix - mZstarFix*(3*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) - mZ*(5*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(mZ - mZstarFix,2)/(2.*Power(Gamma,2))) - Gamma*t3loA*sqrt(TMath::Pi()/2.)*TMath::Erf((mZ - mZstarFix)/(sqrt(2)*Gamma))))/Power(mZZ,4);
	integralT3Lo/=sqrt(2*TMath::Pi()*Gamma*Gamma);
	mZstarFix=mZZ-mZ;
	double t3hiA = (3*Power(Gamma,4) + 4*Power(Gamma,2)*Power(mZ,2) - 2*(Power(Gamma,2) + 2*Power(mZ,2))*Power(mZZ,2) + Power(mZZ,4));
	double integralT3Hi=(f*((Power(Gamma,2)*(Power(mZ,3) + Power(mZ,2)*mZstarFix - mZstarFix*(3*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) - mZ*(5*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(mZ - mZstarFix,2)/(2.*Power(Gamma,2))) - Gamma*t3hiA*sqrt(TMath::Pi()/2.)*TMath::Erf((mZ - mZstarFix)/(sqrt(2)*Gamma))))/Power(mZZ,4);
	integralT3Hi/=sqrt(2*TMath::Pi()*Gamma*Gamma);
	
	integral = integralT1Hi + integralT2Hi + integralT3Hi - integralT1Lo - integralT2Lo - integralT3Lo;
	
	if (integral <= 0.) {
		std::cout << "integral is less than zero...." << std::endl;
		//integral = 100000000.;
	}
	
	return (final/integral)*ZZshape;
	//*/
	//return final;
}


// 2D ggzz

ClassImp(Roo4lMasses2D_BkgGGZZ) 

Roo4lMasses2D_BkgGGZZ::Roo4lMasses2D_BkgGGZZ(){}

Roo4lMasses2D_BkgGGZZ::Roo4lMasses2D_BkgGGZZ(const char *name, const char *title, 
									 RooAbsReal& _mZstar,	     
									 RooAbsReal& _mZZ,	     
									 RooAbsReal& _channelVal
									 ) :
RooAbsPdf(name,title), 
mZstar("mZstar","mZstar",this,_mZstar),
mZZ("mZZ","mZZ",this,_mZZ),
channelVal("channelVal","channelVal",this,_channelVal)
{ 
} 


Roo4lMasses2D_BkgGGZZ::Roo4lMasses2D_BkgGGZZ(const Roo4lMasses2D_BkgGGZZ& other, const char* name) :  
RooAbsPdf(other,name), 
mZstar("mZstar",this,other.mZstar),
mZZ("mZZ",this,other.mZZ),
channelVal("channelVal",this,other.channelVal)
{ 
} 


double Roo4lMasses2D_BkgGGZZ::UnitStep(double arg) const
{
	if(arg<0.0){
		return 0.0;
	}
	else{
		return 1.0;
	}
}


double Roo4lMasses2D_BkgGGZZ::evaluate() const 
{ 
	
	double Gamma=0, mZ=0, a0=0, a1=0, a2=0, a3=0, Gamma0=0, m0=0, f=0, f0=0;
	double a0_bkgd=0, a1_bkgd=0, a2_bkgd=0, a3_bkgd=0, a4_bkgd=0, a5_bkgd=0, a6_bkgd=0, a7_bkgd=0;
	double a8_bkgd=0, a9_bkgd=0;//, a10_bkgd=0, a11_bkgd=0, a12_bkgd=0, a13_bkgd=0;
	
	a2 = 0.; a3 = 0.;
	if (channelVal == 1.){
		mZ = 91.2;
		Gamma0 = 24.3063 + -0.169814*mZZ + 0.00274393*(mZZ - 108.86)*(mZZ - 108.86);
		m0 = -273.264 + 6.81815*mZZ + -0.0520652*TMath::Power(mZZ,2) + 0.000129716*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = -34.7502 + 1.15173*mZZ +-0.0152287*TMath::Power(mZZ,2) + 0.000100538*TMath::Power(mZZ,3)+ -3.317e-07*TMath::Power(mZZ,4)+ 4.37829e-10*TMath::Power(mZZ,5);
		a2 = 0.0710702 + -0.00232026*mZZ + 3.00042e-05*TMath::Power(mZZ,2) + -1.92831e-07*TMath::Power(mZZ,3)+ 6.18009e-10*TMath::Power(mZZ,4)+ -7.92751e-13*TMath::Power(mZZ,5);
		f = 4.81681e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);

		a0_bkgd = 138.962;
		a1_bkgd = 30.6644;
		a2_bkgd = 102.09;
		a3_bkgd = 0.0416272;
		a4_bkgd = 178.464;
		a5_bkgd = 12.6718;
		a6_bkgd = 38.1437;
		a7_bkgd = 0.13019;
		a8_bkgd = 38.1146;
		a9_bkgd = 0.0277349;
	}
	else if (channelVal == 2.){
		mZ = 91.2;
		Gamma0 = 64.5791 + -0.465027*mZZ + 0.00680112*(mZZ - 109.083)*(mZZ - 109.083);
		m0 = -4528.09 + 98.4925*mZZ + -0.70193*TMath::Power(mZZ,2) + 0.00164952*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = 0.0549555 + -0.000367875*mZZ + 6.16188e-10*(mZZ - 106.028)*(mZZ - 106.028)*(mZZ - 106.028)*(mZZ - 106.028);
		a2 = -0.000222052 + 1.51606e-06*mZZ + -9.57278e-08*(mZZ - 144.667)*(mZZ - 144.667);
		f = 4.82988e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);	
		
		a0_bkgd = 138.962;
		a1_bkgd = 30.6644;
		a2_bkgd = 102.09;
		a3_bkgd = 0.0416272;
		a4_bkgd = 178.464;
		a5_bkgd = 12.6718;
		a6_bkgd = 38.1437;
		a7_bkgd = 0.13019;
		a8_bkgd = 38.1146;
		a9_bkgd = 0.0277349;		
	}
	else if (channelVal == 3.){
		mZ = 91.2;
		Gamma0 = 36.9986 + -0.250136*mZZ + 0.00474031*(mZZ - 108.732)*(mZZ - 108.732);
		m0 = -3266.9 + 73.5969*mZZ + -0.544448*TMath::Power(mZZ,2) + 0.00133127*TMath::Power(mZZ,3);
		a0 = 0.00205023 + 1.49685e-05*mZZ + 2.45743e-07*(mZZ - 110.553)*(mZZ - 110.553);
		a1 = 0.0206819 + -0.000168968*mZZ + 4.26934e-11*(mZZ - 40.4456)*(mZZ - 40.4456)*(mZZ - 40.4456)*(mZZ - 40.4456);
		a2 = -1.271e-05 + 8.742e-08*mZZ + -6.34276e-08*(mZZ - 144.032)*(mZZ - 144.032);
		f = 4.82988e-01*(TMath::Erf( (mZZ-1.47620e+02)/8.12397e-01 ) + 1.);
		Gamma =  252.406 + -5.11437*mZZ + 0.0345238*TMath::Power(mZZ,2) + -7.44284e-05*TMath::Power(mZZ,3);
		f0 = 3.68735e-01*(TMath::Erf( (mZZ-1.17687e+02)/2.61381e-01 ) + 1.);	

		a0_bkgd = 138.962;
		a1_bkgd = 30.6644;
		a2_bkgd = 102.09;
		a3_bkgd = 0.0416272;
		a4_bkgd = 178.464;
		a5_bkgd = 12.6718;
		a6_bkgd = 38.1437;
		a7_bkgd = 0.13019;
		a8_bkgd = 38.1146;
		a9_bkgd = 0.0277349;
	}
	else{
		std::cout << "Invalid parameters for this PDF!" << std::endl;
	}
	
	
	//  ZZ shape
	double ZZshape = (.5+.5*TMath::Erf((mZZ-a0_bkgd)/a1_bkgd))*(a3_bkgd/(1+exp((mZZ-a0_bkgd)/a2_bkgd)))+(.5+.5*TMath::Erf((mZZ-a4_bkgd)/a5_bkgd))*(a7_bkgd/(1+exp((mZZ-a4_bkgd)/a6_bkgd))+a9_bkgd/(1+exp((mZZ-a4_bkgd)/a8_bkgd)));
	
	
	double mZDistribution, acceptance;
	
	//================= mZ Distribution =================
	double beta= (1-(mZstar-mZ)*(mZstar-mZ)/(mZZ*mZZ))*(1-(mZstar+mZ)*(mZstar+mZ)/(mZZ*mZZ));
	if(beta<0) return 0.0000001;
	else mZDistribution = beta;
	
	//================= acceptance (chebychev polynomials) 
	acceptance = a0+a1*mZstar+a2*(2*mZstar*mZstar-1)+a3*(4*mZstar*mZstar*mZstar-3*mZstar);
	
	//================= on-shell Z contamination ========
	
	double onShellZ = exp(-(mZstar-mZ)*(mZstar-mZ)/(2*Gamma*Gamma))/sqrt(2*3.1415*Gamma*Gamma);
	//double errorFunc = TMath::Erf((mZstar - m0)/Gamma);
	double errorFunc = exp(-(mZstar-m0)*(mZstar-m0)/(2*Gamma0*Gamma0))/sqrt(2*3.1415*Gamma0*Gamma0);
	
	//std::cout << mZZDistribution << " : " << mZDistribution << " : " << acceptance << std::endl;
	double final =  mZDistribution*((1-f)*(f0*acceptance+(1-f0)*errorFunc)+f*onShellZ);
	if (final <= 0) final = 1e-12;
	
	///*
	double integral;
	double E = 2.71828183;
	
	// ========================================================
	// integral
	///*
	double mZstarFix=12.;
	double integralT1Lo=((1 - f)*f0*((a1*Power(mZstarFix,6))/6. + (2*a2*Power(mZstarFix,7))/7. + (a0 - a2)*mZstarFix*Power(Power(mZ,2) - Power(mZZ,2),2) + (a1*Power(mZstarFix,2)*Power(Power(mZ,2) - Power(mZZ,2),2))/2. - (a1*Power(mZstarFix,4)*(Power(mZ,2) + Power(mZZ,2)))/2. + (Power(mZstarFix,5)*(a0 - a2*(1 + 4*Power(mZ,2) + 4*Power(mZZ,2))))/5. + (2*Power(mZstarFix,3)*(-(a0*(Power(mZ,2) + Power(mZZ,2))) + a2*(Power(mZ,4) + Power(mZZ,2) + Power(mZZ,4) + Power(mZ,2)*(1 - 2*Power(mZZ,2)))))/3.))/Power(mZZ,4);
	mZstarFix=mZZ-mZ;
	double integralT1Hi=((1 - f)*f0*((a1*Power(mZstarFix,6))/6. + (2*a2*Power(mZstarFix,7))/7. + (a0 - a2)*mZstarFix*Power(Power(mZ,2) - Power(mZZ,2),2) + (a1*Power(mZstarFix,2)*Power(Power(mZ,2) - Power(mZZ,2),2))/2. - (a1*Power(mZstarFix,4)*(Power(mZ,2) + Power(mZZ,2)))/2. + (Power(mZstarFix,5)*(a0 - a2*(1 + 4*Power(mZ,2) + 4*Power(mZZ,2))))/5. + (2*Power(mZstarFix,3)*(-(a0*(Power(mZ,2) + Power(mZZ,2))) + a2*(Power(mZ,4) + Power(mZZ,2) + Power(mZZ,4) + Power(mZ,2)*(1 - 2*Power(mZZ,2)))))/3.))/Power(mZZ,4);
	
	//std::cout << "IntegralHI: " << integralT1Hi << ", IntegralLO: " << integralT1Lo << std::endl;
	
	mZstarFix=12.;
	double t2loA = (3*Power(Gamma0,4) + Power(m0,4) + Power(mZ,4) - 2*Power(Gamma0,2)*Power(mZZ,2) + Power(mZZ,4) + Power(m0,2)*(6*Power(Gamma0,2) - 2*Power(mZ,2) - 2*Power(mZZ,2)) - 2*Power(mZ,2)*(Power(Gamma0,2) + Power(mZZ,2)));
	double integralT2Lo=((-1 + f)*(-1 + f0)*(-1.*((Power(Gamma0,2)*(Power(m0,3) + Power(m0,2)*mZstarFix + mZstarFix*(3*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) + m0*(5*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(m0 - mZstarFix,2)/(2.*Power(Gamma0,2)))) - Gamma0*t2loA*sqrt(TMath::Pi()/2.)*TMath::Erf((m0 - mZstarFix)/(sqrt(2)*Gamma0))))/Power(mZZ,4);
	integralT2Lo/=sqrt(2*TMath::Pi()*Gamma0*Gamma0);
	mZstarFix=mZZ-mZ;
	double t2hiA = (3*Power(Gamma0,4) + Power(m0,4) + Power(mZ,4) - 2*Power(Gamma0,2)*Power(mZZ,2) + Power(mZZ,4) + Power(m0,2)*(6*Power(Gamma0,2) - 2*Power(mZ,2) - 2*Power(mZZ,2)) - 2*Power(mZ,2)*(Power(Gamma0,2) + Power(mZZ,2)));
	double integralT2Hi=((-1 + f)*(-1 + f0)*(-1.*((Power(Gamma0,2)*(Power(m0,3) + Power(m0,2)*mZstarFix + mZstarFix*(3*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) + m0*(5*Power(Gamma0,2) - 2*Power(mZ,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(m0 - mZstarFix,2)/(2.*Power(Gamma0,2)))) - Gamma0*t2hiA*sqrt(TMath::Pi()/2.)*TMath::Erf((m0 - mZstarFix)/(sqrt(2)*Gamma0))))/Power(mZZ,4);
	integralT2Hi/=sqrt(2*TMath::Pi()*Gamma0*Gamma0);
	
	mZstarFix=12.;
	double t3loA = (3*Power(Gamma,4) + 4*Power(Gamma,2)*Power(mZ,2) - 2*(Power(Gamma,2) + 2*Power(mZ,2))*Power(mZZ,2) + Power(mZZ,4));
	double integralT3Lo=(f*((Power(Gamma,2)*(Power(mZ,3) + Power(mZ,2)*mZstarFix - mZstarFix*(3*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) - mZ*(5*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(mZ - mZstarFix,2)/(2.*Power(Gamma,2))) - Gamma*t3loA*sqrt(TMath::Pi()/2.)*TMath::Erf((mZ - mZstarFix)/(sqrt(2)*Gamma))))/Power(mZZ,4);
	integralT3Lo/=sqrt(2*TMath::Pi()*Gamma*Gamma);
	mZstarFix=mZZ-mZ;
	double t3hiA = (3*Power(Gamma,4) + 4*Power(Gamma,2)*Power(mZ,2) - 2*(Power(Gamma,2) + 2*Power(mZ,2))*Power(mZZ,2) + Power(mZZ,4));
	double integralT3Hi=(f*((Power(Gamma,2)*(Power(mZ,3) + Power(mZ,2)*mZstarFix - mZstarFix*(3*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2)) - mZ*(5*Power(Gamma,2) + Power(mZstarFix,2) - 2*Power(mZZ,2))))/Power(E,Power(mZ - mZstarFix,2)/(2.*Power(Gamma,2))) - Gamma*t3hiA*sqrt(TMath::Pi()/2.)*TMath::Erf((mZ - mZstarFix)/(sqrt(2)*Gamma))))/Power(mZZ,4);
	integralT3Hi/=sqrt(2*TMath::Pi()*Gamma*Gamma);
	
	integral = integralT1Hi + integralT2Hi + integralT3Hi - integralT1Lo - integralT2Lo - integralT3Lo;
	
	if (integral <= 0.) {
		std::cout << "integral is less than zero...." << std::endl;
		//integral = 100000000.;
	}
	
	return (final/integral)*ZZshape;
	//*/
	//return final;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// backgrounds above
// --------------------------------------------------------------------
// --------------------------------------------------------------------

// --------------------------------------------------------------------
// 2D signal

ClassImp(Roo4lMasses2D) 

Roo4lMasses2D::Roo4lMasses2D(){}

Roo4lMasses2D::Roo4lMasses2D(const char *name, const char *title, 
							 RooAbsReal& _mZstar,
							 RooAbsReal& _mZ,
							 RooAbsReal& _mZZ,
							 RooAbsReal& _Gamma,
							 RooAbsReal& _p0,
							 RooAbsReal& _p1,
							 RooAbsReal& _p2,
							 RooAbsReal& _CBmean,
							 RooAbsReal& _CBwidth,
							 RooAbsReal& _CBalpha,
							 RooAbsReal& _CBn
							 
							 ) :
RooAbsPdf(name,title), 
mZstar("mZstar","mZstar",this,_mZstar),
mZ("mZ","mZ",this,_mZ),
mZZ("mZZ","mZZ",this,_mZZ),
Gamma("Gamma","Gamma",this,_Gamma),
p0("p0","p0",this,_p0),
p1("p1","p1",this,_p1),
p2("p2","p2",this,_p2),
CBmean("CBmean","CBmean",this,_CBmean),
CBwidth("CBwidth","CBwidth",this,_CBwidth),
CBalpha("CBalpha","CBalpha",this,_CBalpha),
CBn("CBn","CBn",this,_CBn)

{ 
} 

Roo4lMasses2D::Roo4lMasses2D(const Roo4lMasses2D& other, const char* name) :  
RooAbsPdf(other,name), 
mZstar("mZstar",this,other.mZstar),
mZ("mZ",this,other.mZ),
mZZ("mZZ",this,other.mZZ),
Gamma("Gamma",this,other.Gamma),
p0("p0",this,other.p0),
p1("p1",this,other.p1),
p2("p2",this,other.p2),
CBmean("CBmean"  ,this,other.CBmean),
CBwidth("CBwidth",this,other.CBwidth),
CBalpha("CBalpha",this,other.CBalpha),
CBn("CBn",this,other.CBn)

{ 
} 

double Roo4lMasses2D::evaluate() const 
{ 
	double mZstarDistribution;
	double Numerator,Denominator,beta;
	beta = (1-(mZstar-mZ)*(mZstar-mZ)/(mZZ*mZZ))*(1-(mZstar+mZ)*(mZstar+mZ)/(mZZ*mZZ));
	
	if(beta<0.0) return 1e-12;
	
	///*
	//Numerator = mZstar*mZstar*sqrt(beta)*(2+(mZstar*mZstar-mZZ*mZZ+mZ*mZ)*(mZstar*mZstar-mZZ*mZZ+mZ*mZ)/(4*mZstar*mZstar*mZ*mZ));
	Numerator = mZstar*mZstar*sqrt(beta)*(2+(mZstar*mZstar-mZZ*mZZ+mZ*mZ)*(mZstar*mZstar-mZZ*mZZ+mZ*mZ)/(4*mZstar*mZstar*mZ*mZ));
	Denominator = 2*(mZstar*mZstar-mZ*mZ)*(mZstar*mZstar-mZ*mZ)-2*mZ*mZ*Gamma*Gamma;
	
	//mZstarDistribution = Numerator/Denominator;
	mZstarDistribution = Numerator/Denominator;
	//*/
	//mZstarDistribution = 1;
	
	double acceptance = p0+p1*mZstar+p2*mZstar*mZstar;
	//double acceptance = 1;
	
	
	double mZZDistribution;
	/*
	 Double_t arg= mZZ - CBmean;  
	 Double_t ret =exp(-0.5*arg*arg/(CBwidth*CBwidth)) ;
	 //cout << "gauss x = " << x << " mean = " << mean << " sigma = " << sigma << endl ;
	 mZZDistribution = ret ;
	 //*/
	///*
	//// CB Shape
	double t = (mZZ-CBmean)/CBwidth;
	//std:cout << "CBwidth: " << CBwidth << " CBmean: " << CBmean << " mZZ: " << mZZ << std::endl;
	if (CBalpha < 0) t = -t;
	
	double absAlpha = fabs((double)CBalpha);
	
	//std::cout << "t = " << t << std::endl;
	if (t >= -absAlpha) {
		
		mZZDistribution = exp(-0.5*t*t);
		//std::cout << "Gaussian!! " << -.5*t*t << " : " << mZZDistribution << std::endl;
		
	}
	else {
		double a =  TMath::Power(CBn/absAlpha,CBn)*exp(-0.5*absAlpha*absAlpha);
		double b = CBn/absAlpha - absAlpha;
		
		//std::cout << "power!!! a: " << a << " b: " << b << endl;
		
		mZZDistribution = a/TMath::Power(b - t, CBn);
		
	}
	//*/
	
	
	//std::cout << "t: " << t << " absAlpha: " << absAlpha << std::endl;
	
	//double final = mZstarDistribution*acceptance;
	double final = mZstarDistribution*acceptance*mZZDistribution;
	if (final <= 0) final = 1e-12;
	
	return final;
}



// --------------------------------------------------------------------


/************RooFourMuMassShapePdf2***********/

ClassImp(RooFourMuMassShapePdf2)


RooFourMuMassShapePdf2::RooFourMuMassShapePdf2(const char *name, const char *title,
											   RooAbsReal& _m4l,
											   RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}


RooFourMuMassShapePdf2::RooFourMuMassShapePdf2(const RooFourMuMassShapePdf2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}

/*
 void RooFourMuMassShapePdf2::readFile() const{
 
 
 ifstream file;
 
 file.open("txtFiles/HiggsBR_7TeV_Official.txt");
 for(int k = 0; k < 217; k++){
 
 file >> scratchMass >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> BR[8][k] >> BR[9][k]
 >> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k] >> BR[18][k] >> BR[19][k] >> BR[20][k]
 >> BR[21][k] >> BR[22][k] >> BR[23][k] >> BR[24][k] >> BR[25][k];
 
 
 }
 file.close();
 
 }
 */

Double_t RooFourMuMassShapePdf2::evaluate() const
{
	
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t x = m4l;
	
	//min and max for integration                                                                                                                            
	//Start the tail at x/.998 to control area around ratio = 1                                                                                              
	Double_t min = x/.998;
	Double_t max = 0;
	Double_t FracEventsNoBrem = FracEventsNoBrem_4mu;
	
	Double_t higgsTotWidth = HiggsWidth(0,mHreq);
	//std::cout << "higgsTotWidth: " << higgsTotWidth << std::endl;
	
	if( 1000 < (mH + 5*higgsTotWidth) ){ max = 1000;}
	else{ max = mH + 5*higgsTotWidth;}
	
	Double_t dmTil = (max-min)/1000;
	
	Double_t sum = 0;
	Double_t pdf_1_NoBrem = FracEventsNoBrem*pdf1(x,mHreq);
	
	for(Double_t mTil = min; mTil < max; mTil += dmTil)
    {
		
		Double_t TC = (dmTil/mTil) * rho(x/mTil,"4mu") * pdf1(mTil,mHreq);
		sum = sum + TC;
		
    }
	
	Double_t pdf_2 = sum + pdf_1_NoBrem;
	
	
	return pdf_2;
	
}


/************RooFourEMassShapePdf2***********/

ClassImp(RooFourEMassShapePdf2)


RooFourEMassShapePdf2::RooFourEMassShapePdf2(const char *name, const char *title,
											 RooAbsReal& _m4l,
											 RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}

RooFourEMassShapePdf2::RooFourEMassShapePdf2(const RooFourEMassShapePdf2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}
/*
 void RooFourEMassShapePdf2::readFile() const{
 
 
 ifstream file;
 
 file.open("txtFiles/HiggsBR_7TeV_Official.txt");
 for(int k = 0; k < 217; k++){
 
 file >> scratchMass >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> BR[8][k] >> BR[9][k]
 >> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k] >> BR[18][k] >> BR[19][k] >> BR[20][k]
 >> BR[21][k] >> BR[22][k] >> BR[23][k] >> BR[24][k] >> BR[25][k];
 
 
 }
 file.close();
 
 }
 */

Double_t RooFourEMassShapePdf2::evaluate() const
{
	using namespace RooFit;
	
 	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t x = m4l;
	
	//min and max for integration                                                                                                                            
	//Start the tail at x/.998 to control area around ratio = 1                                                                      
	
    Double_t min = x/.998;
    Double_t max = 0;
    Double_t FracEventsNoBrem = FracEventsNoBrem_4e;
	
    Double_t higgsTotWidth = HiggsWidth(0,mHreq);
	
    if( 1000 < (mH + 5*higgsTotWidth) ){ max = 1000;}
    else{ max = mH + 5*higgsTotWidth;}
	
    Double_t dmTil = (max-min)/1000;
	
    Double_t sum = 0;
    Double_t pdf_1_NoBrem = FracEventsNoBrem*pdf1(x,mHreq);
	
    for(Double_t mTil = min; mTil < max; mTil += dmTil)
	{
		
		Double_t TC = (dmTil/mTil) * rho(x/mTil,"4e") * pdf1(mTil,mHreq);
		sum = sum + TC;
		
	}
	
    Double_t pdf_2 = sum + pdf_1_NoBrem;
	
    return pdf_2;
	
}


/************RooTwoETwoMuMassShapePdf2*************/

ClassImp(RooTwoETwoMuMassShapePdf2)


RooTwoETwoMuMassShapePdf2::RooTwoETwoMuMassShapePdf2(const char *name, const char *title,
													 RooAbsReal& _m4l,
													 RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}

RooTwoETwoMuMassShapePdf2::RooTwoETwoMuMassShapePdf2(const RooTwoETwoMuMassShapePdf2& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}

/*
 void RooTwoETwoMuMassShapePdf2::readFile() const{
 
 
 ifstream file;
 
 file.open("txtFiles/HiggsBR_7TeV_Official.txt");
 for(int k = 0; k < 217; k++){
 
 file >> scratchMass >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> BR[8][k] >> BR[9][k]
 >> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k] >> BR[18][k] >> BR[19][k] >> BR[20][k]
 >> BR[21][k] >> BR[22][k] >> BR[23][k] >> BR[24][k] >> BR[25][k];
 
 
 }
 file.close();
 
 }
 */

Double_t RooTwoETwoMuMassShapePdf2::evaluate() const
{
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t x = m4l;
	
	//min and max for integration                                                                                                                            
	//Start the tail at x/.998 to control area around ratio = 1                                                                                              
	Double_t min = x/.998;
	Double_t max = 0;
	Double_t FracEventsNoBrem = FracEventsNoBrem_2e2mu;
	
	Double_t higgsTotWidth = HiggsWidth(0,mHreq);
	
	if( 1000 < (mH + 5*higgsTotWidth) ){ max = 1000;}
	else{ max = mH + 5*higgsTotWidth;}
	
	Double_t dmTil = (max-min)/1000;
	
	Double_t sum = 0;
	Double_t pdf_1_NoBrem = FracEventsNoBrem*pdf1(x,mHreq);
	
	for(Double_t mTil = min; mTil < max; mTil += dmTil)
    {
		
		Double_t TC = (dmTil/mTil) * rho(x/mTil,"2e2mu") * pdf1(mTil,mHreq);
		sum = sum + TC;
		
    }
	
	Double_t pdf_2 = sum + pdf_1_NoBrem;
	
	return pdf_2;
	
}

/************RooRelBWUF*************/

ClassImp(RooRelBWUF)


RooRelBWUF::RooRelBWUF(const char *name, const char *title,
													 RooAbsReal& _m4l,
													 RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}

RooRelBWUF::RooRelBWUF(const RooRelBWUF& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}

Double_t RooRelBWUF::evaluate() const
{
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t x = m4l;
	
	/*
	Double_t Gamma_gg = HiggsWidth(7,x);
	Double_t Gamma_ZZ = HiggsWidth(11,x);
	Double_t Gamma_TOT = HiggsWidth(0,x);
		
	Double_t pdf = Gamma_gg*Gamma_ZZ/( pow(pow(x,2)-pow(mHreq,2),2 )+pow(x,2)*pow(Gamma_TOT,2) );
	
	return pdf;
	*/
	///*
	Double_t pdf_1_NoBrem = pdf1(x,mHreq);
	return pdf_1_NoBrem;
	//*/
}



/************RooRelBWUF_SM4*************/

ClassImp(RooRelBWUF_SM4)


RooRelBWUF_SM4::RooRelBWUF_SM4(const char *name, const char *title,
			       RooAbsReal& _m4l,
			       RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}

RooRelBWUF_SM4::RooRelBWUF_SM4(const RooRelBWUF_SM4& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}

Double_t RooRelBWUF_SM4::evaluate() const
{
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t x = m4l;
	

	Double_t pdf_1_NoBrem = pdf1SM4(x,mHreq);
	return pdf_1_NoBrem;

}

/************RooRelBWUFParamWidth*************/

ClassImp(RooRelBWUFParamWidth)


RooRelBWUFParamWidth::RooRelBWUFParamWidth(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _mH,
					   RooAbsReal& _width) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH),
width("width","width",this,_width)
{
}

RooRelBWUFParamWidth::RooRelBWUFParamWidth(const RooRelBWUFParamWidth& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH),
width("width",this,other.width)
{
}

Double_t RooRelBWUFParamWidth::evaluate() const
{
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t mStar = m4l;
	Double_t x = width;
	
	// use non-relativistic Breit Wigner 
	Double_t dm = mStar - mHreq;
	Double_t pdf = 1./(dm*dm + 0.25*x*x);  // x is full width 
	return pdf;
}


/************RooRelBWUFParam*************/

ClassImp(RooRelBWUFParam)


RooRelBWUFParam::RooRelBWUFParam(const char *name, const char *title,
				 RooAbsReal& _m4l,
				 RooAbsReal& _mH,
				 RooAbsReal& _scaleParam,
				 Double_t _widthSF) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH),
scaleParam("scaleParam","scaleParam",this,_scaleParam),
widthSF(_widthSF)
{
}

RooRelBWUFParam::RooRelBWUFParam(const RooRelBWUFParam& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH),
scaleParam("scaleParam",this,other.scaleParam),
widthSF(other.widthSF)
{
}

Double_t RooRelBWUFParam::evaluate() const
{
	using namespace RooFit;
	
	if( BR[0][0] == 0 ){ readFile(); }
	
	
	Double_t mHreq = mH;
	Double_t mStar = m4l;
	Double_t x = scaleParam;
	
	
	/*
	 Double_t Gamma_gg = HiggsWidth(7,x);
	 Double_t Gamma_ZZ = HiggsWidth(11,x);
	 Double_t Gamma_TOT = HiggsWidth(0,x);
	 
	 Double_t pdf = Gamma_gg*Gamma_ZZ/( pow(pow(x,2)-pow(mHreq,2),2 )+pow(x,2)*pow(Gamma_TOT,2) );
	 
	 return pdf;
	 */
	///*
	Double_t pdf_1_NoBrem = pdf1Param(mStar,mHreq,x,widthSF);
	return pdf_1_NoBrem;
	//*/
}

/************RooRelBWHighMass*************/

ClassImp(RooRelBWHighMass)


RooRelBWHighMass::RooRelBWHighMass(const char *name, const char *title,
					   RooAbsReal& _m4l,
					   RooAbsReal& _mH,
					   RooAbsReal& _gamma) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH),
gamma("gamma","gamma",this,_gamma)
{
}

RooRelBWHighMass::RooRelBWHighMass(const RooRelBWHighMass& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH),
gamma("gamma",this,other.gamma)
{
}

Double_t RooRelBWHighMass::evaluate() const
{
	using namespace RooFit;

	Double_t pdf = m4l/( pow( pow(m4l,2) - pow(mH,2),2) + pow(m4l,2)*pow(gamma,2)  );
	return pdf;
}




/***********************RooFourMuMassRes*************/

ClassImp(RooFourMuMassRes)

RooFourMuMassRes::RooFourMuMassRes(const char *name, const char *title,
								   RooAbsReal& _m4l,
								   RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}


RooFourMuMassRes::RooFourMuMassRes(const RooFourMuMassRes& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}


Double_t RooFourMuMassRes::evaluate() const
{
	
	Double_t x = m4l;
	Double_t mHreq = mH;
	Double_t mean_ = RooFit::mean(mHreq,"4mu");
	Double_t alpha1 = -1.3;
	Double_t alpha2 = 1.3;
	Double_t sigma_ = RooFit::sigma_CB(mHreq,"4mu");
	Double_t n1 = 60;
	Double_t n2 = 60;
	Double_t N = RooFit::N(mHreq,"4mu");
	
	Double_t t = (x-mean_)/sigma_;
	if (alpha1 < 0) t = -t;
	
	Double_t absAlpha = TMath::Abs(alpha1);
	Double_t absAlpha2 = TMath::Abs(alpha2);
	
	if (t >= alpha1 && t <= alpha2) {
		return N*exp(-0.5*t*t);
	}
	else if ( t > alpha2) {
		Double_t a1 =  pow(n2/absAlpha2,n2)*exp(-0.5*absAlpha2*absAlpha2);
		Double_t b1 = n2/absAlpha2 - absAlpha2;
		
		return N*(a1/pow(b1 + t, n2));
		
	}
	else {
		Double_t a2 =  pow(n1/absAlpha,n1)*exp(-0.5*absAlpha*absAlpha);
		Double_t b2 = n1/absAlpha - absAlpha;
		
		return N*(a2/pow(b2 - t, n1));
		
	}
	
}


/***********RooFourEMassRes**************/

ClassImp(RooFourEMassRes)

RooFourEMassRes::RooFourEMassRes(const char *name, const char *title,
								 RooAbsReal& _m4l,
								 RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}


RooFourEMassRes::RooFourEMassRes(const RooFourEMassRes& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}

Double_t RooFourEMassRes::evaluate() const
{
	
	Double_t x = m4l;
	Double_t mHreq = mH;
	Double_t mean_ = RooFit::mean(mHreq,"4e");
	Double_t alpha1 = -2.8;
	Double_t alpha2 = 0.7;
	Double_t sigma_ = RooFit::sigma_CB(mHreq,"4e");
	Double_t n1 = 10;
	Double_t n2 = 132;
	Double_t N = RooFit::N(mHreq,"4e");
	
	Double_t t = (x-mean_)/sigma_;
	if (alpha1 < 0) t = -t;
	
	Double_t absAlpha = TMath::Abs(alpha1);
	Double_t absAlpha2 = TMath::Abs(alpha2);
	
	if (t >= alpha1 && t <= alpha2) {
		return N*exp(-0.5*t*t);
	}
	else if ( t > alpha2) {
		Double_t a1 =  pow(n2/absAlpha2,n2)*exp(-0.5*absAlpha2*absAlpha2);
		Double_t b1 = n2/absAlpha2 - absAlpha2;
		
		return N*(a1/pow(b1 + t, n2));
		
	}
	else {
		Double_t a2 =  pow(n1/absAlpha,n1)*exp(-0.5*absAlpha*absAlpha);
		Double_t b2 = n1/absAlpha - absAlpha;
		
		return N*(a2/pow(b2 - t, n1));
	}
	
}


/************RooTwoETwoMuMassRes************/

ClassImp(RooTwoETwoMuMassRes)

RooTwoETwoMuMassRes::RooTwoETwoMuMassRes(const char *name, const char *title,
										 RooAbsReal& _m4l,
										 RooAbsReal& _mH) :
RooAbsPdf(name,title),
m4l("m4l","m4l",this,_m4l),
mH("mH","mH",this,_mH)
{
}


RooTwoETwoMuMassRes::RooTwoETwoMuMassRes(const RooTwoETwoMuMassRes& other, const char* name) :
RooAbsPdf(other,name),
m4l("m4l",this,other.m4l),
mH("mH",this,other.mH)
{
}


Double_t RooTwoETwoMuMassRes::evaluate() const
{
	
	Double_t x = m4l;
	Double_t mHreq = mH;
	Double_t mean_ = RooFit::mean(mHreq,"2e2mu");
	Double_t alpha1 = -1.4;
	Double_t alpha2 = 0.85;
	Double_t sigma_ = RooFit::sigma_CB(mHreq,"2e2mu");
	Double_t n1 = 70;
	Double_t n2 = 15;
	Double_t N = RooFit::N(mHreq,"2e2mu");
	
	Double_t t = (x-mean_)/sigma_;
	if (alpha1 < 0) t = -t;
	
	Double_t absAlpha = TMath::Abs(alpha1);
	Double_t absAlpha2 = TMath::Abs(alpha2);
	
	if (t >= alpha1 && t <= alpha2) {
		return N*exp(-0.5*t*t);
	}
	else if ( t > alpha2) {
		Double_t a1 =  pow(n2/absAlpha2,n2)*exp(-0.5*absAlpha2*absAlpha2);
		Double_t b1 = n2/absAlpha2 - absAlpha2;
		
		return N*(a1/pow(b1 + t, n2));
		
	}
	else {
		Double_t a2 =  pow(n1/absAlpha,n1)*exp(-0.5*absAlpha*absAlpha);
		Double_t b2 = n1/absAlpha - absAlpha;
		
		return N*(a2/pow(b2 - t, n1));
		
	}
	
}




/*************RooRelBW1***********/

ClassImp(RooRelBW1)

RooRelBW1::RooRelBW1(const char *name, const char *title,
					 RooAbsReal& _m,
					 RooAbsReal& _mean,
					 RooAbsReal& _gamma) :
RooAbsPdf(name,title),
m("m","m",this,_m),
mean("mean","mean",this,_mean),
gamma("gamma","gamma",this,_gamma)
{
}


RooRelBW1::RooRelBW1(const RooRelBW1& other, const char* name) :
RooAbsPdf(other,name),
m("m",this,other.m),
mean("mean",this,other.mean),
gamma("gamma",this,other.gamma)
{
}


Double_t RooRelBW1::evaluate() const
{
	
	Double_t arg= m*m - mean*mean;
	Double_t relBW =  (mean*gamma) / (arg*arg + mean*mean*gamma*gamma);
	return relBW;
}



/****************RooTsallis******************/

ClassImp(RooTsallis) 

RooTsallis::RooTsallis() {}

 RooTsallis::RooTsallis(const char *name, const char *title, 
				  RooAbsReal& _x,
				  RooAbsReal& _m,
				  RooAbsReal& _n,
                                  RooAbsReal& _n2,
                                  RooAbsReal& _bb,
			          RooAbsReal& _bb2,
			          RooAbsReal& _T,
                                  RooAbsReal& _fexp):
   
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   m("m","m",this,_m),
   n("n","n",this,_n),
   n2("n2","n2",this,_n2),
   bb("bb","bb",this,_bb),
   bb2("bb2","bb2",this,_bb2),
   T("T","T",this,_T),
   fexp("fexp","fexp",this,_fexp)
 { 
 } 


 RooTsallis::RooTsallis(const RooTsallis& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   m("m",this,other.m),
   n("n",this,other.n),
   n2("n2",this,other.n2),
   bb("bb",this,other.bb),
   bb2("bb2",this,other.bb2),
   T("T",this,other.T),
   fexp("fexp",this,other.fexp)
 {
 } 



 double RooTsallis::evaluate() const 
 { 
   // cout<<"In rooTsallis3::evaluate()"<<endl;
   return pow(x,n2)*exp(-bb*x)*pow(1 + (sqrt(x*x + m*m) - m)/(n*T),-n) + fexp*exp(-bb2*x);
 } 



/*************RooaDoubleCBxBW***********/

ClassImp(RooaDoubleCBxBW)

RooaDoubleCBxBW::RooaDoubleCBxBW(){}

RooaDoubleCBxBW::RooaDoubleCBxBW(const char *name, const char *title,
         RooAbsReal& _x,
         RooAbsReal& _shift,
         RooAbsReal& _sigma,
         RooAbsReal& _alphaL,
         RooAbsReal& _alphaR,
         RooAbsReal& _mean,
         RooAbsReal& _width,
         unsigned _nL,
         unsigned _nR,
         RooAbsReal& _thetaL,
         RooAbsReal& _thetaR,
         bool computeActualCB_=true
         ) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  shift("shift","shift",this,_shift),
  sigma("sigma","sigma",this,_sigma),
  alphaL("alphaL","alphaL",this,_alphaL),
  alphaR("alphaR","alphaR",this,_alphaR),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  nL(_nL),
  nR(_nR),
  thetaL("thetaL","thetaL",this,_thetaL),
  thetaR("thetaR","thetaR",this,_thetaR),
  computeActualCB(computeActualCB_)
{
}

RooaDoubleCBxBW::RooaDoubleCBxBW(const RooaDoubleCBxBW& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  shift("shift",this,other.shift),
  sigma("sigma",this,other.sigma),
  alphaL("alphaL",this,other.alphaL),
  alphaR("alphaR",this,other.alphaR),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  nL(other.nL),
  nR(other.nR),
  thetaL("thetaL",this,other.thetaL),
  thetaR("thetaR",this,other.thetaR),
  computeActualCB(other.computeActualCB)
{
}

double RooaDoubleCBxBW::evaluateDoubleCB() const
{

  double al = alphaL;
  double absaL = fabs((double)al);
  double AL = TMath::Power(nL/absaL,(int)nL)*exp(-0.5*absaL*absaL);
  double BL = nL/absaL-absaL;
  double t0L = -absaL;
  double tsL = (-absaL-thetaL);
  double c1L = ( (AL/TMath::Power(BL-tsL,(int)nL+1))*nL*(tsL-t0L) - AL/TMath::Power(BL-tsL,(int)nL) ) / ((tsL-t0L)*(tsL-t0L));
  double c2L = ( AL/TMath::Power(BL-tsL,(int)nL) - c1L*(tsL-t0L)*(tsL-t0L) ) / (tsL-t0L);

  double ar = alphaR;
  double absaR = fabs((double)ar);
  double AR = TMath::Power(nR/absaR,(int)nR)*exp(-0.5*absaR*absaR);
  double BR = nR/absaR-absaR;
  double t0R = absaR;
  double tsR = (absaR+thetaR);
  double c1R = ( -(AR/TMath::Power(BR+tsR,(int)nR+1))*nR*(tsR-t0R) - AR/TMath::Power(BR+tsR,(int)nR) ) / ((tsR-t0R)*(tsR-t0R));
  double c2R = ( AR/TMath::Power(BR+tsR,(int)nR) - c1R*(tsR-t0R)*(tsR-t0R) ) / (tsR-t0R);

  double t = (x-(mean+shift))/sigma;

  double gaussval = exp(-0.5*t*t);
  double polval = 0.0;

  if (!computeActualCB) {
    if(t >= -absaL && t <= absaR) polval = 0;
    else if (t >= tsL && t < -absaL) polval = c1L*(t-t0L)*(t-t0L) + c2L*(t-t0L);
    else if (t <  tsR && t >= absaR) polval = c1R*(t-t0R)*(t-t0R) + c2R*(t-t0R);
    else if (t <  tsL) polval = AL/TMath::Power(BL-t,(int)nL);
    else if (t >  tsR) polval = AR/TMath::Power(BR+t,(int)nR);
  }
  else {
    if(t >= -absaL && t <= absaR) polval = 0;
    else if (t < -absaL) polval = AL/TMath::Power(BL-t,(int)nL) - exp(-0.5*t*t);
    else polval = AR/TMath::Power(BR+t,(int)nR) - exp(-0.5*t*t);
  }

  return gaussval + polval;
}

double RooaDoubleCBxBW::evaluatePowerLaw(double lim, unsigned power, bool isLeft) const
{
  double al = alphaL;
  double absaL = fabs((double)al);
  double BL = nL/absaL-absaL;

  double ar = alphaR;
  double absaR = fabs((double)ar);
  double BR = nR/absaR-absaR;

  double d1 = -mean/width;
  double d2 = (isLeft ? BL*sigma + shift - x : -BR*sigma + shift - x)/width;
  double expr1 = (d1-d2)*(d1-d2)  +  1.0;

  if      (power ==  0) return -(atan2(d1 - lim, 1.));

  else if (power ==  1) return ((2*d2*atan2((d1 - lim),1.) + 2*d1*atan2((-d1 + lim),1.) + 2*log(fabs(-d2 + lim)) - log(d1*d1 + 1.0 - 2*d1*lim + lim*lim))/ (2*expr1));

  else return ( -evaluatePowerLaw(lim, power-2, isLeft) + 2*(d1-d2)*evaluatePowerLaw(lim, power-1, isLeft) - 1.0/((power-1)*pow(lim-d2,power-1)) )  /  expr1;
}

double RooaDoubleCBxBW::evaluateVoigtian() const
{
  double arg = x - mean - shift;
  double c = 1./(sqrt(2.)*sigma);
  double a = 0.5*c*width;
  double u = c*arg;
  std::complex<double> z(u,a);
  std::complex<double> v = RooMath::erfc_fast(z);
  return atan2(0.,-1.)*v.real()/width;
}

double RooaDoubleCBxBW::evaluateQuadratic(double lim1, double lim2, bool isLeft) const
{

  double al = alphaL;
  double absaL = fabs((double)al);
  double AL = TMath::Power(nL/absaL,(int)nL)*exp(-0.5*absaL*absaL);
  double BL = nL/absaL-absaL;
  double t0L = -absaL;
  double tsL = (-absaL-thetaL);
  double c1L = ( (AL/TMath::Power(BL-tsL,(int)nL+1))*nL*(tsL-t0L) - AL/TMath::Power(BL-tsL,(int)nL) ) / ((tsL-t0L)*(tsL-t0L));
  double c2L = ( AL/TMath::Power(BL-tsL,(int)nL) - c1L*(tsL-t0L)*(tsL-t0L) ) / (tsL-t0L);

  double ar = alphaR;
  double absaR = fabs((double)ar);
  double AR = TMath::Power(nR/absaR,(int)nR)*exp(-0.5*absaR*absaR);
  double BR = nR/absaR-absaR;
  double t0R = absaR;
  double tsR = (absaR+thetaR);
  double c1R = ( -(AR/TMath::Power(BR+tsR,(int)nR+1))*nR*(tsR-t0R) - AR/TMath::Power(BR+tsR,(int)nR) ) / ((tsR-t0R)*(tsR-t0R));
  double c2R = ( AR/TMath::Power(BR+tsR,(int)nR) - c1R*(tsR-t0R)*(tsR-t0R) ) / (tsR-t0R);

  double c1 = c1R;
  double c2 = c2R;
  double t0 = t0R;
  if (isLeft) {
    c1 = c1L;
    c2 = c2L;
    t0 = t0L;
  }

  double sigma2 = sigma*sigma;
  double mean2 = mean*mean;
  double shift2 = shift*shift;
  double width2 = width*width;
  double t02 = t0*t0;
  double x2 = x*x;

  double val1 = (c1*(lim1 + mean) + ((-(c2*sigma*(shift + sigma*t0 + mean - x)) + c1*(shift2 + sigma2*t02 - width2 + mean2 + 2*sigma*t0*(mean - x) + 2*shift*(sigma*t0 + mean - x) - 2*mean*x + x2))*atan2((lim1 + mean),width))/width + ((c2*sigma - 2*c1*(shift + sigma*t0 + mean - x))* log(width2 + (lim1 + mean)*(lim1 + mean)))/2)/sigma2;

  double val2 = (c1*(lim2 + mean) + ((-(c2*sigma*(shift + sigma*t0 + mean - x)) + c1*(shift2 + sigma2*t02 - width2 + mean2 + 2*sigma*t0*(mean - x) + 2*shift*(sigma*t0 + mean - x) - 2*mean*x + x2))*atan2((lim2 + mean),width))/width + ((c2*sigma - 2*c1*(shift + sigma*t0 + mean - x))* log(width2 + (lim2 + mean)*(lim2 + mean)))/2)/sigma2;


  return val2 - val1;
}

double RooaDoubleCBxBW::evaluate() const
{
  double al = alphaL;
  double absaL = fabs((double)al);
  double AL = TMath::Power(nL/absaL,(int)nL)*exp(-0.5*absaL*absaL);
  double tsL = (-absaL-thetaL);

  double ar = alphaR;
  double absaR = fabs((double)ar);
  double AR = TMath::Power(nR/absaR,(int)nR)*exp(-0.5*absaR*absaR);
  double tsR = (absaR+thetaR);

  double inf = 1000000.0;

  double point1 = tsL*(sigma)+shift-x;
  double point2 = (-absaL)*(sigma)+shift-x;
  double point3 = ( absaR)*(sigma)+shift-x;
  double point4 = tsR*(sigma)+shift-x;

  if (width == 0.0) return evaluateDoubleCB();
  else return std::max(0., (AL*pow(-sigma/width,nL)/width)*(evaluatePowerLaw(point1/width, nL, true) - evaluatePowerLaw(-inf, nL, true))) +  evaluateQuadratic(point1, point2, true) + evaluateVoigtian() + evaluateQuadratic(point3, point4, false) + std::max(0., (AR*pow(sigma/width, nR)/width)*(evaluatePowerLaw(inf, nR, false) - evaluatePowerLaw(point4/width, nR, false)));

}



/*************RooCPSHighMassGGH***********/
ClassImp(RooCPSHighMassGGH) 

 RooCPSHighMassGGH::RooCPSHighMassGGH(){}

 RooCPSHighMassGGH::RooCPSHighMassGGH(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _KPrime,
			RooAbsReal& _BRnew,
			RooAbsReal& _IntStr,
			RooAbsReal& _WidthScl,
			Bool_t is8TeV_b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   KPrime("KPrime","KPrime",this,_KPrime),
   BRnew("BRnew","BRnew",this,_BRnew),
   IntStr("IntStr","IntStr",this,_IntStr),
   WidthScl("WidthScl","WidthScl",this,_WidthScl),
   is8TeV(is8TeV_b)
 {    
   if (KPrime/(1-BRnew)<0.1) std::cout << "[RooCPSHighMassGGH] WARNING: parameters of the shape will be extrapolated outside the grid, (KPrime/(1-BRnew)=" << KPrime/(1-BRnew) << std::endl;
   initMatrices();
 } 


 RooCPSHighMassGGH::RooCPSHighMassGGH(const RooCPSHighMassGGH& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   KPrime("KPrime",this,other.KPrime),
   BRnew("BRnew",this,other.BRnew),
   IntStr("IntStr",this,other.IntStr),
   WidthScl("WidthScl",this,other.WidthScl),
   is8TeV(other.is8TeV)
 { 
   
   for (Int_t i=0; i<7; ++i){
     for (Int_t j=0; j<5; ++j){
       a_width[i][j] = other.a_width[i][j];
       a_delta[i][j] = other.a_delta[i][j];
       a_alpha[i][j] = other.a_alpha[i][j];
       a_r[i][j] = other.a_r[i][j];
       a_beta[i][j] = other.a_beta[i][j];
     }
   }

 } 



 Double_t RooCPSHighMassGGH::evaluate() const 
 {    

   Double_t effKPrime = KPrime/(1-BRnew);

   Double_t width;
   Double_t delta;
   Double_t alpha;
   Double_t r;
   Double_t beta ;

   if (KPrime/(1-BRnew)<0.1){
     width = WidthScl*getWidth(mH,0.1);
     delta = getDelta(mH,0.1);
     alpha = getAlpha(mH,0.1);
     r = getR(mH,0.1);
     beta = getBeta(mH,0.1);
   } else {
     width = WidthScl*getWidth(mH,effKPrime);
     delta = getDelta(mH,effKPrime);
     alpha = getAlpha(mH,effKPrime);
     r = getR(mH,effKPrime);
     beta = getBeta(mH,effKPrime);
   }

   // Cast the "x" proxy to double only once to avoid overhead and ambiguous
   // overloads the TComplex operators.
   const double xVal = x;

   Double_t bwHM = xVal / ( TMath::Power( TMath::Power(xVal,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(xVal,2)*TMath::Power(effKPrime*width,2) ); 

   Double_t splineFactor;
   if (xVal<1850) splineFactor = Spline(xVal);
   else splineFactor = Spline(1850);

   Double_t signal = bwHM*splineFactor;
     
   Double_t k=0.25;
   Double_t mH_eff = mH*TMath::Sqrt(1-k*(effKPrime*width/mH)*(effKPrime*width/mH));
   TComplex MSquared(mH_eff*mH_eff,-effKPrime*width*mH_eff);
   TComplex M = MSquared.Sqrt(MSquared);
   TComplex Exp1 = MSquared.Exp((TComplex)alpha);
   TComplex Exp2 = MSquared.Exp(-(TComplex)alpha);
   TComplex Exp3 = MSquared.Exp(-xVal*(TComplex)alpha/M);

   double interference = -r*(1-TMath::Exp(-beta*(xVal-150.)/mH_eff))*((Exp1/(xVal-M)-Exp2/(xVal+M))*Exp3).Re();

   Double_t fValue = signal + IntStr*interference/( TMath::Sqrt(1-BRnew) );
   if (fValue > 0) return fValue;
   else return 0.;   
   
 } 




void RooCPSHighMassGGH::initMatrices() {


  if (is8TeV){

    a_delta[0][0] = 0.015081577003;
    a_delta[0][1] = 0.0576512292027;
    a_delta[0][2] = 0.131081297994;
    a_delta[0][3] = 0.243380963802;
    a_delta[0][4] = 0.39500105381;
    a_delta[1][0] = 0.0972077324986;
    a_delta[1][1] = 0.434791147709;
    a_delta[1][2] = 0.861800074577;
    a_delta[1][3] = 1.58613061905;
    a_delta[1][4] = 2.41115784645;
    a_delta[2][0] = 0.242783889174;
    a_delta[2][1] = 0.827510118484;
    a_delta[2][2] = 1.78883588314;
    a_delta[2][3] = 3.17804527283;
    a_delta[2][4] = 4.9168419838;
    a_delta[3][0] = 0.368439763784;
    a_delta[3][1] = 1.45783913136;
    a_delta[3][2] = 3.27468585968;
    a_delta[3][3] = 5.93577003479;
    a_delta[3][4] = 9.49385356903;
    a_delta[4][0] = 0.679009974003;
    a_delta[4][1] = 2.96808123589;
    a_delta[4][2] = 6.60964488983;
    a_delta[4][3] = 11.7470684052;
    a_delta[4][4] = 18.1083869934;
    a_delta[5][0] = 1.51261091232;
    a_delta[5][1] = 5.85586261749;
    a_delta[5][2] = 12.3102693558;
    a_delta[5][3] = 21.393032074;
    a_delta[5][4] = 31.7420043945;
    a_delta[6][0] = 1.7427560091;
    a_delta[6][1] = 7.33418560028;
    a_delta[6][2] = 15.4683561325;
    a_delta[6][3] = 26.5165843964;
    a_delta[6][4] = 37.535987854;


    a_alpha[0][0] = 0.598146140575;
    a_alpha[0][1] = 0.578208088875;
    a_alpha[0][2] = 0.535880565643;
    a_alpha[0][3] = 0.552722096443;
    a_alpha[0][4] = 0.578202664852;
    a_alpha[1][0] = 1.35998690128;
    a_alpha[1][1] = 1.33517253399;
    a_alpha[1][2] = 1.28430593014;
    a_alpha[1][3] = 1.24044024944;
    a_alpha[1][4] = 1.20355796814;
    a_alpha[2][0] = 1.85280632973;
    a_alpha[2][1] = 1.8602206707;
    a_alpha[2][2] = 1.78636014462;
    a_alpha[2][3] = 1.7336935997;
    a_alpha[2][4] = 1.67925310135;
    a_alpha[3][0] = 2.18563699722;
    a_alpha[3][1] = 2.19838356972;
    a_alpha[3][2] = 2.13869380951;
    a_alpha[3][3] = 2.07461500168;
    a_alpha[3][4] = 1.96054816246;
    a_alpha[4][0] = 2.66806054115;
    a_alpha[4][1] = 2.46829032898;
    a_alpha[4][2] = 2.36234903336;
    a_alpha[4][3] = 2.26285719872;
    a_alpha[4][4] = 2.12015795708;
    a_alpha[5][0] = 2.88549137115;
    a_alpha[5][1] = 2.63377022743;
    a_alpha[5][2] = 2.46883964539;
    a_alpha[5][3] = 2.30901432037;
    a_alpha[5][4] = 2.18072772026;
    a_alpha[6][0] = 4.01096820831;
    a_alpha[6][1] = 3.22500920296;
    a_alpha[6][2] = 3.07853531837;
    a_alpha[6][3] = 2.73314976692;
    a_alpha[6][4] = 2.48712253571;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 2.42344376602e-05;
    a_r[0][1] = 1.1623291357e-05;
    a_r[0][2] = 7.65204003983e-06;
    a_r[0][3] = 5.83182281844e-06;
    a_r[0][4] = 4.40988515038e-06;
    a_r[1][0] = 1.79161961569e-05;
    a_r[1][1] = 9.16974022402e-06;
    a_r[1][2] = 5.99601571594e-06;
    a_r[1][3] = 4.60673072666e-06;
    a_r[1][4] = 3.48981484422e-06;
    a_r[2][0] = 1.22421797641e-05;
    a_r[2][1] = 6.29486339676e-06;
    a_r[2][2] = 4.21144795837e-06;
    a_r[2][3] = 3.10019527205e-06;
    a_r[2][4] = 2.44980742536e-06;
    a_r[3][0] = 8.66732079885e-06;
    a_r[3][1] = 4.37205198978e-06;
    a_r[3][2] = 2.93707080345e-06;
    a_r[3][3] = 2.25974667956e-06;
    a_r[3][4] = 1.80137590178e-06;
    a_r[4][0] = 5.94665698372e-06;
    a_r[4][1] = 2.9800298762e-06;
    a_r[4][2] = 2.08772257793e-06;
    a_r[4][3] = 1.58348893819e-06;
    a_r[4][4] = 1.34879110192e-06;
    a_r[5][0] = 4.0092877498e-06;
    a_r[5][1] = 2.09302174881e-06;
    a_r[5][2] = 1.51016160999e-06;
    a_r[5][3] = 1.20712650187e-06;
    a_r[5][4] = 9.9290059552e-07;
    a_r[6][0] = 2.70367399935e-06;
    a_r[6][1] = 1.51722326791e-06;
    a_r[6][2] = 1.14653687433e-06;
    a_r[6][3] = 9.66354377852e-07;
    a_r[6][4] = 8.75749890383e-07;


    a_width[0][0] = 26.5976848602;
    a_width[0][1] = 26.5976848602;
    a_width[0][2] = 26.7212715149;
    a_width[0][3] = 26.5424861908;
    a_width[0][4] = 26.4101047516;
    a_width[1][0] = 58.7029762268;
    a_width[1][1] = 58.7029762268;
    a_width[1][2] = 58.7458114624;
    a_width[1][3] = 58.1389541626;
    a_width[1][4] = 58.612991333;
    a_width[2][0] = 103.93271637;
    a_width[2][1] = 103.93271637;
    a_width[2][2] = 103.803878784;
    a_width[2][3] = 103.800476074;
    a_width[2][4] = 103.74822998;
    a_width[3][0] = 162.97102356;
    a_width[3][1] = 162.97102356;
    a_width[3][2] = 162.841094971;
    a_width[3][3] = 162.802963257;
    a_width[3][4] = 162.688156128;
    a_width[4][0] = 235.570587158;
    a_width[4][1] = 235.570587158;
    a_width[4][2] = 235.355682373;
    a_width[4][3] = 235.014694214;
    a_width[4][4] = 235.070922852;
    a_width[5][0] = 320.552429199;
    a_width[5][1] = 320.552429199;
    a_width[5][2] = 320.081390381;
    a_width[5][3] = 319.721862793;
    a_width[5][4] = 319.82611084;
    a_width[6][0] = 416.118835449;
    a_width[6][1] = 416.118835449;
    a_width[6][2] = 430.821716309;
    a_width[6][3] = 416.076812744;
    a_width[6][4] = 415.9637146;

  }
  else {

    a_delta[0][0] = 0.015081577003;
    a_delta[0][1] = 0.0576512292027;
    a_delta[0][2] = 0.131081297994;
    a_delta[0][3] = 0.243380963802;
    a_delta[0][4] = 0.39500105381;
    a_delta[1][0] = 0.0972077324986;
    a_delta[1][1] = 0.434791147709;
    a_delta[1][2] = 0.861800074577;
    a_delta[1][3] = 1.58613061905;
    a_delta[1][4] = 2.41115784645;
    a_delta[2][0] = 0.242783889174;
    a_delta[2][1] = 0.827510118484;
    a_delta[2][2] = 1.78883588314;
    a_delta[2][3] = 3.17804527283;
    a_delta[2][4] = 4.9168419838;
    a_delta[3][0] = 0.368439763784;
    a_delta[3][1] = 1.45783913136;
    a_delta[3][2] = 3.27468585968;
    a_delta[3][3] = 5.93577003479;
    a_delta[3][4] = 9.49385356903;
    a_delta[4][0] = 0.679009974003;
    a_delta[4][1] = 2.96808123589;
    a_delta[4][2] = 6.60964488983;
    a_delta[4][3] = 11.7470684052;
    a_delta[4][4] = 18.1083869934;
    a_delta[5][0] = 1.51261091232;
    a_delta[5][1] = 5.85586261749;
    a_delta[5][2] = 12.3102693558;
    a_delta[5][3] = 21.393032074;
    a_delta[5][4] = 31.7420043945;
    a_delta[6][0] = 1.7427560091;
    a_delta[6][1] = 7.33418560028;
    a_delta[6][2] = 15.4683561325;
    a_delta[6][3] = 26.5165843964;
    a_delta[6][4] = 37.535987854;


    a_alpha[0][0] = 0.693464934826;
    a_alpha[0][1] = 0.764990329742;
    a_alpha[0][2] = 0.725805997849;
    a_alpha[0][3] = 0.772053778172;
    a_alpha[0][4] = 0.747024357319;
    a_alpha[1][0] = 1.3640486002;
    a_alpha[1][1] = 1.3855907917;
    a_alpha[1][2] = 1.32408082485;
    a_alpha[1][3] = 1.28046572208;
    a_alpha[1][4] = 1.25549399853;
    a_alpha[2][0] = 1.75767993927;
    a_alpha[2][1] = 1.73263728619;
    a_alpha[2][2] = 1.71535384655;
    a_alpha[2][3] = 1.61282861233;
    a_alpha[2][4] = 1.5456840992;
    a_alpha[3][0] = 1.97163009644;
    a_alpha[3][1] = 1.95817959309;
    a_alpha[3][2] = 1.93145811558;
    a_alpha[3][3] = 1.84361171722;
    a_alpha[3][4] = 1.74427652359;
    a_alpha[4][0] = 2.38409137726;
    a_alpha[4][1] = 2.25097537041;
    a_alpha[4][2] = 2.11196112633;
    a_alpha[4][3] = 2.01256418228;
    a_alpha[4][4] = 1.83331000805;
    a_alpha[5][0] = 2.666233778;
    a_alpha[5][1] = 2.43228340149;
    a_alpha[5][2] = 2.22717642784;
    a_alpha[5][3] = 2.09035038948;
    a_alpha[5][4] = 1.91368079185;
    a_alpha[6][0] = 2.98137307167;
    a_alpha[6][1] = 2.68790650368;
    a_alpha[6][2] = 2.48271894455;
    a_alpha[6][3] = 2.61428976059;
    a_alpha[6][4] = 2.04934287071;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 2.2784628527e-05;
    a_r[0][1] = 1.20308704936e-05;
    a_r[0][2] = 8.08533604868e-06;
    a_r[0][3] = 5.96977088207e-06;
    a_r[0][4] = 4.75232536701e-06;
    a_r[1][0] = 1.79455546458e-05;
    a_r[1][1] = 9.10651215236e-06;
    a_r[1][2] = 5.98334918323e-06;
    a_r[1][3] = 4.73274121759e-06;
    a_r[1][4] = 3.70004795514e-06;
    a_r[2][0] = 1.22285728139e-05;
    a_r[2][1] = 6.1476503106e-06;
    a_r[2][2] = 4.39912128058e-06;
    a_r[2][3] = 3.30409216076e-06;
    a_r[2][4] = 2.65621315521e-06;
    a_r[3][0] = 8.09348966868e-06;
    a_r[3][1] = 4.27752138421e-06;
    a_r[3][2] = 2.96758730656e-06;
    a_r[3][3] = 2.29305896937e-06;
    a_r[3][4] = 1.8912643327e-06;
    a_r[4][0] = 5.57275279789e-06;
    a_r[4][1] = 3.00000010611e-06;
    a_r[4][2] = 2.09141035157e-06;
    a_r[4][3] = 1.65903634297e-06;
    a_r[4][4] = 1.38937730298e-06;
    a_r[5][0] = 3.85146495319e-06;
    a_r[5][1] = 1.99999999495e-06;
    a_r[5][2] = 1.47301966535e-06;
    a_r[5][3] = 1.27019291085e-06;
    a_r[5][4] = 1.09771178813e-06; 
    a_r[6][0] = 2.4772159577e-06;
    a_r[6][1] = 1.43584577472e-06;
    a_r[6][2] = 1.08408983124e-06;
    a_r[6][3] = 1.03147738173e-06;
    a_r[6][4] = 8.84041753579e-07;


    a_width[0][0] = 26.5976848602;
    a_width[0][1] = 26.5976848602;
    a_width[0][2] = 26.6827316284;
    a_width[0][3] = 26.6214389801;
    a_width[0][4] = 26.3873176575;
    a_width[1][0] = 58.7029762268;
    a_width[1][1] = 58.7029762268;
    a_width[1][2] = 58.8415145874;
    a_width[1][3] = 58.4030990601;
    a_width[1][4] = 58.6177330017;
    a_width[2][0] = 103.93271637;
    a_width[2][1] = 103.93271637;
    a_width[2][2] = 104.202682495;
    a_width[2][3] = 103.806724548;
    a_width[2][4] = 103.740852356;
    a_width[3][0] = 162.97102356;
    a_width[3][1] = 162.97102356;
    a_width[3][2] = 162.97102356;
    a_width[3][3] = 162.682937622;
    a_width[3][4] = 162.612564087;
    a_width[4][0] = 235.570587158;
    a_width[4][1] = 235.570587158;
    a_width[4][2] = 235.356552124;
    a_width[4][3] = 234.996078491;
    a_width[4][4] = 235.060302734;
    a_width[5][0] = 320.552429199;
    a_width[5][1] = 320.552429199;
    a_width[5][2] = 344.464599609;
    a_width[5][3] = 320.296142578;
    a_width[5][4] = 320.367614746;
    a_width[6][0] = 416.118835449;
    a_width[6][1] = 416.118835449;
    a_width[6][2] = 464.375823975;
    a_width[6][3] = 455.52545166;
    a_width[6][4] = 415.974060059;
    
  }

  
}



Double_t RooCPSHighMassGGH::interpolateMatrix(const Double_t matrix[][5] , const Double_t& x, const Double_t& y) const
{

  // Dummy copy of ROOT TH2::Interpolate needed because there the function is not constant

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring static const nMass, nCPrime" << std::endl; 
  static const Int_t nMass   = 7;
  static const Int_t nCPrime = 5;

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring arrays masses, cprimes" << std::endl; 
  Double_t masses[nMass]    = {400.,500.,600.,700.,800.,900.,1000.};
  Double_t cprimes[nCPrime] = {0.2,0.4,0.6,0.8,1.0};

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring stepX, stepY, lowX, lowY" << std::endl; 
  Double_t stepX = masses[1]-masses[0];
  Double_t stepY = cprimes[1]-cprimes[0];
  Double_t lowX = masses[0]-stepX/2;  //Double_t highX = masses[nMass-1]+stepX/2;
  Double_t lowY = cprimes[0]-stepY/2; //Double_t highY = cprimes[nMass-1]+stepY/2;
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: stepX=" << stepX << ", stepY=" << stepY << ", lowX=" << lowX << ", lowY=" << lowY << std::endl; 
  Int_t nX = nMass;
  Int_t nY = nCPrime;

  // Get the indexes to identify matrix element
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring i_x, i_y" << std::endl; 
  Int_t i_x = (Int_t)((x-lowX)/stepX);
  Int_t i_y = (Int_t)((y-lowY)/stepY);
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: i_x=" << i_x << ", i_y=" << i_y << std::endl; 

  // Check if we are in the domain for the interpolation
  if(i_x<0 || i_x>(nX-1) || i_y<0 || i_y>(nY-1)){
    std::cout << "[RooCPSHighMassGGH::interpolateMatrix] Cannot interpolate outside array domain, returning 0." << std::endl;
    return 0;
  }

  // Find out the quadrant in the bin
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring quadrant, dx, dy" << std::endl; 
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  dx = masses[i_x]+stepX/2-x;
  dy = cprimes[i_y]+stepY/2-y;
  if (dx<=stepX/2 && dy<=stepY/2)
    quadrant = 1; // upper right
  if (dx>stepX/2 && dy<=stepY/2)
    quadrant = 2; // upper left
  if (dx>stepX/2 && dy>stepY/2)
    quadrant = 3; // lower left
  if (dx<=stepX/2 && dy>stepY/2)
    quadrant = 4; // lower right
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: quadrant=" << quadrant << ", dx=" << dx << ", dy=" << dy << std::endl; 

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring x_1, y1, x2, y2" << std::endl; 
  switch(quadrant) {
  case 1:
    x1 = masses[i_x];
    y1 = cprimes[i_y];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 2:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    y1 = cprimes[i_y];
    x2 = masses[i_x];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 3:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    x2 = masses[i_x];
    y2 = cprimes[i_y];
    break;
  case 4:
    x1 = masses[i_x];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    y2 = cprimes[i_y];
    break;
  }
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: x1=" << x1 << ", y1=" << y1 << ", x2=" << x2 << ", y2=" << y2 << std::endl; 

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring i_x1, i_x2, i_y1, i_y2" << std::endl; 
  Int_t i_x1 = (Int_t)((x1-lowX)/stepX);
  if(i_x1<0) i_x1=0;
  Int_t i_x2 = (Int_t)((x2-lowX)/stepX);
  if(i_x2>nMass-1) i_x2=nMass-1;
  Int_t i_y1 = (Int_t)((y1-lowY)/stepY);
  if(i_y1<0) i_y1=0;
  Int_t i_y2 = (Int_t)((y2-lowY)/stepY);
  if(i_y2>nCPrime-1) i_y2=nCPrime-1;
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: i_x1=" << i_x1 << ", i_x2=" << i_x2 << ", i_y1=" << i_y1 << ", i_y2=" << i_y2 << std::endl; 

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring q11, q12, q21, q22" << std::endl; 
  Double_t q11 = matrix[i_x1][i_y1];
  Double_t q12 = matrix[i_x1][i_y2];
  Double_t q21 = matrix[i_x2][i_y1];
  Double_t q22 = matrix[i_x2][i_y2];
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: q11=" << q11 << ", q12=" << q12 << ", q21=" << q21 << ", q22=" << q22 << std::endl; 

  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: declaring d, f" << std::endl;   
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  // std::cout << "[RooCPSHighMassGGH::interpolateMatrix] DEBUG: d=" << d << ", f=" << f << std::endl;
  
  return f;
}


Double_t RooCPSHighMassGGH::getWidth(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_width,mass,cprime);}

Double_t RooCPSHighMassGGH::getDelta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_delta,mass,cprime);}

Double_t RooCPSHighMassGGH::getAlpha(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_alpha,mass,cprime);}

Double_t RooCPSHighMassGGH::getR(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_r,mass,cprime);}

Double_t RooCPSHighMassGGH::getBeta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_beta,mass,cprime);}




Double_t RooCPSHighMassGGH::Spline(Double_t xx) const{
   const int fNp = 20, fKstep = 0;
   const double fDelta = -1, fXmin = 200, fXmax = 2150;
   const double fX[20] = { 200, 250, 300, 350, 400,
                        450, 500, 550, 600, 650,
                        750, 850, 950, 1050, 1150,
                        1350, 1550, 1750, 1950, 2150 };
   const double fY[20] = { 0.139988, 0.17804, 0.352505, 0.801982, 1.35134,
                        1.75174, 1.84163, 1.86637, 1.78052, 1.7297,
                        1.53036, 1.31241, 1.11751, 0.946671, 0.812524,
                        0.640316, 0.490843, 0.396263, 0.315344, 0.210836 };
   const double fB[20] = { 0.000266891, 0.00147333, 0.00628998, 0.0108012, 0.0105333,
                        0.00444666, 0.000737389, -0.000594766, -0.00180599, -0.000721876,
                        -0.00269494, -0.00186614, -0.00192752, -0.00150674, -0.00115716,
                        -0.000743197, -0.000649875, -0.000369596, -0.000461366, -0.000583113 };
   const double fC[20] = { 8.42921e-06, 2.29701e-05, 6.55918e-05, 1.67598e-05, -2.50909e-05,
                        -7.98341e-05, 4.46636e-06, -3.60991e-05, 1.78298e-05, -7.22121e-06,
                        3.25258e-06, -4.32435e-07, 1.81479e-06, 1.62255e-06, 2.04841e-06,
                        -2.33522e-07, 1.01254e-06, 5.76421e-08, -3.11836e-07, -3.01879e-07 };
   const double fD[20] = { 0, 2.90818e-07, -2.00191e-08, -3.34947e-07, -4.53847e-07,
                        2.18457e-07, 1.71166e-07, -2.84202e-08, 2.09784e-07, -2.33159e-07,
                        8.20808e-08, -2.72765e-08, 1.26473e-08, -2.76285e-09, 7.40335e-10,
                        -1.80875e-09, 1.31533e-09, -1.99675e-09, 4.97867e-11, 0 };
   const double fE[20] = { 0, 2.90818e-09, -6.01655e-09, 2.86727e-09, -4.05627e-09,
                        1.07793e-08, -1.12522e-08, 9.25635e-09, -6.87431e-09, 2.44487e-09,
                        -8.68676e-10, 3.2189e-10, -1.22271e-10, 4.52199e-11, -2.77039e-11,
                        2.13312e-11, -1.3521e-11, 5.24081e-12, -1.24467e-13, 0 };
   const double fF[20] = { 1.16327e-11, -3.56989e-11, 3.55353e-11, -2.76942e-11, 5.93423e-11,
                        -8.81261e-11, 8.20343e-11, -6.45226e-11, 3.72767e-11, -6.6271e-12,
                        2.38113e-12, -8.8832e-13, 3.34981e-13, -1.45848e-13, 4.90352e-14,
                        -3.48522e-14, 1.87618e-14, -5.36528e-15, 1.24467e-16, 0 };
   int klow=0;
   if(xx<=fXmin) klow=0;
   else if(xx>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((xx-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=xx-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*(fD[klow]+dx*(fE[klow]+dx*fF[klow])))));
}




/*************RooBWHighMassGGH***********/
ClassImp(RooBWHighMassGGH) 

 RooBWHighMassGGH::RooBWHighMassGGH(){}

 RooBWHighMassGGH::RooBWHighMassGGH(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _KPrime,
			RooAbsReal& _BRnew,
			RooAbsReal& _IntStr,
			Bool_t is8TeV_b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   KPrime("KPrime","KPrime",this,_KPrime),
   BRnew("BRnew","BRnew",this,_BRnew),
   IntStr("IntStr","IntStr",this,_IntStr),
   is8TeV(is8TeV_b)
 {    
   initMatrices();
 } 


 RooBWHighMassGGH::RooBWHighMassGGH(const RooBWHighMassGGH& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   KPrime("KPrime",this,other.KPrime),
   BRnew("BRnew",this,other.BRnew),
   IntStr("IntStr",this,other.IntStr),
   is8TeV(other.is8TeV)
 { 
   
   for (Int_t i=0; i<7; ++i){
     for (Int_t j=0; j<5; ++j){
       a_width[i][j] = other.a_width[i][j];
       a_delta[i][j] = other.a_delta[i][j];
       a_alpha[i][j] = other.a_alpha[i][j];
       a_r[i][j] = other.a_r[i][j];
       a_beta[i][j] = other.a_beta[i][j];
     }
   }

 } 



 Double_t RooBWHighMassGGH::evaluate() const 
 {    

   Double_t effKPrime = KPrime/(1-BRnew);

   Double_t width = getWidth(mH,effKPrime);
   Double_t delta = getDelta(mH,effKPrime);
   Double_t alpha = getAlpha(mH,effKPrime);
   Double_t r = getR(mH,effKPrime);
   Double_t beta = getBeta(mH,effKPrime);

   // Cast the "x" proxy to double only once to avoid overhead and ambiguous
   // overloads the TComplex operators.
   const double xVal = x;

   Double_t bwHM = xVal / ( TMath::Power( TMath::Power(xVal,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(xVal,2)*TMath::Power(effKPrime*width,2) ); 
   Double_t signal = bwHM;
     
   Double_t k=0.25;
   Double_t mH_eff = mH*TMath::Sqrt(1-k*(effKPrime*width/mH)*(effKPrime*width/mH));
   TComplex MSquared(mH_eff*mH_eff,-effKPrime*width*mH_eff);
   TComplex M = MSquared.Sqrt(MSquared);
   TComplex Exp1 = MSquared.Exp((TComplex)alpha);
   TComplex Exp2 = MSquared.Exp(-(TComplex)alpha);
   TComplex Exp3 = MSquared.Exp(-x*(TComplex)alpha/M);

   double interference = -r*(1-TMath::Exp(-beta*(xVal-150.)/mH_eff))*((Exp1/(xVal-M)-Exp2/(xVal+M))*Exp3).Re();

   Double_t fValue = signal + IntStr*interference/( TMath::Sqrt(1-BRnew) );
   if (fValue > 0) return fValue;
   else return 0.;   
   
 } 




void RooBWHighMassGGH::initMatrices() {


  if (is8TeV){

    a_delta[0][0] = -0.0356573872268;
    a_delta[0][1] = 0.24223549664;
    a_delta[0][2] = 0.61275935173;
    a_delta[0][3] = 0.383159637451;
    a_delta[0][4] = 0.210959866643;
    a_delta[1][0] = -0.12527577579;
    a_delta[1][1] = 0.481088936329;
    a_delta[1][2] = 1.28957521915;
    a_delta[1][3] = 1.21549904346;
    a_delta[1][4] = 1.1599419117;
    a_delta[2][0] = -0.238580748439;
    a_delta[2][1] = 1.60643100739;
    a_delta[2][2] = 4.06644678116;
    a_delta[2][3] = 4.52600479126;
    a_delta[2][4] = 4.87067365646;
    a_delta[3][0] = -0.577017605305;
    a_delta[3][1] = 1.585511446;
    a_delta[3][2] = 4.4688835144;
    a_delta[3][3] = 4.78664445877;
    a_delta[3][4] = 5.02496528625;
    a_delta[4][0] = -1.91204202175;
    a_delta[4][1] = 2.39801764488;
    a_delta[4][2] = 8.14476394653;
    a_delta[4][3] = 7.69154310226;
    a_delta[4][4] = 7.35162782669;
    a_delta[5][0] = -2.2465929985;
    a_delta[5][1] = 3.01269602776;
    a_delta[5][2] = 10.0250816345;
    a_delta[5][3] = 13.1820316315;
    a_delta[5][4] = 15.5497436523;
    a_delta[6][0] = -7.72830677032;
    a_delta[6][1] = -9.27864837646;
    a_delta[6][2] = -11.3457698822;
    a_delta[6][3] = 0.881492912769;
    a_delta[6][4] = 10.0519399643;


    a_alpha[0][0] = 0.623292386532;
    a_alpha[0][1] = 0.640495181084;
    a_alpha[0][2] = 0.663432240486;
    a_alpha[0][3] = 0.694594204426;
    a_alpha[0][4] = 0.717965662479;
    a_alpha[1][0] = 0.606636106968;
    a_alpha[1][1] = 0.853205323219;
    a_alpha[1][2] = 1.18196427822;
    a_alpha[1][3] = 0.890077710152;
    a_alpha[1][4] = 0.6711627841;
    a_alpha[2][0] = 0.653616249561;
    a_alpha[2][1] = 0.913539767265;
    a_alpha[2][2] = 1.2601044178;
    a_alpha[2][3] = 1.01195764542;
    a_alpha[2][4] = 0.825847506523;
    a_alpha[3][0] = 0.792090892792;
    a_alpha[3][1] = 0.818761646748;
    a_alpha[3][2] = 0.85432267189;
    a_alpha[3][3] = 0.832669138908;
    a_alpha[3][4] = 0.816429018974;
    a_alpha[4][0] = 0.646219670773;
    a_alpha[4][1] = 0.859919428825;
    a_alpha[4][2] = 1.14485239983;
    a_alpha[4][3] = 0.955607414246;
    a_alpha[4][4] = 0.81367367506;
    a_alpha[5][0] = 0.77941018343;
    a_alpha[5][1] = 0.968531548977;
    a_alpha[5][2] = 1.22069334984;
    a_alpha[5][3] = 1.1067135334;
    a_alpha[5][4] = 1.02122867107;
    a_alpha[6][0] = 0.578494369984;
    a_alpha[6][1] = 0.760364711285;
    a_alpha[6][2] = 1.00285851955;
    a_alpha[6][3] = 0.962706923485;
    a_alpha[6][4] = 0.932593226433;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 1.48761073433e-05;
    a_r[0][1] = 1.13308524305e-05;
    a_r[0][2] = 6.60384648654e-06;
    a_r[0][3] = 3.95643291995e-06;
    a_r[0][4] = 1.97087274501e-06;
    a_r[1][0] = 8.98416146811e-06;
    a_r[1][1] = 6.72903252052e-06;
    a_r[1][2] = 3.7221943785e-06;
    a_r[1][3] = 2.58310797108e-06;
    a_r[1][4] = 1.72879299498e-06;
    a_r[2][0] = 6.09106336924e-06;
    a_r[2][1] = 4.64774257125e-06;
    a_r[2][2] = 2.72331431006e-06;
    a_r[2][3] = 2.00821682483e-06;
    a_r[2][4] = 1.4718938246e-06;
    a_r[3][0] = 4.51792948297e-06;
    a_r[3][1] = 3.430521474e-06;
    a_r[3][2] = 1.98064412871e-06;
    a_r[3][3] = 1.46860270434e-06;
    a_r[3][4] = 1.08457163606e-06;
    a_r[4][0] = 3.64733136848e-06;
    a_r[4][1] = 2.79037953987e-06;
    a_r[4][2] = 1.64777725331e-06;
    a_r[4][3] = 1.22761582588e-06;
    a_r[4][4] = 9.12494670047e-07;
    a_r[5][0] = 2.97835867968e-06;
    a_r[5][1] = 2.2536275992e-06;
    a_r[5][2] = 1.28731971927e-06;
    a_r[5][3] = 1.00490672139e-06;
    a_r[5][4] = 7.93097001406e-07;
    a_r[6][0] = 2.38152711063e-06;
    a_r[6][1] = 1.79617472895e-06;
    a_r[6][2] = 1.01570481092e-06;
    a_r[6][3] = 8.09386847322e-07;
    a_r[6][4] = 6.54648374621e-07;


    a_width[0][0] = 29.3698158264;
    a_width[0][1] = 29.3002681732;
    a_width[0][2] = 29.2075366974;
    a_width[0][3] = 29.2782478333;
    a_width[0][4] = 29.3312797546;
    a_width[1][0] = 68.235206604;
    a_width[1][1] = 67.8260879517;
    a_width[1][2] = 67.2806015015;
    a_width[1][3] = 67.083480835;
    a_width[1][4] = 66.9356384277;
    a_width[2][0] = 124.218963623;
    a_width[2][1] = 121.670410156;
    a_width[2][2] = 118.272354126;
    a_width[2][3] = 117.718635559;
    a_width[2][4] = 117.303344727;
    a_width[3][0] = 199.596191406;
    a_width[3][1] = 191.995681763;
    a_width[3][2] = 181.861663818;
    a_width[3][3] = 183.690856934;
    a_width[3][4] = 185.062744141;
    a_width[4][0] = 310.622528076;
    a_width[4][1] = 297.16607666;
    a_width[4][2] = 279.224121094;
    a_width[4][3] = 276.43447876;
    a_width[4][4] = 274.342224121;
    a_width[5][0] = 461.420013428;
    a_width[5][1] = 449.248168945;
    a_width[5][2] = 433.019042969;
    a_width[5][3] = 417.291137695;
    a_width[5][4] = 405.495178223;
    a_width[6][0] = 678.944702148;
    a_width[6][1] = 678.29888916;
    a_width[6][2] = 677.437805176;
    a_width[6][3] = 635.253601074;
    a_width[6][4] = 603.615478516;

  }
  else {

    a_delta[0][0] = -0.0387836955488;
    a_delta[0][1] = 0.206786364317;
    a_delta[0][2] = 0.534213125706;
    a_delta[0][3] = 0.356079429388;
    a_delta[0][4] = 0.2224791646;
    a_delta[1][0] = -0.142410874367;
    a_delta[1][1] = 0.410860836506;
    a_delta[1][2] = 1.14855647087;
    a_delta[1][3] = 0.902282416821;
    a_delta[1][4] = 0.717576861382;
    a_delta[2][0] = -0.349563479424;
    a_delta[2][1] = 1.12769818306;
    a_delta[2][2] = 3.0973803997;
    a_delta[2][3] = 3.09016609192;
    a_delta[2][4] = 3.08475542068;
    a_delta[3][0] = -0.729484498501;
    a_delta[3][1] = 1.12449371815;
    a_delta[3][2] = 3.59646463394;
    a_delta[3][3] = 3.09707832336;
    a_delta[3][4] = 2.72253847122;
    a_delta[4][0] = -2.05952262878;
    a_delta[4][1] = 0.973166584969;
    a_delta[4][2] = 5.01675224304;
    a_delta[4][3] = 4.60968160629;
    a_delta[4][4] = 4.30437898636;
    a_delta[5][0] = -1.13842725754;
    a_delta[5][1] = 2.01005077362;
    a_delta[5][2] = 6.20802164078;
    a_delta[5][3] = 8.58551692963;
    a_delta[5][4] = 10.3686380386;
    a_delta[6][0] = -16.7097110748;
    a_delta[6][1] = -17.208190918;
    a_delta[6][2] = -17.8728294373;
    a_delta[6][3] = -11.6410417557;
    a_delta[6][4] = -6.96720075607;


    a_alpha[0][0] = 0.657900691032;
    a_alpha[0][1] = 1.35420918465;
    a_alpha[0][2] = 2.28262042999;
    a_alpha[0][3] = 2.36626386642;
    a_alpha[0][4] = 2.42899632454;
    a_alpha[1][0] = 0.66292309761;
    a_alpha[1][1] = 0.851784825325;
    a_alpha[1][2] = 1.10360050201;
    a_alpha[1][3] = 0.840485334396;
    a_alpha[1][4] = 0.643148958683;
    a_alpha[2][0] = 0.710033237934;
    a_alpha[2][1] = 0.847059130669;
    a_alpha[2][2] = 1.02976036072;
    a_alpha[2][3] = 0.860089302063;
    a_alpha[2][4] = 0.732836008072;
    a_alpha[3][0] = 1.42068183422;
    a_alpha[3][1] = 1.22611534595;
    a_alpha[3][2] = 0.966693341732;
    a_alpha[3][3] = 0.86537194252;
    a_alpha[3][4] = 0.789380848408;
    a_alpha[4][0] = 1.44478917122;
    a_alpha[4][1] = 1.32054758072;
    a_alpha[4][2] = 1.15489220619;
    a_alpha[4][3] = 0.964189589024;
    a_alpha[4][4] = 0.821162641048;
    a_alpha[5][0] = 0.878096818924;
    a_alpha[5][1] = 1.01728570461;
    a_alpha[5][2] = 1.202870965;
    a_alpha[5][3] = 1.09377551079;
    a_alpha[5][4] = 1.01195383072;
    a_alpha[6][0] = 2.99999547005;
    a_alpha[6][1] = 2.05522418022;
    a_alpha[6][2] = 0.795529186726;
    a_alpha[6][3] = 0.774422287941;
    a_alpha[6][4] = 0.758592069149;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 1.49143634189e-05;
    a_r[0][1] = 1.14036974992e-05;
    a_r[0][2] = 6.72280884828e-06;
    a_r[0][3] = 4.13136740463e-06;
    a_r[0][4] = 2.18778609451e-06;
    a_r[1][0] = 8.97418158274e-06;
    a_r[1][1] = 6.73627982906e-06;
    a_r[1][2] = 3.75241097572e-06;
    a_r[1][3] = 2.63809306489e-06;
    a_r[1][4] = 1.80235451808e-06;
    a_r[2][0] = 6.20700529907e-06;
    a_r[2][1] = 4.73968111692e-06;
    a_r[2][2] = 2.78324932879e-06;
    a_r[2][3] = 2.06628055821e-06;
    a_r[2][4] = 1.52855409397e-06;
    a_r[3][0] = 4.77096682516e-06;
    a_r[3][1] = 3.59284967999e-06;
    a_r[3][2] = 2.02202704713e-06;
    a_r[3][3] = 1.51344966071e-06;
    a_r[3][4] = 1.13201656404e-06;
    a_r[4][0] = 3.66865469914e-06;
    a_r[4][1] = 2.83416306956e-06;
    a_r[4][2] = 1.72150771505e-06;
    a_r[4][3] = 1.28625231355e-06;
    a_r[4][4] = 9.59810677159e-07;
    a_r[5][0] = 2.97442375086e-06;
    a_r[5][1] = 2.27610507864e-06;
    a_r[5][2] = 1.34501340199e-06;
    a_r[5][3] = 1.05685342078e-06;
    a_r[5][4] = 8.40733491714e-07;
    a_r[6][0] = 3.23288031723e-06;
    a_r[6][1] = 2.31977242038e-06;
    a_r[6][2] = 1.10229530037e-06;
    a_r[6][3] = 8.74006730101e-07;
    a_r[6][4] = 7.02790316609e-07;


    a_width[0][0] = 29.3776245117;
    a_width[0][1] = 29.2529735565;
    a_width[0][2] = 29.0867729187;
    a_width[0][3] = 29.2772274017;
    a_width[0][4] = 29.4200706482;
    a_width[1][0] = 68.4645004272;
    a_width[1][1] = 67.8752212524;
    a_width[1][2] = 67.089515686;
    a_width[1][3] = 66.9777069092;
    a_width[1][4] = 66.8938446045;
    a_width[2][0] = 123.957496643;
    a_width[2][1] = 121.468582153;
    a_width[2][2] = 118.150032043;
    a_width[2][3] = 117.622924805;
    a_width[2][4] = 117.227584839;
    a_width[3][0] = 202.743499756;
    a_width[3][1] = 196.633117676;
    a_width[3][2] = 188.485946655;
    a_width[3][3] = 189.171157837;
    a_width[3][4] = 189.685073853;
    a_width[4][0] = 319.84085083;
    a_width[4][1] = 312.385772705;
    a_width[4][2] = 302.445678711;
    a_width[4][3] = 295.498382568;
    a_width[4][4] = 290.287902832;
    a_width[5][0] = 453.523803711;
    a_width[5][1] = 451.083343506;
    a_width[5][2] = 447.829406738;
    a_width[5][3] = 437.375549316;
    a_width[5][4] = 429.53515625;
    a_width[6][0] = 692.998046875;
    a_width[6][1] = 689.624084473;
    a_width[6][2] = 685.125427246;
    a_width[6][3] = 640.699401855;
    a_width[6][4] = 607.379943848;
    
  }

  
}



Double_t RooBWHighMassGGH::interpolateMatrix(const Double_t matrix[][5] , const Double_t& x, const Double_t& y) const
{

  // Dummy copy of ROOT TH2::Interpolate needed because there the function is not constant

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring static const nMass, nCPrime" << std::endl; 
  static const Int_t nMass   = 7;
  static const Int_t nCPrime = 5;

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring arrays masses, cprimes" << std::endl; 
  Double_t masses[nMass]    = {400.,500.,600.,700.,800.,900.,1000.};
  Double_t cprimes[nCPrime] = {0.2,0.4,0.6,0.8,1.0};

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring stepX, stepY, lowX, lowY" << std::endl; 
  Double_t stepX = masses[1]-masses[0];
  Double_t stepY = cprimes[1]-cprimes[0];
  Double_t lowX = masses[0]-stepX/2;  //Double_t highX = masses[nMass-1]+stepX/2;
  Double_t lowY = cprimes[0]-stepY/2; //Double_t highY = cprimes[nMass-1]+stepY/2;
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: stepX=" << stepX << ", stepY=" << stepY << ", lowX=" << lowX << ", lowY=" << lowY << std::endl; 
  Int_t nX = nMass;
  Int_t nY = nCPrime;

  // Get the indexes to identify matrix element
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring i_x, i_y" << std::endl; 
  Int_t i_x = (Int_t)((x-lowX)/stepX);
  Int_t i_y = (Int_t)((y-lowY)/stepY);
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: i_x=" << i_x << ", i_y=" << i_y << std::endl; 

  // Check if we are in the domain for the interpolation
  if(i_x<0 || i_x>(nX-1) || i_y<0 || i_y>(nY-1)){
    std::cout << "[RooBWHighMassGGH::interpolateMatrix] Cannot interpolate outside array domain, returning 0." << std::endl;
    return 0;
  }

  // Find out the quadrant in the bin
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring quadrant, dx, dy" << std::endl; 
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  dx = masses[i_x]+stepX/2-x;
  dy = cprimes[i_y]+stepY/2-y;
  if (dx<=stepX/2 && dy<=stepY/2)
    quadrant = 1; // upper right
  if (dx>stepX/2 && dy<=stepY/2)
    quadrant = 2; // upper left
  if (dx>stepX/2 && dy>stepY/2)
    quadrant = 3; // lower left
  if (dx<=stepX/2 && dy>stepY/2)
    quadrant = 4; // lower right
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: quadrant=" << quadrant << ", dx=" << dx << ", dy=" << dy << std::endl; 

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring x_1, y1, x2, y2" << std::endl; 
  switch(quadrant) {
  case 1:
    x1 = masses[i_x];
    y1 = cprimes[i_y];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 2:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    y1 = cprimes[i_y];
    x2 = masses[i_x];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 3:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    x2 = masses[i_x];
    y2 = cprimes[i_y];
    break;
  case 4:
    x1 = masses[i_x];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    y2 = cprimes[i_y];
    break;
  }
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: x1=" << x1 << ", y1=" << y1 << ", x2=" << x2 << ", y2=" << y2 << std::endl; 

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring i_x1, i_x2, i_y1, i_y2" << std::endl; 
  Int_t i_x1 = (Int_t)((x1-lowX)/stepX);
  if(i_x1<0) i_x1=0;
  Int_t i_x2 = (Int_t)((x2-lowX)/stepX);
  if(i_x2>nMass-1) i_x2=nMass-1;
  Int_t i_y1 = (Int_t)((y1-lowY)/stepY);
  if(i_y1<0) i_y1=0;
  Int_t i_y2 = (Int_t)((y2-lowY)/stepY);
  if(i_y2>nCPrime-1) i_y2=nCPrime-1;
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: i_x1=" << i_x1 << ", i_x2=" << i_x2 << ", i_y1=" << i_y1 << ", i_y2=" << i_y2 << std::endl; 

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring q11, q12, q21, q22" << std::endl; 
  Double_t q11 = matrix[i_x1][i_y1];
  Double_t q12 = matrix[i_x1][i_y2];
  Double_t q21 = matrix[i_x2][i_y1];
  Double_t q22 = matrix[i_x2][i_y2];
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: q11=" << q11 << ", q12=" << q12 << ", q21=" << q21 << ", q22=" << q22 << std::endl; 

  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: declaring d, f" << std::endl;   
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  // std::cout << "[RooBWHighMassGGH::interpolateMatrix] DEBUG: d=" << d << ", f=" << f << std::endl;
  
  return f;
}


Double_t RooBWHighMassGGH::getWidth(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_width,mass,cprime);}

Double_t RooBWHighMassGGH::getDelta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_delta,mass,cprime);}

Double_t RooBWHighMassGGH::getAlpha(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_alpha,mass,cprime);}

Double_t RooBWHighMassGGH::getR(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_r,mass,cprime);}

Double_t RooBWHighMassGGH::getBeta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_beta,mass,cprime);}


Double_t RooBWHighMassGGH::H_width(Double_t mass) const
{
  // SM Higgs widths calculated with HDECAY [hep-ph/9704448]
  if      ((int)mass == 100) return  0.2788E-02;
  else if ((int)mass == 101) return  0.2820E-02;
  else if ((int)mass == 102) return  0.2853E-02;
  else if ((int)mass == 103) return  0.2887E-02;
  else if ((int)mass == 104) return  0.2923E-02;
  else if ((int)mass == 105) return  0.2960E-02;
  else if ((int)mass == 106) return  0.2999E-02;
  else if ((int)mass == 107) return  0.3039E-02;
  else if ((int)mass == 108) return  0.3082E-02;
  else if ((int)mass == 109) return  0.3127E-02;
  else if ((int)mass == 110) return  0.3175E-02;
  else if ((int)mass == 111) return  0.3225E-02;
  else if ((int)mass == 112) return  0.3278E-02;
  else if ((int)mass == 113) return  0.3335E-02;
  else if ((int)mass == 114) return  0.3396E-02;
  else if ((int)mass == 115) return  0.3460E-02;
  else if ((int)mass == 116) return  0.3529E-02;
  else if ((int)mass == 117) return  0.3603E-02;
  else if ((int)mass == 118) return  0.3682E-02;
  else if ((int)mass == 119) return  0.3767E-02;
  else if ((int)mass == 120) return  0.3859E-02;
  else if ((int)mass == 121) return  0.3957E-02;
  else if ((int)mass == 122) return  0.4063E-02;
  else if ((int)mass == 123) return  0.4177E-02;
  else if ((int)mass == 124) return  0.4301E-02;
  else if ((int)mass == 125) return  0.4434E-02;
  else if ((int)mass == 126) return  0.4578E-02;
  else if ((int)mass == 127) return  0.4734E-02;
  else if ((int)mass == 128) return  0.4903E-02;
  else if ((int)mass == 129) return  0.5086E-02;
  else if ((int)mass == 130) return  0.5284E-02;
  else if ((int)mass == 131) return  0.5500E-02;
  else if ((int)mass == 132) return  0.5734E-02;
  else if ((int)mass == 133) return  0.5989E-02;
  else if ((int)mass == 134) return  0.6266E-02;
  else if ((int)mass == 135) return  0.6569E-02;
  else if ((int)mass == 136) return  0.6900E-02;
  else if ((int)mass == 137) return  0.7261E-02;
  else if ((int)mass == 138) return  0.7657E-02;
  else if ((int)mass == 139) return  0.8091E-02;
  else if ((int)mass == 140) return  0.8569E-02;
  else if ((int)mass == 141) return  0.9096E-02;
  else if ((int)mass == 142) return  0.9677E-02;
  else if ((int)mass == 143) return  0.1032E-01;
  else if ((int)mass == 144) return  0.1104E-01;
  else if ((int)mass == 145) return  0.1184E-01;
  else if ((int)mass == 146) return  0.1274E-01;
  else if ((int)mass == 147) return  0.1375E-01;
  else if ((int)mass == 148) return  0.1489E-01;
  else if ((int)mass == 149) return  0.1620E-01;
  else if ((int)mass == 150) return  0.1770E-01;
  else if ((int)mass == 151) return  0.1945E-01;
  else if ((int)mass == 152) return  0.2150E-01;
  else if ((int)mass == 153) return  0.2396E-01;
  else if ((int)mass == 154) return  0.2693E-01;
  else if ((int)mass == 155) return  0.3063E-01;
  else if ((int)mass == 156) return  0.3535E-01;
  else if ((int)mass == 157) return  0.4162E-01;
  else if ((int)mass == 158) return  0.5035E-01;
  else if ((int)mass == 159) return  0.6316E-01;
  else if ((int)mass == 160) return  0.8272E-01;
  else if ((int)mass == 161) return  0.1111    ;
  else if ((int)mass == 162) return  0.1455    ;
  else if ((int)mass == 163) return  0.1804    ;
  else if ((int)mass == 164) return  0.2134    ;
  else if ((int)mass == 165) return  0.2442    ;
  else if ((int)mass == 166) return  0.2730    ;
  else if ((int)mass == 167) return  0.3004    ;
  else if ((int)mass == 168) return  0.3267    ;
  else if ((int)mass == 169) return  0.3521    ;
  else if ((int)mass == 170) return  0.3769    ;
  else if ((int)mass == 171) return  0.4013    ;
  else if ((int)mass == 172) return  0.4254    ;
  else if ((int)mass == 173) return  0.4493    ;
  else if ((int)mass == 174) return  0.4733    ;
  else if ((int)mass == 175) return  0.4974    ;
  else if ((int)mass == 176) return  0.5220    ;
  else if ((int)mass == 177) return  0.5471    ;
  else if ((int)mass == 178) return  0.5732    ;
  else if ((int)mass == 179) return  0.6006    ;
  else if ((int)mass == 180) return  0.6302    ;
  else if ((int)mass == 181) return  0.6630    ;
  else if ((int)mass == 182) return  0.7001    ;
  else if ((int)mass == 183) return  0.7417    ;
  else if ((int)mass == 184) return  0.7858    ;
  else if ((int)mass == 185) return  0.8300    ;
  else if ((int)mass == 186) return  0.8732    ;
  else if ((int)mass == 187) return  0.9154    ;
  else if ((int)mass == 188) return  0.9566    ;
  else if ((int)mass == 189) return  0.9971    ;
  else if ((int)mass == 190) return   1.037    ;
  else if ((int)mass == 191) return   1.076    ;
  else if ((int)mass == 192) return   1.116    ;
  else if ((int)mass == 193) return   1.154    ;
  else if ((int)mass == 194) return   1.193    ;
  else if ((int)mass == 195) return   1.232    ;
  else if ((int)mass == 196) return   1.271    ;
  else if ((int)mass == 197) return   1.310    ;
  else if ((int)mass == 198) return   1.349    ;
  else if ((int)mass == 199) return   1.388    ;
  else if ((int)mass == 200) return   1.428    ;
  else if ((int)mass == 201) return   1.468    ;
  else if ((int)mass == 202) return   1.508    ;
  else if ((int)mass == 203) return   1.548    ;
  else if ((int)mass == 204) return   1.589    ;
  else if ((int)mass == 205) return   1.630    ;
  else if ((int)mass == 206) return   1.672    ;
  else if ((int)mass == 207) return   1.714    ;
  else if ((int)mass == 208) return   1.756    ;
  else if ((int)mass == 209) return   1.799    ;
  else if ((int)mass == 210) return   1.842    ;
  else if ((int)mass == 211) return   1.886    ;
  else if ((int)mass == 212) return   1.930    ;
  else if ((int)mass == 213) return   1.975    ;
  else if ((int)mass == 214) return   2.020    ;
  else if ((int)mass == 215) return   2.066    ;
  else if ((int)mass == 216) return   2.112    ;
  else if ((int)mass == 217) return   2.159    ;
  else if ((int)mass == 218) return   2.206    ;
  else if ((int)mass == 219) return   2.254    ;
  else if ((int)mass == 220) return   2.303    ;
  else if ((int)mass == 221) return   2.352    ;
  else if ((int)mass == 222) return   2.401    ;
  else if ((int)mass == 223) return   2.451    ;
  else if ((int)mass == 224) return   2.502    ;
  else if ((int)mass == 225) return   2.553    ;
  else if ((int)mass == 226) return   2.605    ;
  else if ((int)mass == 227) return   2.658    ;
  else if ((int)mass == 228) return   2.711    ;
  else if ((int)mass == 229) return   2.765    ;
  else if ((int)mass == 230) return   2.819    ;
  else if ((int)mass == 231) return   2.874    ;
  else if ((int)mass == 232) return   2.930    ;
  else if ((int)mass == 233) return   2.986    ;
  else if ((int)mass == 234) return   3.043    ;
  else if ((int)mass == 235) return   3.101    ;
  else if ((int)mass == 236) return   3.159    ;
  else if ((int)mass == 237) return   3.218    ;
  else if ((int)mass == 238) return   3.278    ;
  else if ((int)mass == 239) return   3.338    ;
  else if ((int)mass == 240) return   3.399    ;
  else if ((int)mass == 241) return   3.461    ;
  else if ((int)mass == 242) return   3.523    ;
  else if ((int)mass == 243) return   3.586    ;
  else if ((int)mass == 244) return   3.650    ;
  else if ((int)mass == 245) return   3.714    ;
  else if ((int)mass == 246) return   3.780    ;
  else if ((int)mass == 247) return   3.845    ;
  else if ((int)mass == 248) return   3.912    ;
  else if ((int)mass == 249) return   3.979    ;
  else if ((int)mass == 250) return   4.048    ;
  else if ((int)mass == 251) return   4.116    ;
  else if ((int)mass == 252) return   4.186    ;
  else if ((int)mass == 253) return   4.256    ;
  else if ((int)mass == 254) return   4.327    ;
  else if ((int)mass == 255) return   4.399    ;
  else if ((int)mass == 256) return   4.472    ;
  else if ((int)mass == 257) return   4.545    ;
  else if ((int)mass == 258) return   4.619    ;
  else if ((int)mass == 259) return   4.694    ;
  else if ((int)mass == 260) return   4.770    ;
  else if ((int)mass == 261) return   4.846    ;
  else if ((int)mass == 262) return   4.923    ;
  else if ((int)mass == 263) return   5.001    ;
  else if ((int)mass == 264) return   5.080    ;
  else if ((int)mass == 265) return   5.159    ;
  else if ((int)mass == 266) return   5.239    ;
  else if ((int)mass == 267) return   5.320    ;
  else if ((int)mass == 268) return   5.402    ;
  else if ((int)mass == 269) return   5.485    ;
  else if ((int)mass == 270) return   5.568    ;
  else if ((int)mass == 271) return   5.652    ;
  else if ((int)mass == 272) return   5.737    ;
  else if ((int)mass == 273) return   5.823    ;
  else if ((int)mass == 274) return   5.910    ;
  else if ((int)mass == 275) return   5.997    ;
  else if ((int)mass == 276) return   6.085    ;
  else if ((int)mass == 277) return   6.175    ;
  else if ((int)mass == 278) return   6.264    ;
  else if ((int)mass == 279) return   6.355    ;
  else if ((int)mass == 280) return   6.447    ;
  else if ((int)mass == 281) return   6.539    ;
  else if ((int)mass == 282) return   6.632    ;
  else if ((int)mass == 283) return   6.726    ;
  else if ((int)mass == 284) return   6.821    ;
  else if ((int)mass == 285) return   6.917    ;
  else if ((int)mass == 286) return   7.013    ;
  else if ((int)mass == 287) return   7.111    ;
  else if ((int)mass == 288) return   7.209    ;
  else if ((int)mass == 289) return   7.308    ;
  else if ((int)mass == 290) return   7.408    ;
  else if ((int)mass == 291) return   7.509    ;
  else if ((int)mass == 292) return   7.611    ;
  else if ((int)mass == 293) return   7.714    ;
  else if ((int)mass == 294) return   7.817    ;
  else if ((int)mass == 295) return   7.921    ;
  else if ((int)mass == 296) return   8.027    ;
  else if ((int)mass == 297) return   8.133    ;
  else if ((int)mass == 298) return   8.240    ;
  else if ((int)mass == 299) return   8.348    ;
  else if ((int)mass == 300) return   8.456    ;
  else if ((int)mass == 301) return   8.566    ;
  else if ((int)mass == 302) return   8.676    ;
  else if ((int)mass == 303) return   8.788    ;
  else if ((int)mass == 304) return   8.900    ;
  else if ((int)mass == 305) return   9.013    ;
  else if ((int)mass == 306) return   9.127    ;
  else if ((int)mass == 307) return   9.242    ;
  else if ((int)mass == 308) return   9.358    ;
  else if ((int)mass == 309) return   9.474    ;
  else if ((int)mass == 310) return   9.592    ;
  else if ((int)mass == 311) return   9.710    ;
  else if ((int)mass == 312) return   9.829    ;
  else if ((int)mass == 313) return   9.950    ;
  else if ((int)mass == 314) return   10.07    ;
  else if ((int)mass == 315) return   10.19    ;
  else if ((int)mass == 316) return   10.32    ;
  else if ((int)mass == 317) return   10.44    ;
  else if ((int)mass == 318) return   10.56    ;
  else if ((int)mass == 319) return   10.69    ;
  else if ((int)mass == 320) return   10.82    ;
  else if ((int)mass == 321) return   10.94    ;
  else if ((int)mass == 322) return   11.07    ;
  else if ((int)mass == 323) return   11.20    ;
  else if ((int)mass == 324) return   11.33    ;
  else if ((int)mass == 325) return   11.46    ;
  else if ((int)mass == 326) return   11.59    ;
  else if ((int)mass == 327) return   11.72    ;
  else if ((int)mass == 328) return   11.86    ;
  else if ((int)mass == 329) return   11.99    ;
  else if ((int)mass == 330) return   12.13    ;
  else if ((int)mass == 331) return   12.26    ;
  else if ((int)mass == 332) return   12.40    ;
  else if ((int)mass == 333) return   12.53    ;
  else if ((int)mass == 334) return   12.67    ;
  else if ((int)mass == 335) return   12.81    ;
  else if ((int)mass == 336) return   12.95    ;
  else if ((int)mass == 337) return   13.09    ;
  else if ((int)mass == 338) return   13.22    ;
  else if ((int)mass == 339) return   13.36    ;
  else if ((int)mass == 340) return   13.50    ;
  else if ((int)mass == 341) return   13.64    ;
  else if ((int)mass == 342) return   13.78    ;
  else if ((int)mass == 343) return   13.92    ;
  else if ((int)mass == 344) return   14.06    ;
  else if ((int)mass == 345) return   14.20    ;
  else if ((int)mass == 346) return   14.37    ;
  else if ((int)mass == 347) return   14.54    ;
  else if ((int)mass == 348) return   14.74    ;
  else if ((int)mass == 349) return   14.95    ;
  else if ((int)mass == 350) return   15.17    ;
  else if ((int)mass == 351) return   15.40    ;
  else if ((int)mass == 352) return   15.63    ;
  else if ((int)mass == 353) return   15.86    ;
  else if ((int)mass == 354) return   16.09    ;
  else if ((int)mass == 355) return   16.33    ;
  else if ((int)mass == 356) return   16.58    ;
  else if ((int)mass == 357) return   16.82    ;
  else if ((int)mass == 358) return   17.07    ;
  else if ((int)mass == 359) return   17.32    ;
  else if ((int)mass == 360) return   17.57    ;
  else if ((int)mass == 361) return   17.83    ;
  else if ((int)mass == 362) return   18.08    ;
  else if ((int)mass == 363) return   18.34    ;
  else if ((int)mass == 364) return   18.60    ;
  else if ((int)mass == 365) return   18.87    ;
  else if ((int)mass == 366) return   19.13    ;
  else if ((int)mass == 367) return   19.40    ;
  else if ((int)mass == 368) return   19.67    ;
  else if ((int)mass == 369) return   19.94    ;
  else if ((int)mass == 370) return   20.21    ;
  else if ((int)mass == 371) return   20.49    ;
  else if ((int)mass == 372) return   20.77    ;
  else if ((int)mass == 373) return   21.04    ;
  else if ((int)mass == 374) return   21.33    ;
  else if ((int)mass == 375) return   21.61    ;
  else if ((int)mass == 376) return   21.89    ;
  else if ((int)mass == 377) return   22.18    ;
  else if ((int)mass == 378) return   22.46    ;
  else if ((int)mass == 379) return   22.75    ;
  else if ((int)mass == 380) return   23.04    ;
  else if ((int)mass == 381) return   23.33    ;
  else if ((int)mass == 382) return   23.63    ;
  else if ((int)mass == 383) return   23.92    ;
  else if ((int)mass == 384) return   24.22    ;
  else if ((int)mass == 385) return   24.52    ;
  else if ((int)mass == 386) return   24.82    ;
  else if ((int)mass == 387) return   25.12    ;
  else if ((int)mass == 388) return   25.42    ;
  else if ((int)mass == 389) return   25.72    ;
  else if ((int)mass == 390) return   26.03    ;
  else if ((int)mass == 391) return   26.34    ;
  else if ((int)mass == 392) return   26.64    ;
  else if ((int)mass == 393) return   26.95    ;
  else if ((int)mass == 394) return   27.26    ;
  else if ((int)mass == 395) return   27.58    ;
  else if ((int)mass == 396) return   27.89    ;
  else if ((int)mass == 397) return   28.21    ;
  else if ((int)mass == 398) return   28.52    ;
  else if ((int)mass == 399) return   28.84    ;
  else if ((int)mass == 400) return   29.16    ;
  else if ((int)mass == 401) return   29.48    ;
  else if ((int)mass == 402) return   29.80    ;
  else if ((int)mass == 403) return   30.13    ;
  else if ((int)mass == 404) return   30.45    ;
  else if ((int)mass == 405) return   30.78    ;
  else if ((int)mass == 406) return   31.10    ;
  else if ((int)mass == 407) return   31.43    ;
  else if ((int)mass == 408) return   31.76    ;
  else if ((int)mass == 409) return   32.09    ;
  else if ((int)mass == 410) return   32.43    ;
  else if ((int)mass == 411) return   32.76    ;
  else if ((int)mass == 412) return   33.10    ;
  else if ((int)mass == 413) return   33.43    ;
  else if ((int)mass == 414) return   33.77    ;
  else if ((int)mass == 415) return   34.11    ;
  else if ((int)mass == 416) return   34.45    ;
  else if ((int)mass == 417) return   34.79    ;
  else if ((int)mass == 418) return   35.14    ;
  else if ((int)mass == 419) return   35.48    ;
  else if ((int)mass == 420) return   35.83    ;
  else if ((int)mass == 421) return   36.18    ;
  else if ((int)mass == 422) return   36.52    ;
  else if ((int)mass == 423) return   36.87    ;
  else if ((int)mass == 424) return   37.23    ;
  else if ((int)mass == 425) return   37.58    ;
  else if ((int)mass == 426) return   37.93    ;
  else if ((int)mass == 427) return   38.29    ;
  else if ((int)mass == 428) return   38.64    ;
  else if ((int)mass == 429) return   39.00    ;
  else if ((int)mass == 430) return   39.36    ;
  else if ((int)mass == 431) return   39.72    ;
  else if ((int)mass == 432) return   40.08    ;
  else if ((int)mass == 433) return   40.45    ;
  else if ((int)mass == 434) return   40.81    ;
  else if ((int)mass == 435) return   41.18    ;
  else if ((int)mass == 436) return   41.55    ;
  else if ((int)mass == 437) return   41.91    ;
  else if ((int)mass == 438) return   42.28    ;
  else if ((int)mass == 439) return   42.66    ;
  else if ((int)mass == 440) return   43.03    ;
  else if ((int)mass == 441) return   43.40    ;
  else if ((int)mass == 442) return   43.78    ;
  else if ((int)mass == 443) return   44.15    ;
  else if ((int)mass == 444) return   44.53    ;
  else if ((int)mass == 445) return   44.91    ;
  else if ((int)mass == 446) return   45.29    ;
  else if ((int)mass == 447) return   45.67    ;
  else if ((int)mass == 448) return   46.06    ;
  else if ((int)mass == 449) return   46.44    ;
  else if ((int)mass == 450) return   46.83    ;
  else if ((int)mass == 451) return   47.22    ;
  else if ((int)mass == 452) return   47.61    ;
  else if ((int)mass == 453) return   48.00    ;
  else if ((int)mass == 454) return   48.39    ;
  else if ((int)mass == 455) return   48.78    ;
  else if ((int)mass == 456) return   49.17    ;
  else if ((int)mass == 457) return   49.57    ;
  else if ((int)mass == 458) return   49.97    ;
  else if ((int)mass == 459) return   50.37    ;
  else if ((int)mass == 460) return   50.77    ;
  else if ((int)mass == 461) return   51.17    ;
  else if ((int)mass == 462) return   51.57    ;
  else if ((int)mass == 463) return   51.97    ;
  else if ((int)mass == 464) return   52.38    ;
  else if ((int)mass == 465) return   52.79    ;
  else if ((int)mass == 466) return   53.20    ;
  else if ((int)mass == 467) return   53.61    ;
  else if ((int)mass == 468) return   54.02    ;
  else if ((int)mass == 469) return   54.43    ;
  else if ((int)mass == 470) return   54.84    ;
  else if ((int)mass == 471) return   55.26    ;
  else if ((int)mass == 472) return   55.68    ;
  else if ((int)mass == 473) return   56.09    ;
  else if ((int)mass == 474) return   56.51    ;
  else if ((int)mass == 475) return   56.94    ;
  else if ((int)mass == 476) return   57.36    ;
  else if ((int)mass == 477) return   57.78    ;
  else if ((int)mass == 478) return   58.21    ;
  else if ((int)mass == 479) return   58.63    ;
  else if ((int)mass == 480) return   59.06    ;
  else if ((int)mass == 481) return   59.49    ;
  else if ((int)mass == 482) return   59.92    ;
  else if ((int)mass == 483) return   60.36    ;
  else if ((int)mass == 484) return   60.79    ;
  else if ((int)mass == 485) return   61.23    ;
  else if ((int)mass == 486) return   61.67    ;
  else if ((int)mass == 487) return   62.10    ;
  else if ((int)mass == 488) return   62.55    ;
  else if ((int)mass == 489) return   62.99    ;
  else if ((int)mass == 490) return   63.43    ;
  else if ((int)mass == 491) return   63.88    ;
  else if ((int)mass == 492) return   64.32    ;
  else if ((int)mass == 493) return   64.77    ;
  else if ((int)mass == 494) return   65.22    ;
  else if ((int)mass == 495) return   65.67    ;
  else if ((int)mass == 496) return   66.12    ;
  else if ((int)mass == 497) return   66.58    ;
  else if ((int)mass == 498) return   67.03    ;
  else if ((int)mass == 499) return   67.49    ;
  else if ((int)mass == 500) return   67.95    ;
  else if ((int)mass == 501) return   68.41    ;
  else if ((int)mass == 502) return   68.87    ;
  else if ((int)mass == 503) return   69.33    ;
  else if ((int)mass == 504) return   69.80    ;
  else if ((int)mass == 505) return   70.27    ;
  else if ((int)mass == 506) return   70.73    ;
  else if ((int)mass == 507) return   71.20    ;
  else if ((int)mass == 508) return   71.68    ;
  else if ((int)mass == 509) return   72.15    ;
  else if ((int)mass == 510) return   72.62    ;
  else if ((int)mass == 511) return   73.10    ;
  else if ((int)mass == 512) return   73.58    ;
  else if ((int)mass == 513) return   74.06    ;
  else if ((int)mass == 514) return   74.54    ;
  else if ((int)mass == 515) return   75.02    ;
  else if ((int)mass == 516) return   75.50    ;
  else if ((int)mass == 517) return   75.99    ;
  else if ((int)mass == 518) return   76.48    ;
  else if ((int)mass == 519) return   76.97    ;
  else if ((int)mass == 520) return   77.46    ;
  else if ((int)mass == 521) return   77.95    ;
  else if ((int)mass == 522) return   78.44    ;
  else if ((int)mass == 523) return   78.94    ;
  else if ((int)mass == 524) return   79.44    ;
  else if ((int)mass == 525) return   79.94    ;
  else if ((int)mass == 526) return   80.44    ;
  else if ((int)mass == 527) return   80.94    ;
  else if ((int)mass == 528) return   81.44    ;
  else if ((int)mass == 529) return   81.95    ;
  else if ((int)mass == 530) return   82.46    ;
  else if ((int)mass == 531) return   82.97    ;
  else if ((int)mass == 532) return   83.48    ;
  else if ((int)mass == 533) return   83.99    ;
  else if ((int)mass == 534) return   84.50    ;
  else if ((int)mass == 535) return   85.02    ;
  else if ((int)mass == 536) return   85.54    ;
  else if ((int)mass == 537) return   86.06    ;
  else if ((int)mass == 538) return   86.58    ;
  else if ((int)mass == 539) return   87.10    ;
  else if ((int)mass == 540) return   87.63    ;
  else if ((int)mass == 541) return   88.15    ;
  else if ((int)mass == 542) return   88.68    ;
  else if ((int)mass == 543) return   89.21    ;
  else if ((int)mass == 544) return   89.74    ;
  else if ((int)mass == 545) return   90.28    ;
  else if ((int)mass == 546) return   90.81    ;
  else if ((int)mass == 547) return   91.35    ;
  else if ((int)mass == 548) return   91.89    ;
  else if ((int)mass == 549) return   92.43    ;
  else if ((int)mass == 550) return   92.97    ;
  else if ((int)mass == 551) return   93.51    ;
  else if ((int)mass == 552) return   94.06    ;
  else if ((int)mass == 553) return   94.61    ;
  else if ((int)mass == 554) return   95.16    ;
  else if ((int)mass == 555) return   95.71    ;
  else if ((int)mass == 556) return   96.26    ;
  else if ((int)mass == 557) return   96.82    ;
  else if ((int)mass == 558) return   97.38    ;
  else if ((int)mass == 559) return   97.93    ;
  else if ((int)mass == 560) return   98.50    ;
  else if ((int)mass == 561) return   99.06    ;
  else if ((int)mass == 562) return   99.62    ;
  else if ((int)mass == 563) return   100.2    ;
  else if ((int)mass == 564) return   100.8    ;
  else if ((int)mass == 565) return   101.3    ;
  else if ((int)mass == 566) return   101.9    ;
  else if ((int)mass == 567) return   102.5    ;
  else if ((int)mass == 568) return   103.0    ;
  else if ((int)mass == 569) return   103.6    ;
  else if ((int)mass == 570) return   104.2    ;
  else if ((int)mass == 571) return   104.8    ;
  else if ((int)mass == 572) return   105.4    ;
  else if ((int)mass == 573) return   106.0    ;
  else if ((int)mass == 574) return   106.5    ;
  else if ((int)mass == 575) return   107.1    ;
  else if ((int)mass == 576) return   107.7    ;
  else if ((int)mass == 577) return   108.3    ;
  else if ((int)mass == 578) return   108.9    ;
  else if ((int)mass == 579) return   109.5    ;
  else if ((int)mass == 580) return   110.1    ;
  else if ((int)mass == 581) return   110.7    ;
  else if ((int)mass == 582) return   111.3    ;
  else if ((int)mass == 583) return   111.9    ;
  else if ((int)mass == 584) return   112.5    ;
  else if ((int)mass == 585) return   113.1    ;
  else if ((int)mass == 586) return   113.7    ;
  else if ((int)mass == 587) return   114.4    ;
  else if ((int)mass == 588) return   115.0    ;
  else if ((int)mass == 589) return   115.6    ;
  else if ((int)mass == 590) return   116.2    ;
  else if ((int)mass == 591) return   116.8    ;
  else if ((int)mass == 592) return   117.5    ;
  else if ((int)mass == 593) return   118.1    ;
  else if ((int)mass == 594) return   118.7    ;
  else if ((int)mass == 595) return   119.3    ;
  else if ((int)mass == 596) return   120.0    ;
  else if ((int)mass == 597) return   120.6    ;
  else if ((int)mass == 598) return   121.2    ;
  else if ((int)mass == 599) return   121.9    ;
  else if ((int)mass == 600) return   122.5    ;
  else if ((int)mass == 601) return   123.2    ;
  else if ((int)mass == 602) return   123.8    ;
  else if ((int)mass == 603) return   124.4    ;
  else if ((int)mass == 604) return   125.1    ;
  else if ((int)mass == 605) return   125.7    ;
  else if ((int)mass == 606) return   126.4    ;
  else if ((int)mass == 607) return   127.1    ;
  else if ((int)mass == 608) return   127.7    ;
  else if ((int)mass == 609) return   128.4    ;
  else if ((int)mass == 610) return   129.0    ;
  else if ((int)mass == 611) return   129.7    ;
  else if ((int)mass == 612) return   130.4    ;
  else if ((int)mass == 613) return   131.0    ;
  else if ((int)mass == 614) return   131.7    ;
  else if ((int)mass == 615) return   132.4    ;
  else if ((int)mass == 616) return   133.0    ;
  else if ((int)mass == 617) return   133.7    ;
  else if ((int)mass == 618) return   134.4    ;
  else if ((int)mass == 619) return   135.1    ;
  else if ((int)mass == 620) return   135.8    ;
  else if ((int)mass == 621) return   136.4    ;
  else if ((int)mass == 622) return   137.1    ;
  else if ((int)mass == 623) return   137.8    ;
  else if ((int)mass == 624) return   138.5    ;
  else if ((int)mass == 625) return   139.2    ;
  else if ((int)mass == 626) return   139.9    ;
  else if ((int)mass == 627) return   140.6    ;
  else if ((int)mass == 628) return   141.3    ;
  else if ((int)mass == 629) return   142.0    ;
  else if ((int)mass == 630) return   142.7    ;
  else if ((int)mass == 631) return   143.4    ;
  else if ((int)mass == 632) return   144.1    ;
  else if ((int)mass == 633) return   144.8    ;
  else if ((int)mass == 634) return   145.6    ;
  else if ((int)mass == 635) return   146.3    ;
  else if ((int)mass == 636) return   147.0    ;
  else if ((int)mass == 637) return   147.7    ;
  else if ((int)mass == 638) return   148.4    ;
  else if ((int)mass == 639) return   149.2    ;
  else if ((int)mass == 640) return   149.9    ;
  else if ((int)mass == 641) return   150.6    ;
  else if ((int)mass == 642) return   151.4    ;
  else if ((int)mass == 643) return   152.1    ;
  else if ((int)mass == 644) return   152.8    ;
  else if ((int)mass == 645) return   153.6    ;
  else if ((int)mass == 646) return   154.3    ;
  else if ((int)mass == 647) return   155.1    ;
  else if ((int)mass == 648) return   155.8    ;
  else if ((int)mass == 649) return   156.6    ;
  else if ((int)mass == 650) return   157.3    ;
  else if ((int)mass == 651) return   158.1    ;
  else if ((int)mass == 652) return   158.8    ;
  else if ((int)mass == 653) return   159.6    ;
  else if ((int)mass == 654) return   160.3    ;
  else if ((int)mass == 655) return   161.1    ;
  else if ((int)mass == 656) return   161.9    ;
  else if ((int)mass == 657) return   162.6    ;
  else if ((int)mass == 658) return   163.4    ;
  else if ((int)mass == 659) return   164.2    ;
  else if ((int)mass == 660) return   165.0    ;
  else if ((int)mass == 661) return   165.7    ;
  else if ((int)mass == 662) return   166.5    ;
  else if ((int)mass == 663) return   167.3    ;
  else if ((int)mass == 664) return   168.1    ;
  else if ((int)mass == 665) return   168.9    ;
  else if ((int)mass == 666) return   169.7    ;
  else if ((int)mass == 667) return   170.5    ;
  else if ((int)mass == 668) return   171.3    ;
  else if ((int)mass == 669) return   172.1    ;
  else if ((int)mass == 670) return   172.9    ;
  else if ((int)mass == 671) return   173.7    ;
  else if ((int)mass == 672) return   174.5    ;
  else if ((int)mass == 673) return   175.3    ;
  else if ((int)mass == 674) return   176.1    ;
  else if ((int)mass == 675) return   176.9    ;
  else if ((int)mass == 676) return   177.7    ;
  else if ((int)mass == 677) return   178.6    ;
  else if ((int)mass == 678) return   179.4    ;
  else if ((int)mass == 679) return   180.2    ;
  else if ((int)mass == 680) return   181.0    ;
  else if ((int)mass == 681) return   181.9    ;
  else if ((int)mass == 682) return   182.7    ;
  else if ((int)mass == 683) return   183.5    ;
  else if ((int)mass == 684) return   184.4    ;
  else if ((int)mass == 685) return   185.2    ;
  else if ((int)mass == 686) return   186.1    ;
  else if ((int)mass == 687) return   186.9    ;
  else if ((int)mass == 688) return   187.8    ;
  else if ((int)mass == 689) return   188.6    ;
  else if ((int)mass == 690) return   189.5    ;
  else if ((int)mass == 691) return   190.3    ;
  else if ((int)mass == 692) return   191.2    ;
  else if ((int)mass == 693) return   192.1    ;
  else if ((int)mass == 694) return   192.9    ;
  else if ((int)mass == 695) return   193.8    ;
  else if ((int)mass == 696) return   194.7    ;
  else if ((int)mass == 697) return   195.5    ;
  else if ((int)mass == 698) return   196.4    ;
  else if ((int)mass == 699) return   197.3    ;
  else if ((int)mass == 700) return   198.2    ;
  else if ((int)mass == 701) return   199.1    ;
  else if ((int)mass == 702) return   200.0    ;
  else if ((int)mass == 703) return   200.8    ;
  else if ((int)mass == 704) return   201.7    ;
  else if ((int)mass == 705) return   202.6    ;
  else if ((int)mass == 706) return   203.5    ;
  else if ((int)mass == 707) return   204.4    ;
  else if ((int)mass == 708) return   205.3    ;
  else if ((int)mass == 709) return   206.3    ;
  else if ((int)mass == 710) return   207.2    ;
  else if ((int)mass == 711) return   208.1    ;
  else if ((int)mass == 712) return   209.0    ;
  else if ((int)mass == 713) return   209.9    ;
  else if ((int)mass == 714) return   210.8    ;
  else if ((int)mass == 715) return   211.8    ;
  else if ((int)mass == 716) return   212.7    ;
  else if ((int)mass == 717) return   213.6    ;
  else if ((int)mass == 718) return   214.6    ;
  else if ((int)mass == 719) return   215.5    ;
  else if ((int)mass == 720) return   216.4    ;
  else if ((int)mass == 721) return   217.4    ;
  else if ((int)mass == 722) return   218.3    ;
  else if ((int)mass == 723) return   219.3    ;
  else if ((int)mass == 724) return   220.2    ;
  else if ((int)mass == 725) return   221.2    ;
  else if ((int)mass == 726) return   222.2    ;
  else if ((int)mass == 727) return   223.1    ;
  else if ((int)mass == 728) return   224.1    ;
  else if ((int)mass == 729) return   225.0    ;
  else if ((int)mass == 730) return   226.0    ;
  else if ((int)mass == 731) return   227.0    ;
  else if ((int)mass == 732) return   228.0    ;
  else if ((int)mass == 733) return   229.0    ;
  else if ((int)mass == 734) return   229.9    ;
  else if ((int)mass == 735) return   230.9    ;
  else if ((int)mass == 736) return   231.9    ;
  else if ((int)mass == 737) return   232.9    ;
  else if ((int)mass == 738) return   233.9    ;
  else if ((int)mass == 739) return   234.9    ;
  else if ((int)mass == 740) return   235.9    ;
  else if ((int)mass == 741) return   236.9    ;
  else if ((int)mass == 742) return   237.9    ;
  else if ((int)mass == 743) return   238.9    ;
  else if ((int)mass == 744) return   239.9    ;
  else if ((int)mass == 745) return   241.0    ;
  else if ((int)mass == 746) return   242.0    ;
  else if ((int)mass == 747) return   243.0    ;
  else if ((int)mass == 748) return   244.0    ;
  else if ((int)mass == 749) return   245.1    ;
  else if ((int)mass == 750) return   246.1    ;
  else if ((int)mass == 751) return   247.1    ;
  else if ((int)mass == 752) return   248.2    ;
  else if ((int)mass == 753) return   249.2    ;
  else if ((int)mass == 754) return   250.3    ;
  else if ((int)mass == 755) return   251.3    ;
  else if ((int)mass == 756) return   252.4    ;
  else if ((int)mass == 757) return   253.4    ;
  else if ((int)mass == 758) return   254.5    ;
  else if ((int)mass == 759) return   255.6    ;
  else if ((int)mass == 760) return   256.6    ;
  else if ((int)mass == 761) return   257.7    ;
  else if ((int)mass == 762) return   258.8    ;
  else if ((int)mass == 763) return   259.9    ;
  else if ((int)mass == 764) return   260.9    ;
  else if ((int)mass == 765) return   262.0    ;
  else if ((int)mass == 766) return   263.1    ;
  else if ((int)mass == 767) return   264.2    ;
  else if ((int)mass == 768) return   265.3    ;
  else if ((int)mass == 769) return   266.4    ;
  else if ((int)mass == 770) return   267.5    ;
  else if ((int)mass == 771) return   268.6    ;
  else if ((int)mass == 772) return   269.7    ;
  else if ((int)mass == 773) return   270.8    ;
  else if ((int)mass == 774) return   272.0    ;
  else if ((int)mass == 775) return   273.1    ;
  else if ((int)mass == 776) return   274.2    ;
  else if ((int)mass == 777) return   275.3    ;
  else if ((int)mass == 778) return   276.5    ;
  else if ((int)mass == 779) return   277.6    ;
  else if ((int)mass == 780) return   278.7    ;
  else if ((int)mass == 781) return   279.9    ;
  else if ((int)mass == 782) return   281.0    ;
  else if ((int)mass == 783) return   282.2    ;
  else if ((int)mass == 784) return   283.3    ;
  else if ((int)mass == 785) return   284.5    ;
  else if ((int)mass == 786) return   285.6    ;
  else if ((int)mass == 787) return   286.8    ;
  else if ((int)mass == 788) return   288.0    ;
  else if ((int)mass == 789) return   289.1    ;
  else if ((int)mass == 790) return   290.3    ;
  else if ((int)mass == 791) return   291.5    ;
  else if ((int)mass == 792) return   292.7    ;
  else if ((int)mass == 793) return   293.8    ;
  else if ((int)mass == 794) return   295.0    ;
  else if ((int)mass == 795) return   296.2    ;
  else if ((int)mass == 796) return   297.4    ;
  else if ((int)mass == 797) return   298.6    ;
  else if ((int)mass == 798) return   299.8    ;
  else if ((int)mass == 799) return   301.0    ;
  else if ((int)mass == 800) return   302.2    ;
  else if ((int)mass == 801) return   303.5    ;
  else if ((int)mass == 802) return   304.7    ;
  else if ((int)mass == 803) return   305.9    ;
  else if ((int)mass == 804) return   307.1    ;
  else if ((int)mass == 805) return   308.4    ;
  else if ((int)mass == 806) return   309.6    ;
  else if ((int)mass == 807) return   310.8    ;
  else if ((int)mass == 808) return   312.1    ;
  else if ((int)mass == 809) return   313.3    ;
  else if ((int)mass == 810) return   314.6    ;
  else if ((int)mass == 811) return   315.8    ;
  else if ((int)mass == 812) return   317.1    ;
  else if ((int)mass == 813) return   318.3    ;
  else if ((int)mass == 814) return   319.6    ;
  else if ((int)mass == 815) return   320.9    ;
  else if ((int)mass == 816) return   322.1    ;
  else if ((int)mass == 817) return   323.4    ;
  else if ((int)mass == 818) return   324.7    ;
  else if ((int)mass == 819) return   326.0    ;
  else if ((int)mass == 820) return   327.3    ;
  else if ((int)mass == 821) return   328.6    ;
  else if ((int)mass == 822) return   329.9    ;
  else if ((int)mass == 823) return   331.2    ;
  else if ((int)mass == 824) return   332.5    ;
  else if ((int)mass == 825) return   333.8    ;
  else if ((int)mass == 826) return   335.1    ;
  else if ((int)mass == 827) return   336.4    ;
  else if ((int)mass == 828) return   337.7    ;
  else if ((int)mass == 829) return   339.1    ;
  else if ((int)mass == 830) return   340.4    ;
  else if ((int)mass == 831) return   341.7    ;
  else if ((int)mass == 832) return   343.1    ;
  else if ((int)mass == 833) return   344.4    ;
  else if ((int)mass == 834) return   345.8    ;
  else if ((int)mass == 835) return   347.1    ;
  else if ((int)mass == 836) return   348.5    ;
  else if ((int)mass == 837) return   349.8    ;
  else if ((int)mass == 838) return   351.2    ;
  else if ((int)mass == 839) return   352.6    ;
  else if ((int)mass == 840) return   353.9    ;
  else if ((int)mass == 841) return   355.3    ;
  else if ((int)mass == 842) return   356.7    ;
  else if ((int)mass == 843) return   358.1    ;
  else if ((int)mass == 844) return   359.5    ;
  else if ((int)mass == 845) return   360.9    ;
  else if ((int)mass == 846) return   362.3    ;
  else if ((int)mass == 847) return   363.7    ;
  else if ((int)mass == 848) return   365.1    ;
  else if ((int)mass == 849) return   366.5    ;
  else if ((int)mass == 850) return   367.9    ;
  else if ((int)mass == 851) return   369.3    ;
  else if ((int)mass == 852) return   370.8    ;
  else if ((int)mass == 853) return   372.2    ;
  else if ((int)mass == 854) return   373.6    ;
  else if ((int)mass == 855) return   375.1    ;
  else if ((int)mass == 856) return   376.5    ;
  else if ((int)mass == 857) return   378.0    ;
  else if ((int)mass == 858) return   379.4    ;
  else if ((int)mass == 859) return   380.9    ;
  else if ((int)mass == 860) return   382.3    ;
  else if ((int)mass == 861) return   383.8    ;
  else if ((int)mass == 862) return   385.3    ;
  else if ((int)mass == 863) return   386.7    ;
  else if ((int)mass == 864) return   388.2    ;
  else if ((int)mass == 865) return   389.7    ;
  else if ((int)mass == 866) return   391.2    ;
  else if ((int)mass == 867) return   392.7    ;
  else if ((int)mass == 868) return   394.2    ;
  else if ((int)mass == 869) return   395.7    ;
  else if ((int)mass == 870) return   397.2    ;
  else if ((int)mass == 871) return   398.7    ;
  else if ((int)mass == 872) return   400.2    ;
  else if ((int)mass == 873) return   401.8    ;
  else if ((int)mass == 874) return   403.3    ;
  else if ((int)mass == 875) return   404.8    ;
  else if ((int)mass == 876) return   406.4    ;
  else if ((int)mass == 877) return   407.9    ;
  else if ((int)mass == 878) return   409.4    ;
  else if ((int)mass == 879) return   411.0    ;
  else if ((int)mass == 880) return   412.5    ;
  else if ((int)mass == 881) return   414.1    ;
  else if ((int)mass == 882) return   415.7    ;
  else if ((int)mass == 883) return   417.2    ;
  else if ((int)mass == 884) return   418.8    ;
  else if ((int)mass == 885) return   420.4    ;
  else if ((int)mass == 886) return   422.0    ;
  else if ((int)mass == 887) return   423.6    ;
  else if ((int)mass == 888) return   425.2    ;
  else if ((int)mass == 889) return   426.8    ;
  else if ((int)mass == 890) return   428.4    ;
  else if ((int)mass == 891) return   430.0    ;
  else if ((int)mass == 892) return   431.6    ;
  else if ((int)mass == 893) return   433.2    ;
  else if ((int)mass == 894) return   434.9    ;
  else if ((int)mass == 895) return   436.5    ;
  else if ((int)mass == 896) return   438.1    ;
  else if ((int)mass == 897) return   439.8    ;
  else if ((int)mass == 898) return   441.4    ;
  else if ((int)mass == 899) return   443.1    ;
  else if ((int)mass == 900) return   444.7    ;
  else if ((int)mass == 901) return   446.4    ;
  else if ((int)mass == 902) return   448.1    ;
  else if ((int)mass == 903) return   449.7    ;
  else if ((int)mass == 904) return   451.4    ;
  else if ((int)mass == 905) return   453.1    ;
  else if ((int)mass == 906) return   454.8    ;
  else if ((int)mass == 907) return   456.5    ;
  else if ((int)mass == 908) return   458.2    ;
  else if ((int)mass == 909) return   459.9    ;
  else if ((int)mass == 910) return   461.6    ;
  else if ((int)mass == 911) return   463.3    ;
  else if ((int)mass == 912) return   465.0    ;
  else if ((int)mass == 913) return   466.7    ;
  else if ((int)mass == 914) return   468.5    ;
  else if ((int)mass == 915) return   470.2    ;
  else if ((int)mass == 916) return   472.0    ;
  else if ((int)mass == 917) return   473.7    ;
  else if ((int)mass == 918) return   475.5    ;
  else if ((int)mass == 919) return   477.2    ;
  else if ((int)mass == 920) return   479.0    ;
  else if ((int)mass == 921) return   480.7    ;
  else if ((int)mass == 922) return   482.5    ;
  else if ((int)mass == 923) return   484.3    ;
  else if ((int)mass == 924) return   486.1    ;
  else if ((int)mass == 925) return   487.9    ;
  else if ((int)mass == 926) return   489.7    ;
  else if ((int)mass == 927) return   491.5    ;
  else if ((int)mass == 928) return   493.3    ;
  else if ((int)mass == 929) return   495.1    ;
  else if ((int)mass == 930) return   496.9    ;
  else if ((int)mass == 931) return   498.8    ;
  else if ((int)mass == 932) return   500.6    ;
  else if ((int)mass == 933) return   502.4    ;
  else if ((int)mass == 934) return   504.3    ;
  else if ((int)mass == 935) return   506.1    ;
  else if ((int)mass == 936) return   508.0    ;
  else if ((int)mass == 937) return   509.8    ;
  else if ((int)mass == 938) return   511.7    ;
  else if ((int)mass == 939) return   513.6    ;
  else if ((int)mass == 940) return   515.4    ;
  else if ((int)mass == 941) return   517.3    ;
  else if ((int)mass == 942) return   519.2    ;
  else if ((int)mass == 943) return   521.1    ;
  else if ((int)mass == 944) return   523.0    ;
  else if ((int)mass == 945) return   524.9    ;
  else if ((int)mass == 946) return   526.8    ;
  else if ((int)mass == 947) return   528.8    ;
  else if ((int)mass == 948) return   530.7    ;
  else if ((int)mass == 949) return   532.6    ;
  else if ((int)mass == 950) return   534.6    ;
  else if ((int)mass == 951) return   536.5    ;
  else if ((int)mass == 952) return   538.5    ;
  else if ((int)mass == 953) return   540.4    ;
  else if ((int)mass == 954) return   542.4    ;
  else if ((int)mass == 955) return   544.3    ;
  else if ((int)mass == 956) return   546.3    ;
  else if ((int)mass == 957) return   548.3    ;
  else if ((int)mass == 958) return   550.3    ;
  else if ((int)mass == 959) return   552.3    ;
  else if ((int)mass == 960) return   554.3    ;
  else if ((int)mass == 961) return   556.3    ;
  else if ((int)mass == 962) return   558.3    ;
  else if ((int)mass == 963) return   560.3    ;
  else if ((int)mass == 964) return   562.3    ;
  else if ((int)mass == 965) return   564.4    ;
  else if ((int)mass == 966) return   566.4    ;
  else if ((int)mass == 967) return   568.5    ;
  else if ((int)mass == 968) return   570.5    ;
  else if ((int)mass == 969) return   572.6    ;
  else if ((int)mass == 970) return   574.6    ;
  else if ((int)mass == 971) return   576.7    ;
  else if ((int)mass == 972) return   578.8    ;
  else if ((int)mass == 973) return   580.9    ;
  else if ((int)mass == 974) return   582.9    ;
  else if ((int)mass == 975) return   585.0    ;
  else if ((int)mass == 976) return   587.1    ;
  else if ((int)mass == 977) return   589.3    ;
  else if ((int)mass == 978) return   591.4    ;
  else if ((int)mass == 979) return   593.5    ;
  else if ((int)mass == 980) return   595.6    ;
  else if ((int)mass == 981) return   597.8    ;
  else if ((int)mass == 982) return   599.9    ;
  else if ((int)mass == 983) return   602.1    ;
  else if ((int)mass == 984) return   604.2    ;
  else if ((int)mass == 985) return   606.4    ;
  else if ((int)mass == 986) return   608.5    ;
  else if ((int)mass == 987) return   610.7    ;
  else if ((int)mass == 988) return   612.9    ;
  else if ((int)mass == 989) return   615.1    ;
  else if ((int)mass == 990) return   617.3    ;
  else if ((int)mass == 991) return   619.5    ;
  else if ((int)mass == 992) return   621.7    ;
  else if ((int)mass == 993) return   623.9    ;
  else if ((int)mass == 994) return   626.2    ;
  else if ((int)mass == 995) return   628.4    ;
  else if ((int)mass == 996) return   630.6    ;
  else if ((int)mass == 997) return   632.9    ;
  else if ((int)mass == 998) return   635.1    ;
  else if ((int)mass == 999) return   637.4    ;
  else if ((int)mass == 1000) return  639.7    ;
  else {
    std::cout << "[RooBWHighMassGGH::H_width]: mass is out of range" << std::endl;
      return 0.;
  }
}





/*************RooCPSHighMassGGHNoInterf***********/
ClassImp(RooCPSHighMassGGHNoInterf) 

 RooCPSHighMassGGHNoInterf::RooCPSHighMassGGHNoInterf(){}

 RooCPSHighMassGGHNoInterf::RooCPSHighMassGGHNoInterf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _KPrime,
			RooAbsReal& _BRnew,
			Bool_t is8TeV_b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   KPrime("KPrime","KPrime",this,_KPrime),
   BRnew("BRnew","BRnew",this,_BRnew),
   is8TeV(is8TeV_b)
 {    
   initMatrices();
 } 


 RooCPSHighMassGGHNoInterf::RooCPSHighMassGGHNoInterf(const RooCPSHighMassGGHNoInterf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   KPrime("KPrime",this,other.KPrime),
   BRnew("BRnew",this,other.BRnew),
   is8TeV(other.is8TeV)
 { 
   
   for (Int_t i=0; i<7; ++i){
     for (Int_t j=0; j<5; ++j){
       a_width[i][j] = other.a_width[i][j];
       a_delta[i][j] = other.a_delta[i][j];
     }
   }

 } 



 Double_t RooCPSHighMassGGHNoInterf::evaluate() const 
 {    

   Double_t effKPrime = KPrime/(1-BRnew);

   Double_t width = getWidth(mH,effKPrime);
   Double_t delta = getDelta(mH,effKPrime);

   Double_t bwHM = x / ( TMath::Power( TMath::Power(x,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(x,2)*TMath::Power(effKPrime*width,2) ); 

   Double_t splineFactor;
   if (x<1850) splineFactor = Spline(x);
   else splineFactor = Spline(1850);

   Double_t signal = bwHM*splineFactor;
     
   Double_t fValue = signal;
   if (fValue > 0) return fValue;
   else return 0.;   
   
 } 




void RooCPSHighMassGGHNoInterf::initMatrices() {


  if (is8TeV){

    a_delta[0][0] = 0.015081577003;
    a_delta[0][1] = 0.0576512292027;
    a_delta[0][2] = 0.131081297994;
    a_delta[0][3] = 0.243380963802;
    a_delta[0][4] = 0.39500105381;
    a_delta[1][0] = 0.0972077324986;
    a_delta[1][1] = 0.434791147709;
    a_delta[1][2] = 0.861800074577;
    a_delta[1][3] = 1.58613061905;
    a_delta[1][4] = 2.41115784645;
    a_delta[2][0] = 0.242783889174;
    a_delta[2][1] = 0.827510118484;
    a_delta[2][2] = 1.78883588314;
    a_delta[2][3] = 3.17804527283;
    a_delta[2][4] = 4.9168419838;
    a_delta[3][0] = 0.368439763784;
    a_delta[3][1] = 1.45783913136;
    a_delta[3][2] = 3.27468585968;
    a_delta[3][3] = 5.93577003479;
    a_delta[3][4] = 9.49385356903;
    a_delta[4][0] = 0.679009974003;
    a_delta[4][1] = 2.96808123589;
    a_delta[4][2] = 6.60964488983;
    a_delta[4][3] = 11.7470684052;
    a_delta[4][4] = 18.1083869934;
    a_delta[5][0] = 1.51261091232;
    a_delta[5][1] = 5.85586261749;
    a_delta[5][2] = 12.3102693558;
    a_delta[5][3] = 21.393032074;
    a_delta[5][4] = 31.7420043945;
    a_delta[6][0] = 1.7427560091;
    a_delta[6][1] = 7.33418560028;
    a_delta[6][2] = 15.4683561325;
    a_delta[6][3] = 26.5165843964;
    a_delta[6][4] = 37.535987854;

    a_width[0][0] = 26.5976848602;
    a_width[0][1] = 26.5976848602;
    a_width[0][2] = 26.7212715149;
    a_width[0][3] = 26.5424861908;
    a_width[0][4] = 26.4101047516;
    a_width[1][0] = 58.7029762268;
    a_width[1][1] = 58.7029762268;
    a_width[1][2] = 58.7458114624;
    a_width[1][3] = 58.1389541626;
    a_width[1][4] = 58.612991333;
    a_width[2][0] = 103.93271637;
    a_width[2][1] = 103.93271637;
    a_width[2][2] = 103.803878784;
    a_width[2][3] = 103.800476074;
    a_width[2][4] = 103.74822998;
    a_width[3][0] = 162.97102356;
    a_width[3][1] = 162.97102356;
    a_width[3][2] = 162.841094971;
    a_width[3][3] = 162.802963257;
    a_width[3][4] = 162.688156128;
    a_width[4][0] = 235.570587158;
    a_width[4][1] = 235.570587158;
    a_width[4][2] = 235.355682373;
    a_width[4][3] = 235.014694214;
    a_width[4][4] = 235.070922852;
    a_width[5][0] = 320.552429199;
    a_width[5][1] = 320.552429199;
    a_width[5][2] = 320.081390381;
    a_width[5][3] = 319.721862793;
    a_width[5][4] = 319.82611084;
    a_width[6][0] = 416.118835449;
    a_width[6][1] = 416.118835449;
    a_width[6][2] = 430.821716309;
    a_width[6][3] = 416.076812744;
    a_width[6][4] = 415.9637146;

  }
  else {

    a_delta[0][0] = 0.015081577003;
    a_delta[0][1] = 0.0576512292027;
    a_delta[0][2] = 0.131081297994;
    a_delta[0][3] = 0.243380963802;
    a_delta[0][4] = 0.39500105381;
    a_delta[1][0] = 0.0972077324986;
    a_delta[1][1] = 0.434791147709;
    a_delta[1][2] = 0.861800074577;
    a_delta[1][3] = 1.58613061905;
    a_delta[1][4] = 2.41115784645;
    a_delta[2][0] = 0.242783889174;
    a_delta[2][1] = 0.827510118484;
    a_delta[2][2] = 1.78883588314;
    a_delta[2][3] = 3.17804527283;
    a_delta[2][4] = 4.9168419838;
    a_delta[3][0] = 0.368439763784;
    a_delta[3][1] = 1.45783913136;
    a_delta[3][2] = 3.27468585968;
    a_delta[3][3] = 5.93577003479;
    a_delta[3][4] = 9.49385356903;
    a_delta[4][0] = 0.679009974003;
    a_delta[4][1] = 2.96808123589;
    a_delta[4][2] = 6.60964488983;
    a_delta[4][3] = 11.7470684052;
    a_delta[4][4] = 18.1083869934;
    a_delta[5][0] = 1.51261091232;
    a_delta[5][1] = 5.85586261749;
    a_delta[5][2] = 12.3102693558;
    a_delta[5][3] = 21.393032074;
    a_delta[5][4] = 31.7420043945;
    a_delta[6][0] = 1.7427560091;
    a_delta[6][1] = 7.33418560028;
    a_delta[6][2] = 15.4683561325;
    a_delta[6][3] = 26.5165843964;
    a_delta[6][4] = 37.535987854;


    a_width[0][0] = 26.5976848602;
    a_width[0][1] = 26.5976848602;
    a_width[0][2] = 26.6827316284;
    a_width[0][3] = 26.6214389801;
    a_width[0][4] = 26.3873176575;
    a_width[1][0] = 58.7029762268;
    a_width[1][1] = 58.7029762268;
    a_width[1][2] = 58.8415145874;
    a_width[1][3] = 58.4030990601;
    a_width[1][4] = 58.6177330017;
    a_width[2][0] = 103.93271637;
    a_width[2][1] = 103.93271637;
    a_width[2][2] = 104.202682495;
    a_width[2][3] = 103.806724548;
    a_width[2][4] = 103.740852356;
    a_width[3][0] = 162.97102356;
    a_width[3][1] = 162.97102356;
    a_width[3][2] = 162.97102356;
    a_width[3][3] = 162.682937622;
    a_width[3][4] = 162.612564087;
    a_width[4][0] = 235.570587158;
    a_width[4][1] = 235.570587158;
    a_width[4][2] = 235.356552124;
    a_width[4][3] = 234.996078491;
    a_width[4][4] = 235.060302734;
    a_width[5][0] = 320.552429199;
    a_width[5][1] = 320.552429199;
    a_width[5][2] = 344.464599609;
    a_width[5][3] = 320.296142578;
    a_width[5][4] = 320.367614746;
    a_width[6][0] = 416.118835449;
    a_width[6][1] = 416.118835449;
    a_width[6][2] = 464.375823975;
    a_width[6][3] = 455.52545166;
    a_width[6][4] = 415.974060059;
    
  }

  
}



Double_t RooCPSHighMassGGHNoInterf::interpolateMatrix(const Double_t matrix[][5] , const Double_t& x, const Double_t& y) const
{

  // Dummy copy of ROOT TH2::Interpolate needed because there the function is not constant

  static const Int_t nMass   = 7;
  static const Int_t nCPrime = 5;

  Double_t masses[nMass]    = {400.,500.,600.,700.,800.,900.,1000.};
  Double_t cprimes[nCPrime] = {0.2,0.4,0.6,0.8,1.0};

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  Double_t stepX = masses[1]-masses[0];
  Double_t stepY = cprimes[1]-cprimes[0];
  Double_t lowX = masses[0]-stepX/2;  //Double_t highX = masses[nMass-1]+stepX/2;
  Double_t lowY = cprimes[0]-stepY/2; //Double_t highY = cprimes[nMass-1]+stepY/2;
  Int_t nX = nMass;
  Int_t nY = nCPrime;

  // Get the indexes to identify matrix element
  Int_t i_x = (Int_t)((x-lowX)/stepX);
  Int_t i_y = (Int_t)((y-lowY)/stepY);

  // Check if we are in the domain for the interpolation
  if(i_x<0 || i_x>(nX-1) || i_y<0 || i_y>(nY-1)){
    // std::cout << "[RooCPSHighMassGGHNoInterf::interpolateMatrix] Cannot interpolate outside array domain, returning 0." << std::endl;
    return 0;
  }

  // Find out the quadrant in the bin
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  dx = masses[i_x]+stepX/2-x;
  dy = cprimes[i_y]+stepY/2-y;
  if (dx<=stepX/2 && dy<=stepY/2)
    quadrant = 1; // upper right
  if (dx>stepX/2 && dy<=stepY/2)
    quadrant = 2; // upper left
  if (dx>stepX/2 && dy>stepY/2)
    quadrant = 3; // lower left
  if (dx<=stepX/2 && dy>stepY/2)
    quadrant = 4; // lower right

  switch(quadrant) {
  case 1:
    x1 = masses[i_x];
    y1 = cprimes[i_y];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 2:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    y1 = cprimes[i_y];
    x2 = masses[i_x];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 3:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    x2 = masses[i_x];
    y2 = cprimes[i_y];
    break;
  case 4:
    x1 = masses[i_x];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    y2 = cprimes[i_y];
    break;
  }

  Int_t i_x1 = (Int_t)((x1-lowX)/stepX);
  if(i_x1<0) i_x1=0;
  Int_t i_x2 = (Int_t)((x2-lowX)/stepX);
  if(i_x2>nMass-1) i_x2=nMass-1;
  Int_t i_y1 = (Int_t)((y1-lowY)/stepY);
  if(i_y1<0) i_y1=0;
  Int_t i_y2 = (Int_t)((y2-lowY)/stepY);
  if(i_y2>nCPrime-1) i_y2=nCPrime-1;
  Double_t q11 = matrix[i_x1][i_y1];
  Double_t q12 = matrix[i_x1][i_y2];
  Double_t q21 = matrix[i_x2][i_y1];
  Double_t q22 = matrix[i_x2][i_y2];
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  
  // std::cout << "[RooCPSHighMassGGHNoInterf::interpolateMatrix] DEBUG: " << x << " - " << y << "->" << q11 << " / " << q12 << " / " << q21 << " / " << q22 << " / " << f << std::endl;
  
  return f;
}


Double_t RooCPSHighMassGGHNoInterf::getWidth(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_width,mass,cprime);}

Double_t RooCPSHighMassGGHNoInterf::getDelta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_delta,mass,cprime);}




Double_t RooCPSHighMassGGHNoInterf::Spline(Double_t xx) const{
   const int fNp = 20, fKstep = 0;
   const double fDelta = -1, fXmin = 200, fXmax = 2150;
   const double fX[20] = { 200, 250, 300, 350, 400,
                        450, 500, 550, 600, 650,
                        750, 850, 950, 1050, 1150,
                        1350, 1550, 1750, 1950, 2150 };
   const double fY[20] = { 0.139988, 0.17804, 0.352505, 0.801982, 1.35134,
                        1.75174, 1.84163, 1.86637, 1.78052, 1.7297,
                        1.53036, 1.31241, 1.11751, 0.946671, 0.812524,
                        0.640316, 0.490843, 0.396263, 0.315344, 0.210836 };
   const double fB[20] = { 0.000266891, 0.00147333, 0.00628998, 0.0108012, 0.0105333,
                        0.00444666, 0.000737389, -0.000594766, -0.00180599, -0.000721876,
                        -0.00269494, -0.00186614, -0.00192752, -0.00150674, -0.00115716,
                        -0.000743197, -0.000649875, -0.000369596, -0.000461366, -0.000583113 };
   const double fC[20] = { 8.42921e-06, 2.29701e-05, 6.55918e-05, 1.67598e-05, -2.50909e-05,
                        -7.98341e-05, 4.46636e-06, -3.60991e-05, 1.78298e-05, -7.22121e-06,
                        3.25258e-06, -4.32435e-07, 1.81479e-06, 1.62255e-06, 2.04841e-06,
                        -2.33522e-07, 1.01254e-06, 5.76421e-08, -3.11836e-07, -3.01879e-07 };
   const double fD[20] = { 0, 2.90818e-07, -2.00191e-08, -3.34947e-07, -4.53847e-07,
                        2.18457e-07, 1.71166e-07, -2.84202e-08, 2.09784e-07, -2.33159e-07,
                        8.20808e-08, -2.72765e-08, 1.26473e-08, -2.76285e-09, 7.40335e-10,
                        -1.80875e-09, 1.31533e-09, -1.99675e-09, 4.97867e-11, 0 };
   const double fE[20] = { 0, 2.90818e-09, -6.01655e-09, 2.86727e-09, -4.05627e-09,
                        1.07793e-08, -1.12522e-08, 9.25635e-09, -6.87431e-09, 2.44487e-09,
                        -8.68676e-10, 3.2189e-10, -1.22271e-10, 4.52199e-11, -2.77039e-11,
                        2.13312e-11, -1.3521e-11, 5.24081e-12, -1.24467e-13, 0 };
   const double fF[20] = { 1.16327e-11, -3.56989e-11, 3.55353e-11, -2.76942e-11, 5.93423e-11,
                        -8.81261e-11, 8.20343e-11, -6.45226e-11, 3.72767e-11, -6.6271e-12,
                        2.38113e-12, -8.8832e-13, 3.34981e-13, -1.45848e-13, 4.90352e-14,
                        -3.48522e-14, 1.87618e-14, -5.36528e-15, 1.24467e-16, 0 };
   int klow=0;
   if(xx<=fXmin) klow=0;
   else if(xx>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((xx-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=xx-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*(fD[klow]+dx*(fE[klow]+dx*fF[klow])))));
}





/*************RooCPSHighMassVBF***********/
ClassImp(RooCPSHighMassVBF) 

 RooCPSHighMassVBF::RooCPSHighMassVBF(){}

 RooCPSHighMassVBF::RooCPSHighMassVBF(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _KPrime,
			RooAbsReal& _BRnew,
			RooAbsReal& _IntStr,
			RooAbsReal& _WidthScl,
			Bool_t is8TeV_b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   KPrime("KPrime","KPrime",this,_KPrime),
   BRnew("BRnew","BRnew",this,_BRnew),
   IntStr("IntStr","IntStr",this,_IntStr),
   WidthScl("WidthScl","WidthScl",this,_WidthScl),
   is8TeV(is8TeV_b)
 { 
   if (KPrime/(1-BRnew)<0.1) std::cout << "[RooCPSHighMassVBF] WARNING: parameters of the shape will be extrapolated outside the grid, (KPrime/(1-BRnew)=" << KPrime/(1-BRnew) << std::endl;
   initMatrices();
 } 


 RooCPSHighMassVBF::RooCPSHighMassVBF(const RooCPSHighMassVBF& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   KPrime("KPrime",this,other.KPrime),
   BRnew("BRnew",this,other.BRnew),
   IntStr("IntStr",this,other.IntStr),
   WidthScl("WidthScl",this,other.WidthScl),
   is8TeV(other.is8TeV)
 { 

   for (Int_t i=0; i<7; ++i){
     for (Int_t j=0; j<5; ++j){
       a_width[i][j] = other.a_width[i][j];
       a_delta[i][j] = other.a_delta[i][j];
       a_alpha[i][j] = other.a_alpha[i][j];
       a_r[i][j] = other.a_r[i][j];
       a_beta[i][j] = other.a_beta[i][j];
     }
   }

 } 



 Double_t RooCPSHighMassVBF::evaluate() const 
 {    

   Double_t effKPrime = KPrime/(1-BRnew);

   Double_t width;
   Double_t delta;
   Double_t alpha;
   Double_t r;
   Double_t beta ;

   if (KPrime/(1-BRnew)<0.1){
     width = WidthScl*getWidth(mH,0.1);
     delta = getDelta(mH,0.1);
     alpha = getAlpha(mH,0.1);
     r = getR(mH,0.1);
     beta = getBeta(mH,0.1);
   } else {
     width = WidthScl*getWidth(mH,effKPrime);
     delta = getDelta(mH,effKPrime);
     alpha = getAlpha(mH,effKPrime);
     r = getR(mH,effKPrime);
     beta = getBeta(mH,effKPrime);
   }

   // Cast the "x" proxy to double only once to avoid overhead and ambiguous
   // overloads the TComplex operators.
   const double xVal = x;

   Double_t bwHM = xVal / ( TMath::Power( TMath::Power(xVal,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(xVal,2)*TMath::Power(effKPrime*width,2) ); 

   Double_t splineFactor;
   if (xVal<1850) splineFactor = Spline(xVal);
   else splineFactor = Spline(1850);

   Double_t signal = bwHM*splineFactor;

   Double_t k=0.25;      
   Double_t mH_eff = mH*TMath::Sqrt(1-k*(effKPrime*width/mH)*(effKPrime*width/mH));
   TComplex MSquared(mH_eff*mH_eff,-effKPrime*width*mH_eff);
   TComplex M = MSquared.Sqrt(MSquared);
   TComplex Exp1 = MSquared.Exp((TComplex)alpha);
   TComplex Exp2 = MSquared.Exp(-(TComplex)alpha);
   TComplex Exp3 = MSquared.Exp(-xVal*(TComplex)alpha/M);

   double interference = -r*(1-TMath::Exp(-beta*(xVal-150.)/mH_eff))*((Exp1/(xVal-M)-Exp2/(xVal+M))*Exp3).Re();

   Double_t fValue = signal*(1+IntStr*interference/(bwHM*( TMath::Sqrt(1-BRnew) ))); 
   if (fValue > 0) return fValue;
   else return 0.;
   
   
 } 


void RooCPSHighMassVBF::initMatrices() {

  if (is8TeV){

    a_delta[0][0] = -0.0185139477253;
    a_delta[0][1] = 0.0201042648405;
    a_delta[0][2] = 0.156271114945;
    a_delta[0][3] = 0.310714632273;
    a_delta[0][4] = 0.622391521931;
    a_delta[1][0] = 0.0428073592484;
    a_delta[1][1] = 0.15564674139;
    a_delta[1][2] = 0.44218480587;
    a_delta[1][3] = 0.845812380314;
    a_delta[1][4] = 1.2718859911;
    a_delta[2][0] = 0.1351634413;
    a_delta[2][1] = 0.623583912849;
    a_delta[2][2] = 1.33247685432;
    a_delta[2][3] = 2.8496837616;
    a_delta[2][4] = 4.2275352478;
    a_delta[3][0] = 0.407202005386;
    a_delta[3][1] = 1.78870689869;
    a_delta[3][2] = 4.4582362175;
    a_delta[3][3] = 8.72150039673;
    a_delta[3][4] = 12.9632072449;
    a_delta[4][0] = 1.10362958908;
    a_delta[4][1] = 1.65279102325;
    a_delta[4][2] = 2.63262009621;
    a_delta[4][3] = 5.56618690491;
    a_delta[4][4] = 8.73703479767;
    a_delta[5][0] = 2.47303795815;
    a_delta[5][1] = 3.7929110527;
    a_delta[5][2] = 0.907296419144;
    a_delta[5][3] = -3.98578429222;
    a_delta[5][4] = -1.68104398251;
    a_delta[6][0] = 3.5346019268;
    a_delta[6][1] = 14.8740081787;
    a_delta[6][2] = 22.4615650177;
    a_delta[6][3] = 26.1576747894;
    a_delta[6][4] = 25.2779407501;


    a_alpha[0][0] = 7.81095170975;
    a_alpha[0][1] = 0.0101348012686;
    a_alpha[0][2] = 0.0100047122687;
    a_alpha[0][3] = 0.0101482300088;
    a_alpha[0][4] = 0.0108778737485;
    a_alpha[1][0] = 5.38699150085;
    a_alpha[1][1] = 1.85048913956;
    a_alpha[1][2] = 0.77756267786;
    a_alpha[1][3] = 0.0100484453142;
    a_alpha[1][4] = 0.0100035723299;
    a_alpha[2][0] = 4.18369340897;
    a_alpha[2][1] = 1.64176881313;
    a_alpha[2][2] = 1.22510075569;
    a_alpha[2][3] = 1.17058205605;
    a_alpha[2][4] = 0.258560866117;
    a_alpha[3][0] = 0.0100004831329;
    a_alpha[3][1] = 0.375370591879;
    a_alpha[3][2] = 0.327079296112;
    a_alpha[3][3] = 0.578940987587;
    a_alpha[3][4] = 0.543993651867;
    a_alpha[4][0] = 0.737967014313;
    a_alpha[4][1] = 1.17760062218;
    a_alpha[4][2] = 0.787738025188;
    a_alpha[4][3] = 0.742407977581;
    a_alpha[4][4] = 0.761043250561;
    a_alpha[5][0] = 1.90185630322;
    a_alpha[5][1] = 1.78479862213;
    a_alpha[5][2] = 1.24611675739;
    a_alpha[5][3] = 0.713018894196;
    a_alpha[5][4] = 0.814603388309;
    a_alpha[6][0] = 2.93534183502;
    a_alpha[6][1] = 2.2874276638;
    a_alpha[6][2] = 1.45080220699;
    a_alpha[6][3] = 0.476683199406;
    a_alpha[6][4] = 0.556163012981;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 1.22830260807e-05;
    a_r[0][1] = 8.98443158803e-06;
    a_r[0][2] = 4.95080030305e-06;
    a_r[0][3] = 2.42228702518e-06;
    a_r[0][4] = 8.00000350409e-07;
    a_r[1][0] = 9.54647111939e-06;
    a_r[1][1] = 6.83476037011e-06;
    a_r[1][2] = 3.71724195247e-06;
    a_r[1][3] = 2.55762461165e-06;
    a_r[1][4] = 1.65574124367e-06;
    a_r[2][0] = 6.41826864012e-06;
    a_r[2][1] = 4.66322171633e-06;
    a_r[2][2] = 2.57499277723e-06;
    a_r[2][3] = 1.95199231712e-06;
    a_r[2][4] = 1.3401167962e-06;
    a_r[3][0] = 4.33439527114e-06;
    a_r[3][1] = 3.28254350279e-06;
    a_r[3][2] = 1.93596497411e-06;
    a_r[3][3] = 1.59018213708e-06;
    a_r[3][4] = 1.25049041344e-06;
    a_r[4][0] = 4.21222694058e-06;
    a_r[4][1] = 2.9631075904e-06;
    a_r[4][2] = 1.67966402387e-06;
    a_r[4][3] = 1.31089291244e-06;
    a_r[4][4] = 1.01920193174e-06;
    a_r[5][0] = 3.77969490728e-06;
    a_r[5][1] = 2.66709412244e-06;
    a_r[5][2] = 1.36012192797e-06;
    a_r[5][3] = 9.07517403448e-07;
    a_r[5][4] = 7.24386666207e-07;
    a_r[6][0] = 3.1886186207e-06;
    a_r[6][1] = 2.66552501671e-06;
    a_r[6][2] = 1.34276660901e-06;
    a_r[6][3] = 8.55406483424e-07;
    a_r[6][4] = 6.22862899036e-07;


    a_width[0][0] = 26.2755947113;
    a_width[0][1] = 26.5309486389;
    a_width[0][2] = 26.4451007843;
    a_width[0][3] = 26.2696399689;
    a_width[0][4] = 26.1615753174;
    a_width[1][0] = 56.4464111328;
    a_width[1][1] = 58.1373405457;
    a_width[1][2] = 59.370475769;
    a_width[1][3] = 58.9584884644;
    a_width[1][4] = 57.898525238;
    a_width[2][0] = 98.7043075562;
    a_width[2][1] = 104.071395874;
    a_width[2][2] = 107.650276184;
    a_width[2][3] = 107.190132141;
    a_width[2][4] = 108.819702148;
    a_width[3][0] = 159.671554565;
    a_width[3][1] = 160.143630981;
    a_width[3][2] = 160.779037476;
    a_width[3][3] = 156.512420654;
    a_width[3][4] = 153.055923462;
    a_width[4][0] = 226.030212402;
    a_width[4][1] = 230.492599487;
    a_width[4][2] = 236.527374268;
    a_width[4][3] = 225.275009155;
    a_width[4][4] = 215.201370239;
    a_width[5][0] = 317.596405029;
    a_width[5][1] = 326.611419678;
    a_width[5][2] = 344.662902832;
    a_width[5][3] = 349.542297363;
    a_width[5][4] = 336.28414917;
    a_width[6][0] = 361.693115234;
    a_width[6][1] = 353.281433105;
    a_width[6][2] = 403.341369629;
    a_width[6][3] = 466.035919189;
    a_width[6][4] = 470.777160645;

  }
  else {

    a_delta[0][0] = -0.00602485053241;
    a_delta[0][1] = 0.0497528798878;
    a_delta[0][2] = 0.138858646154;
    a_delta[0][3] = 0.138038516045;
    a_delta[0][4] = 0.382582753897;
    a_delta[1][0] = 0.0145467305556;
    a_delta[1][1] = 0.116387352347;
    a_delta[1][2] = 0.48050069809;
    a_delta[1][3] = 0.990342915058;
    a_delta[1][4] = 1.39260792732;
    a_delta[2][0] = 0.0406219810247;
    a_delta[2][1] = 0.561077475548;
    a_delta[2][2] = 1.47835934162;
    a_delta[2][3] = 3.00171518326;
    a_delta[2][4] = 4.36777114868;
    a_delta[3][0] = 0.662394285202;
    a_delta[3][1] = 1.9265152216;
    a_delta[3][2] = 3.76084780693;
    a_delta[3][3] = 6.87533044815;
    a_delta[3][4] = 10.7636861801;
    a_delta[4][0] = 0.698600530624;
    a_delta[4][1] = 4.1973285675;
    a_delta[4][2] = 9.02881240845;
    a_delta[4][3] = 13.0142173767;
    a_delta[4][4] = 16.7301330566;
    a_delta[5][0] = 1.05567407608;
    a_delta[5][1] = 1.65154206753;
    a_delta[5][2] = 0.957142710686;
    a_delta[5][3] = 1.80052042007;
    a_delta[5][4] = 1.28974568844;
    a_delta[6][0] = 3.60177016258;
    a_delta[6][1] = 12.6570005417;
    a_delta[6][2] = 13.6016483307;
    a_delta[6][3] = 16.3001346588;
    a_delta[6][4] = 16.1479625702;


    a_alpha[0][0] = 7.99999427795;
    a_alpha[0][1] = 0.0100004309788;
    a_alpha[0][2] = 3.2126493454;
    a_alpha[0][3] = 2.80285358429;
    a_alpha[0][4] = 4.83957672119;
    a_alpha[1][0] = 0.0100059742108;
    a_alpha[1][1] = 0.0100066978484;
    a_alpha[1][2] = 0.744444966316;
    a_alpha[1][3] = 0.327008664608;
    a_alpha[1][4] = 0.0100099137053;
    a_alpha[2][0] = 1.97709751129;
    a_alpha[2][1] = 0.951394557953;
    a_alpha[2][2] = 0.730503380299;
    a_alpha[2][3] = 0.668072938919;
    a_alpha[2][4] = 0.0101733189076;
    a_alpha[3][0] = 1.12670946121;
    a_alpha[3][1] = 1.19023668766;
    a_alpha[3][2] = 0.592940509319;
    a_alpha[3][3] = 0.360803157091;
    a_alpha[3][4] = 0.0100000798702;
    a_alpha[4][0] = 0.0125287985429;
    a_alpha[4][1] = 1.37840366364;
    a_alpha[4][2] = 1.3054047823;
    a_alpha[4][3] = 1.01578176022;
    a_alpha[4][4] = 0.909670114517;
    a_alpha[5][0] = 1.48734974861;
    a_alpha[5][1] = 1.68027830124;
    a_alpha[5][2] = 1.29963386059;
    a_alpha[5][3] = 1.13843405247;
    a_alpha[5][4] = 0.938163638115;
    a_alpha[6][0] = 2.49235057831;
    a_alpha[6][1] = 2.76869463921;
    a_alpha[6][2] = 0.648731172085;
    a_alpha[6][3] = 0.47552099824;
    a_alpha[6][4] = 0.0100004505366;


    a_beta[0][0] = 5.0;
    a_beta[0][1] = 5.0;
    a_beta[0][2] = 5.0;
    a_beta[0][3] = 5.0;
    a_beta[0][4] = 5.0;
    a_beta[1][0] = 5.0;
    a_beta[1][1] = 5.0;
    a_beta[1][2] = 5.0;
    a_beta[1][3] = 5.0;
    a_beta[1][4] = 5.0;
    a_beta[2][0] = 5.0;
    a_beta[2][1] = 5.0;
    a_beta[2][2] = 5.0;
    a_beta[2][3] = 5.0;
    a_beta[2][4] = 5.0;
    a_beta[3][0] = 5.0;
    a_beta[3][1] = 5.0;
    a_beta[3][2] = 5.0;
    a_beta[3][3] = 5.0;
    a_beta[3][4] = 5.0;
    a_beta[4][0] = 5.0;
    a_beta[4][1] = 5.0;
    a_beta[4][2] = 5.0;
    a_beta[4][3] = 5.0;
    a_beta[4][4] = 5.0;
    a_beta[5][0] = 5.0;
    a_beta[5][1] = 5.0;
    a_beta[5][2] = 5.0;
    a_beta[5][3] = 5.0;
    a_beta[5][4] = 5.0;
    a_beta[6][0] = 5.0;
    a_beta[6][1] = 5.0;
    a_beta[6][2] = 5.0;
    a_beta[6][3] = 5.0;
    a_beta[6][4] = 5.0;


    a_r[0][0] = 1.31995675474e-05;
    a_r[0][1] = 1.00640645542e-05;
    a_r[0][2] = 5.54175812795e-06;
    a_r[0][3] = 2.26255224334e-06;
    a_r[0][4] = 8.00000009349e-07;
    a_r[1][0] = 8.69973064255e-06;
    a_r[1][1] = 6.27477720627e-06;
    a_r[1][2] = 3.45648322764e-06;
    a_r[1][3] = 2.40017675424e-06;
    a_r[1][4] = 1.46912748278e-06;
    a_r[2][0] = 6.0946326812e-06;
    a_r[2][1] = 4.60652609036e-06;
    a_r[2][2] = 2.65052358372e-06;
    a_r[2][3] = 2.03129661713e-06;
    a_r[2][4] = 1.45131468798e-06;
    a_r[3][0] = 4.87701117891e-06;
    a_r[3][1] = 3.61279262506e-06;
    a_r[3][2] = 1.97546864911e-06;
    a_r[3][3] = 1.51487427047e-06;
    a_r[3][4] = 1.17223669349e-06;
    a_r[4][0] = 3.69253280041e-06;
    a_r[4][1] = 3.07117966258e-06;
    a_r[4][2] = 2.04037019103e-06;
    a_r[4][3] = 1.55809368607e-06;
    a_r[4][4] = 1.19911203456e-06;
    a_r[5][0] = 3.46195929524e-06;
    a_r[5][1] = 2.55810459748e-06;
    a_r[5][2] = 1.41951807109e-06;
    a_r[5][3] = 1.08365168217e-06;
    a_r[5][4] = 7.98482290065e-07;
    a_r[6][0] = 3.63989897778e-06;
    a_r[6][1] = 2.6686595902e-06;
    a_r[6][2] = 1.2225913224e-06;
    a_r[6][3] = 8.71984127571e-07;
    a_r[6][4] = 6.06827427418e-07;


    a_width[0][0] = 26.0809555054;
    a_width[0][1] = 26.3398857117;
    a_width[0][2] = 26.1370944977;
    a_width[0][3] = 26.1412086487;
    a_width[0][4] = 26.0180282593;
    a_width[1][0] = 59.4263000488;
    a_width[1][1] = 59.2995262146;
    a_width[1][2] = 59.1552963257;
    a_width[1][3] = 58.8251228333;
    a_width[1][4] = 58.127784729;
    a_width[2][0] = 100.632019043;
    a_width[2][1] = 103.334915161;
    a_width[2][2] = 105.333137512;
    a_width[2][3] = 104.352851868;
    a_width[2][4] = 103.870521545;
    a_width[3][0] = 162.842437744;
    a_width[3][1] = 163.783569336;
    a_width[3][2] = 167.428451538;
    a_width[3][3] = 165.722122192;
    a_width[3][4] = 162.353118896;
    a_width[4][0] = 231.11781311;
    a_width[4][1] = 226.819869995;
    a_width[4][2] = 222.608901978;
    a_width[4][3] = 216.191223145;
    a_width[4][4] = 209.027923584;
    a_width[5][0] = 318.513092041;
    a_width[5][1] = 318.359008789;
    a_width[5][2] = 325.207244873;
    a_width[5][3] = 320.036865234;
    a_width[5][4] = 317.293792725;
    a_width[6][0] = 385.060577393;
    a_width[6][1] = 373.778320312;
    a_width[6][2] = 435.293792725;
    a_width[6][3] = 450.927459717;
    a_width[6][4] = 462.96295166;
    
  }


}


Double_t RooCPSHighMassVBF::interpolateMatrix(const Double_t matrix[][5] , const Double_t& x, const Double_t& y) const
{

  // Dummy copy of ROOT TH2::Interpolate needed because there the function is not constant

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring static const nMass, nCPrime" << std::endl; 
  static const Int_t nMass   = 7;
  static const Int_t nCPrime = 5;

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring arrays masses, cprimes" << std::endl; 
  Double_t masses[nMass]    = {400.,500.,600.,700.,800.,900.,1000.};
  Double_t cprimes[nCPrime] = {0.2,0.4,0.6,0.8,1.0};

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring stepX, stepY, lowX, lowY" << std::endl; 
  Double_t stepX = masses[1]-masses[0];
  Double_t stepY = cprimes[1]-cprimes[0];
  Double_t lowX = masses[0]-stepX/2;  //Double_t highX = masses[nMass-1]+stepX/2;
  Double_t lowY = cprimes[0]-stepY/2; //Double_t highY = cprimes[nMass-1]+stepY/2;
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: stepX=" << stepX << ", stepY=" << stepY << ", lowX=" << lowX << ", lowY=" << lowY << std::endl; 
  Int_t nX = nMass;
  Int_t nY = nCPrime;

  // Get the indexes to identify matrix element
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring i_x, i_y" << std::endl; 
  Int_t i_x = (Int_t)((x-lowX)/stepX);
  Int_t i_y = (Int_t)((y-lowY)/stepY);
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: i_x=" << i_x << ", i_y=" << i_y << std::endl; 

  // Check if we are in the domain for the interpolation
  if(i_x<0 || i_x>(nX-1) || i_y<0 || i_y>(nY-1)){
    std::cout << "[RooCPSHighMassVBF::interpolateMatrix] Cannot interpolate outside array domain, returning 0." << std::endl;
    return 0;
  }

  // Find out the quadrant in the bin
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring quadrant, dx, dy" << std::endl; 
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  dx = masses[i_x]+stepX/2-x;
  dy = cprimes[i_y]+stepY/2-y;
  if (dx<=stepX/2 && dy<=stepY/2)
    quadrant = 1; // upper right
  if (dx>stepX/2 && dy<=stepY/2)
    quadrant = 2; // upper left
  if (dx>stepX/2 && dy>stepY/2)
    quadrant = 3; // lower left
  if (dx<=stepX/2 && dy>stepY/2)
    quadrant = 4; // lower right
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: quadrant=" << quadrant << ", dx=" << dx << ", dy=" << dy << std::endl; 

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring x_1, y1, x2, y2" << std::endl; 
  switch(quadrant) {
  case 1:
    x1 = masses[i_x];
    y1 = cprimes[i_y];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 2:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    y1 = cprimes[i_y];
    x2 = masses[i_x];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 3:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    x2 = masses[i_x];
    y2 = cprimes[i_y];
    break;
  case 4:
    x1 = masses[i_x];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    y2 = cprimes[i_y];
    break;
  }

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: x1=" << x1 << ", y1=" << y1 << ", x2=" << x2 << ", y2=" << y2 << std::endl; 

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring i_x1, i_x2, i_y1, i_y2" << std::endl; 
  Int_t i_x1 = (Int_t)((x1-lowX)/stepX);
  if(i_x1<0) i_x1=0;
  Int_t i_x2 = (Int_t)((x2-lowX)/stepX);
  if(i_x2>nMass-1) i_x2=nMass-1;
  Int_t i_y1 = (Int_t)((y1-lowY)/stepY);
  if(i_y1<0) i_y1=0;
  Int_t i_y2 = (Int_t)((y2-lowY)/stepY);
  if(i_y2>nCPrime-1) i_y2=nCPrime-1;
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: i_x1=" << i_x1 << ", i_x2=" << i_x2 << ", i_y1=" << i_y1 << ", i_y2=" << i_y2 << std::endl; 

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring q11, q12, q21, q22" << std::endl; 
  Double_t q11 = matrix[i_x1][i_y1];
  Double_t q12 = matrix[i_x1][i_y2];
  Double_t q21 = matrix[i_x2][i_y1];
  Double_t q22 = matrix[i_x2][i_y2];
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: q11=" << q11 << ", q12=" << q12 << ", q21=" << q21 << ", q22=" << q22 << std::endl; 

  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: declaring d, f" << std::endl;   
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  // std::cout << "[RooCPSHighMassVBF::interpolateMatrix] DEBUG: d=" << d << ", f=" << f << std::endl;
  
  return f;
}


Double_t RooCPSHighMassVBF::getWidth(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_width,mass,cprime);}

Double_t RooCPSHighMassVBF::getDelta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_delta,mass,cprime);}

Double_t RooCPSHighMassVBF::getAlpha(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_alpha,mass,cprime);}

Double_t RooCPSHighMassVBF::getR(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_r,mass,cprime);}

Double_t RooCPSHighMassVBF::getBeta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_beta,mass,cprime);}




Double_t RooCPSHighMassVBF::Spline(Double_t xx) const{
  

  if (is8TeV){

    const int fNp = 17, fKstep = 0;
    const double fDelta = -1, fXmin = 225, fXmax = 1925;
    const double fX[17] = { 225, 425, 525, 625, 725,
			    825, 925, 1025, 1125, 1225,
			    1325, 1425, 1525, 1625, 1725,
			    1825, 1925 };
    const double fY[17] = { 0.0306502, 0.326692, 0.487167, 0.642114, 0.696507,
			    0.842851, 1.09799, 1.29813, 1.40643, 1.48621,
			    1.35416, 1.29446, 1.31814, 1.14274, 1.01259,
			    0.804478, 0.674614 };
    const double fB[17] = { 0.0016145, 0.0014908, 0.00175492, 0.00095217, 0.000716595,
			    0.00220355, 0.00251384, 0.0013996, 0.00114068, -0.000319938,
			    -0.00142889, 0.000282881, -0.000783237, -0.00170154, -0.00157723,
			    -0.0021374, -1.23446e-05 };
    const double fC[17] = { -1.39581e-06, 7.77332e-07, 1.86391e-06, -9.89143e-06, 7.53568e-06,
			    7.33385e-06, -4.23088e-06, -6.91158e-06, 4.32242e-06, -1.89286e-05,
			    7.83905e-06, 9.2787e-06, -1.99399e-05, 1.07569e-05, -9.51384e-06,
			    3.91223e-06, 100 };
    const double fD[17] = { 3.62191e-09, 3.62191e-09, -3.91845e-08, 5.80904e-08, -6.72782e-10,
			    -3.85491e-08, -8.93568e-09, 3.74467e-08, -7.75034e-08, 8.92255e-08,
			    4.79882e-09, -9.73953e-08, 1.02323e-07, -6.7569e-08, 4.47535e-08,
			    4.47535e-08, 46.4102 };
    int klow=0;
    if(xx<=fXmin) klow=0;
    else if(xx>=fXmax) klow=fNp-1;
    else {
      if(fKstep) {
	// Equidistant knots, use histogramming
	klow = int((xx-fXmin)/fDelta);
	if (klow < fNp-1) klow = fNp-1;
      } else {
	int khig=fNp-1, khalf;
	// Non equidistant knots, binary search
	while(khig-klow>1)
	  if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
	  else khig=khalf;
      }
    }
    // Evaluate now
    double dx=xx-fX[klow];

    return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));

  } 
  else{

    const int fNp = 14, fKstep = 0;
    const double fDelta = -1, fXmin = 250, fXmax = 1950;
    const double fX[14] = { 250, 450, 550, 650, 750,
			    850, 950, 1050, 1150, 1250,
			    1450, 1650, 1850, 1950 };
    const double fY[14] = { 0.0863593, 0.465588, 0.619189, 0.736075, 0.833095,
			    0.913433, 1.17228, 1.33452, 1.2992, 1.25595,
			    1.07763, 0.810053, 0.552352, 0.490991 };
    const double fB[14] = { 0.00194059, 0.00172127, 0.00131814, 0.00112081, 0.000615821,
			    0.00173663, 0.00261322, 0.000443189, -0.000578349, -0.000487042,
			    -0.00119102, -0.00143731, -0.000938952, -0.000238068 };
    const double fC[14] = { 4.2989e-07, -1.52654e-06, -2.50475e-06, 5.31463e-07, -5.58133e-06,
			    1.67894e-05, -8.02358e-06, -1.36767e-05, 3.46131e-06, -2.54824e-06,
			    -9.71639e-07, -2.59803e-07, 2.75157e-06, 100 };
    const double fD[14] = { -3.26071e-09, -3.26071e-09, 1.01207e-08, -2.0376e-08, 7.45692e-08,
			    -8.271e-08, -1.88437e-08, 5.71267e-08, -2.00318e-08, 2.62766e-09,
			    1.18639e-09, 5.01896e-09, 5.01896e-09, 95.3332 };
    int klow=0;
    if(xx<=fXmin) klow=0;
    else if(xx>=fXmax) klow=fNp-1;
    else {
      if(fKstep) {
	// Equidistant knots, use histogramming
	klow = int((xx-fXmin)/fDelta);
	if (klow < fNp-1) klow = fNp-1;
      } else {
	int khig=fNp-1, khalf;
	// Non equidistant knots, binary search
	while(khig-klow>1)
	  if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
	  else khig=khalf;
      }
    }
    // Evaluate now
    double dx=xx-fX[klow];
    return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));

  }


}



/*************RooCPSHighMassVBFNoInterf***********/
ClassImp(RooCPSHighMassVBFNoInterf) 

 RooCPSHighMassVBFNoInterf::RooCPSHighMassVBFNoInterf(){}

 RooCPSHighMassVBFNoInterf::RooCPSHighMassVBFNoInterf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _KPrime,
			RooAbsReal& _BRnew,
			Bool_t is8TeV_b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   KPrime("KPrime","KPrime",this,_KPrime),
   BRnew("BRnew","BRnew",this,_BRnew),
   is8TeV(is8TeV_b)
 { 
   initMatrices();
 } 


 RooCPSHighMassVBFNoInterf::RooCPSHighMassVBFNoInterf(const RooCPSHighMassVBFNoInterf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   KPrime("KPrime",this,other.KPrime),
   BRnew("BRnew",this,other.BRnew),
   is8TeV(other.is8TeV)
 { 

   for (Int_t i=0; i<7; ++i){
     for (Int_t j=0; j<5; ++j){
       a_width[i][j] = other.a_width[i][j];
       a_delta[i][j] = other.a_delta[i][j];
     }
   }

 } 



 Double_t RooCPSHighMassVBFNoInterf::evaluate() const 
 {    

   Double_t effKPrime = KPrime/(1-BRnew); 
   Double_t width = getWidth(mH,effKPrime);
   Double_t delta = getDelta(mH,effKPrime);

   Double_t bwHM = x / ( TMath::Power( TMath::Power(x,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(x,2)*TMath::Power(effKPrime*width,2) ); 

   Double_t splineFactor;
   if (x<1850) splineFactor = Spline(x);
   else splineFactor = Spline(1850);

   Double_t signal = bwHM*splineFactor;

   Double_t fValue = signal;
   if (fValue > 0) return fValue;
   else return 0.;
   
   
 } 


void RooCPSHighMassVBFNoInterf::initMatrices() {

  if (is8TeV){

    a_delta[0][0] = -0.0185139477253;
    a_delta[0][1] = 0.0201042648405;
    a_delta[0][2] = 0.156271114945;
    a_delta[0][3] = 0.310714632273;
    a_delta[0][4] = 0.622391521931;
    a_delta[1][0] = 0.0428073592484;
    a_delta[1][1] = 0.15564674139;
    a_delta[1][2] = 0.44218480587;
    a_delta[1][3] = 0.845812380314;
    a_delta[1][4] = 1.2718859911;
    a_delta[2][0] = 0.1351634413;
    a_delta[2][1] = 0.623583912849;
    a_delta[2][2] = 1.33247685432;
    a_delta[2][3] = 2.8496837616;
    a_delta[2][4] = 4.2275352478;
    a_delta[3][0] = 0.407202005386;
    a_delta[3][1] = 1.78870689869;
    a_delta[3][2] = 4.4582362175;
    a_delta[3][3] = 8.72150039673;
    a_delta[3][4] = 12.9632072449;
    a_delta[4][0] = 1.10362958908;
    a_delta[4][1] = 1.65279102325;
    a_delta[4][2] = 2.63262009621;
    a_delta[4][3] = 5.56618690491;
    a_delta[4][4] = 8.73703479767;
    a_delta[5][0] = 2.47303795815;
    a_delta[5][1] = 3.7929110527;
    a_delta[5][2] = 0.907296419144;
    a_delta[5][3] = -3.98578429222;
    a_delta[5][4] = -1.68104398251;
    a_delta[6][0] = 3.5346019268;
    a_delta[6][1] = 14.8740081787;
    a_delta[6][2] = 22.4615650177;
    a_delta[6][3] = 26.1576747894;
    a_delta[6][4] = 25.2779407501;


    a_width[0][0] = 26.2755947113;
    a_width[0][1] = 26.5309486389;
    a_width[0][2] = 26.4451007843;
    a_width[0][3] = 26.2696399689;
    a_width[0][4] = 26.1615753174;
    a_width[1][0] = 56.4464111328;
    a_width[1][1] = 58.1373405457;
    a_width[1][2] = 59.370475769;
    a_width[1][3] = 58.9584884644;
    a_width[1][4] = 57.898525238;
    a_width[2][0] = 98.7043075562;
    a_width[2][1] = 104.071395874;
    a_width[2][2] = 107.650276184;
    a_width[2][3] = 107.190132141;
    a_width[2][4] = 108.819702148;
    a_width[3][0] = 159.671554565;
    a_width[3][1] = 160.143630981;
    a_width[3][2] = 160.779037476;
    a_width[3][3] = 156.512420654;
    a_width[3][4] = 153.055923462;
    a_width[4][0] = 226.030212402;
    a_width[4][1] = 230.492599487;
    a_width[4][2] = 236.527374268;
    a_width[4][3] = 225.275009155;
    a_width[4][4] = 215.201370239;
    a_width[5][0] = 317.596405029;
    a_width[5][1] = 326.611419678;
    a_width[5][2] = 344.662902832;
    a_width[5][3] = 349.542297363;
    a_width[5][4] = 336.28414917;
    a_width[6][0] = 361.693115234;
    a_width[6][1] = 353.281433105;
    a_width[6][2] = 403.341369629;
    a_width[6][3] = 466.035919189;
    a_width[6][4] = 470.777160645;

  }
  else {

    a_delta[0][0] = -0.00602485053241;
    a_delta[0][1] = 0.0497528798878;
    a_delta[0][2] = 0.138858646154;
    a_delta[0][3] = 0.138038516045;
    a_delta[0][4] = 0.382582753897;
    a_delta[1][0] = 0.0145467305556;
    a_delta[1][1] = 0.116387352347;
    a_delta[1][2] = 0.48050069809;
    a_delta[1][3] = 0.990342915058;
    a_delta[1][4] = 1.39260792732;
    a_delta[2][0] = 0.0406219810247;
    a_delta[2][1] = 0.561077475548;
    a_delta[2][2] = 1.47835934162;
    a_delta[2][3] = 3.00171518326;
    a_delta[2][4] = 4.36777114868;
    a_delta[3][0] = 0.662394285202;
    a_delta[3][1] = 1.9265152216;
    a_delta[3][2] = 3.76084780693;
    a_delta[3][3] = 6.87533044815;
    a_delta[3][4] = 10.7636861801;
    a_delta[4][0] = 0.698600530624;
    a_delta[4][1] = 4.1973285675;
    a_delta[4][2] = 9.02881240845;
    a_delta[4][3] = 13.0142173767;
    a_delta[4][4] = 16.7301330566;
    a_delta[5][0] = 1.05567407608;
    a_delta[5][1] = 1.65154206753;
    a_delta[5][2] = 0.957142710686;
    a_delta[5][3] = 1.80052042007;
    a_delta[5][4] = 1.28974568844;
    a_delta[6][0] = 3.60177016258;
    a_delta[6][1] = 12.6570005417;
    a_delta[6][2] = 13.6016483307;
    a_delta[6][3] = 16.3001346588;
    a_delta[6][4] = 16.1479625702;


    a_width[0][0] = 26.0809555054;
    a_width[0][1] = 26.3398857117;
    a_width[0][2] = 26.1370944977;
    a_width[0][3] = 26.1412086487;
    a_width[0][4] = 26.0180282593;
    a_width[1][0] = 59.4263000488;
    a_width[1][1] = 59.2995262146;
    a_width[1][2] = 59.1552963257;
    a_width[1][3] = 58.8251228333;
    a_width[1][4] = 58.127784729;
    a_width[2][0] = 100.632019043;
    a_width[2][1] = 103.334915161;
    a_width[2][2] = 105.333137512;
    a_width[2][3] = 104.352851868;
    a_width[2][4] = 103.870521545;
    a_width[3][0] = 162.842437744;
    a_width[3][1] = 163.783569336;
    a_width[3][2] = 167.428451538;
    a_width[3][3] = 165.722122192;
    a_width[3][4] = 162.353118896;
    a_width[4][0] = 231.11781311;
    a_width[4][1] = 226.819869995;
    a_width[4][2] = 222.608901978;
    a_width[4][3] = 216.191223145;
    a_width[4][4] = 209.027923584;
    a_width[5][0] = 318.513092041;
    a_width[5][1] = 318.359008789;
    a_width[5][2] = 325.207244873;
    a_width[5][3] = 320.036865234;
    a_width[5][4] = 317.293792725;
    a_width[6][0] = 385.060577393;
    a_width[6][1] = 373.778320312;
    a_width[6][2] = 435.293792725;
    a_width[6][3] = 450.927459717;
    a_width[6][4] = 462.96295166;
    
  }


}



Double_t RooCPSHighMassVBFNoInterf::interpolateMatrix(const Double_t matrix[][5] , const Double_t& x, const Double_t& y) const
{

  // Dummy copy of ROOT TH2::Interpolate needed because there the function is not constant

  static const Int_t nMass   = 7;
  static const Int_t nCPrime = 5;

  Double_t masses[nMass]    = {400.,500.,600.,700.,800.,900.,1000.};
  Double_t cprimes[nCPrime] = {0.2,0.4,0.6,0.8,1.0};

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  Double_t stepX = masses[1]-masses[0];
  Double_t stepY = cprimes[1]-cprimes[0];
  Double_t lowX = masses[0]-stepX/2;  //Double_t highX = masses[nMass-1]+stepX/2;
  Double_t lowY = cprimes[0]-stepY/2; //Double_t highY = cprimes[nMass-1]+stepY/2;
  Int_t nX = nMass;
  Int_t nY = nCPrime;

  // Get the indexes to identify matrix element
  Int_t i_x = (Int_t)((x-lowX)/stepX);
  Int_t i_y = (Int_t)((y-lowY)/stepY);

  // Check if we are in the domain for the interpolation
  if(i_x<0 || i_x>(nX-1) || i_y<0 || i_y>(nY-1)){
    // std::cout << "[RooCPSHighMassVBFNoInterf::interpolateMatrix] Cannot interpolate outside array domain, returning 0." << std::endl;
    return 0;
  }

  // Find out the quadrant in the bin
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  dx = masses[i_x]+stepX/2-x;
  dy = cprimes[i_y]+stepY/2-y;
  if (dx<=stepX/2 && dy<=stepY/2)
    quadrant = 1; // upper right
  if (dx>stepX/2 && dy<=stepY/2)
    quadrant = 2; // upper left
  if (dx>stepX/2 && dy>stepY/2)
    quadrant = 3; // lower left
  if (dx<=stepX/2 && dy>stepY/2)
    quadrant = 4; // lower right

  switch(quadrant) {
  case 1:
    x1 = masses[i_x];
    y1 = cprimes[i_y];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 2:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    y1 = cprimes[i_y];
    x2 = masses[i_x];
    (i_y+1) > nCPrime ? y2 = cprimes[i_y] + stepY : y2 = cprimes[i_y+1];
    break;
  case 3:
    (i_x-1) < 0 ? x1 = masses[i_x] - stepX : x1 = masses[i_x-1];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    x2 = masses[i_x];
    y2 = cprimes[i_y];
    break;
  case 4:
    x1 = masses[i_x];
    (i_y-1) < 0 ? y1 = cprimes[i_y] - stepY : y1 = cprimes[i_y-1];
    (i_x+1) > nMass-1 ? x2 = masses[i_x] + stepX : x2 = masses[i_x+1];
    y2 = cprimes[i_y];
    break;
  }

  Int_t i_x1 = (Int_t)((x1-lowX)/stepX);
  if(i_x1<0) i_x1=0;
  Int_t i_x2 = (Int_t)((x2-lowX)/stepX);
  if(i_x2>nMass-1) i_x2=nMass-1;
  Int_t i_y1 = (Int_t)((y1-lowY)/stepY);
  if(i_y1<0) i_y1=0;
  Int_t i_y2 = (Int_t)((y2-lowY)/stepY);
  if(i_y2>nCPrime-1) i_y2=nCPrime-1;
  Double_t q11 = matrix[i_x1][i_y1];
  Double_t q12 = matrix[i_x1][i_y2];
  Double_t q21 = matrix[i_x2][i_y1];
  Double_t q22 = matrix[i_x2][i_y2];
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  
  // std::cout << "[RooCPSHighMassVBFNoInterf::interpolateMatrix] DEBUG: " << x << " - " << y << "->" << q11 << " / " << q12 << " / " << q21 << " / " << q22 << " / " << f << std::endl;
  
  return f;
}


Double_t RooCPSHighMassVBFNoInterf::getWidth(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_width,mass,cprime);}

Double_t RooCPSHighMassVBFNoInterf::getDelta(const Double_t& mass, const Double_t& cprime) const {return interpolateMatrix(a_delta,mass,cprime);}



Double_t RooCPSHighMassVBFNoInterf::Spline(Double_t xx) const{
  

  if (is8TeV){

    const int fNp = 17, fKstep = 0;
    const double fDelta = -1, fXmin = 225, fXmax = 1925;
    const double fX[17] = { 225, 425, 525, 625, 725,
			    825, 925, 1025, 1125, 1225,
			    1325, 1425, 1525, 1625, 1725,
			    1825, 1925 };
    const double fY[17] = { 0.0306502, 0.326692, 0.487167, 0.642114, 0.696507,
			    0.842851, 1.09799, 1.29813, 1.40643, 1.48621,
			    1.35416, 1.29446, 1.31814, 1.14274, 1.01259,
			    0.804478, 0.674614 };
    const double fB[17] = { 0.0016145, 0.0014908, 0.00175492, 0.00095217, 0.000716595,
			    0.00220355, 0.00251384, 0.0013996, 0.00114068, -0.000319938,
			    -0.00142889, 0.000282881, -0.000783237, -0.00170154, -0.00157723,
			    -0.0021374, -1.23446e-05 };
    const double fC[17] = { -1.39581e-06, 7.77332e-07, 1.86391e-06, -9.89143e-06, 7.53568e-06,
			    7.33385e-06, -4.23088e-06, -6.91158e-06, 4.32242e-06, -1.89286e-05,
			    7.83905e-06, 9.2787e-06, -1.99399e-05, 1.07569e-05, -9.51384e-06,
			    3.91223e-06, 100 };
    const double fD[17] = { 3.62191e-09, 3.62191e-09, -3.91845e-08, 5.80904e-08, -6.72782e-10,
			    -3.85491e-08, -8.93568e-09, 3.74467e-08, -7.75034e-08, 8.92255e-08,
			    4.79882e-09, -9.73953e-08, 1.02323e-07, -6.7569e-08, 4.47535e-08,
			    4.47535e-08, 46.4102 };
    int klow=0;
    if(xx<=fXmin) klow=0;
    else if(xx>=fXmax) klow=fNp-1;
    else {
      if(fKstep) {
	// Equidistant knots, use histogramming
	klow = int((xx-fXmin)/fDelta);
	if (klow < fNp-1) klow = fNp-1;
      } else {
	int khig=fNp-1, khalf;
	// Non equidistant knots, binary search
	while(khig-klow>1)
	  if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
	  else khig=khalf;
      }
    }
    // Evaluate now
    double dx=xx-fX[klow];

    return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));

  } 
  else{

    const int fNp = 14, fKstep = 0;
    const double fDelta = -1, fXmin = 250, fXmax = 1950;
    const double fX[14] = { 250, 450, 550, 650, 750,
			    850, 950, 1050, 1150, 1250,
			    1450, 1650, 1850, 1950 };
    const double fY[14] = { 0.0863593, 0.465588, 0.619189, 0.736075, 0.833095,
			    0.913433, 1.17228, 1.33452, 1.2992, 1.25595,
			    1.07763, 0.810053, 0.552352, 0.490991 };
    const double fB[14] = { 0.00194059, 0.00172127, 0.00131814, 0.00112081, 0.000615821,
			    0.00173663, 0.00261322, 0.000443189, -0.000578349, -0.000487042,
			    -0.00119102, -0.00143731, -0.000938952, -0.000238068 };
    const double fC[14] = { 4.2989e-07, -1.52654e-06, -2.50475e-06, 5.31463e-07, -5.58133e-06,
			    1.67894e-05, -8.02358e-06, -1.36767e-05, 3.46131e-06, -2.54824e-06,
			    -9.71639e-07, -2.59803e-07, 2.75157e-06, 100 };
    const double fD[14] = { -3.26071e-09, -3.26071e-09, 1.01207e-08, -2.0376e-08, 7.45692e-08,
			    -8.271e-08, -1.88437e-08, 5.71267e-08, -2.00318e-08, 2.62766e-09,
			    1.18639e-09, 5.01896e-09, 5.01896e-09, 95.3332 };
    int klow=0;
    if(xx<=fXmin) klow=0;
    else if(xx>=fXmax) klow=fNp-1;
    else {
      if(fKstep) {
	// Equidistant knots, use histogramming
	klow = int((xx-fXmin)/fDelta);
	if (klow < fNp-1) klow = fNp-1;
      } else {
	int khig=fNp-1, khalf;
	// Non equidistant knots, binary search
	while(khig-klow>1)
	  if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
	  else khig=khalf;
      }
    }
    // Evaluate now
    double dx=xx-fX[klow];
    return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));

  }


}



/*************RooSigPlusInt***********/
ClassImp(RooSigPlusInt) 

 RooSigPlusInt::RooSigPlusInt(){}

 RooSigPlusInt::RooSigPlusInt(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mH,
			RooAbsReal& _delta,
                        RooAbsReal& _width,
                        RooAbsReal& _k,
                        RooAbsReal& _CSquared,
                        RooAbsReal& _BRnew,
                        RooAbsReal& _alpha,
                        RooAbsReal& _beta,
                        RooAbsReal& _r) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mH("mH","mH",this,_mH),
   delta("delta","delta",this,_delta),
   width("width","width",this,_width),
   k("k","k",this,_k),
   CSquared("CSquared","CSquared",this,_CSquared),
   BRnew("BRnew","BRnew",this,_BRnew),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   r("r","r",this,_r)
 { 
 } 


 RooSigPlusInt::RooSigPlusInt(const RooSigPlusInt& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mH("mH",this,other.mH),
   delta("delta",this,other.delta),
   width("width",this,other.width),
   k("k",this,other.k),
   CSquared("CSquared",this,other.CSquared),
   BRnew("BRnew",this,other.BRnew),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   r("r",this,other.r)
 { 
 } 



 Double_t RooSigPlusInt::evaluate() const 
 { 
   // Cast the "x" proxy to double only once to avoid overhead and ambiguous
   // overloads the TComplex operators.
   const double xVal = x;

   Double_t totWidthSF = CSquared/(1-BRnew);

   Double_t mH_eff = mH*TMath::Sqrt(1-k*(totWidthSF*width/mH)*(totWidthSF*width/mH));
 
   Double_t splineFactor = Spline(x);
   Double_t bwHM = xVal / ( TMath::Power( TMath::Power(xVal,2) - TMath::Power(mH+delta,2) , 2 ) + TMath::Power(xVal,2)*TMath::Power(totWidthSF*width,2) ); 

   Double_t signal = splineFactor*bwHM;
      
   TComplex MSquared(mH_eff*mH_eff,-totWidthSF*width*mH_eff);
   TComplex M = MSquared.Sqrt(MSquared);
   TComplex Exp1 = MSquared.Exp((TComplex)alpha);
   TComplex Exp2 = MSquared.Exp(-(TComplex)alpha);
   TComplex Exp3 = MSquared.Exp(-xVal*(TComplex)alpha/M);

   double interference = -r*sqrt(1/1-BRnew)*(1-TMath::Exp(-beta*(xVal-150.)/mH_eff))*((Exp1/(xVal-M)-Exp2/(xVal+M))*Exp3).Re();

   return signal + interference;   
   
 } 


Double_t RooSigPlusInt::Spline(Double_t xx) const{
   const int fNp = 20, fKstep = 0;
   const double fDelta = -1, fXmin = 200, fXmax = 2150;
   const double fX[20] = { 200, 250, 300, 350, 400,
                        450, 500, 550, 600, 650,
                        750, 850, 950, 1050, 1150,
                        1350, 1550, 1750, 1950, 2150 };
   const double fY[20] = { 0.139988, 0.17804, 0.352505, 0.801982, 1.35134,
                        1.75174, 1.84163, 1.86637, 1.78052, 1.7297,
                        1.53036, 1.31241, 1.11751, 0.946671, 0.812524,
                        0.640316, 0.490843, 0.396263, 0.315344, 0.210836 };
   const double fB[20] = { 0.000266891, 0.00147333, 0.00628998, 0.0108012, 0.0105333,
                        0.00444666, 0.000737389, -0.000594766, -0.00180599, -0.000721876,
                        -0.00269494, -0.00186614, -0.00192752, -0.00150674, -0.00115716,
                        -0.000743197, -0.000649875, -0.000369596, -0.000461366, -0.000583113 };
   const double fC[20] = { 8.42921e-06, 2.29701e-05, 6.55918e-05, 1.67598e-05, -2.50909e-05,
                        -7.98341e-05, 4.46636e-06, -3.60991e-05, 1.78298e-05, -7.22121e-06,
                        3.25258e-06, -4.32435e-07, 1.81479e-06, 1.62255e-06, 2.04841e-06,
                        -2.33522e-07, 1.01254e-06, 5.76421e-08, -3.11836e-07, -3.01879e-07 };
   const double fD[20] = { 0, 2.90818e-07, -2.00191e-08, -3.34947e-07, -4.53847e-07,
                        2.18457e-07, 1.71166e-07, -2.84202e-08, 2.09784e-07, -2.33159e-07,
                        8.20808e-08, -2.72765e-08, 1.26473e-08, -2.76285e-09, 7.40335e-10,
                        -1.80875e-09, 1.31533e-09, -1.99675e-09, 4.97867e-11, 0 };
   const double fE[20] = { 0, 2.90818e-09, -6.01655e-09, 2.86727e-09, -4.05627e-09,
                        1.07793e-08, -1.12522e-08, 9.25635e-09, -6.87431e-09, 2.44487e-09,
                        -8.68676e-10, 3.2189e-10, -1.22271e-10, 4.52199e-11, -2.77039e-11,
                        2.13312e-11, -1.3521e-11, 5.24081e-12, -1.24467e-13, 0 };
   const double fF[20] = { 1.16327e-11, -3.56989e-11, 3.55353e-11, -2.76942e-11, 5.93423e-11,
                        -8.81261e-11, 8.20343e-11, -6.45226e-11, 3.72767e-11, -6.6271e-12,
                        2.38113e-12, -8.8832e-13, 3.34981e-13, -1.45848e-13, 4.90352e-14,
                        -3.48522e-14, 1.87618e-14, -5.36528e-15, 1.24467e-16, 0 };
   int klow=0;
   if(xx<=fXmin) klow=0;
   else if(xx>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((xx-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(xx>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=xx-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*(fD[klow]+dx*(fE[klow]+dx*fF[klow])))));
}

