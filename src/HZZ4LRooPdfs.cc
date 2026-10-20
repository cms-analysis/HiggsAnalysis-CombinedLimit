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
#include "RtypesCore.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>

using namespace TMath;

namespace {
  using HighMassMatrix = std::array<std::array<Double_t, 5>, 7>;

  constexpr std::array<std::array<double, 20>, 7> gghSplineCoefficients = {{
#include "../data/hzz4l/SplineGGH.inc"
  }};

  constexpr std::array<std::array<double, 17>, 5> vbf8TeVSplineCoefficients = {{
#include "../data/hzz4l/SplineVBF8TeV.inc"
  }};

  constexpr std::array<std::array<double, 14>, 5> vbf7TeVSplineCoefficients = {{
#include "../data/hzz4l/SplineVBF7TeV.inc"
  }};
}  // namespace

// Retain the persisted member array types while copying compiled C++17 data.
template <typename Destination>
static void copyHighMassMatrix(const HighMassMatrix& source, Destination& destination) {
  std::size_t row = 0;
  for (auto& values : destination) {
    const auto& sourceRow = source.at(row++);
    std::copy_n(sourceRow.begin(), sourceRow.size(), std::begin(values));
  }
}

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
  static constexpr std::array<std::array<HighMassMatrix, 5>, 2> matrices = {{
#include "../data/hzz4l/CPSHighMassGGHMatrices.inc"
  }};
  const auto& parameters = matrices.at(is8TeV ? 1 : 0);
  copyHighMassMatrix(parameters.at(0), a_width);
  copyHighMassMatrix(parameters.at(1), a_delta);
  copyHighMassMatrix(parameters.at(2), a_alpha);
  copyHighMassMatrix(parameters.at(3), a_r);
  copyHighMassMatrix(parameters.at(4), a_beta);
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
   const auto& coefficients = gghSplineCoefficients;
   const auto& fX = coefficients.at(0);
   const auto& fY = coefficients.at(1);
   const auto& fB = coefficients.at(2);
   const auto& fC = coefficients.at(3);
   const auto& fD = coefficients.at(4);
   const auto& fE = coefficients.at(5);
   const auto& fF = coefficients.at(6);
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
  static constexpr std::array<std::array<HighMassMatrix, 5>, 2> matrices = {{
#include "../data/hzz4l/BWHighMassGGHMatrices.inc"
  }};
  const auto& parameters = matrices.at(is8TeV ? 1 : 0);
  copyHighMassMatrix(parameters.at(0), a_width);
  copyHighMassMatrix(parameters.at(1), a_delta);
  copyHighMassMatrix(parameters.at(2), a_alpha);
  copyHighMassMatrix(parameters.at(3), a_r);
  copyHighMassMatrix(parameters.at(4), a_beta);
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
  static constexpr std::array<Double_t, 901> widths = {{
#include "../data/hzz4l/HiggsWidthHDECAY.inc"
  }};
  const int integerMass = static_cast<int>(mass);
  if (integerMass < 100 || integerMass > 1000) {
    std::cout << "[RooBWHighMassGGH::H_width]: mass is out of range" << std::endl;
    return 0.;
  }
  return widths.at(static_cast<std::size_t>(integerMass - 100));
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
  static constexpr std::array<std::array<HighMassMatrix, 2>, 2> matrices = {{
#include "../data/hzz4l/CPSHighMassGGHNoInterfMatrices.inc"
  }};
  const auto& parameters = matrices.at(is8TeV ? 1 : 0);
  copyHighMassMatrix(parameters.at(0), a_width);
  copyHighMassMatrix(parameters.at(1), a_delta);
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
   const auto& coefficients = gghSplineCoefficients;
   const auto& fX = coefficients.at(0);
   const auto& fY = coefficients.at(1);
   const auto& fB = coefficients.at(2);
   const auto& fC = coefficients.at(3);
   const auto& fD = coefficients.at(4);
   const auto& fE = coefficients.at(5);
   const auto& fF = coefficients.at(6);
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
  static constexpr std::array<std::array<HighMassMatrix, 5>, 2> matrices = {{
#include "../data/hzz4l/CPSHighMassVBFMatrices.inc"
  }};
  const auto& parameters = matrices.at(is8TeV ? 1 : 0);
  copyHighMassMatrix(parameters.at(0), a_width);
  copyHighMassMatrix(parameters.at(1), a_delta);
  copyHighMassMatrix(parameters.at(2), a_alpha);
  copyHighMassMatrix(parameters.at(3), a_r);
  copyHighMassMatrix(parameters.at(4), a_beta);
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
    const auto& coefficients = vbf8TeVSplineCoefficients;
    const auto& fX = coefficients.at(0);
    const auto& fY = coefficients.at(1);
    const auto& fB = coefficients.at(2);
    const auto& fC = coefficients.at(3);
    const auto& fD = coefficients.at(4);
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
    const auto& coefficients = vbf7TeVSplineCoefficients;
    const auto& fX = coefficients.at(0);
    const auto& fY = coefficients.at(1);
    const auto& fB = coefficients.at(2);
    const auto& fC = coefficients.at(3);
    const auto& fD = coefficients.at(4);
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
  static constexpr std::array<std::array<HighMassMatrix, 2>, 2> matrices = {{
#include "../data/hzz4l/CPSHighMassVBFNoInterfMatrices.inc"
  }};
  const auto& parameters = matrices.at(is8TeV ? 1 : 0);
  copyHighMassMatrix(parameters.at(0), a_width);
  copyHighMassMatrix(parameters.at(1), a_delta);
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
    const auto& coefficients = vbf8TeVSplineCoefficients;
    const auto& fX = coefficients.at(0);
    const auto& fY = coefficients.at(1);
    const auto& fB = coefficients.at(2);
    const auto& fC = coefficients.at(3);
    const auto& fD = coefficients.at(4);
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
    const auto& coefficients = vbf7TeVSplineCoefficients;
    const auto& fX = coefficients.at(0);
    const auto& fY = coefficients.at(1);
    const auto& fB = coefficients.at(2);
    const auto& fC = coefficients.at(3);
    const auto& fD = coefficients.at(4);
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
   const auto& coefficients = gghSplineCoefficients;
   const auto& fX = coefficients.at(0);
   const auto& fY = coefficients.at(1);
   const auto& fB = coefficients.at(2);
   const auto& fC = coefficients.at(3);
   const auto& fD = coefficients.at(4);
   const auto& fE = coefficients.at(5);
   const auto& fF = coefficients.at(6);
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

