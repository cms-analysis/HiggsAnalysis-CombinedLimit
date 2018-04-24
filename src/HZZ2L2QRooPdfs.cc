#include <iostream>
#include <math.h>
#include "TMath.h"

//#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCB.h"
//#include "HiggsAnalysis/CombinedLimit/interface/RooFermi.h"
//#include "HiggsAnalysis/CombinedLimit/interface/RooRelBW.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;
using namespace std; 

ClassImp(RooCB)

  RooCB::RooCB(){}

RooCB::RooCB(const char *name, const char *title,
	     RooAbsReal& _x,
	     RooAbsReal& _mean,
	     RooAbsReal& _width,
	     RooAbsReal& _alpha,
	     RooAbsReal& _n,
              RooAbsReal& _theta
	     ) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  alpha("alpha","alpha",this,_alpha),
  n("n","n",this,_n),
  theta("theta","theta",this,_theta)
{
}

RooCB::RooCB(const RooCB& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  alpha("alpha",this,other.alpha),
  n("n",this,other.n),
  theta("theta",this,other.theta)
{
}

double RooCB::evaluate() const
{
  double a = cos(theta)*alpha - sin(theta)*width;
  double w = sin(theta)*alpha + cos(theta)*width;

  double t = (x-mean)/w;
  if(a<0) t = -t;

  double absa = fabs((double)a);

  double A = TMath::Power(n/absa,n)*exp(-0.5*absa*absa);
  double B = n/absa-absa;

  if(t >= -absa){
    return exp(-0.5*t*t);
  }else{
    return A/TMath::Power(B-t,n);
  }
}


 ClassImp(RooDoubleCB) 

 RooDoubleCB::RooDoubleCB() :
   RooAbsPdf(),
   x("x", "x", this),
   xp("xp", "xp", this),
   mean("mean", "mean", this),
   width("width", "width", this),
   alpha1("alpha1", "alpha1", this),
   n1("n1", "n1", this),
   alpha2("alpha2", "alpha2", this),
   n2("n2", "n2", this)
 {}

 RooDoubleCB::RooDoubleCB(const char *name, const char *title, 
		    RooAbsReal& _x,
		    RooAbsReal& _mean,
		    RooAbsReal& _width,
		    RooAbsReal& _alpha1,
		    RooAbsReal& _n1,
		    RooAbsReal& _alpha2,
		    RooAbsReal& _n2
		    ) :
   RooAbsPdf(name,title), 
   x("x", "x", this, _x),
   xp("xp", "xp", this),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   alpha1("alpha1","alpha1",this,_alpha1),
   n1("n1","n1",this,_n1),
   alpha2("alpha2","alpha2",this,_alpha2),
   n2("n2","n2",this,_n2)
 { 
 } 


 RooDoubleCB::RooDoubleCB(const char *name, const char *title,
        RooAbsReal& _x,
        RooAbsReal& _xp,
        RooAbsReal& _mean,
        RooAbsReal& _width,
        RooAbsReal& _alpha1,
        RooAbsReal& _n1,
        RooAbsReal& _alpha2,
        RooAbsReal& _n2
 ) :
   RooAbsPdf(name, title),
   x("x", "x", this, _x),
   xp("xp", "xp", this, _xp),
   mean("mean", "mean", this, _mean),
   width("width", "width", this, _width),
   alpha1("alpha1", "alpha1", this, _alpha1),
   n1("n1", "n1", this, _n1),
   alpha2("alpha2", "alpha2", this, _alpha2),
   n2("n2", "n2", this, _n2)
 {
 }

 RooDoubleCB::RooDoubleCB(const RooDoubleCB& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x", this, other.x),
   xp("xp", this, other.xp),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   alpha1("alpha1",this,other.alpha1),
   n1("n1",this,other.n1),
   alpha2("alpha2",this,other.alpha2),
   n2("n2",this,other.n2)

 { 
 } 

 double RooDoubleCB::evaluate() const 
 { 
   double dx=0;
   if (x.absArg()) dx += x;
   if (xp.absArg()) dx -= xp;
   double t = (dx-mean)/width;
   if(t>=-alpha1 && t<=alpha2){
     return exp(-0.5*t*t);
   }else if(t<-alpha1){
     double A1 = pow(n1/fabs(alpha1),n1)*exp(-alpha1*alpha1/2);
     double B1 = n1/fabs(alpha1)-fabs(alpha1);
     return A1*pow(B1-t,-n1);
   }else if(t>alpha2){
     double A2 = pow(n2/fabs(alpha2),n2)*exp(-alpha2*alpha2/2);
     double B2 = n2/fabs(alpha2)-fabs(alpha2);
     return A2*pow(B2+t,-n2);
   }else{
     cout << "ERROR evaluating range..." << endl;
     return 99;
   }
   
 } 

 Int_t RooDoubleCB::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const 
 {
   if (_forceNumInt) return 0;

   Int_t code=0;

   RooArgSet deps_x, deps_xp, deps_mean, deps_width, deps_alpha1, deps_alpha2, deps_n1, deps_n2;
   RooRealVar* rrv_x = dynamic_cast<RooRealVar*>(x.absArg()); if (!rrv_x && x.absArg()) x.absArg()->leafNodeServerList(&deps_x, 0, true);
   RooRealVar* rrv_xp = dynamic_cast<RooRealVar*>(xp.absArg()); if (!rrv_xp && xp.absArg()) xp.absArg()->leafNodeServerList(&deps_xp, 0, true);
   RooRealVar* rrv_mean = dynamic_cast<RooRealVar*>(mean.absArg()); if (!rrv_mean) mean.absArg()->leafNodeServerList(&deps_mean, 0, true);
   RooRealVar* rrv_width = dynamic_cast<RooRealVar*>(width.absArg()); if (!rrv_width) width.absArg()->leafNodeServerList(&deps_width, 0, true);
   RooRealVar* rrv_alpha1 = dynamic_cast<RooRealVar*>(alpha1.absArg()); if (!rrv_alpha1) alpha1.absArg()->leafNodeServerList(&deps_alpha1, 0, true);
   RooRealVar* rrv_alpha2 = dynamic_cast<RooRealVar*>(alpha2.absArg()); if (!rrv_alpha2) alpha2.absArg()->leafNodeServerList(&deps_alpha2, 0, true);
   RooRealVar* rrv_n1 = dynamic_cast<RooRealVar*>(n1.absArg()); if (!rrv_n1) n1.absArg()->leafNodeServerList(&deps_n1, 0, true);
   RooRealVar* rrv_n2 = dynamic_cast<RooRealVar*>(n2.absArg()); if (!rrv_n2) n2.absArg()->leafNodeServerList(&deps_n2, 0, true);

   if (rrv_x){
     if (
       (deps_xp.find(*rrv_x)==nullptr || rrv_xp)
       ||
       (deps_mean.find(*rrv_x)==nullptr || rrv_mean)
       ||
       (deps_width.find(*rrv_x)==nullptr || rrv_width)
       ||
       (deps_alpha1.find(*rrv_x)==nullptr || rrv_alpha1)
       ||
       (deps_alpha2.find(*rrv_x)==nullptr || rrv_alpha2)
       ||
       (deps_n1.find(*rrv_x)==nullptr || rrv_n1)
       ||
       (deps_n2.find(*rrv_x)==nullptr || rrv_n2)
       ){
       if (matchArgs(allVars, analVars, x)) code+=1;
     }
   }
   if (rrv_xp){
     if (
       (deps_x.find(*rrv_xp)==nullptr || rrv_x)
       ||
       (deps_mean.find(*rrv_xp)==nullptr || rrv_mean)
       ||
       (deps_width.find(*rrv_xp)==nullptr || rrv_width)
       ||
       (deps_alpha1.find(*rrv_xp)==nullptr || rrv_alpha1)
       ||
       (deps_alpha2.find(*rrv_xp)==nullptr || rrv_alpha2)
       ||
       (deps_n1.find(*rrv_xp)==nullptr || rrv_n1)
       ||
       (deps_n2.find(*rrv_xp)==nullptr || rrv_n2)
       ){
       if (matchArgs(allVars, analVars, xp)) code+=2;
     }
   }
   code=code%3;
   return code;
 }

 Double_t RooDoubleCB::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   static const Double_t root2 = sqrt(2.);
   static const Double_t rootPiBy2 = sqrt(atan2(0., -1.)/2.);

   assert(code>0 && code<3);

   double central=0;
   double left=0;
   double right=0;

   switch (code){
   case 1:
   {
     double xmin = x.min(rangeName);
     double xmax = x.max(rangeName);
     if (xp.absArg()){
       xmin -= xp;
       xmax -= xp;
     }

     Double_t xscale = root2*width;

     //compute gaussian contribution
     double central_low =max(xmin, mean - alpha1*width);
     double central_high=min(xmax, mean + alpha2*width);
     if (central_low < central_high) central = rootPiBy2*width*(TMath::Erf((central_high-mean)/xscale)-TMath::Erf((central_low-mean)/xscale));

     //compute left tail;
     double A1 = pow(n1/fabs(alpha1), n1)*exp(-alpha1*alpha1/2);
     double B1 = n1/fabs(alpha1)-fabs(alpha1);

     double left_low=xmin;
     double left_high=min(xmax, mean - alpha1*width);
     if (left_low < left_high){ //is the left tail in range?
       if (fabs(n1-1.)>1.e-5) left = A1/(-n1+1.)*width*(pow(B1-(left_low-mean)/width, -n1+1.)-pow(B1-(left_high-mean)/width, -n1+1.));
       else left = A1*width*(log(B1-(left_low-mean)/width) - log(B1-(left_high-mean)/width));
     }

     //compute right tail;
     double A2 = pow(n2/fabs(alpha2), n2)*exp(-alpha2*alpha2/2);
     double B2 = n2/fabs(alpha2)-fabs(alpha2);

     double right_low=max(xmin, mean + alpha2*width);
     double right_high=xmax;
     if (right_low < right_high){ //is the right tail in range?
       if (fabs(n2-1.)>1.e-5) right = A2/(-n2+1.)*width*(pow(B2+(right_high-mean)/width, -n2+1.)-pow(B2+(right_low-mean)/width, -n2+1.));
       else right = A2*width*(log(B2+(right_high-mean)/width) - log(B2+(right_low-mean)/width));
     }

     break;
   }
   case 2:
   {
     double xmin = -xp.max(rangeName);
     double xmax = -xp.min(rangeName);
     if (x.absArg()){
       xmin += x;
       xmax += x;
     }

     Double_t xscale = root2*width;

     //compute gaussian contribution
     double central_low =max(xmin, mean - alpha1*width);
     double central_high=min(xmax, mean + alpha2*width);
     if (central_low < central_high) central = rootPiBy2*width*(TMath::Erf((central_high-mean)/xscale)-TMath::Erf((central_low-mean)/xscale));

     //compute left tail;
     double A1 = pow(n1/fabs(alpha1), n1)*exp(-alpha1*alpha1/2);
     double B1 = n1/fabs(alpha1)-fabs(alpha1);

     double left_low=xmin;
     double left_high=min(xmax, mean - alpha1*width);
     if (left_low < left_high){ //is the left tail in range?
       if (fabs(n1-1.)>1.e-5) left = A1/(-n1+1.)*width*(pow(B1-(left_low-mean)/width, -n1+1.)-pow(B1-(left_high-mean)/width, -n1+1.));
       else left = A1*width*(log(B1-(left_low-mean)/width) - log(B1-(left_high-mean)/width));
     }

     //compute right tail;
     double A2 = pow(n2/fabs(alpha2), n2)*exp(-alpha2*alpha2/2);
     double B2 = n2/fabs(alpha2)-fabs(alpha2);

     double right_low=max(xmin, mean + alpha2*width);
     double right_high=xmax;
     if (right_low < right_high){ //is the right tail in range?
       if (fabs(n2-1.)>1.e-5) right = A2/(-n2+1.)*width*(pow(B2+(right_high-mean)/width, -n2+1.)-pow(B2+(right_low-mean)/width, -n2+1.));
       else right = A2*width*(log(B2+(right_high-mean)/width) - log(B2+(right_low-mean)/width));
     }

     break;
   }
   }

   return left+central+right;

 }

 ClassImp(RooFermi) 

 RooFermi::RooFermi(){}

 RooFermi::RooFermi(const char *name, const char *title, 
		      RooAbsReal& _x,
		      RooAbsReal& _cutOff,
		    RooAbsReal& _beta
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   cutOff("cutOff","cutOff",this,_cutOff),
   beta("beta","beta",this,_beta)
 { 
 } 


 RooFermi::RooFermi(const RooFermi& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   cutOff("cutOff",this,other.cutOff),
   beta("beta",this,other.beta)

 { 
 } 



 double RooFermi::evaluate() const 
 { 
   return 1.0/(exp((cutOff-x)/beta)+1);
 } 

 ClassImp(RooRelBW) 

 RooRelBW::RooRelBW(){}

 RooRelBW::RooRelBW(const char *name, const char *title, 
		    RooAbsReal& _x,
		    RooAbsReal& _mean,
		    RooAbsReal& _width,
		    RooAbsReal& _n
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   n("n","n",this,_n)
 { 
 } 


 RooRelBW::RooRelBW(const RooRelBW& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   n("n",this,other.n)

 { 
 } 



 double RooRelBW::evaluate() const 
 { 
   return pow(x*x,n)/((x*x-mean*mean)*(x*x-mean*mean)+pow(x*x/(mean*mean),2*n)*mean*mean*width*width);
 } 


ClassImp(Triangle)

  Triangle::Triangle(){}

Triangle::Triangle(const char *name, const char *title,                
		   RooAbsReal& _m,
		   RooAbsReal& _start,
		   RooAbsReal& _turn,
		   RooAbsReal& _stop
		   ):
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  start("start","start",this,_start),
  turn("turn","turn",this,_turn),
  stop("stop","stop",this,_stop)
{
}

Triangle::Triangle(const Triangle& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m),start("start", this, other.start), turn("turn", this, other.turn), stop("stop", this, other.stop)
{
}

Double_t Triangle::evaluate() const 
{
  //std::cout << m << " "<<1.+(start-m)/turn << " " << 1+(turn-m)/stop << std::endl;
  if(m<turn  && m > turn+start)
    return 1.+(turn-m)/start;
  if(m>=turn && m < turn+stop)
    return 1.+(turn-m)/stop;
  
  return 0;
}


Int_t Triangle::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const 
{
  if (matchArgs(allVars,analVars,m)) return 1;
  return 0;
}

Double_t Triangle::analyticalIntegral(Int_t code, const char* rangeName) const 
{

  // WARNING, ASSSUMES TURN TO BE IN INTERVAL
  assert(code==1) ;
  //whole triangle
  Double_t sumleft = sqrt(1+ (turn+start)*(turn+start) ) ;
  Double_t sumright= sqrt(1+ (turn+stop)*(turn+stop) );


  if(m.min() < turn+start)// correct for left missing bit
    sumleft -= sumleft*(m.min()-(turn+start))/fabs(start);


  if(m.max() > turn+stop)// correct for right missing bit
    sumright -= sumright*(turn+stop -m.max())/fabs(stop);

  

  return sumleft+sumright;    
}




ClassImp(RooLevelledExp)

  RooLevelledExp::RooLevelledExp(){}

RooLevelledExp::RooLevelledExp(const char *name, const char *title,
			       RooAbsReal& _x,
			       RooAbsReal& _sigma, 
			       RooAbsReal& _alpha,
			       RooAbsReal& _m,
			       RooAbsReal& _theta):
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  sigma("sigma","sigma",this,_sigma),
  alpha("alpha","alpha",this,_alpha),
  m("m","m",this,_m),
  //  k("k","k",this,_k),
  theta("theta","theta",this,_theta)
{
}

RooLevelledExp::RooLevelledExp(const RooLevelledExp& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  sigma("sigma",this,other.sigma),
  alpha("alpha",this,other.alpha),
  m("m",this,other.m),
  theta("theta",this,other.theta)
{
}

double RooLevelledExp::evaluate() const
{
  double res=0.0;
  double s = cos(theta)*sigma - sin(theta)*alpha;
  double a = sin(theta)*sigma + cos(theta)*alpha;
    
  //original
  double t = fabs(x-m);
  double den = (s + a*t);
  res=exp(-1.0*t/den);
  

  return res;
}
