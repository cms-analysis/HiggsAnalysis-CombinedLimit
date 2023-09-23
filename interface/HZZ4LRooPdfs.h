#ifndef HZZ4LROOPDFS
#define HZZ4LROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "TDirectory.h"
#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <string>

/*
namespace RooFit{
	
	void readFile();
	
	const Double_t FracEventsNoBrem_4mu = 0.703618;
	const Double_t FracEventsNoBrem_4e = 0.583196;
	const Double_t FracEventsNoBrem_2e2mu = 0.641297;
	
	Double_t Lgg_7(Double_t mHstar);
	Double_t HiggsWidth(Int_t ID,Double_t mHrequested);
	Double_t pdf1(double mHstar,double mHreq);
	Double_t rho(double r, TString proc);
	Double_t Sigma(double mHreq, TString proc);
	
	Double_t N(Double_t mH, TString proc);
	Double_t sigma_CB(Double_t mH, TString proc);
	Double_t mean(Double_t mH, TString proc);
	
	Double_t scratchMass;
	Double_t BR[26][217];
	Double_t CS[6][197];
	
}
*/
class RooqqZZPdf : public RooAbsPdf {
public:
	RooqqZZPdf() {} ;
	RooqqZZPdf(const char *name, const char *title,
			   RooAbsReal& _m4l,
			   RooAbsReal& _a1,
			   RooAbsReal& _a2,
			   RooAbsReal& _a3,
			   RooAbsReal& _b1,
			   RooAbsReal& _b2,
			   RooAbsReal& _b3,
			   RooAbsReal& _frac);
	RooqqZZPdf(const RooqqZZPdf& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooqqZZPdf(*this,newname); }
	inline virtual ~RooqqZZPdf() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a1 ;
	RooRealProxy a2 ;
	RooRealProxy a3 ;
	RooRealProxy b1 ;
	RooRealProxy b2 ;
	RooRealProxy b3 ;
	RooRealProxy frac ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooqqZZPdf,1) // Your description goes here...                                                                                                   
};



class RooggZZPdf : public RooAbsPdf {
public:
	RooggZZPdf() {} ;
	RooggZZPdf(const char *name, const char *title,
			   RooAbsReal& _m4l,
			   RooAbsReal& _a1,
			   RooAbsReal& _a2,
			   RooAbsReal& _a3,
			   RooAbsReal& _b1,
			   RooAbsReal& _b2,
			   RooAbsReal& _b3,
			   RooAbsReal& _frac);
	RooggZZPdf(const RooggZZPdf& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooggZZPdf(*this,newname); }
	inline virtual ~RooggZZPdf() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a1 ;
	RooRealProxy a2 ;
	RooRealProxy a3 ;
	RooRealProxy b1 ;
	RooRealProxy b2 ;
	RooRealProxy b3 ;
	RooRealProxy frac ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooggZZPdf,1) // Your description goes here...                                                                                                   
};

// ------- RooqqZZPdf_v2 -------


class RooqqZZPdf_v2 : public RooAbsPdf {
public:
	RooqqZZPdf_v2() {} ;
	RooqqZZPdf_v2(const char *name, const char *title,
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
			   
			   );
	RooqqZZPdf_v2(const RooqqZZPdf_v2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooqqZZPdf_v2(*this,newname); }
	inline virtual ~RooqqZZPdf_v2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a0 ;
	RooRealProxy a1 ;
	RooRealProxy a2 ;
	RooRealProxy a3 ;
	RooRealProxy a4 ;
	RooRealProxy a5 ;
	RooRealProxy a6 ;
	RooRealProxy a7 ;
	RooRealProxy a8 ;
	RooRealProxy a9 ;
	RooRealProxy a10 ;
	RooRealProxy a11 ;
	RooRealProxy a12 ;
	RooRealProxy a13 ;
	
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooqqZZPdf_v2,1) // Your description goes here...                                                                                                   
};


// ------- RooVBFZZPdf -------


class RooVBFZZPdf : public RooAbsPdf {
public:
	RooVBFZZPdf() {} ;
	RooVBFZZPdf(const char *name, const char *title,
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
			   );
	RooVBFZZPdf(const RooVBFZZPdf& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooVBFZZPdf(*this,newname); }
	inline virtual ~RooVBFZZPdf() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a0 ;
	RooRealProxy a1 ;
	RooRealProxy a2 ;
	RooRealProxy a3 ;
	RooRealProxy a4 ;
	RooRealProxy a5 ;
	RooRealProxy a6 ;
	RooRealProxy a7 ;
	RooRealProxy a8 ;
	RooRealProxy a9 ;
	RooRealProxy a10 ;
	RooRealProxy a11 ;
	RooRealProxy a12 ;
	RooRealProxy a13 ;
	RooRealProxy a14 ;
	RooRealProxy a15 ;
	RooRealProxy a16 ;
	
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooVBFZZPdf,1) // Your description goes here...                                                                                                   
};



// ------- RooVBFZZPdf_v2 -------


class RooVBFZZPdf_v2 : public RooAbsPdf {
public:
	RooVBFZZPdf_v2() {} ;
	RooVBFZZPdf_v2(const char *name, const char *title,
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
			   );
	RooVBFZZPdf_v2(const RooVBFZZPdf_v2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooVBFZZPdf_v2(*this,newname); }
	inline virtual ~RooVBFZZPdf_v2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a4 ;
	RooRealProxy a5 ;
	RooRealProxy a6 ;
	RooRealProxy a7 ;
	RooRealProxy a8 ;
	RooRealProxy a9 ;
	RooRealProxy a10 ;
	RooRealProxy a11 ;
	RooRealProxy a12 ;
	RooRealProxy a13 ;
	RooRealProxy a14 ;
	RooRealProxy a15 ;
	RooRealProxy a16 ;
	
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooVBFZZPdf_v2,1) // Your description goes here...                                                                                                   
};




class RooggZZPdf_v2 : public RooAbsPdf {
public:
	RooggZZPdf_v2() {} ;
	RooggZZPdf_v2(const char *name, const char *title,
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
			   
			   );
	RooggZZPdf_v2(const RooggZZPdf_v2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooggZZPdf_v2(*this,newname); }
	inline virtual ~RooggZZPdf_v2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy a0 ;
	RooRealProxy a1 ;
	RooRealProxy a2 ;
	RooRealProxy a3 ;
	RooRealProxy a4 ;
	RooRealProxy a5 ;
	RooRealProxy a6 ;
	RooRealProxy a7 ;
	RooRealProxy a8 ;
	RooRealProxy a9 ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooggZZPdf_v2,1) // Your description goes here...                                                                                                   
};

class RooBetaFunc_v2 : public RooAbsPdf {
public:
	RooBetaFunc_v2();
	RooBetaFunc_v2(const char *name, const char *title,
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
				   );
	RooBetaFunc_v2(const RooBetaFunc_v2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooBetaFunc_v2(*this,newname); }
	inline virtual ~RooBetaFunc_v2() { }
	
protected:
	
    RooRealProxy mZstar;	     
    RooRealProxy mZ;     
	RooRealProxy m0;     
    RooRealProxy mZZ;     
    RooRealProxy Gamma;
	RooRealProxy Gamma0;
    RooRealProxy a0;  // mZZ distribution vars
    RooRealProxy a1;  // mZZ distribution vars
    RooRealProxy a2;
    RooRealProxy a3;
    RooRealProxy f;
	RooRealProxy f0;
	
    Double_t evaluate() const ;
	
private:
	
	ClassDef(RooBetaFunc_v2,1) 
};

class Roo4lMasses2D_Bkg : public RooAbsPdf {
public:
	Roo4lMasses2D_Bkg();
	Roo4lMasses2D_Bkg(const char *name, const char *title,
					  RooAbsReal& _mZstar,	       
					  RooAbsReal& _mZZ,	     
					  RooAbsReal& _channelVal	     
					  );
	Roo4lMasses2D_Bkg(const Roo4lMasses2D_Bkg& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new Roo4lMasses2D_Bkg(*this,newname); }
	inline virtual ~Roo4lMasses2D_Bkg() { }
	
protected:
	
    RooRealProxy mZstar;	     
    RooRealProxy mZZ;     
    RooRealProxy channelVal;     
	
    Double_t evaluate() const ;
    Double_t UnitStep(double arg) const;
private:
	
	ClassDef(Roo4lMasses2D_Bkg,1) 
};

//------------------------

class Roo4lMasses2D_BkgGGZZ : public RooAbsPdf {
public:
	Roo4lMasses2D_BkgGGZZ();
	Roo4lMasses2D_BkgGGZZ(const char *name, const char *title,
					  RooAbsReal& _mZstar,	       
					  RooAbsReal& _mZZ,	     
					  RooAbsReal& _channelVal	     
					  );
	Roo4lMasses2D_BkgGGZZ(const Roo4lMasses2D_BkgGGZZ& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new Roo4lMasses2D_BkgGGZZ(*this,newname); }
	inline virtual ~Roo4lMasses2D_BkgGGZZ() { }
	
protected:
	
    RooRealProxy mZstar;	     
    RooRealProxy mZZ;     
    RooRealProxy channelVal;     
	
    Double_t evaluate() const ;
    Double_t UnitStep(double arg) const;
private:
	
	ClassDef(Roo4lMasses2D_BkgGGZZ,1) 
};

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// backgrounds above
// --------------------------------------------------------------------
// --------------------------------------------------------------------

// --------------------------------------
// 2D signal
class Roo4lMasses2D : public RooAbsPdf {
public:
	Roo4lMasses2D();
	Roo4lMasses2D(const char *name, const char *title,
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
				  );
	Roo4lMasses2D(const Roo4lMasses2D& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new Roo4lMasses2D(*this,newname); }
	inline virtual ~Roo4lMasses2D() { }
	
protected:
	
    RooRealProxy mZstar;             
    RooRealProxy mZ;     
    RooRealProxy mZZ;     
    RooRealProxy Gamma;     
    RooRealProxy p0;     
    RooRealProxy p1;     
    RooRealProxy p2;     
    RooRealProxy CBmean;    
    RooRealProxy CBwidth;    
    RooRealProxy CBalpha;            
    RooRealProxy CBn;             
	
    Double_t evaluate() const ;
	
private:
	
	ClassDef(Roo4lMasses2D,1) 
};

// --------------------------------------

class RooFourMuMassShapePdf2 : public RooAbsPdf {
public:
	RooFourMuMassShapePdf2() {} ;
	RooFourMuMassShapePdf2(const char *name, const char *title,
						   RooAbsReal& _m4l,
						   RooAbsReal& _mH);
	RooFourMuMassShapePdf2(const RooFourMuMassShapePdf2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooFourMuMassShapePdf2(*this,newname); }
	inline virtual ~RooFourMuMassShapePdf2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	
	Double_t evaluate() const ;
	//void readFile() const ;
	
private:
	
	ClassDef(RooFourMuMassShapePdf2,2) // Your description goes here...                                                                                       
};


class RooFourEMassShapePdf2 : public RooAbsPdf {
public:
	RooFourEMassShapePdf2() {} ;
	RooFourEMassShapePdf2(const char *name, const char *title,
						  RooAbsReal& _m4l,
						  RooAbsReal& _mH);
	RooFourEMassShapePdf2(const RooFourEMassShapePdf2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooFourEMassShapePdf2(*this,newname); }
	inline virtual ~RooFourEMassShapePdf2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	
	Double_t evaluate() const ;
	//void readFile() const ;
	
private:
	
	ClassDef(RooFourEMassShapePdf2,2) // Your description goes here...                                                                                        
};



class RooTwoETwoMuMassShapePdf2 : public RooAbsPdf {
public:
	RooTwoETwoMuMassShapePdf2() {} ;
	RooTwoETwoMuMassShapePdf2(const char *name, const char *title,
							  RooAbsReal& _m4l,
							  RooAbsReal& _mH);
	RooTwoETwoMuMassShapePdf2(const RooTwoETwoMuMassShapePdf2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooTwoETwoMuMassShapePdf2(*this,newname); }
	inline virtual ~RooTwoETwoMuMassShapePdf2() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	
	Double_t evaluate() const ;
	//void readFile() const ;
	
private:	
	
	ClassDef(RooTwoETwoMuMassShapePdf2,2) // Your description goes here...                                                                                    
};


class RooFourMuMassRes : public RooAbsPdf {
public:
	RooFourMuMassRes() {} ;
	RooFourMuMassRes(const char *name, const char *title,
					 RooAbsReal& _m4l,
					 RooAbsReal& _mH);
	RooFourMuMassRes(const RooFourMuMassRes& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooFourMuMassRes(*this,newname); }
	inline virtual ~RooFourMuMassRes() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH  ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooFourMuMassRes,1) // Your description goes here...                                                                                             
};

class RooFourEMassRes : public RooAbsPdf {
public:
	RooFourEMassRes() {} ;
	RooFourEMassRes(const char *name, const char *title,
					RooAbsReal& _m4l,
					RooAbsReal& _mH);
	RooFourEMassRes(const RooFourEMassRes& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooFourEMassRes(*this,newname); }
	inline virtual ~RooFourEMassRes() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH  ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooFourEMassRes,1) // Your description goes here...                                                                                              
};


class RooTwoETwoMuMassRes : public RooAbsPdf {
public:
	RooTwoETwoMuMassRes() {} ;
	RooTwoETwoMuMassRes(const char *name, const char *title,
						RooAbsReal& _m4l,
						RooAbsReal& _mH);
	RooTwoETwoMuMassRes(const RooTwoETwoMuMassRes& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooTwoETwoMuMassRes(*this,newname); }
	inline virtual ~RooTwoETwoMuMassRes() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH  ;
	
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooTwoETwoMuMassRes,1) // Your description goes here...                                                                                          
};

class RooRelBW1 : public RooAbsPdf {
public:
	RooRelBW1() {} ;
	RooRelBW1(const char *name, const char *title,
			  RooAbsReal& _m,
			  RooAbsReal& _mean,
			  RooAbsReal& _gamma);
	RooRelBW1(const RooRelBW1& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBW1(*this,newname); }
	inline virtual ~RooRelBW1() { }
	
protected:
	
	RooRealProxy m ;
	RooRealProxy mean ;
	RooRealProxy gamma ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBW1,1) // Your description goes here...                                                                                                    
};

//////////////////////////////////////////////

class RooRelBWUF : public RooAbsPdf {
public:
	RooRelBWUF() {} ;
	RooRelBWUF(const char *name, const char *title,
			  RooAbsReal& _m4l,
			  RooAbsReal& _mH);
	RooRelBWUF(const RooRelBWUF& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBWUF(*this,newname); }
	inline virtual ~RooRelBWUF() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBWUF,2) // Your description goes here...                                                                                                    
};


//////////////////////////////////////////////

class RooRelBWUF_SM4 : public RooAbsPdf {
public:
	RooRelBWUF_SM4() {} ;
	RooRelBWUF_SM4(const char *name, const char *title,
			  RooAbsReal& _m4l,
			  RooAbsReal& _mH);
	RooRelBWUF_SM4(const RooRelBWUF_SM4& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBWUF_SM4(*this,newname); }
	inline virtual ~RooRelBWUF_SM4() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBWUF_SM4,2) // Your description goes here...                                                                                                    
};


//////////////////////////////////////////////

class RooRelBWUFParamWidth : public RooAbsPdf {
public:
	RooRelBWUFParamWidth() {} ;
	RooRelBWUFParamWidth(const char *name, const char *title,
					RooAbsReal& _m4l,
					RooAbsReal& _mH,
					RooAbsReal& _width);
	RooRelBWUFParamWidth(const RooRelBWUFParamWidth& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBWUFParamWidth(*this,newname); }
	inline virtual ~RooRelBWUFParamWidth() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	RooRealProxy width;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBWUFParamWidth,2) // Your description goes here...                                                                                                    
};


class RooRelBWUFParam : public RooAbsPdf {
public:
	RooRelBWUFParam() {} ;
	RooRelBWUFParam(const char *name, const char *title,
			RooAbsReal& _m4l,
			RooAbsReal& _mH,
			RooAbsReal& _scaleParam,
			Double_t _widthSF=1.
			);
	RooRelBWUFParam(const RooRelBWUFParam& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBWUFParam(*this,newname); }
	inline virtual ~RooRelBWUFParam() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	RooRealProxy scaleParam ;
	Double_t widthSF ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBWUFParam,3) // Your description goes here...                                                                                                    
};

//////////////////////////////////////////////

class RooRelBWHighMass : public RooAbsPdf {
public:
	RooRelBWHighMass() {} ;
	RooRelBWHighMass(const char *name, const char *title,
					RooAbsReal& _m4l,
					RooAbsReal& _mH,
					RooAbsReal& _gamma);
	RooRelBWHighMass(const RooRelBWHighMass& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRelBWHighMass(*this,newname); }
	inline virtual ~RooRelBWHighMass() { }
	
protected:
	
	RooRealProxy m4l ;
	RooRealProxy mH ;
	RooRealProxy gamma ;
	
	Double_t evaluate() const ;
	
private:
	
	ClassDef(RooRelBWHighMass,2) // Your description goes here...                                                                                                    
};

///////////////////////////////////////////////////

class RooTsallis : public RooAbsPdf {
public:
  RooTsallis();
  RooTsallis(const char *name, const char *title,
	          RooAbsReal& _x,
        	  RooAbsReal& _m,
              	  RooAbsReal& _n,
	          RooAbsReal& _n2,
                  RooAbsReal& _bb,
	          RooAbsReal& _bb2,
	          RooAbsReal& _T,
	          RooAbsReal& _fexp);

  RooTsallis(const RooTsallis& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTsallis(*this,newname); }
  inline virtual ~RooTsallis() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x ;
  RooRealProxy m ;
  RooRealProxy n ;
  RooRealProxy n2 ;
  RooRealProxy bb ;
  RooRealProxy bb2 ;
  RooRealProxy T ;
  RooRealProxy fexp ;

  Double_t evaluate() const ;

private:

 ClassDef(RooTsallis,1) // Your description goes here...
};



/////////////////////////////////////////////////////////////

class RooaDoubleCBxBW : public RooAbsPdf {
 public:
  RooaDoubleCBxBW();
  RooaDoubleCBxBW(const char *name, const char *title,
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
        bool _computeActualCB
    );
  RooaDoubleCBxBW(const RooaDoubleCBxBW& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooaDoubleCBxBW(*this,newname); }
  inline virtual ~RooaDoubleCBxBW() { }

 protected:

  RooRealProxy x ;
  RooRealProxy shift;
  RooRealProxy sigma;
  RooRealProxy alphaL;
  RooRealProxy alphaR;
  RooRealProxy mean;
  RooRealProxy width;
  unsigned nL;
  unsigned nR;
  RooRealProxy thetaL;
  RooRealProxy thetaR;
  bool computeActualCB;

  Double_t evaluate() const ;
  Double_t evaluateDoubleCB() const ;
  Double_t evaluatePowerLaw(double lim, unsigned power, bool isLeft) const ;
  Double_t evaluateQuadratic(double lim1, double lim2, bool isLeft) const ;
  Double_t evaluateVoigtian() const ;

 private:

  ClassDef(RooaDoubleCBxBW,1)
    };



///////////////////////////////////////////////////
class RooCPSHighMassGGH : public RooAbsPdf {
public:
  RooCPSHighMassGGH(); 
  RooCPSHighMassGGH(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _KPrime,
	      RooAbsReal& _BRnew,
	      RooAbsReal& _IntStr,
	      RooAbsReal& _WidthScl,
	      Bool_t is8TeV);
  RooCPSHighMassGGH(const RooCPSHighMassGGH& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCPSHighMassGGH(*this,newname); }
  Double_t Spline(Double_t xx) const;
  inline virtual ~RooCPSHighMassGGH() {  }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy KPrime ;
  RooRealProxy BRnew ;
  RooRealProxy IntStr ;
  RooRealProxy WidthScl ;

  Bool_t is8TeV;

  Double_t a_width[7][5];
  Double_t a_delta[7][5];
  Double_t a_alpha[7][5];
  Double_t a_r[7][5];
  Double_t a_beta[7][5];
  
  Double_t evaluate() const ;
  void initMatrices() ;
  Double_t interpolateMatrix(const Double_t matrix[][5], const Double_t& x, const Double_t& y) const;
  Double_t getWidth(const Double_t& mass, const Double_t& cprime) const;
  Double_t getDelta(const Double_t& mass, const Double_t& cprime) const;
  Double_t getAlpha(const Double_t& mass, const Double_t& cprime) const;
  Double_t getR(const Double_t& mass, const Double_t& cprime) const;
  Double_t getBeta(const Double_t& mass, const Double_t& cprime) const;
  


private:

  ClassDef(RooCPSHighMassGGH,5)
};



///////////////////////////////////////////////////
class RooBWHighMassGGH : public RooAbsPdf {
public:
  RooBWHighMassGGH(); 
  RooBWHighMassGGH(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _KPrime,
	      RooAbsReal& _BRnew,
	      RooAbsReal& _IntStr,
	      Bool_t is8TeV);
  RooBWHighMassGGH(const RooBWHighMassGGH& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBWHighMassGGH(*this,newname); }
  inline virtual ~RooBWHighMassGGH() {  }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy KPrime ;
  RooRealProxy BRnew ;
  RooRealProxy IntStr ;

  Bool_t is8TeV;

  Double_t a_width[7][5];
  Double_t a_delta[7][5];
  Double_t a_alpha[7][5];
  Double_t a_r[7][5];
  Double_t a_beta[7][5];
  
  Double_t evaluate() const ;
  void initMatrices() ;
  Double_t interpolateMatrix(const Double_t matrix[][5], const Double_t& x, const Double_t& y) const;
  Double_t getWidth(const Double_t& mass, const Double_t& cprime) const;
  Double_t getDelta(const Double_t& mass, const Double_t& cprime) const;
  Double_t getAlpha(const Double_t& mass, const Double_t& cprime) const;
  Double_t getR(const Double_t& mass, const Double_t& cprime) const;
  Double_t getBeta(const Double_t& mass, const Double_t& cprime) const;

  Double_t H_width(Double_t mass) const;
  


private:

  ClassDef(RooBWHighMassGGH,4)
};


///////////////////////////////////////////////////
class RooCPSHighMassGGHNoInterf : public RooAbsPdf {
public:
  RooCPSHighMassGGHNoInterf(); 
  RooCPSHighMassGGHNoInterf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _KPrime,
	      RooAbsReal& _BRnew,
	      Bool_t is8TeV);
  RooCPSHighMassGGHNoInterf(const RooCPSHighMassGGHNoInterf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCPSHighMassGGHNoInterf(*this,newname); }
  Double_t Spline(Double_t xx) const;
  inline virtual ~RooCPSHighMassGGHNoInterf() {  }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy KPrime ;
  RooRealProxy BRnew ;

  Bool_t is8TeV;

  Double_t a_width[7][5];
  Double_t a_delta[7][5];
  
  Double_t evaluate() const ;
  void initMatrices() ;
  Double_t interpolateMatrix(const Double_t matrix[][5], const Double_t& x, const Double_t& y) const;
  Double_t getWidth(const Double_t& mass, const Double_t& cprime) const;
  Double_t getDelta(const Double_t& mass, const Double_t& cprime) const;
  


private:

  ClassDef(RooCPSHighMassGGHNoInterf,4)
};




///////////////////////////////////////////////////
class RooCPSHighMassVBF : public RooAbsPdf {
public:
  RooCPSHighMassVBF(); 
  RooCPSHighMassVBF(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _KPrime,
	      RooAbsReal& _BRnew,
	      RooAbsReal& _IntStr,
	      RooAbsReal& _WidthScl,
	      Bool_t is8TeV);
  RooCPSHighMassVBF(const RooCPSHighMassVBF& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCPSHighMassVBF(*this,newname); }
  Double_t Spline(Double_t xx) const;
  inline virtual ~RooCPSHighMassVBF() {  }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy KPrime ;
  RooRealProxy BRnew ;
  RooRealProxy IntStr ;
  RooRealProxy WidthScl ;

  Bool_t is8TeV;

  Double_t a_width[7][5];
  Double_t a_delta[7][5];
  Double_t a_alpha[7][5];
  Double_t a_r[7][5];
  Double_t a_beta[7][5];
  
  Double_t evaluate() const ;
  void initMatrices() ;
  Double_t interpolateMatrix(const Double_t matrix[][5], const Double_t& x, const Double_t& y) const;
  Double_t getWidth(const Double_t& mass, const Double_t& cprime) const;
  Double_t getDelta(const Double_t& mass, const Double_t& cprime) const;
  Double_t getAlpha(const Double_t& mass, const Double_t& cprime) const;
  Double_t getR(const Double_t& mass, const Double_t& cprime) const;
  Double_t getBeta(const Double_t& mass, const Double_t& cprime) const;



private:

  ClassDef(RooCPSHighMassVBF,4)
};


///////////////////////////////////////////////////
class RooCPSHighMassVBFNoInterf : public RooAbsPdf {
public:
  RooCPSHighMassVBFNoInterf(); 
  RooCPSHighMassVBFNoInterf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _KPrime,
	      RooAbsReal& _BRnew,
	      Bool_t is8TeV);
  RooCPSHighMassVBFNoInterf(const RooCPSHighMassVBFNoInterf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCPSHighMassVBFNoInterf(*this,newname); }
  Double_t Spline(Double_t xx) const;
  inline virtual ~RooCPSHighMassVBFNoInterf() {  }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy KPrime ;
  RooRealProxy BRnew ;

  Bool_t is8TeV;

  Double_t a_width[7][5];
  Double_t a_delta[7][5];
  
  Double_t evaluate() const ;
  void initMatrices() ;
  Double_t interpolateMatrix(const Double_t matrix[][5], const Double_t& x, const Double_t& y) const;
  Double_t getWidth(const Double_t& mass, const Double_t& cprime) const;
  Double_t getDelta(const Double_t& mass, const Double_t& cprime) const;



private:

  ClassDef(RooCPSHighMassVBFNoInterf,4)
};


///////////////////////////////////////////////////
class RooSigPlusInt : public RooAbsPdf {
public:
  RooSigPlusInt(); 
  RooSigPlusInt(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mH,
	      RooAbsReal& _delta,
	      RooAbsReal& _width,
	      RooAbsReal& _k,
	      RooAbsReal& _CSquared,
	      RooAbsReal& _BRnew,
	      RooAbsReal& _alpha,
	      RooAbsReal& _beta,
	      RooAbsReal& _r);
  RooSigPlusInt(const RooSigPlusInt& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooSigPlusInt(*this,newname); }
  inline virtual ~RooSigPlusInt() { }
  Double_t Spline(Double_t xx) const;

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy delta ;
  RooRealProxy width ;
  RooRealProxy k ;
  RooRealProxy CSquared ;
  RooRealProxy BRnew ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  RooRealProxy r ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooSigPlusInt,1) // Your description goes here...
};




#endif
