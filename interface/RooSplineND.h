#ifndef HiggsAnalysis_CombinedLimit_RooSplineND_h
#define HiggsAnalysis_CombinedLimit_RooSplineND_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooListProxy.h>
#include "TMath.h"
#include "TMatrixTSym.h"
#include "TMatrixTSparse.h"
#include "TMatrix.h"
#include "TMatrixF.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TDecompLU.h"
#include "TDecompBK.h"
#include "TDecompQRH.h"
#include "TDecompSparse.h"
#include "TVectorD.h"
#include "TTree.h"
#include "TEventList.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "RooRealVar.h"

#include <map>
#include <vector>
#include <string>
 
/**********************************************************************
Original Author -- Nicholas Wardle

BEGIN_HTML
<p>
Use of radial basis functions for interpolation 
between points in multidimensional space.

Produces N ->1 function from TTree

Branch to be considered as F(x) should be passed in fName, 
otherwise it is assumed to be called, "f"

TODO : 
1) Add additional Radial basis function to choose from (via enums?)

2) Make use of multiple trees in root file to produce N->M mapping 
it should be possible to keep the decomposition once produced to solve 
for different f-vectors

3) Make the decomposition a persistent member. Can be used to apply to 
different functions after decomposition is performed.

4) Any better linear algebra packages from / rather than ROOT

</p>
END_HTML
************************************************************************/

class RooSplineND : public RooAbsReal {

   public:
      RooSplineND() : ndim_(0),M_(0),eps_(3.) {}
      RooSplineND(const char *name, const char *title, RooArgList &vars, TTree *tree, const char* fName="f", double eps=3., std::string cutstring="" ) ;
      RooSplineND(const RooSplineND& other, const char *name) ; 
      RooSplineND(const char *name, const char *title, const RooListProxy &vars, int ndim, int M, double eps, std::vector<double> &w, std::map<int,std::vector<double> > &map, std::map<int,std::pair<double,double> > & ,double,double) ;
      ~RooSplineND() ;

      TObject * clone(const char *newname) const ;

      TGraph* getGraph(const char *xvar, double step) ;

    protected:
        Double_t evaluate() const;

    private:
        RooListProxy vars_;
 
	mutable std::vector<double> w_;
	mutable std::map<int,std::vector<double> > v_map;
	mutable std::map<int,std::pair<double,double> > r_map;
	
	int ndim_;
	int M_;
	double eps_;
  	double axis_pts_;
	
	double w_mean, w_rms;

	void calculateWeights(std::vector<double> &);
	double getDistSquare(int i, int j);
	double getDistFromSquare(int i) const;
	double radialFunc(double d2, double eps) const;
	

  ClassDef(RooSplineND,1) 
};

#endif
