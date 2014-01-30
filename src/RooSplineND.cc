#include "../interface/RooSplineND.h"

RooSplineND::RooSplineND(const char *name, const char *title, RooArgList &vars, TTree *tree, const char *fName, double eps) :
  RooAbsReal(name,title),
  vars_("vars","Variables", this)
{
  ndim_ = vars.getSize();
  M_    = tree->GetEntries();

  std::cout << "RooSplineND -- Num Dimensions == " << ndim_ <<std::endl;
  std::cout << "RooSplineND -- Num Samples    == " << M_ << std::endl;
  double *b_map = new double(ndim_);

  RooAbsReal *rIt;	
  TIterator *iter = vars.createIterator(); int it_c=0;
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
    std::vector<double >tmpv(M_,0); 
    v_map.insert(std::pair<int, std::vector<double> >(it_c,tmpv));
    r_map.insert(std::pair<int, std::pair<double,double> >(it_c,std::pair<double,double>(1.e6,-1e6))); 
    tree->SetBranchAddress(rIt->GetName(),&b_map[it_c]);
    it_c++;
  }
  // Assume the function val (yi) branch is f
  double F;
  tree->SetBranchAddress(fName,&F);
  std::vector<double> F_vec;

  // Run through tree and store points 
  for (int i=0;i<M_;i++){
    tree->GetEntry(i);
    for (int k=0;k<ndim_;k++){
      double cval = b_map[k];
      if (cval < r_map[k].first) r_map[k].first=cval;
      if (cval > r_map[k].second) r_map[k].second=cval;
      v_map[k][i] = cval;
    }
    F_vec.push_back(F);
  }

  // Try to re-scale axis to even out dimensions.
  axis_pts_ = TMath::Power(M_,1./ndim_);
  eps_= eps;
  calculateWeights(F_vec); 
  delete b_map;	
}

//_____________________________________________________________________________
// Copy Constructor
RooSplineND::RooSplineND(const RooSplineND& other, const char *name) :
 RooAbsReal(other, name),vars_("vars",this,RooListProxy())
{
  ndim_ = other.ndim_;
  M_    = other.M_;
  eps_  = other.eps_;
  axis_pts_ = other.axis_pts_;

  RooAbsReal *rIt;	
  TIterator *iter = other.vars_.createIterator();
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
  }

  w_mean = other.w_mean;
  w_rms  = other.w_rms;

  // STL copy constructors
  w_    = other.w_;
  v_map = other.v_map; 
  r_map = other.r_map;

}
//_____________________________________________________________________________
// Clone Constructor
RooSplineND::RooSplineND(const char *name, const char *title, const RooListProxy &vars, 
 int ndim, int M, double eps, std::vector<double> &w, std::map<int,std::vector<double> > &map, std::map<int,std::pair<double,double> > &rmap,double wmean, double wrms) :
 RooAbsReal(name, title),vars_("vars",this,RooListProxy()) 
{

  RooAbsReal *rIt;	
  TIterator *iter = vars.createIterator();
  while( (rIt = (RooAbsReal*) iter->Next()) ){ 
    vars_.add(*rIt);
  }

  ndim_ = ndim;
  M_    = M;
  eps_  = eps;
  axis_pts_ = TMath::Power(M_,1./ndim_);

  w_    = w;
  v_map = map;
  r_map = rmap;
	
  w_rms  = wrms;
  w_mean = wrms;
  
}

//_____________________________________________________________________________
TObject *RooSplineND::clone(const char *newname) const 
{
    return new RooSplineND(newname, this->GetTitle(), 
	vars_,ndim_,M_,eps_,w_,v_map,r_map,w_mean,w_rms);
}
//_____________________________________________________________________________
RooSplineND::~RooSplineND() 
{
}
//_____________________________________________________________________________
TGraph * RooSplineND::getGraph(const char *xvar, double step){

  TGraph *gr = new TGraph();
  gr->SetLineWidth(2);
  RooRealVar* v = (RooRealVar*) vars_.find(xvar);
  gr->SetTitle(v->GetTitle());
  gr->GetXaxis()->SetTitle(xvar);
  double vorig = v->getVal();
  int cp=0;

  for (double xv=v->getMin();xv<=v->getMax();xv+=step){
    v->setVal(xv);
    gr->SetPoint(cp,xv,evaluate());
    cp++;
  }

  v->setVal(vorig);
  return gr;
}
//_____________________________________________________________________________
void RooSplineND::calculateWeights(std::vector<double> &f){

  if (M_==0) {
	w_mean = 0;
	w_rms = 1;
	return;
  }
  // Solve system of Linear equations for weights vector 
  TMatrixTSym<double> fMatrix(M_);
  // Fill the Matrix
  for (int i=0;i<M_;i++){
    fMatrix(i,i)=1.;
    for (int j=i+1;j<M_;j++){
        double d2  = getDistSquare(i,j);
	double rad = radialFunc(d2,eps_);
        fMatrix(i,j) =  rad;
	fMatrix(j,i) =  rad; // it is symmetric	
    }
  }

  TVectorD weights(M_);
  for (int i=0;i<M_;i++){
    weights[i]=(double)f[i];
  }

  TDecompQRH decomp(fMatrix);
  std::cout << "RooSplineND -- Solving for Weights" << std::endl;
  
  decomp.Solve(weights); // Solution now in weights
  std::cout << "RooSplineND -- ........ Done" << std::endl;

  w_mean = 0.;
  for (int i=0;i<M_;i++){
    double tw = weights[i];
    w_.push_back(tw);
    w_mean+=(1./M_)*TMath::Abs(tw);
    w_rms+=(1./M_)*(tw*tw);
  }
  w_rms -= (w_mean*w_mean);
  w_rms = TMath::Sqrt(w_rms);
}
//_____________________________________________________________________________
double RooSplineND::getDistSquare(int i, int j){
  double D = 0.; 
  for (int k=0;k<ndim_;k++){
    double v_i = v_map[k][i];
    double v_j = v_map[k][j];
    double dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    //std::cout << "dimension - " << k << ", Values at pts " << i <<"," << j << " are "<<v_i<< ","<<v_j<< " Distance " <<dk<<std::endl;
    D += dk*dk;
  }
  return D; // only ever use square of distance!
}
//_____________________________________________________________________________
double RooSplineND::getDistFromSquare(int i) const{
  // Read parameters distance from point i in the sample
  double D = 0.; 
  for (int k=0;k<ndim_;k++){
    double v_i = v_map[k][i];
    RooAbsReal *v = (RooAbsReal*)vars_.at(k);
    double v_j = v->getVal();
    double dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    D += dk*dk;
  }
  return D; // only ever use square of distance!
  
}
//_____________________________________________________________________________
double RooSplineND::radialFunc(double d2, double eps) const{
  double expo = (d2/(eps*eps));
  double retval = 1./(1+(TMath::Power(expo,1.5)));
  return retval;
}
//_____________________________________________________________________________
Double_t RooSplineND::evaluate() const {
 double ret = 0;
 for (int i=0;i<M_;i++){
 //  std::cout << "EVAL == "<< i << " " << w_[i] << " " << getDistFromSquare(i) << std::endl;
   double w = w_[i];
   if (w==0) continue;
   ret+=((w/w_mean)*radialFunc(getDistFromSquare(i),eps_));
 }
 ret*=w_mean;
 return ret;
}
//_____________________________________________________________________________

ClassImp(RooSplineND)
