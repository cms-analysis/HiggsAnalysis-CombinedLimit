#include "../interface/RooSplineND.h"
//#include </afs/cern.ch/work/n/nckw/combine-versions/102x/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/cpStudies/eigen/Eigen/Dense>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

RooSplineND::RooSplineND(const char *name, const char *title, RooArgList &vars, TTree *tree, const char *fName, double eps, bool rescale, std::string cutstring) :
  RooAbsReal(name,title),
  vars_("vars","Variables", this)
{
  rescaleAxis = rescale;
  ndim_ = vars.getSize();
  int nentries  = tree->GetEntries();

  // Get Selection list
  tree->Draw(">>cutlist",cutstring.c_str());
  TEventList *keep_points = (TEventList*)gDirectory->Get("cutlist");
  M_	= keep_points->GetN(); 

  std::cout << "RooSplineND -- Num Dimensions == " << ndim_ <<std::endl;
  std::cout << "RooSplineND -- Num Samples    == " << M_ << std::endl;

  float *b_map = new float(ndim_);

  int it_c=0;
  for (RooAbsArg *rIt : vars) {
    vars_.add(*rIt);
    std::vector<double >tmpv(M_,0); 
    v_map.insert(std::pair<int, std::vector<double> >(it_c,tmpv));
    r_map.insert(std::pair<int, std::pair<double,double> >(it_c,std::pair<double,double>(1.e6,-1e6))); 
    tree->SetBranchAddress(rIt->GetName(),&b_map[it_c]);
    it_c++;
  }

  float F;
  tree->SetBranchAddress(fName,&F);
  std::vector<double> F_vec;

  // Run through tree and store points if selected
  int nselect=0;
  for (int i=0;i<nentries;i++){
    if ( !(keep_points->Contains(i)) ) continue;
    tree->GetEntry(i);
    //std::cout <<"Adding point " << i << ", ";
    for (int k=0;k<ndim_;k++){
      float cval = b_map[k];
      if (cval < r_map[k].first) r_map[k].first=cval;
      if (cval > r_map[k].second) r_map[k].second=cval;
      v_map[k][nselect] = (double)cval;
      //std::cout << "x" << k << "=" << cval << ", ";
    }
    //std::cout << "F=" << F << std::endl;
    F_vec.push_back(F);
    nselect++;
  }
  //std::cout << "... N(selected) = " << nselect << std::endl;
  keep_points->Reset();
  tree->SetEventList(0);
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

  vars_.add(other.vars_);

  w_mean = other.w_mean;
  w_rms  = other.w_rms;

  // STL copy constructors not so helpful :/
  std::vector<double>::const_iterator 		     w_it; 
  std::map<int,std::vector<double> >::const_iterator v_it;
  std::map<int,std::pair<double,double> >::const_iterator  r_it;
  
  for (w_it = other.w_.begin(); w_it != other.w_.end(); w_it++){
    w_.push_back(*w_it); 
  }
  
  for (v_it = other.v_map.begin(); v_it != other.v_map.end(); v_it++){
    std::vector<double>::const_iterator it; 
    std::vector<double> this_vec;
    for (it = (v_it->second).begin(); it!=(v_it->second).end(); it++) this_vec.push_back(*it);
    v_map.insert(std::pair<int, std::vector<double> >(v_it->first,this_vec)); 
  }
  
  for (r_it = other.r_map.begin(); r_it != other.r_map.end(); r_it++){
    r_map.insert(std::pair<int, std::pair<double, double> >(r_it->first,std::pair<double,double>((r_it->second).first,(r_it->second).second))); 
  }
  
  rescaleAxis=other.rescaleAxis;
}
//_____________________________________________________________________________
// Clone Constructor

RooSplineND::RooSplineND(const char *name, const char *title, const RooListProxy &vars, 
 int ndim, int M, double eps, bool rescale, std::vector<double> &w, std::map<int,std::vector<double> > &map, std::map<int,std::pair<double,double> > &rmap,double wmean, double wrms) :
 RooAbsReal(name, title),vars_("vars",this,RooListProxy()) 
{
  vars_.add(vars);
  ndim_ = ndim;
  M_    = M;
  eps_  = eps;
  axis_pts_ = TMath::Power(M_,1./ndim_);

  // STL copy constructors not so helpful :/
  std::vector<double>::const_iterator 		     w_it; 
  std::map<int,std::vector<double> >::const_iterator v_it;
  std::map<int,std::pair<double,double> >::const_iterator  r_it;
  
  for (w_it = w.begin(); w_it != w.end(); w_it++){
    w_.push_back(*w_it); 
  }
  
  for (v_it = map.begin(); v_it != map.end(); v_it++){
    std::vector<double>::const_iterator it; 
    std::vector<double> this_vec;
    for (it = (v_it->second).begin(); it!=(v_it->second).end(); it++) this_vec.push_back(*it);
    v_map.insert(std::pair<int, std::vector<double> >(v_it->first,this_vec)); 
  }
  
  for (r_it = rmap.begin(); r_it != rmap.end(); r_it++){
    r_map.insert(std::pair<int, std::pair<double, double> >(r_it->first,std::pair<double,double>((r_it->second).first,(r_it->second).second))); 
  }
	
  w_rms  = wrms;
  w_mean = wmean;
  
  rescaleAxis = rescale;
}

//_____________________________________________________________________________
TObject *RooSplineND::clone(const char *newname) const 
{
    return new RooSplineND(*this, newname);
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
void RooSplineND::printPoint(int i) const{

  std::cout  << " point - " << i ;
  for (int k=0;k<ndim_;k++){
    double v_i = v_map[k][i];
    std::cout << ", x"<<k<<"="<<v_i; 
  }
  std::cout  << std::endl;
  
}
//_____________________________________________________________________________
void RooSplineND::calculateWeights(std::vector<double> &f){

  std::cout << "RooSplineND -- Solving for Weights" << std::endl;
  if (M_==0) {
	w_mean = 0;
	w_rms = 1;
	return;
  }
  
  MatrixXd fMatrix(M_,M_);
  for (int i=0;i<M_;i++){
    fMatrix(i,i)=1.;
    for (int j=i+1;j<M_;j++){
        double d2  = getDistSquare(i,j);
	if (d2 < 0.0001) {
		std::cout << " ERROR  - points likely duplicated, which will lead to errors in solving for weights. \
		The distance^2 is smaller than 0.0001 for points "<< i << " and " << j << " ... " <<  std::endl;
		printPoint(i);
		printPoint(j);
	}
	double rad = radialFunc(d2,eps_);
        fMatrix(i,j) =  rad;
	fMatrix(j,i) =  rad; // it is symmetric	
    }
  }
  
  VectorXd weights(M_);
  for (int i=0;i<M_;i++) weights(i)=f[i];
  
  VectorXd x = fMatrix.colPivHouseholderQr().solve(weights);

  std::cout << "RooSplineND -- ........ Done" << std::endl;

  w_mean = 0.;
  for (int i=0;i<M_;i++){
    //double tw = weights[i];
    double tw = x(i);
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
    double dk; 
    if (rescaleAxis) dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    else  dk = (v_i-v_j);
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
    double dk; 
    if (rescaleAxis) dk = axis_pts_*(v_i-v_j)/(r_map[k].second-r_map[k].first);
    else dk = (v_i-v_j);
    D += dk*dk;
  }
  return D; // only ever use square of distance!
  
}
//_____________________________________________________________________________
double RooSplineND::radialFunc(double d2, double eps, double cutoff) const{
  double expo = (d2/(eps*eps));
  //double retval = 1./(1+(TMath::Power(expo,1.5)));
  //if (cutoff > 0){
    //if ( TMath::Sqrt(d2) > cutoff*eps ) return 0.;
  //}
  double retval = TMath::Exp(-1*expo);
  return retval;
}
//_____________________________________________________________________________
Double_t RooSplineND::evaluate() const {
 double ret = 0;
 for (int i=0;i<M_;i++){
   //std::cout << "EVAL == "<< i << " " << w_[i] << " " << getDistFromSquare(i) << std::endl;
   double w = w_[i];
   if (w==0) continue;
   //if ( TMath::Abs(w)< 0.01*TMath::Abs(w_mean) )  continue;  
   ret+=((w)*radialFunc(getDistFromSquare(i),eps_));
 }
 //ret*=w_mean;
 return ret;
}
//_____________________________________________________________________________

ClassImp(RooSplineND)
