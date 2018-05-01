#include "HiggsAnalysis/CombinedLimit/interface/RooPiecewisePolynomial.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TIterator.h"


using namespace std;


RooPiecewisePolynomial::RooPiecewisePolynomial(const int nfcn_, const int polyndof_) :
  RooAbsReal(),
  xvar("xvar", "xvar", this),
  parList("parList", "parList", this),
  nfcn(nfcn_), polyndof(polyndof_),
  nnodes(nfcn-1), // How many nodes are in between
  ndof_endfcn(polyndof-1), // 1 degree for the function value
  ndof_middlefcn(polyndof-2) // +1 for slope of the other node
{
  assert((nfcn>2 && polyndof>=2) || (nfcn==2 && polyndof>=1) || (nfcn==1 && polyndof>=0));
}
RooPiecewisePolynomial::RooPiecewisePolynomial(const char* name, const char* title, const int nfcn_, const int polyndof_) :
  RooAbsReal(name, title),
  xvar("xvar", "xvar", this),
  parList("parList", "parList", this),
  nfcn(nfcn_), polyndof(polyndof_),
  nnodes(nfcn-1), // How many nodes are in between
  ndof_endfcn(polyndof-1), // 1 degree for the function value
  ndof_middlefcn(polyndof-2) // +1 for slope of the other node
{
  assert((nfcn>2 && polyndof>=2) || (nfcn==2 && polyndof>=1) || (nfcn==1 && polyndof>=0));
}
RooPiecewisePolynomial::RooPiecewisePolynomial(const char* name, const char* title, RooAbsReal& xvar_, RooArgList const& parList_, const int nfcn_, const int polyndof_) :
  RooAbsReal(name, title),
  xvar("xvar", "xvar", this, xvar_),
  parList("parList", "parList", this),
  nfcn(nfcn_), polyndof(polyndof_),
  nnodes(nfcn-1), // How many nodes are in between
  ndof_endfcn(polyndof-1), // 1 degree for the function value
  ndof_middlefcn(polyndof-2) // +1 for slope of the other node
{
  assert((nfcn>2 && polyndof>=2) || (nfcn==2 && polyndof>=1) || (nfcn==1 && polyndof>=0));

  TIterator* coefIter = parList_.createIterator();
  RooAbsArg* func;
  while ((func = (RooAbsArg*) coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      cerr << "RooPiecewisePolynomial::RooPiecewisePolynomial(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    parList.add(*func);
  }
  delete coefIter;
}
RooPiecewisePolynomial::RooPiecewisePolynomial(RooPiecewisePolynomial const& other, const char* name) :
  RooAbsReal(other, name),
  xvar("xvar", this, other.xvar),
  parList("parList", this, other.parList),
  nfcn(other.nfcn), polyndof(other.polyndof),
  nnodes(other.nnodes),
  ndof_endfcn(other.ndof_endfcn),
  ndof_middlefcn(other.ndof_middlefcn)
{}

double RooPiecewisePolynomial::eval(double x, std::vector<double> const& par)const{
  // If we say the form of the polynomial is [0] + [1]*x + [2]*x2 + [3]*x3 + [4]*x4...,
  // use the highest two orders for matching at the nodes and free the rest.
  const double d_epsilon = 0;
  if (
    (nfcn>=2 && (int) par.size()!=2*ndof_endfcn+(nfcn-2)*ndof_middlefcn+nnodes)
    ||
    (nfcn==1 && (int) par.size()!=polyndof)
    ||
    nfcn<=0
    ) assert(0);

  if (nfcn==1){
    double res = 0;
    for (int ip=0; ip<polyndof; ip++) res += par[ip]*pow(x, ip);
    return res;
  }

  int npars_reduced[nfcn];
  for (int index=0; index<nfcn; index++){
    if (index==0 || index==nfcn-1) npars_reduced[index] = ndof_endfcn;
    else npars_reduced[index] = ndof_middlefcn;
  }

  vector<double> node(nnodes, 0); // First [0,...,nnodes-1] parameters are nodes
  vector<vector<double>> pars_full(nfcn, vector<double>(polyndof, 0)); // The full coefficients array

  for (int ip=0; ip<nnodes; ip++) node[ip] = par[ip];
  // Check for the nodes to be consecutive
  for (int ip=0; ip<nnodes; ip++){
    for (int ip2=ip+1; ip2<nnodes; ip2++){
      if (node[ip]>node[ip2]) return d_epsilon;
    }
  }
  int pos_ctr = nnodes;
  for (int index=0; index<nfcn; index++){
    for (int ipar=0; ipar<npars_reduced[index]; ipar++){
      if (!(index==(nfcn-1) && ipar==(npars_reduced[index]-1))) pars_full[index][ipar] = par[pos_ctr];
      else pars_full[index][ipar+1] = par[pos_ctr]; // Special case to avoid singular matrix. This corresponds to having the x^n contribution free instead of x^(n-1)
      pos_ctr++;
    }
  }

  vector<vector<double>> xton(nnodes, vector<double>(polyndof, 0)); // Array of node^power
  vector<vector<double>> nxtom(nnodes, vector<double>(polyndof, 0)); // Array of power*node^(power-1)
  for (int inode=0; inode<nnodes; inode++){
    for (int ipow=0; ipow<polyndof; ipow++){
      if (ipow==0) xton[inode][ipow]=1; // nxtom==0
      else if (ipow==1){
        xton[inode][ipow]=node[inode];
        nxtom[inode][ipow]=1;
      }
      else{
        xton[inode][ipow]=pow(node[inode], ipow);
        nxtom[inode][ipow]=((double) ipow)*pow(node[inode], ipow-1);
      }
    }
  }

  vector<double> ysbar_nodes(2*nnodes, 0);
  vector<vector<double>> coeff_ysbar(2*nnodes, vector<double>(2*nnodes, 0));
  int cstart=-1;
  for (int inode=0; inode<nnodes; inode++){
    int i=inode;
    int j=i+1;
    double sign_i = 1, sign_j=-1;
    for (int ip=0; ip<npars_reduced[i]; ip++){
      ysbar_nodes[inode] += sign_i*pars_full[i][ip]*xton[inode][ip];
      ysbar_nodes[nnodes+inode] += sign_i*pars_full[i][ip]*nxtom[inode][ip];
    }
    for (int ip=0; ip<npars_reduced[j]; ip++){
      if (!(j==(nfcn-1) && ip==(npars_reduced[j]-1))){
        ysbar_nodes[inode] += sign_j*pars_full[j][ip]*xton[inode][ip];
        ysbar_nodes[nnodes+inode] += sign_j*pars_full[j][ip]*nxtom[inode][ip];
      }
      else{
        ysbar_nodes[inode] += sign_j*pars_full[j][ip+1]*xton[inode][ip+1];
        ysbar_nodes[nnodes+inode] += sign_j*pars_full[j][ip+1]*nxtom[inode][ip+1];
      }
    }

    if (cstart>=0){
      coeff_ysbar[inode][cstart] = -sign_i*xton[inode][polyndof-2];
      coeff_ysbar[nnodes + inode][cstart] = -sign_i*nxtom[inode][polyndof-2];
    }
    coeff_ysbar[inode][cstart+1] = -sign_i*xton[inode][polyndof-1];
    coeff_ysbar[nnodes + inode][cstart+1] = -sign_i*nxtom[inode][polyndof-1];
    coeff_ysbar[inode][cstart+2] = -sign_j*xton[inode][polyndof-2];
    coeff_ysbar[nnodes + inode][cstart+2] = -sign_j*nxtom[inode][polyndof-2];
    if ((cstart+3)<2*nnodes){
      coeff_ysbar[inode][cstart+3] = -sign_j*xton[inode][polyndof-1];
      coeff_ysbar[nnodes + inode][cstart+3] = -sign_j*nxtom[inode][polyndof-1];
    }
    cstart+=2;
  }

  TVectorD polyvec(2*nnodes, ysbar_nodes.data());
  TMatrixD polycoeff(2*nnodes, 2*nnodes);
  for (int i=0; i<2*nnodes; i++){ for (int j=0; j<2*nnodes; j++) polycoeff(i, j)=coeff_ysbar[i][j]; }
  double testdet=0;
  TMatrixD polycoeff_inv = polycoeff.Invert(&testdet);
  if (testdet!=0){
    TVectorD unknowncoeffs = polycoeff_inv*polyvec;
    pos_ctr=0;
    for (int index=0; index<nfcn; index++){
      for (int ip=npars_reduced[index]; ip<polyndof; ip++){
        if (!(index==(nfcn-1) && ip==npars_reduced[index])) pars_full[index][ip] = unknowncoeffs[pos_ctr];
        else pars_full[index][ip-1] = unknowncoeffs[pos_ctr];
        pos_ctr++;
      }
    }

    int index_chosen=0;
    for (int index=0; index<nnodes-1; index++){
      if (x>=node[index] && x<node[index+1]){
        index_chosen = index+1;
        break;
      }
    }
    if (x>=node[nnodes-1]) index_chosen = nfcn-1;

    double res = 0;
    for (int ip=0; ip<polyndof; ip++) res += pars_full[index_chosen][ip]*pow(x, ip);
    return res;
  }
  else{
    cerr << "RooPiecewisePolynomial::eval: Something went wrong, and the determinant is 0!" << endl;
    return d_epsilon;
  }
}

double RooPiecewisePolynomial::evaluate()const{
  if (parList.getSize()==0) return 0;
  std::vector<double> par; par.reserve(parList.getSize());
  for (int ip=0; ip<parList.getSize(); ip++) par.push_back((dynamic_cast<RooAbsReal*>(parList.at(ip)))->getVal());
  return eval(xvar, par);
}
double RooPiecewisePolynomial::evaluate(double* x, double* p)const{ // For calling in a TF1 object
  // No size check is done, so be careful!
  unsigned int psize=2*ndof_endfcn+(nfcn-2)*ndof_middlefcn+nnodes;
  //cout << "RooPiecewisePolynomial::evaluate: N parameters = " << psize << endl;
  std::vector<double> par; par.reserve(psize);
  for (unsigned int ip=0; ip<psize; ip++) par.push_back(p[ip]);
  //cout << "RooPiecewisePolynomial::evaluate: Calling eval at x=" << x[0] << endl;
  return eval(x[0], par);
}


ClassImp(RooPiecewisePolynomial)
