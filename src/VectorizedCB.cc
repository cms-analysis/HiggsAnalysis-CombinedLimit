#include "HiggsAnalysis/CombinedLimit/interface/VectorizedCB.h"
#include "RooMath.h"
#include "vectorized.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"
#include <RooRealVar.h>
#include <stdexcept>

VectorizedCBShape::VectorizedCBShape(const RooCBShape &gaus, const RooAbsData &data, bool includeZeroWeights)
{
    RooArgSet obs(*data.get());
    if (obs.getSize() != 1) throw std::invalid_argument("Multi-dimensional dataset?");
    RooRealVar *x = dynamic_cast<RooRealVar*>(obs.first());

    Worker w(gaus);
    if (obs.contains(w.mvar())) {
        x_ = dynamic_cast<const RooRealVar*>(& w.mvar());
    } else {
        throw std::invalid_argument("CBShape observable is not m: if this is intended, set --X-rtd ADDNLL_CBNLL=0 to disable CBShape vectorization in NLL.");
    }
    m0_ = & w.m0var();
    n_  = & w.nvar();
    alpha_ = & w.alphavar();
    sigma_ = & w.sigmavar();

    xvals_.reserve(data.numEntries());
    for (unsigned int i = 0, n = data.numEntries(); i < n; ++i) {
        obs.assignValueOnly(*data.get(i), true);
        if (data.weight() || includeZeroWeights) xvals_.push_back(x->getVal());        
    }
    work1_.resize(xvals_.size());
    work2_.resize(xvals_.size());
    bool sorted = std::is_sorted(xvals_.begin(), xvals_.end());
    int selmode = runtimedef::get("VECTCB_MODE");
    switch (selmode) {
        case -1: // conservative setting (slower) 
            mode_ = (sorted ? Sorted : Plain);
            break;
        case 0:  // default: pre-sort if needed, use faster VDT math (but still in double precision, seems to give identical results as 1)
            mode_ = (sorted ? FastSorted : FastPreSort);
            break;
        case 1:  
            mode_ = Plain;
            break;
        case 2:
            mode_ = FastPlain;
            break;
        case 3:
            mode_ = PreSort;
            break;
        case 4:
            mode_ = FastPreSort;
            break;

    }
    if (sigma_->getVal() < 0) {
        mode_ = Plain; // anything else will be broken in these circumstances
    }
    if (hasPreSort()) {
        xindex_.resize(xvals_.size());
        for (unsigned int i = 0, n = xvals_.size(); i < n; ++i) {
            xindex_[i] = i;
        }
        std::sort(xindex_.begin(), xindex_.end(), [this](int i1, int i2) { return xvals_[i1] < xvals_[i2]; });
        xsorted_.resize(xvals_.size());
        for (unsigned int i = 0, n = xvals_.size(); i < n; ++i) {
            xsorted_[i] = xvals_[xindex_[i]];
        }
    }
}


void VectorizedCBShape::fill(std::vector<Double_t> &out) const {
    double norm = 1.0/getIntegral();

    double mean  = m0_->getVal();
    double alpha1 = alpha_->getVal();
    double n1 = n_->getVal();
    double sigma = sigma_->getVal(), invw = 1/sigma;
    double invn1 = 1.0/n1;

    out.resize(xvals_.size());
    unsigned int n = xvals_.size();

    if (hasSort()) {
        const std::vector<Double_t> & x = (hasPreSort() ? xsorted_ : xvals_);
        // t = (x-mean)*invw
        // and then check if t > -alpha1
        // ---> check if (x-mean)/sigma > -alpha1
        //            if x > -alpha1*sigma + mean
        //      check if (x-mean)/sigma > 10
        //            if x > mean + 10*sigma
        auto cut1 = std::lower_bound(x.begin(), x.end(), mean - alpha1 * sigma);
        auto cut2 = std::lower_bound(cut1,      x.end(), mean +   10   * sigma);
        unsigned int ibegin = cut1 - x.begin();
        unsigned int iend   = cut2 - x.begin();
        for (unsigned int i = 0; i < iend; ++i) {
            work1_[i] = (x[i]-mean)*invw;
        }
        if (cut1 > x.begin()) {
            cbCB(&work1_[0], ibegin, norm, &out[0], &work2_[0]);
        }
        if (cut2 > cut1) {
            cbGauss(&work1_[ibegin], iend-ibegin, norm, &out[ibegin], &work2_[ibegin]);
        }
        std::fill(&out[iend], &*out.end(), 0.);
        if (hasPreSort()) { // put back in place
            for (unsigned int i = 0; i < n; ++i) work2_[xindex_[i]] = out[i];
            std::swap(out, work2_);
        }
    } else {
        double norm2 = norm * std::exp(-0.5*alpha1*alpha1);
        double alpha1invn1 = alpha1*invn1;
        for (unsigned int i = 0; i < n; ++i) {
            work1_[i] = (xvals_[i]-mean)*invw; 
        }
        if (!hasFast()) {
            for (unsigned int i = 0; i < n; ++i) {
                work1_[i] = (xvals_[i]-mean)*invw; 
                if(work1_[i]>-alpha1){
                    out[i] = norm * std::exp(-0.5*std::pow(work1_[i],2));
                } else {
                    out[i] = norm2 * std::exp(-n1 * std::log(1. - alpha1invn1*(alpha1+work1_[i])));
                }
            }
        } else {
            for (unsigned int i = 0; i < n; ++i) {
                if(work1_[i]>-alpha1){
                    out[i] = norm * vdt::fast_exp(-0.5*work1_[i]*work1_[i]);
                } else {
                    out[i] = norm2 * vdt::fast_exp(-n1 * vdt::fast_log(1. - alpha1invn1*(alpha1+work1_[i])));
                }
            }
        }
    }
}

double VectorizedCBShape::getIntegral() const {
    double central=0;
    double left=0;
    double right=0;

    double xmin = x_->getMin();
    double xmax = x_->getMax();

    //printf("xmin = %5e, xmax = %5e\n",xmin,xmax);

    //bool isfullrange = xmin<=-RooNumber::infinity() && xmax>=RooNumber::infinity();

    static const double rootPiBy2 = std::sqrt(std::atan2(0.0,-1.0)/2.0);
    static const double invRoot2 = 1.0/std::sqrt(2);

    double invwidth = 1.0/sigma_->getVal();
    double mean = m0_->getVal();
    double tmin = (xmin-mean)*invwidth;
    double tmax = (xmax-mean)*invwidth;

    bool isfullrange = (tmin<-1000. && tmax>1000.);

    //compute gaussian contribution
    double width = sigma_->getVal();
    double alpha1 = alpha_->getVal();
    double n1 = n_->getVal();

    double central_low =std::max(xmin,mean - alpha1*width );
    double central_high=xmax;
    if(central_low < central_high) {// is the gaussian part in range?
        double tcentral_low = (central_low-mean)*invwidth;
        double tcentral_high = tmax;
        central = rootPiBy2*width*(TMath::Erf(tcentral_high*invRoot2)-TMath::Erf(tcentral_low*invRoot2));
    }

    //compute left tail;
    if (isfullrange && (n1-1.0)>1.e-5) {
        left = width*std::exp(-0.5*alpha1*alpha1)*n1*vdt::inv(alpha1*(n1-1.));
    } else {

        double left_low=xmin;
        double left_high=std::min(xmax,mean - alpha1*width);
        double thigh = (left_high-mean)*invwidth;

        if(left_low < left_high){ //is the left tail in range?
            double n1invalpha1 = n1*vdt::inv(fabs(alpha1));
            if(fabs(n1-1.0)>1.e-5) {
                double invn1m1 = vdt::inv(n1-1.);
                double leftpow = std::pow(n1invalpha1,-n1*invn1m1);
                double left0 = width*std::exp(-0.5*alpha1*alpha1)*invn1m1;
                double left1, left2;

                if (xmax>(mean-alpha1*width)) left1 = n1invalpha1;
                else left1 = std::pow( leftpow*(n1invalpha1 - alpha1 - thigh), 1.-n1);

                if (tmin<-1000.) left2 = 0.;
                else left2 = std::pow( leftpow*(n1invalpha1 - alpha1 - tmin ), 1.-n1);

                left = left0*(left1-left2);

                //left = width*std::exp(-0.5*alpha1*alpha1)*invn1m1*(n1invalpha1 - my_pow( gbrmath::fast_pow(n1invalpha1,-n1*invn1m1)*(n1invalpha1 - alpha1 - tmin), 1.-n1)) ;
                //left = width*std::exp(-0.5*alpha1*alpha1)*invn1m1*(n1invalpha1 - my_pow( gbrmath::fast_pow(n1invalpha1,-n1*invn1m1)*(n1invalpha1 - alpha1 - tmin), 1.-n1)) ;
                //left = A1*vdt::inv(-n1+1.0)*width*(gbrmath::fast_pow(B1-(left_low-mean)*invwidth,-n1+1.)-gbrmath::fast_pow(B1-(left_high-mean)*invwidth,-n1+1.));
            }
            else {
                double A1 = std::pow(n1invalpha1,n1)*std::exp(-0.5*alpha1*alpha1);
                double B1 = n1invalpha1-fabs(alpha1);   
                left = A1*width*(std::log(B1-(left_low-mean)*invwidth) - std::log(B1-(left_high-mean)*invwidth) );
            }
        }
    }

    return left + central + right;
}

void VectorizedCBShape::cbGauss(double* __restrict__ t, unsigned int n, double norm, double* __restrict__ out,  double* __restrict__ work2) const {
    for (unsigned int i = 0; i < n; ++i) {
        work2[i] = -0.5*t[i]*t[i];
    }
    if (hasFast()) {
        vdt::fast_expv(n, work2, t);
    } else {
        vdt::expv(n, work2, t);
    }
    for (unsigned int i = 0; i < n; ++i) {
        out[i] = t[i]*norm;
    }
}
void VectorizedCBShape::cbCB(double* __restrict__ t, unsigned int n, double norm, double* __restrict__ out,  double* __restrict__ work2) const {
    // val = norm * std::exp(-0.5*alpha1*alpha1) * pow(1. - alpha1/n1 * (alpha1+t), -n1);
    double alpha1 = alpha_->getVal();
    double n1 = n_->getVal(), notn1 = -n1;
    double alpha1invn1 = alpha1/n1;
    double prefactor = norm*std::exp(-0.5*std::pow(alpha1, 2));

    for (unsigned int i = 0; i < n; ++i) {
        work2[i] = 1 - alpha1invn1*(alpha1+t[i]);
    }
    if (hasFast()) {
        vdt::fast_logv(n, work2, t);
    } else {
        vdt::logv(n, work2, t);
    }
    for (unsigned int i = 0; i < n; ++i) {
         t[i] *= notn1;
    }
    if (hasFast()) {
        vdt::fast_expv(n, t, work2);
    } else {
        vdt::expv(n, t, work2);
    }
    for (unsigned int i = 0; i < n; ++i) {
        out[i] = prefactor*work2[i];
    }
}
