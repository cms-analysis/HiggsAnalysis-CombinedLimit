#include "HiggsAnalysis/CombinedLimit/interface/Accumulators.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>


namespace {
    template<typename T> void copyVector(const std::vector<T>& source, std::vector<T>& target, const unsigned int n){
        std::copy(source.begin(), source.begin()+n, target.begin());
    }

    /// need the __restrict__ to make them work 
    template<typename T> void subtract(T * __restrict__ out, unsigned int n, T  const * __restrict__ ref) {
        for (unsigned int i = 0; i < n; ++i) out[i] -= ref[i];
    }
    template<typename T> void logratio(T * __restrict__ out, unsigned int n, T  const * __restrict__ ref) {
        for (unsigned int i = 0; i < n; ++i) {
            out[i] = (out[i] > 0 && ref[i] > 0) ? std::log(out[i]/ref[i]) : 0;
        }
    }
    template<typename T> void sumdiff(T * __restrict__ sum, T * __restrict__ diff,
                 unsigned int n, 
                 const T  * __restrict__ h1, const T  * __restrict__ h2) {
        //printf("sumdiff(sum = %p, diff = %p, n = %d, h1 = %p, h2 = %p\n", (void*)sum, (void*)diff, n, (void*)h1, (void*)h2);
        for (unsigned int i = 0; i < n; ++i) {
            sum[i]  = h1[i] + h2[i];
            diff[i] = h1[i] - h2[i];
            //printf("%3d: sum = %.6f, diff = %.6f, h1 = %.6f, h2 = %.6f\n", i, sum[i], diff[i], h1[i], h2[i]);
        }
    }
    template<typename T> void meld(T * __restrict__ out, unsigned int n, T  const * __restrict__ diff, T  const * __restrict__ sum, T x, T y) {
        for (unsigned int i = 0; i < n; ++i) {
            out[i] += x*(diff[i] + y*sum[i]);
        }
    }
}


template<typename T> void FastTemplate_t<T>::Subtract(const FastTemplate_t & ref) {
    subtract(&(this->values_)[0], this->size_, &ref[0]);
}
template<typename T> void FastTemplate_t<T>::LogRatio(const FastTemplate_t & ref) {
    logratio(&(this->values_)[0], this->size_, &ref[0]);
}
template<typename T> void FastTemplate_t<T>::SumDiff(const FastTemplate_t & h1, const FastTemplate_t & h2,
                           FastTemplate_t & sum, FastTemplate_t & diff) {
    sumdiff(&sum[0], &diff[0], h1.size_, &h1[0], &h2[0]);
}

template<typename T> void FastTemplate_t<T>::Meld(const FastTemplate_t & diff, const FastTemplate_t & sum, T x, T y) {
    meld(&(this->values_)[0], this->size_, &diff[0], &sum[0], x, y);
}

template<typename T> void FastTemplate_t<T>::Log() {
    for (unsigned int i = 0; i < this->size_; ++i) {
        //if ((this->values_)[i] <= 0) printf("WARNING: log(%g) at bin %d of %d bins (%d active bins)\n", (this->values_)[i], i, int((this->values_).size()), this->size_);
        (this->values_)[i] = (this->values_)[i] > 0 ? std::log((this->values_)[i]) : T(-999);
    }
}

template<typename T> void FastTemplate_t<T>::Exp() {
    for (unsigned int i = 0; i < this->size_; ++i) {
        (this->values_)[i] = std::exp((this->values_)[i]);
    }
}

template<typename T> void FastTemplate_t<T>::CropUnderflows(T minimum, bool activebinsonly) {
    for (unsigned int i = 0, n = activebinsonly ? this->size_ : (this->values_).size(); i < n; ++i) {
        if ((this->values_)[i] < minimum) (this->values_)[i] = minimum;
    }
}

template<typename T> T FastTemplate_t<T>::Integral() const {
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < this->size_; ++i) total += (this->values_)[i];
    return total.sum();
}

template<typename T> void FastTemplate_t<T>::Scale(T factor) {
    for (unsigned int i = 0; i < this->size_; ++i) (this->values_)[i] *= factor;
}

template<typename T> void FastTemplate_t<T>::Clear() {
    for (unsigned int i = 0; i < this->size_; ++i) (this->values_)[i] = T(0.);
}

template<typename T> void FastTemplate_t<T>::CopyValues(const FastTemplate_t<T> &other) {
    std::copy(other.values_.begin(), other.values_.begin()+this->size_, (this->values_).begin());
}

template<typename T> void FastTemplate_t<T>::CopyValues(const TH1 &other) {
     for (unsigned int i = 0; i < this->size_; ++i) (this->values_)[i] = other.GetBinContent(i+1);
}

template<typename T> void FastTemplate_t<T>::CopyValues(const TH2 &other) {
    for (unsigned int i = 0, ix = 1, nx = other.GetNbinsX(), ny = other.GetNbinsY(); ix <= nx; ++ix) {
        for (unsigned int iy = 1; iy <= ny; ++iy, ++i) {
            (this->values_)[i] = other.GetBinContent(ix,iy);
            //printf("FastTemplate_t<T>::CopyValues from %s: (ix,iy) = (%d/%d,%d/%d), i = %d/%d, val = %.5f\n", other.GetName(), ix, nx, iy, ny,  i, this->size_, (this->values_)[i]);
        }
    }
}
template<typename T> void FastTemplate_t<T>::CopyValues(const TH3 &other) {
    for (unsigned int i = 0, ix = 1, nx = other.GetNbinsX(), ny = other.GetNbinsY(), nz = other.GetNbinsZ(); ix <= nx; ++ix) {
        for (unsigned int iy = 1; iy <= ny; ++iy) {
        	for (unsigned int iz = 1; iz <= nz; ++iz, ++i) {
	            (this->values_)[i] = other.GetBinContent(ix,iy, iz);
            //printf("FastTemplate_t<T>::CopyValues from %s: (ix,iy) = (%d/%d,%d/%d), i = %d/%d, val = %.5f\n", other.GetName(), ix, nx, iy, ny,  i, this->size_, (this->values_)[i]);
					}
        }
    }
}

template<typename T> void FastTemplate_t<T>::Dump() const {
    printf("--- dumping template with %d bins (%d active bins) (@%p) ---\n", int((this->values_).size()), this->size_, (void*)(&(this->values_)[0]));
    for (unsigned int i = 0; i < (this->values_).size(); ++i) printf(" bin %3d: yval = %9.5f\n", i, (this->values_)[i]);
    printf("\n"); 
}



/***************************************/
/*              FastHisto              */
/***************************************/

template<typename T, typename U> FastHisto_t<T,U>::FastHisto_t(const TH1 &hist, bool normX) :
    FastTemplate_t<T>(hist),
    axis_(*(hist.GetXaxis())),
    normX_(normX)
{}
template<typename T, typename U> FastHisto_t<T,U>::FastHisto_t(const FastHisto_t &other) :
    FastTemplate_t<T>(other),
    axis_(other.axis_),
    normX_(other.normX_)
{}

template<typename T, typename U> T FastHisto_t<T,U>::IntegralWidth(int binmin, int binmax) const {
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < std::min(GetNbinsX(), this->size_); ++i){
      if (binmin>=0 && (int)i<binmin) continue;
      if (binmax>=0 && (int)i>binmax) continue;
      total += GetBinContent(i) * GetBinWidth(i);
    }
    return total.sum();
}
template<typename T, typename U> void FastHisto_t<T,U>::Dump() const {
  printf("--- dumping histo template with %d bins (%d active) in range %.2f - %.2f (@%p)---\n", int((this->values_).size()), this->size_, axis_[0], axis_[this->size_], (void*)&(this->values_)[0]);
  for (unsigned int i = 0; i < (this->values_).size(); ++i) {
    printf(" bin %3d, x = %6.2f: yval = %9.5f, width = %6.3f\n",
      i, 0.5*(axis_[i]+axis_[i+1]), (this->values_)[i], GetBinWidth(i));
  }
  printf("\n");
}
template<typename T, typename U> T FastHisto_t<T,U>::GetMax() const {
    return * std::max((this->values_).begin(), (this->values_).end());
}



/***************************************/
/*             FastHisto2D             */
/***************************************/

template<typename T, typename U> FastHisto2D_t<T,U>::FastHisto2D_t(const TH2 &hist, bool normX, bool normY) :
    FastTemplate_t<T>(hist),
    axisX_(*(hist.GetXaxis())),
    axisY_(*(hist.GetYaxis())),
    normX_(normX),
    normY_(normY)
{}

template<typename T, typename U> FastHisto2D_t<T,U>::FastHisto2D_t(const FastHisto2D_t &other) :
    FastTemplate_t<T>(other),
    axisX_(other.axisX_),
    axisY_(other.axisY_),
    normX_(other.normX_),
    normY_(other.normY_)
{}

template<typename T, typename U> T FastHisto2D_t<T,U>::IntegralWidth(int xbinmin, int xbinmax, int ybinmin, int ybinmax) const {
    DefaultAccumulator<T> total = 0;
    const unsigned int binX_ = GetNbinsX();
    const unsigned int binY_ = GetNbinsY();
    for (unsigned int ix=0; ix<binX_; ix++){
      if (xbinmin>=0 && (int)ix<xbinmin) continue;
      if (xbinmax>=0 && (int)ix>xbinmax) continue;

      for (unsigned int iy=0; iy<binY_; iy++){
        if (ybinmin>=0 && (int)iy<ybinmin) continue;
        if (ybinmax>=0 && (int)iy>ybinmax) continue;

        total += (this->values_)[ix * binY_ +iy] * GetBinWidthX(ix) * GetBinWidthY(iy);
      }
    }
    return total.sum();
}

template<typename T, typename U> void FastHisto2D_t<T,U>::NormalizeXSlices() {
  normX_ = true;
  const unsigned int binX_ = GetNbinsX();
  const unsigned int binY_ = GetNbinsY();
  for (unsigned int ix = 0, offs = 0; ix < binX_; ++ix, offs += binY_) {
    T *values = &(this->values_)[offs];
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < binY_; ++i) total += values[i] * GetBinWidthY(i);
    if (total.sum()!=T(0)) {
      for (unsigned int i = 0; i < binY_; ++i) values[i] /= total.sum();
    }
  }
}

template<typename T, typename U> void FastHisto2D_t<T,U>::Dump() const {
  const unsigned int binX_ = GetNbinsX();
  const unsigned int binY_ = GetNbinsY();
  printf("--- dumping histo template with %d x %d bins (@%p)---\n", binX_, binY_, (void*)(&(this->values_)[0]));
  for (unsigned int ix=0; ix<binX_; ix++){
    for (unsigned int iy=0; iy<binY_; iy++){
      int i = ix * binY_ +iy;
      printf(" bin %3d, x = %6.2f, y = %6.2f: yval = %9.5f, width = %6.3f\n",
        i, 0.5*(axisX_[ix]+axisX_[ix+1]),
        0.5*(axisY_[iy]+axisY_[iy+1]),
        (this->values_)[i], GetBinWidthX(ix) * GetBinWidthY(iy));
    }
  }
  printf("\n");
}

template<typename T, typename U> T FastHisto2D_t<T,U>::GetMaxOnXY() const {
    return *std::max((this->values_).begin(), (this->values_).end());
}

template<typename T, typename U> T FastHisto2D_t<T,U>::GetMaxOnX(const U &y) const {
  const unsigned int binX_ = GetNbinsX();
  const unsigned int binY_ = GetNbinsY();
  int iy = FindBinY(y);
  T ret = 0;
  for (unsigned int i = iy; i < binX_*binY_; i += binY_) {
    if (ret < (this->values_)[i]) ret = (this->values_)[i];
  }
  return ret;
}

template<typename T, typename U> T FastHisto2D_t<T,U>::GetMaxOnY(const U &x) const {
  const unsigned int binY_ = GetNbinsY();
  int ix = FindBinX(x);
  return *std::max( &(this->values_)[ix * binY_], &(this->values_)[(ix+1) * binY_] );
}


/***************************************/
/*             FastHisto3D             */
/***************************************/

template<typename T, typename U> FastHisto3D_t<T,U>::FastHisto3D_t(const TH3 &hist, bool normX, bool normY, bool normZ) :
    FastTemplate_t<T>(hist),
    axisX_(*(hist.GetXaxis())),
    axisY_(*(hist.GetYaxis())),
    axisZ_(*(hist.GetZaxis())),
    normX_(normX),
    normY_(normY),
    normZ_(normZ)
{}

template<typename T, typename U> FastHisto3D_t<T,U>::FastHisto3D_t(const FastHisto3D_t &other) :
    FastTemplate_t<T>(other),
    axisX_(other.axisX_),
    axisY_(other.axisY_),
    axisZ_(other.axisZ_),
    normX_(other.normX_),
    normY_(other.normY_),
    normZ_(other.normZ_)
{}

template<typename T, typename U> T FastHisto3D_t<T,U>::IntegralWidth(int xbinmin, int xbinmax, int ybinmin, int ybinmax, int zbinmin, int zbinmax) const {
    DefaultAccumulator<T> total = 0;
    const unsigned int binX_ = GetNbinsX();
    const unsigned int binY_ = GetNbinsY();
    const unsigned int binZ_ = GetNbinsZ();
    for (unsigned int ix=0; ix<binX_; ix++){
      if (xbinmin>=0 && (int)ix<xbinmin) continue;
      if (xbinmax>=0 && (int)ix>xbinmax) continue;

      for (unsigned int iy=0; iy<binY_; iy++){
        if (ybinmin>=0 && (int)iy<ybinmin) continue;
        if (ybinmax>=0 && (int)iy>ybinmax) continue;

        for (unsigned int iz=0; iz<binZ_; iz++){
          if (zbinmin>=0 && (int)iz<zbinmin) continue;
          if (zbinmax>=0 && (int)iz>zbinmax) continue;

          total += (this->values_)[(ix * binY_ +iy)*binZ_ + iz] * GetBinWidthX(ix) * GetBinWidthY(iy) * GetBinWidthZ(iz);
        }
      }
    }
    return total.sum();
}

template<typename T, typename U> void FastHisto3D_t<T,U>::NormalizeXSlices() {
  normX_ = true;
  const unsigned int binX_ = GetNbinsX();
  const unsigned int binY_ = GetNbinsY();
  const unsigned int binZ_ = GetNbinsZ();
  for (unsigned int ix = 0, offs = 0; ix < binX_; ++ix, offs += (binY_*binZ_)) {
    T *values = &(this->values_)[offs];
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < binY_; ++i) {
      T widthY = GetBinWidthY(i);
      for (unsigned int j = 0; j < binZ_; ++j) {
        T widthZ = GetBinWidthZ(j);
        total += values[i*binZ_+j] * widthY * widthZ;
      }
    }
    if (total.sum() != T(0)) {
      for (unsigned int i = 0; i < binY_; ++i) {
        for (unsigned int j = 0; j < binZ_; ++j) {
          values[i*binZ_+j] /= total.sum();
        }
      }
    }
  }
}

template<typename T, typename U> void FastHisto3D_t<T,U>::Dump() const {
  const unsigned int binX_ = GetNbinsX();
  const unsigned int binY_ = GetNbinsY();
  const unsigned int binZ_ = GetNbinsZ();
  printf("--- dumping histo template with %d x %d y %d bins (@%p)---\n", binX_, binY_, binZ_, (void*)(&(this->values_)[0]));
  for (unsigned int ix=0; ix<binX_; ix++){
    for (unsigned int iy=0; iy<binY_; iy++){
      for (unsigned int iz=0; iz<binZ_; iz++){
        int i = (ix * binY_ +iy)*binZ_ + iz;
        printf(" bin %3d, x = %6.2f, y = %6.2f, z = %6.2f : zval = %9.5f, width = %6.3f\n",
          i, 0.5*(axisX_[ix]+axisX_[ix+1]),
          0.5*(axisY_[iy]+axisY_[iy+1]),
          0.5*(axisZ_[iz]+axisZ_[iz+1]),
          (this->values_)[i], GetBinWidthX(ix) * GetBinWidthY(iy) * GetBinWidthZ(iz));
      }
    }
  }
  printf("\n");
}


typedef FastHistoAxis_t<Double_t> FastHistoAxis_d;
typedef FastTemplate_t<Double_t> FastTemplate_d;
typedef FastHisto_t<Double_t> FastHisto_d;
typedef FastHisto2D_t<Double_t> FastHisto2D_d;
typedef FastHisto3D_t<Double_t> FastHisto3D_d;

typedef FastHistoAxis_t<Float_t> FastHistoAxis_f;
typedef FastTemplate_t<Float_t> FastTemplate_f;
typedef FastHisto_t<Float_t> FastHisto_f;
typedef FastHisto2D_t<Float_t> FastHisto2D_f;
typedef FastHisto3D_t<Float_t> FastHisto3D_f;

typedef FastHistoAxis_d FastHistoAxis;
// typedef FastTemplate_d FastTemplate;
// typedef FastHisto_d FastHisto;
// typedef FastHisto2D_d FastHisto2D;
// typedef FastHisto3D_d FastHisto3D;

