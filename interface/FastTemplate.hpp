#include "HiggsAnalysis/CombinedLimit/interface/Accumulators.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

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

template<typename T> FastHisto_t<T>::FastHisto_t(const TH1 &hist) :
    FastTemplate_t<T>(hist),
    binEdges_(this->size()+1),
    binWidths_(this->size())
{
    for (unsigned int i = 0, n = this->size(); i < n; ++i) {
        binEdges_[i] = hist.GetBinLowEdge(i+1);
        binWidths_[i] = hist.GetBinWidth(i+1);
    }
    binEdges_.back() = hist.GetBinLowEdge(this->size()+1);
}

template<typename T> FastHisto_t<T>::FastHisto_t(const FastHisto_t &other) :
    FastTemplate_t<T>(other),
    binEdges_(other.binEdges_),
    binWidths_(other.binWidths_)
{}

template<typename T> int FastHisto_t<T>::FindBin(const T &x) const {
    auto match = std::lower_bound(binEdges_.begin(), binEdges_.end(), x);
    if (match == binEdges_.begin()) return -1;
    if (match == binEdges_.end()) return (this->values_).size();
    return match - binEdges_.begin() - 1;
}


template<typename T> T FastHisto_t<T>::GetAt(const T &x) const {
    auto match = std::lower_bound(binEdges_.begin(), binEdges_.end(), x);
    if (match == binEdges_.begin() || match == binEdges_.end()) return T(0.0);
    return (this->values_)[match - binEdges_.begin() - 1];
}

template<typename T> T FastHisto_t<T>::IntegralWidth() const {
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < this->size_; ++i) total += (this->values_)[i] * binWidths_[i];
    return total.sum();
}

template<typename T> void FastHisto_t<T>::Dump() const {
    printf("--- dumping histo template with %d bins (%d active) in range %.2f - %.2f (@%p)---\n", int((this->values_).size()), this->size_, binEdges_[0], binEdges_[this->size_], (void*)&(this->values_)[0]);
    for (unsigned int i = 0; i < (this->values_).size(); ++i) {
        printf(" bin %3d, x = %6.2f: yval = %9.5f, width = %6.3f\n", 
                    i, 0.5*(binEdges_[i]+binEdges_[i+1]), (this->values_)[i], binWidths_[i]);
    }
    printf("\n"); 
}

template<typename T> T FastHisto_t<T>::GetMax() const {
    return * std::max((this->values_).begin(), (this->values_).end());
}


template<typename T> FastHisto2D_t<T>::FastHisto2D_t(const TH2 &hist, bool normXonly) :
    FastTemplate_t<T>(hist),
    binX_(hist.GetNbinsX()),
    binY_(hist.GetNbinsY()),
    binEdgesX_(hist.GetNbinsX()+1),
    binEdgesY_(hist.GetNbinsY()+1),
    binWidths_(this->size_)
{
    const TAxis *ax = hist.GetXaxis(), *ay = hist.GetYaxis();
    for (unsigned int ix = 0; ix < binX_; ++ix) {
        binEdgesX_[ix] = ax->GetBinLowEdge(ix+1);
    }
    binEdgesX_.back() = ax->GetBinLowEdge(binX_+1);
    for (unsigned int iy = 0; iy < binY_; ++iy) {
        binEdgesY_[iy] = ay->GetBinLowEdge(iy+1);
    }
    binEdgesY_.back() = ay->GetBinLowEdge(binY_+1);
    for (unsigned int ix = 1, i = 0; ix <= binX_; ++ix) {
        for (unsigned int iy = 1; iy <= binY_; ++iy, ++i) {
            binWidths_[i] = (normXonly ? 1 : ax->GetBinWidth(ix))*ay->GetBinWidth(iy);
        }
    }
}

template<typename T> FastHisto2D_t<T>::FastHisto2D_t(const FastHisto2D_t &other) :
    FastTemplate_t<T>(other),
    binX_(other.binX_),
    binY_(other.binY_),
    binEdgesX_(other.binEdgesX_),
    binEdgesY_(other.binEdgesY_),
    binWidths_(other.binWidths_)
{}

template<typename T> T FastHisto2D_t<T>::GetAt(const T &x, const T &y) const {
    auto matchx = std::lower_bound(binEdgesX_.begin(), binEdgesX_.end(), x);
    if (matchx == binEdgesX_.begin() || matchx == binEdgesX_.end()) return T(0.0);
    int ix = (matchx - binEdgesX_.begin() - 1);
    auto matchy = std::lower_bound(binEdgesY_.begin(), binEdgesY_.end(), y);
    if (matchy == binEdgesY_.begin() || matchy == binEdgesY_.end()) return T(0.0);
    int iy = (matchy - binEdgesY_.begin() - 1);
    return (this->values_)[ix * binY_ + iy];
}

template<typename T> T FastHisto2D_t<T>::IntegralWidth() const {
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < this->size_; ++i) total += (this->values_)[i] * binWidths_[i];
    return total.sum();
}

template<typename T> void FastHisto2D_t<T>::NormalizeXSlices() {
    for (unsigned int ix = 0, offs = 0; ix < binX_; ++ix, offs += binY_) {
       T *values = & (this->values_)[offs], *widths = & binWidths_[offs];
       DefaultAccumulator<T> totalc = 0;
       for (unsigned int i = 0; i < binY_; ++i) totalc += values[i] * widths[i];
       double total = totalc.sum();
       if (total > 0) {
            total = T(1.0)/total;
            for (unsigned int i = 0; i < binY_; ++i) values[i] *= total;
       } 
    }
}

template<typename T> void FastHisto2D_t<T>::Dump() const {
    printf("--- dumping histo template with %d x %d bins (@%p)---\n", binX_, binY_, (void*)(&(this->values_)[0]));
    for (unsigned int i = 0; i < this->size_; ++i) {
        printf(" bin %3d, x = %6.2f, y = %6.2f: yval = %9.5f, width = %6.3f\n", 
                    i, 0.5*(binEdgesX_[i/binY_]+binEdgesX_[i/binY_+1]), 
                       0.5*(binEdgesY_[i%binY_]+binEdgesY_[(i%binY_)+1]),
                     (this->values_)[i], binWidths_[i]);
    }
    printf("\n"); 
}

template<typename T> T FastHisto2D_t<T>::GetMaxOnXY() const {
    return *std::max((this->values_).begin(), (this->values_).end());
}

template<typename T> T FastHisto2D_t<T>::GetMaxOnX(const T &y) const {
    auto matchy = std::lower_bound(binEdgesY_.begin(), binEdgesY_.end(), y);
    if (matchy == binEdgesY_.begin() || matchy == binEdgesY_.end()) return T(0.0);
    int iy = (matchy - binEdgesY_.begin() - 1);
    T ret = 0.0;
    for (unsigned int i = iy; i < this->size_; i += binY_) {
        if (ret < (this->values_)[i]) ret = (this->values_)[i];
    }
    return ret;
}

template<typename T> T FastHisto2D_t<T>::GetMaxOnY(const T &x) const {
    auto matchx = std::lower_bound(binEdgesX_.begin(), binEdgesX_.end(), x);
    if (matchx == binEdgesX_.begin() || matchx == binEdgesX_.end()) return T(0.0);
    int ix = (matchx - binEdgesX_.begin() - 1);
    return *std::max( &(this->values_)[ix * binY_], &(this->values_)[(ix+1) * binY_] );
}


template<typename T> FastHisto3D_t<T>::FastHisto3D_t(const TH3 &hist, bool normXonly) :
    FastTemplate_t<T>(hist),
    binX_(hist.GetNbinsX()),
    binY_(hist.GetNbinsY()),
    binZ_(hist.GetNbinsZ()),
    binEdgesX_(hist.GetNbinsX()+1),
    binEdgesY_(hist.GetNbinsY()+1),
    binEdgesZ_(hist.GetNbinsZ()+1),
    binWidths_(this->size_)
{
    const TAxis *ax = hist.GetXaxis(), *ay = hist.GetYaxis(), *az = hist.GetZaxis();
    for (unsigned int ix = 0; ix < binX_; ++ix) {
        binEdgesX_[ix] = ax->GetBinLowEdge(ix+1);
    }
    binEdgesX_.back() = ax->GetBinLowEdge(binX_+1);
    for (unsigned int iy = 0; iy < binY_; ++iy) {
        binEdgesY_[iy] = ay->GetBinLowEdge(iy+1);
    }
    binEdgesY_.back() = ay->GetBinLowEdge(binY_+1);
    for (unsigned int iz = 0; iz < binZ_; ++iz) {
        binEdgesZ_[iz] = az->GetBinLowEdge(iz+1);
    }
    binEdgesZ_.back() = az->GetBinLowEdge(binZ_+1);
    for (unsigned int ix = 1, i = 0; ix <= binX_; ++ix) {
        for (unsigned int iy = 1; iy <= binY_; ++iy ) {
        	for (unsigned int iz = 1; iz <= binZ_; ++iz, ++i) {
              binWidths_[i] = (normXonly ? 1 : ax->GetBinWidth(ix))*ay->GetBinWidth(iy)*az->GetBinWidth(iz);
          }
        }
    }
}

template<typename T> FastHisto3D_t<T>::FastHisto3D_t(const FastHisto3D_t &other) :
    FastTemplate_t<T>(other),
    binX_(other.binX_),
    binY_(other.binY_),
    binZ_(other.binZ_),
    binEdgesX_(other.binEdgesX_),
    binEdgesY_(other.binEdgesY_),
    binEdgesZ_(other.binEdgesZ_),
    binWidths_(other.binWidths_)
{}

template<typename T> T FastHisto3D_t<T>::GetAt(const T &x, const T &y, const T &z) const {
    auto matchx = std::lower_bound(binEdgesX_.begin(), binEdgesX_.end(), x);
    if (matchx == binEdgesX_.begin() || matchx == binEdgesX_.end()) return T(0.0);
    int ix = (matchx - binEdgesX_.begin() - 1);
    auto matchy = std::lower_bound(binEdgesY_.begin(), binEdgesY_.end(), y);
    if (matchy == binEdgesY_.begin() || matchy == binEdgesY_.end()) return T(0.0);
    int iy = (matchy - binEdgesY_.begin() - 1);
    auto matchz = std::lower_bound(binEdgesZ_.begin(), binEdgesZ_.end(), z);
    if (matchz == binEdgesZ_.begin() || matchz == binEdgesZ_.end()) return T(0.0);
    int iz = (matchz - binEdgesZ_.begin() - 1);
    return (this->values_)[ix * binY_ *binZ_ +binZ_*iy + iz];
}

template<typename T> T FastHisto3D_t<T>::IntegralWidth() const {
    DefaultAccumulator<T> total = 0;
    for (unsigned int i = 0; i < this->size_; ++i) total += (this->values_)[i] * binWidths_[i];
    return total.sum();
}

template<typename T> void FastHisto3D_t<T>::NormalizeXSlices() {
    for (unsigned int ix = 0, offs = 0; ix < binX_; ++ix, offs += (binY_*binZ_)) {
       T *values = & (this->values_)[offs], *widths = & binWidths_[offs];
       double total = 0;
       for (unsigned int i = 0; i < binY_; ++i) {
       	for (unsigned int j = 0; j < binZ_; ++j) {
					total += values[i*binZ_+j] * widths[i*binZ_+j];
				}
			 }
       if (total > 0) {
            total = T(1.0)/total;
            for (unsigned int i = 0; i < binY_; ++i) {
            	for (unsigned int j = 0; j < binZ_; ++j) {
								values[i*binZ_+j] *= total;
							}
						}
       } 
    }
}

template<typename T> void FastHisto3D_t<T>::Dump() const {
    printf("--- dumping histo template with %d x %d y %d bins (@%p)---\n", binX_, binY_, binZ_, (void*)(&(this->values_)[0]));
    for (unsigned int i = 0; i < this->size_; ++i) {
        printf(" bin %3d, x = %6.2f, y = %6.2f, z = %6.2f : zval = %9.5f, width = %6.3f\n", 
                    i, 0.5*(binEdgesX_[i/(binY_*binZ_)]+binEdgesX_[i/(binY_*binZ_)+1]), 
                       0.5*(binEdgesY_[i%(binY_*binZ_)/binZ_]+binEdgesY_[(i%(binY_*binZ_)/binZ_)+1]),
                       0.5*(binEdgesZ_[i%binZ_]+binEdgesZ_[(i%binZ_)+1]),
                     (this->values_)[i], binWidths_[i]);
    }
    printf("\n"); 
}


namespace { 
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

typedef FastTemplate_t<Double_t> FastTemplate_d;
typedef FastHisto_t<Double_t> FastHisto_d;
typedef FastHisto2D_t<Double_t> FastHisto2D_d;
typedef FastHisto3D_t<Double_t> FastHisto3D_d;

typedef FastTemplate_t<Float_t> FastTemplate_f;
typedef FastHisto_t<Float_t> FastHisto_f;
typedef FastHisto2D_t<Float_t> FastHisto2D_f;
typedef FastHisto3D_t<Float_t> FastHisto3D_f;

typedef FastTemplate_d FastTemplate;
typedef FastHisto_d FastHisto;
typedef FastHisto2D_d FastHisto2D;
typedef FastHisto3D_d FastHisto3D;

