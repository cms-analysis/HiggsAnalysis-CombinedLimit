#include "../interface/FastTemplate.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

FastTemplate::T FastTemplate::Integral() const {
    T total = 0;
    for (unsigned int i = 0; i < size_; ++i) total += values_[i];
    return total;
}

void FastTemplate::Scale(T factor) {
    for (unsigned int i = 0; i < size_; ++i) values_[i] *= factor;
}

void FastTemplate::Clear() {
    for (unsigned int i = 0; i < size_; ++i) values_[i] = T(0.);
}

void FastTemplate::CopyValues(const FastTemplate &other) {
    std::copy(other.values_.begin(), other.values_.begin()+size_, values_.begin());
}

void FastTemplate::CopyValues(const TH1 &other) {
     for (unsigned int i = 0; i < size_; ++i) values_[i] = other.GetBinContent(i+1);
}

void FastTemplate::CopyValues(const TH2 &other) {
    for (unsigned int i = 0, ix = 1, nx = other.GetNbinsX(), ny = other.GetNbinsY(); ix <= nx; ++ix) {
        for (unsigned int iy = 1; iy <= ny; ++iy, ++i) {
            values_[i] = other.GetBinContent(ix,iy);
            //printf("FastTemplate::CopyValues from %s: (ix,iy) = (%d/%d,%d/%d), i = %d/%d, val = %.5f\n", other.GetName(), ix, nx, iy, ny,  i, size_, values_[i]);
        }
    }
}

void FastTemplate::Dump() const {
    printf("--- dumping template with %d bins (%d active bins) (@%p) ---\n", int(values_.size()), size_, (void*)(&values_[0]));
    for (unsigned int i = 0; i < values_.size(); ++i) printf(" bin %3d: yval = %9.5f\n", i, values_[i]);
    printf("\n"); 
}

FastHisto::FastHisto(const TH1 &hist) :
    FastTemplate(hist),
    binEdges_(size()+1),
    binWidths_(size())
{
    for (unsigned int i = 0, n = size(); i < n; ++i) {
        binEdges_[i] = hist.GetBinLowEdge(i+1);
        binWidths_[i] = hist.GetBinWidth(i+1);
    }
    binEdges_.back() = hist.GetBinLowEdge(size()+1);
}

FastHisto::FastHisto(const FastHisto &other) :
    FastTemplate(other),
    binEdges_(other.binEdges_),
    binWidths_(other.binWidths_)
{
}

int FastHisto::FindBin(const T &x) const {
    auto match = std::lower_bound(binEdges_.begin(), binEdges_.end(), x);
    if (match == binEdges_.begin()) return -1;
    if (match == binEdges_.end()) return values_.size();
    return match - binEdges_.begin() - 1;
}


FastHisto::T FastHisto::GetAt(const T &x) const {
    auto match = std::lower_bound(binEdges_.begin(), binEdges_.end(), x);
    if (match == binEdges_.begin() || match == binEdges_.end()) return T(0.0);
    return values_[match - binEdges_.begin() - 1];
}

FastHisto::T FastHisto::IntegralWidth() const {
    double total = 0;
    for (unsigned int i = 0; i < size_; ++i) total += values_[i] * binWidths_[i];
    return total;
}

void FastHisto::Dump() const {
    printf("--- dumping histo template with %d bins (%d active) in range %.2f - %.2f (@%p)---\n", int(values_.size()), size_, binEdges_[0], binEdges_[size_], (void*)&values_[0]);
    for (unsigned int i = 0; i < values_.size(); ++i) {
        printf(" bin %3d, x = %6.2f: yval = %9.5f, width = %6.3f\n", 
                    i, 0.5*(binEdges_[i]+binEdges_[i+1]), values_[i], binWidths_[i]);
    }
    printf("\n"); 
}

FastHisto2D::FastHisto2D(const TH2 &hist, bool normXonly) :
    FastTemplate(hist),
    binX_(hist.GetNbinsX()),
    binY_(hist.GetNbinsY()),
    binEdgesX_(hist.GetNbinsX()+1),
    binEdgesY_(hist.GetNbinsY()+1),
    binWidths_(size_)
{
    TAxis *ax = hist.GetXaxis(), *ay = hist.GetYaxis();
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

FastHisto2D::FastHisto2D(const FastHisto2D &other) :
    FastTemplate(other),
    binX_(other.binX_),
    binY_(other.binY_),
    binEdgesX_(other.binEdgesX_),
    binEdgesY_(other.binEdgesY_),
    binWidths_(other.binWidths_)
{
}

FastHisto2D::T FastHisto2D::GetAt(const T &x, const T &y) const {
    auto matchx = std::lower_bound(binEdgesX_.begin(), binEdgesX_.end(), x);
    if (matchx == binEdgesX_.begin() || matchx == binEdgesX_.end()) return T(0.0);
    int ix = (matchx - binEdgesX_.begin() - 1);
    auto matchy = std::lower_bound(binEdgesY_.begin(), binEdgesY_.end(), y);
    if (matchy == binEdgesY_.begin() || matchy == binEdgesY_.end()) return T(0.0);
    int iy = (matchy - binEdgesY_.begin() - 1);
    return values_[ix * binY_ + iy];
}

FastHisto2D::T FastHisto2D::IntegralWidth() const {
    double total = 0;
    for (unsigned int i = 0; i < size_; ++i) total += values_[i] * binWidths_[i];
    return total;
}

void FastHisto2D::NormalizeXSlices() {
    for (unsigned int ix = 0, offs = 0; ix < binX_; ++ix, offs += binY_) {
       T *values = & values_[offs], *widths = & binWidths_[offs];
       double total = 0;
       for (unsigned int i = 0; i < binY_; ++i) total += values[i] * widths[i];
       if (total > 0) {
            total = T(1.0)/total;
            for (unsigned int i = 0; i < binY_; ++i) values[i] *= total;
       } 
    }
}

void FastHisto2D::Dump() const {
    printf("--- dumping histo template with %d x %d bins (@%p)---\n", binX_, binY_, (void*)(&values_[0]));
    for (unsigned int i = 0; i < size_; ++i) {
        printf(" bin %3d, x = %6.2f, y = %6.2f: yval = %9.5f, width = %6.3f\n", 
                    i, 0.5*(binEdgesX_[i/binY_]+binEdgesX_[i/binY_+1]), 
                       0.5*(binEdgesY_[i%binY_]+binEdgesY_[(i%binY_)+1]),
                     values_[i], binWidths_[i]);
    }
    printf("\n"); 
}


namespace { 
    /// need the __restrict__ to make them work 
    void subtract(FastTemplate::T * __restrict__ out, unsigned int n, FastTemplate::T  const * __restrict__ ref) {
        for (unsigned int i = 0; i < n; ++i) out[i] -= ref[i];
    }
    void logratio(FastTemplate::T * __restrict__ out, unsigned int n, FastTemplate::T  const * __restrict__ ref) {
        for (unsigned int i = 0; i < n; ++i) {
            out[i] = (out[i] > 0 && ref[i] > 0) ? std::log(out[i]/ref[i]) : 0;
        }
    }
    void sumdiff(FastTemplate::T * __restrict__ sum, FastTemplate::T * __restrict__ diff, 
                 unsigned int n, 
                 const FastTemplate::T  * __restrict__ h1, const FastTemplate::T  * __restrict__ h2) {
        //printf("sumdiff(sum = %p, diff = %p, n = %d, h1 = %p, h2 = %p\n", (void*)sum, (void*)diff, n, (void*)h1, (void*)h2);
        for (unsigned int i = 0; i < n; ++i) {
            sum[i]  = h1[i] + h2[i];
            diff[i] = h1[i] - h2[i];
            //printf("%3d: sum = %.6f, diff = %.6f, h1 = %.6f, h2 = %.6f\n", i, sum[i], diff[i], h1[i], h2[i]);
        }
    }
    void meld(FastTemplate::T * __restrict__ out, unsigned int n, FastTemplate::T  const * __restrict__ diff, FastTemplate::T  const * __restrict__ sum, FastTemplate::T x, FastTemplate::T y) {
        for (unsigned int i = 0; i < n; ++i) {
            out[i] += x*(diff[i] + y*sum[i]);
        }
    }
}

void FastTemplate::Subtract(const FastTemplate & ref) {
    subtract(&values_[0], size_, &ref[0]);
}
void FastTemplate::LogRatio(const FastTemplate & ref) {
    logratio(&values_[0], size_, &ref[0]);
}
void FastTemplate::SumDiff(const FastTemplate & h1, const FastTemplate & h2, 
                           FastTemplate & sum, FastTemplate & diff) {
    sumdiff(&sum[0], &diff[0], h1.size_, &h1[0], &h2[0]);
}

void FastTemplate::Meld(const FastTemplate & diff, const FastTemplate & sum, T x, T y) {
    meld(&values_[0], size_, &diff[0], &sum[0], x, y);
}

void FastTemplate::Log() {
    for (unsigned int i = 0; i < size_; ++i) {
        //if (values_[i] <= 0) printf("WARNING: log(%g) at bin %d of %d bins (%d active bins)\n", values_[i], i, int(values_.size()), size_);
        values_[i] = values_[i] > 0 ? std::log(values_[i]) : T(-999);
    }
}

void FastTemplate::Exp() {
    for (unsigned int i = 0; i < size_; ++i) {
        values_[i] = std::exp(values_[i]);
    }
}

void FastTemplate::CropUnderflows(T minimum, bool activebinsonly) {
    for (unsigned int i = 0, n = activebinsonly ? size_ : values_.size(); i < n; ++i) {
        if (values_[i] < minimum) values_[i] = minimum;
    }
}   
