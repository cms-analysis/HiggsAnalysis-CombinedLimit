#ifndef HiggsAnalysis_CombinedLimit_FastH1
#define HiggsAnalysis_CombinedLimit_FastH1

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <algorithm>
#include <vector>


template <typename U> class FastHistoAxis_t {
protected:
  std::vector<U> binEdges_;

public:
  FastHistoAxis_t() : binEdges_(){}
  FastHistoAxis_t(unsigned int size_) : binEdges_(size_, U(0)){}
  FastHistoAxis_t(const TAxis& axis){
    int nbins = axis.GetNbins();
    for (int ix=0; ix<=nbins; ++ix) binEdges_.push_back(U(axis.GetBinLowEdge(ix+1)));
  }
  FastHistoAxis_t(const FastHistoAxis_t<U>& other) : binEdges_(other.binEdges_){}
  FastHistoAxis_t(const std::vector<U>& other) : binEdges_(other){}
  virtual inline ~FastHistoAxis_t(){}

  unsigned int size() const { return binEdges_.size(); }
  unsigned int GetNbins() const{ int s=size(); return (unsigned int)std::max(s-1, 0); }

  void resize(unsigned int newsize){ if (newsize != size()) binEdges_.resize(newsize); }
  void swap(const FastHistoAxis_t<U>& other){ std::swap(binEdges_, other.binEdges_); }

  U& operator[](unsigned int i) { return binEdges_.at(i); }
  const U& operator[](unsigned int i) const { return binEdges_.at(i); }
  FastHistoAxis_t<U>& operator=(const FastHistoAxis_t<U>& other){ binEdges_ = other.binEdges_; return *this; }
  FastHistoAxis_t<U>& operator=(const TAxis& axis){
    FastHistoAxis_t<U> other(axis);
    swap(other);
    return *this;
  }

  int FindBin(const U& x) const{
    if (x==binEdges_.at(size()-1) && size()>1) return (int)(size()-2);
    auto bbegin = binEdges_.begin();
    auto bend = binEdges_.end();
    auto match = std::upper_bound(bbegin, bend, x);
    return match - bbegin - 1;
  }
  U GetBinWidth(const unsigned int bin) const{
    if (size()==0) return U(0);
    else if (bin>=0 && bin<size()-1) return (binEdges_[bin+1]-binEdges_[bin]);
    else return U(1);
  }
  U GetBinLowEdge(const int bin=-1) const {
    if (size()==0) return U(0);
    else if (bin>=0 && bin<(int)size()) return binEdges_[bin];
    else if (bin>=(int)size()) return binEdges_[size()-1];
    else return binEdges_[0];
  }
  U GetBinUpEdge(const int bin=-1) const {
    if (size()==0) return U(0);
    else if (bin>=0 && bin<(int)(size()-1)) return binEdges_[bin+1];
    else if (bin>=(int)(size()-1)) return binEdges_[size()-1];
    else return binEdges_[0];
  }
};


template <typename T> class FastTemplate_t {
protected:
        unsigned int size_;
        std::vector<T> values_;

public:
        T Integral() const ;
        void Scale(T factor) ;
        void Clear() ; 
        void CopyValues(const FastTemplate_t<T> &other);
        void CopyValues(const TH1 &other) ;
        void CopyValues(const TH2 &other) ;
        void CopyValues(const TH3 &other) ;
        T & operator[](unsigned int i) { return values_[i]; }
        const T & operator[](unsigned int i) const { return values_[i]; }
        /// return the full size of the template
        const unsigned int fullsize() const { return values_.size(); }
        /// return the active size of the template (can be less than the full size if the SetActiveSize
        /// has been used to inform the code that only the first N bins are not empty)
        const unsigned int size() const { return size_; }
        
        /// *this = log(*this) 
        void Log();
        /// *this = exp(*this) 
        void Exp();
        /// *this = *this - reference
        void Subtract(const FastTemplate_t<T> &reference);
        /// *this = log(*this)/(reference)
        void LogRatio(const FastTemplate_t<T> &reference);
        /// assigns sum and diff
        static void SumDiff(const FastTemplate_t<T> &h1, const FastTemplate_t<T> &h2, FastTemplate_t<T> &sum, FastTemplate_t<T> &diff);
        /// Does this += x * (diff + (sum)*y)
        void Meld(const FastTemplate_t<T> & diff, const FastTemplate_t<T> & sum, T x, T y);
        /// protect from underflows (*this = max(*this, minimum));
        void CropUnderflows(T minimum=1e-9, bool activebinsonly=true);

        /// Tell the code that only the first N bins of the template are non-empty,
        /// and so that only those have to be considered when doing operations
        void SetActiveSize(unsigned int size) { size_ = size; }

        virtual void Dump() const ;

        FastTemplate_t() : size_(0), values_() {}
        FastTemplate_t(unsigned int size) : size_(size), values_(size_) {}
        FastTemplate_t(const FastTemplate_t<T> &other) : size_(other.size()), values_(other.values_) {}
        FastTemplate_t(const TH1 &other) : size_(other.GetNbinsX()), values_(size_) { CopyValues(other); }
        FastTemplate_t(const TH2 &other) : size_(other.GetNbinsX()*other.GetNbinsY()), values_(size_) { CopyValues(other); }
        FastTemplate_t(const TH3 &other) : size_(other.GetNbinsX()*other.GetNbinsY()*other.GetNbinsZ()), values_(size_) { CopyValues(other); }
        FastTemplate_t<T> & operator=(const FastTemplate_t<T> &other) {
          if (&other != this) {
            size_ = other.size_;
            values_ = other.values_;
          }
          return *this;
        }
        virtual FastTemplate_t<T> & operator=(const TH1 &other) {
          if ((int)size() != other.GetNbinsX()) {
            size_ = (unsigned int)other.GetNbinsX();
            values_.resize(size_);
          }
          this->CopyValues(other);
          return *this;
        }
        virtual inline ~FastTemplate_t() {}
        void Resize(unsigned int newsize) {
          if (newsize != size()) {
            size_ = newsize;
            values_.resize(size_);
          }
        }
};
template <typename T, typename U=Double_t> class FastHisto_t : public FastTemplate_t<T> {
private:
        FastHistoAxis_t<U> axis_;
        bool normX_;

public:
        void swap(FastHisto_t<T,U> &other) {
            std::swap(this->size_, other.size_);
            std::swap(this->values_, other.values_);
            std::swap(axis_, other.axis_);
            std::swap(normX_, other.normX_);
        }

        int FindBin(const U &x) const { return axis_.FindBin(x); }
        unsigned int GetNbinsX() const { return axis_.GetNbins(); }
        U GetBinWidth(const unsigned int bin) const { if (normX_) return U(1); else return axis_.GetBinWidth(bin); }
        U GetXmin(const int bin=-1) const { return axis_.GetBinLowEdge(bin); }
        U GetXmax(const int bin=-1) const { return axis_.GetBinUpEdge(bin); }
        const T& GetBinContent(const unsigned int bin) const { return (this->values_).at(bin); }
        T& GetBinContent(const unsigned int bin) { return (this->values_).at(bin); }
        T GetAt(const U &x) const {
          int bin = FindBin(x);
          if (bin < 0 || bin >= int((this->values_).size())) {
            return T(0);
          } else {
            return GetBinContent((unsigned int)bin);
          }
        }

        //this should really be least significant of T and U
        //typically that would be T
        T IntegralWidth(int binmin=-1, int binmax=-1) const ;
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum!=T(0)) this->Scale(T(1)/sum);
        }
        T GetMax() const ;

        U GetEdge(unsigned int i) const { return GetXmin(i); }
        U GetWidth(unsigned int i) const { return GetBinWidth(i); }

        void Dump() const ;

        FastHisto_t() : FastTemplate_t<T>(), axis_(), normX_(false) {}
        FastHisto_t(const TH1 &hist, bool normX=false);
        FastHisto_t(const FastHisto_t<T,U> &other);
        FastHisto_t<T,U>& operator=(const FastHisto_t<T,U> &other) {
          normX_ = other.normX_;
          if (this->size() != other.size()) {
            this->size_ = other.size_;
            this->values_ = other.values_;
            axis_ = other.axis_;
          }
          else this->CopyValues(other);
          return *this;
        }
        FastHisto_t<T,U>& operator=(const TH1 &other) {
          if ((int)this->size() != other.GetNbinsX()) {
            FastHisto_t<T,U> fh(other);
            swap(fh);
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto_t<T,U>(){}
};
template <typename T, typename U=Double_t> class FastHisto2D_t : public FastTemplate_t<T> {
private:
        FastHistoAxis_t<U> axisX_;
        FastHistoAxis_t<U> axisY_;
        bool normX_;
        bool normY_;

public:
        void swap(FastHisto2D_t<T,U> &other) {
            std::swap(this->size_, other.size_);
            std::swap(this->values_, other.values_);
            std::swap(axisX_, other.axisX_);
            std::swap(axisY_, other.axisY_);
            std::swap(normX_, other.normX_);
            std::swap(normY_, other.normY_);
        }

        int FindBinX(const U &x) const { return axisX_.FindBin(x); }
        unsigned int GetNbinsX() const { return axisX_.GetNbins(); }
        U GetBinWidthX(const unsigned int bin) const { if (normX_) return U(1); else return axisX_.GetBinWidth(bin); }
        U GetXmin(const int bin=-1) const { return axisX_.GetBinLowEdge(bin); }
        U GetXmax(const int bin=-1) const { return axisX_.GetBinUpEdge(bin); }

        int FindBinY(const U &y) const { return axisY_.FindBin(y); }
        unsigned int GetNbinsY() const { return axisY_.GetNbins(); }
        U GetBinWidthY(const unsigned int bin) const { if (normY_) return U(1); else return axisY_.GetBinWidth(bin); }
        U GetYmin(const int bin=-1) const { return axisY_.GetBinLowEdge(bin); }
        U GetYmax(const int bin=-1) const { return axisY_.GetBinUpEdge(bin); }

        const T& GetBinContent(const unsigned int ix, const unsigned int iy) const {
          const unsigned int bin = ix * GetNbinsY() + iy;
          return (this->values_).at(bin);
        }
        T& GetBinContent(const unsigned int ix, const unsigned int iy){
          const unsigned int bin = ix * GetNbinsY() + iy;
          return (this->values_).at(bin);
        }
        T GetAt(const U &x, const U &y) const {
          int xbin = FindBinX(x);
          int ybin = FindBinY(y);
          if (xbin<0 || ybin<0) return T(0);
          else return GetBinContent((unsigned int)xbin, (unsigned int)ybin);
        }

        T IntegralWidth(int xbinmin=-1, int xbinmax=-1, int ybinmin=-1, int ybinmax=-1) const;
        unsigned int binX() const { return GetNbinsX(); }
        unsigned int binY() const { return GetNbinsY(); }
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum!=T(0)) this->Scale(T(1)/sum);
        }

        // For each X, normalize along Y
        void NormalizeXSlices() ;

        void Dump() const ;

        T GetMaxOnXY() const ;
        T GetMaxOnX(const U &y) const ;
        T GetMaxOnY(const U &x) const ;

        FastHisto2D_t() : FastTemplate_t<T>(), axisX_(), axisY_(), normX_(false), normY_(false) {}
        FastHisto2D_t(const TH2 &hist, bool normX=false, bool normY_=false);
        FastHisto2D_t(const FastHisto2D_t<T,U> &other);
        FastHisto2D_t<T,U> & operator=(const FastHisto2D_t<T,U> &other) {
          normX_ = other.normX_;
          normY_ = other.normY_;
          if (GetNbinsX() != other.GetNbinsX() || GetNbinsY() != other.GetNbinsY()) {
            this->size_      = other.size_;
            this->values_    = other.values_;
            axisX_ = other.axisX_;
            axisY_ = other.axisY_;
          }
          else this->CopyValues(other);
          return *this;
        }
        FastHisto2D_t<T,U>& operator=(const TH2 &other) {
          if (GetNbinsX() != other.GetNbinsX() || GetNbinsY() != other.GetNbinsY()) {
            FastHisto2D_t<T,U> fh(other);
            swap(fh);
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto2D_t(){}
};

template <typename T, typename U=Double_t> class FastHisto3D_t : public FastTemplate_t<T> {
private:
        FastHistoAxis_t<U> axisX_;
        FastHistoAxis_t<U> axisY_;
        FastHistoAxis_t<U> axisZ_;
        bool normX_;
        bool normY_;
        bool normZ_;

public:
        void swap(FastHisto3D_t<T,U> &other) {
            std::swap(this->size_, other.size_);
            std::swap(this->values_, other.values_);
            std::swap(axisX_, other.axisX_);
            std::swap(axisY_, other.axisY_);
            std::swap(axisZ_, other.axisZ_);
            std::swap(normX_, other.normX_);
            std::swap(normY_, other.normY_);
            std::swap(normZ_, other.normZ_);
        }

        int FindBinX(const U &x) const { return axisX_.FindBin(x); }
        unsigned int GetNbinsX() const { return axisX_.GetNbins(); }
        U GetBinWidthX(const unsigned int bin) const { if (normX_) return U(1); else return axisX_.GetBinWidth(bin); }
        U GetXmin(const int bin=-1) const { return axisX_.GetBinLowEdge(bin); }
        U GetXmax(const int bin=-1) const { return axisX_.GetBinUpEdge(bin); }

        int FindBinY(const U &y) const { return axisY_.FindBin(y); }
        unsigned int GetNbinsY() const { return axisY_.GetNbins(); }
        U GetBinWidthY(const unsigned int bin) const { if (normY_) return U(1); else return axisY_.GetBinWidth(bin); }
        U GetYmin(const int bin=-1) const { return axisY_.GetBinLowEdge(bin); }
        U GetYmax(const int bin=-1) const { return axisY_.GetBinUpEdge(bin); }

        int FindBinZ(const U &y) const { return axisZ_.FindBin(y); }
        unsigned int GetNbinsZ() const { return axisZ_.GetNbins(); }
        U GetBinWidthZ(const unsigned int bin) const { if (normZ_) return U(1); else return axisZ_.GetBinWidth(bin); }
        U GetZmin(const int bin=-1) const { return axisZ_.GetBinLowEdge(bin); }
        U GetZmax(const int bin=-1) const { return axisZ_.GetBinUpEdge(bin); }

        const T& GetBinContent(const unsigned int ix, const unsigned int iy, const unsigned int iz) const {
          const unsigned int bin = (ix * GetNbinsY() +iy)*GetNbinsZ() + iz;
          return (this->values_).at(bin);
        }
        T& GetBinContent(const unsigned int ix, const unsigned int iy, const unsigned int iz){
          const unsigned int bin = (ix * GetNbinsY() +iy)*GetNbinsZ() + iz;
          return (this->values_).at(bin);
        }
        T GetAt(const U &x, const U &y, const U &z) const {
          int xbin = FindBinX(x);
          int ybin = FindBinY(y);
          int zbin = FindBinZ(z);
          if (xbin<0 || ybin<0 || zbin<0) return T(0);
          else return GetBinContent((unsigned int)xbin, (unsigned int)ybin, (unsigned int)zbin);
        }

        T IntegralWidth(int xbinmin=-1, int xbinmax=-1, int ybinmin=-1, int ybinmax=-1, int zbinmin=-1, int zbinmax=-1) const;
        unsigned int binX() const { return GetNbinsX(); }
        unsigned int binY() const { return GetNbinsY(); }
        unsigned int binZ() const { return GetNbinsZ(); }
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum!=T(0)) this->Scale(U(1)/sum);
        }
        // For each X, normalize along Y
        void NormalizeXSlices() ;

        void Dump() const ;

        FastHisto3D_t() : FastTemplate_t<T>(), axisX_(), axisY_(), axisZ_(), normX_(false), normY_(false), normZ_(false) {}
        FastHisto3D_t(const TH3 &hist, bool normX=false, bool normY=false, bool normZ=false);
        FastHisto3D_t(const FastHisto3D_t<T,U> &other);
        FastHisto3D_t<T,U>& operator=(const FastHisto3D_t<T,U> &other) {
          normX_ = other.normX_;
          normY_ = other.normY_;
          normZ_ = other.normZ_;
          if (GetNbinsX() != other.GetNbinsX() || GetNbinsY() != other.GetNbinsY() || GetNbinsZ() != other.GetNbinsZ()) {
            this->size_      = other.size_;
            this->values_    = other.values_;
            axisX_ = other.axisX_;
            axisY_ = other.axisY_;
            axisZ_ = other.axisZ_;
          }
          else this->CopyValues(other);
          return *this;
        }
        FastHisto3D_t<T,U>& operator=(const TH3 &other) {
          if (GetNbinsX() != other.GetNbinsX() || GetNbinsY() != other.GetNbinsY() || GetNbinsZ() != other.GetNbinsZ()) {
            FastHisto3D_t<T,U> fh(other);
            swap(fh);
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto3D_t(){}
};

#include "FastTemplate.hpp"
#endif
