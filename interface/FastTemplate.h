#ifndef HiggsAnalysis_CombinedLimit_FastH1
#define HiggsAnalysis_CombinedLimit_FastH1

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <algorithm>
#include <vector>

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

        void Dump() const ;

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
        FastTemplate_t<T> & operator=(const TH1 &other) {
          if (size() != unsigned(other.GetNbinsX())) {
            size_ = other.GetNbinsX();
            values_.resize(size_);
          }
          this->CopyValues(other); return *this;
        }
        ~FastTemplate_t() {}
        void Resize(unsigned int newsize) {
          if (newsize != size()) {
            size_ = newsize;
            values_.resize(size_);
          }
        }
};
template <typename T> class FastHisto_t : public FastTemplate_t<T> {
private:
        std::vector<T> binEdges_;
        bool normX_;

public:
        void swap(FastHisto_t<T> &other) {
            std::swap(this->size_, other.size_);
            std::swap(this->values_, other.values_);
            std::swap(binEdges_, other.binEdges_);
            std::swap(normX_, other.normX_);
        }
        T GetAt(const T &x) const ;
        int FindBin(const T &x) const ;
        const T & GetBinContent(const unsigned int bin) const { return (this->values_)[bin]; }
        T IntegralWidth(int binmin=-1, int binmax=-1) const;
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum > T(0)) this->Scale(1.0f/sum);
        }

        void Dump() const ;

        T GetMax() const ;

        const T & GetEdge(unsigned int i) const { return binEdges_[i]; }
        const T & GetWidth(unsigned int i) const { return binWidths_[i]; }
        T GetBinWidth(const unsigned int bin) const;
        T GetXmin() const { return (this->size_>0 ? binEdges_[0] : T(0)); }
        T GetXmax() const { return (this->size_>0 ? binEdges_[this->size_-1] : T(0)); }

        FastHisto_t() : FastTemplate_t<T>(), binEdges_(), normX_(false) {}
        FastHisto_t(const TH1 &hist, bool normX=false);
        FastHisto_t(const FastHisto_t<T> &other);
        FastHisto_t & operator=(const FastHisto_t<T> &other) {
          normX_ = other.normX_;
          if (this->size() != other.size()) {
            this->size_ = other.size_;
            this->values_    = other.values_;
            binEdges_  = other.binEdges_;
          }
          else this->CopyValues(other);
          return *this;
        }
        FastHisto_t<T> & operator=(const TH1 &other) {
          if (this->size() != unsigned(other.GetNbinsX())) {
            FastHisto_t<T> fh(other);
            this->swap(fh);
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto_t<T>(){}
};
template <typename T> class FastHisto2D_t : public FastTemplate_t<T> {
private:
        unsigned int binX_, binY_;
        std::vector<T> binEdgesX_;
        std::vector<T> binEdgesY_;
        bool normX_;
        bool normY_;

public:
        void swap(FastHisto2D_t<T> &other) {
            std::swap(this->size_, other.size_);
            std::swap(binX_, other.binX_);
            std::swap(binY_, other.binY_);
            std::swap(this->values_, other.values_);
            std::swap(binEdgesX_, other.binEdgesX_);
            std::swap(binEdgesY_, other.binEdgesY_);
            std::swap(normX_, other.normX_);
            std::swap(normY_, other.normY_);
        }
        T GetAt(const T &x, const T &y) const ;
        T IntegralWidth(int xbinmin=-1, int xbinmax=-1, int ybinmin=-1, int ybinmax=-1) const;
        unsigned int binX() const { return binX_; }
        unsigned int binY() const { return binY_; }
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum > T(0)) this->Scale(1.0f/sum);
        }
        /// For each X, normalize along Y
        void NormalizeXSlices() ;

        void Dump() const ;

        T GetMaxOnXY() const ;
        T GetMaxOnX(const T &y) const ;
        T GetMaxOnY(const T &x) const ;

        int FindBinX(const T &t) const;
        int FindBinY(const T &t) const;

        T GetBinWidthX(const unsigned int bin) const;
        T GetBinWidthY(const unsigned int bin) const;

        T GetXmin() const { return (binX_>0 ? binEdgesX_[0] : T(0)); }
        T GetYmin() const { return (binY_>0 ? binEdgesY_[0] : T(0)); }
        T GetXmax() const { return (binX_>0 ? binEdgesX_[binX_-1] : T(0)); }
        T GetYmax() const { return (binY_>0 ? binEdgesY_[binY_-1] : T(0)); }

        FastHisto2D_t() : FastTemplate_t<T>(), binX_(0), binY_(0), binEdgesX_(), binEdgesY_(), normX_(false), normY_(false) {}
        FastHisto2D_t(const TH2 &hist, bool normX=false, bool normY_=false);
        FastHisto2D_t(const FastHisto2D_t<T> &other);
        FastHisto2D_t & operator=(const FastHisto2D_t<T> &other) {
          normX_ = other.normX_;
          normY_ = other.normY_;
          if (binX() != other.binX() || binY() != other.binY()) {
            this->size_      = other.size_;
            this->values_    = other.values_;
            binEdgesX_ = other.binEdgesX_;
            binEdgesY_ = other.binEdgesY_;
            binX_      = other.binX_;
            binY_      = other.binY_;
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto2D_t(){}
};

template <typename T> class FastHisto3D_t : public FastTemplate_t<T> {
private:
        unsigned int binX_, binY_, binZ_;
        std::vector<T> binEdgesX_;
        std::vector<T> binEdgesY_;
        std::vector<T> binEdgesZ_;
        bool normX_;
        bool normY_;
        bool normZ_;

public:
        void swap(FastHisto3D_t<T> &other) {
            std::swap(this->size_, other.size_);
            std::swap(binX_, other.binX_);
            std::swap(binY_, other.binY_);
            std::swap(binZ_, other.binZ_);
            std::swap(this->values_, other.values_);
            std::swap(binEdgesX_, other.binEdgesX_);
            std::swap(binEdgesY_, other.binEdgesY_);
            std::swap(binEdgesZ_, other.binEdgesZ_);
            std::swap(normX_, other.normX_);
            std::swap(normY_, other.normY_);
            std::swap(normZ_, other.normZ_);
        }
        T GetAt(const T &x, const T &y, const T &z) const ;
        T IntegralWidth(int xbinmin=-1, int xbinmax=-1, int ybinmin=-1, int ybinmax=-1, int zbinmin=-1, int zbinmax=-1) const;
        unsigned int binX() const { return binX_; }
        unsigned int binY() const { return binY_; }
        unsigned int binZ() const { return binZ_; }
        void Normalize() {
            T sum = this->IntegralWidth();
            if (sum > T(0)) this->Scale(1.0f/sum);
        }
        /// For each X, normalize along Y
        void NormalizeXSlices() ;

        void Dump() const ;

        int FindBinX(const T &t) const;
        int FindBinY(const T &t) const;
        int FindBinZ(const T &t) const;

        T GetBinWidthX(const unsigned int bin) const;
        T GetBinWidthY(const unsigned int bin) const;
        T GetBinWidthZ(const unsigned int bin) const;
        T GetXmin() const { return (binX_>0 ? binEdgesX_[0] : T(0)); }
        T GetYmin() const { return (binY_>0 ? binEdgesY_[0] : T(0)); }
        T GetZmin() const { return (binZ_>0 ? binEdgesZ_[0] : T(0)); }
        T GetXmax() const { return (binX_>0 ? binEdgesX_[binX_-1] : T(0)); }
        T GetYmax() const { return (binY_>0 ? binEdgesY_[binY_-1] : T(0)); }
        T GetZmax() const { return (binZ_>0 ? binEdgesZ_[binZ_-1] : T(0)); }

        FastHisto3D_t() : FastTemplate_t<T>(), binX_(0), binY_(0), binZ_(0), binEdgesX_(), binEdgesY_(), binEdgesZ_(), normX_(false), normY_(false), normZ_(false) {}
        FastHisto3D_t(const TH3 &hist, bool normX=false, bool normY=false, bool normZ=false);
        FastHisto3D_t(const FastHisto3D_t<T> &other);
        FastHisto3D_t & operator=(const FastHisto3D_t<T> &other) {
          normX_ = other.normX_;
          normY_ = other.normY_;
          normZ_ = other.normZ_;
          if (binX() != other.binX() || binY() != other.binY() || binZ() != other.binZ()) {
            this->size_      = other.size_;
            this->values_    = other.values_;
            binEdgesX_ = other.binEdgesX_;
            binEdgesY_ = other.binEdgesY_;
            binEdgesZ_ = other.binEdgesZ_;
            binX_      = other.binX_;
            binY_      = other.binY_;
            binZ_      = other.binZ_;
          }
          else this->CopyValues(other);
          return *this;
        }
        ~FastHisto3D_t(){}
};

#include "FastTemplate.hpp"
#endif
