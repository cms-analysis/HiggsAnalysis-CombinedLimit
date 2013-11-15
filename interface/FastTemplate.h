#ifndef HiggsAnalysis_CombinedLimit_FastH1
#define HiggsAnalysis_CombinedLimit_FastH1

#include <TH1.h>
#include <TH2.h>
#include <algorithm>
#include <vector>

class FastTemplate {
    public:
        typedef double T;
        typedef std::vector<T> AT;
        FastTemplate() : size_(0), values_() {}
        FastTemplate(unsigned int size) : size_(size), values_(size_) {}
        FastTemplate(const FastTemplate &other) : size_(other.size()), values_(other.values_) {}
        FastTemplate(const TH1 &other) : size_(other.GetNbinsX()), values_(size_) { CopyValues(other); }
        FastTemplate(const TH2 &other) : size_(other.GetNbinsX()*other.GetNbinsY()), values_(size_) { CopyValues(other); }
        FastTemplate & operator=(const FastTemplate &other) { 
            if (&other != this) {
                size_ = other.size_;
                values_ = other.values_;
            }
            return *this; 
        }
        FastTemplate & operator=(const TH1 &other) { 
            if (size() != unsigned(other.GetNbinsX())) { 
                size_ = other.GetNbinsX();
                values_.resize(size_);
            }
            CopyValues(other); return *this;  
        }
        ~FastTemplate() {  }
        void Resize(unsigned int newsize) {
            if (newsize != size()) {
                size_ = newsize;
                values_.resize(size_);
            }
        }
        T Integral() const ;
        void Scale(T factor) ;
        void Clear() ; 
        void CopyValues(const FastTemplate &other) ;
        void CopyValues(const TH1 &other) ;
        void CopyValues(const TH2 &other) ;
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
        void Subtract(const FastTemplate &reference);
        /// *this = log(*this)/(reference)
        void LogRatio(const FastTemplate &reference);
        /// assigns sum and diff
        static void SumDiff(const FastTemplate &h1, const FastTemplate &h2, FastTemplate &sum, FastTemplate &diff);
        /// Does this += x * (diff + (sum)*y)
        void Meld(const FastTemplate & diff, const FastTemplate & sum, T x, T y) ;
        /// protect from underflows (*this = max(*this, minimum));
        void CropUnderflows(T minimum=1e-9, bool activebinsonly=true);

        /// Tell the code that only the first N bins of the template are non-empty,
        /// and so that only those have to be considered when doing operations
        void SetActiveSize(unsigned int size) { size_ = size; }

        void Dump() const ;
    protected:
        unsigned int size_;
        AT values_;
};
class FastHisto : public FastTemplate {
    public:
        FastHisto() : FastTemplate(), binEdges_(), binWidths_() {}
        FastHisto(const TH1 &hist) ;
        FastHisto(const FastHisto &other) ;
        FastHisto & operator=(const FastHisto &other) { 
            if (size() != other.size()) {
                size_ = other.size_;
                values_    = other.values_;
                binWidths_ = other.binWidths_;
                binEdges_  = other.binEdges_;
            } else CopyValues(other); 
            return *this; 
        }
        FastHisto & operator=(const TH1 &other) { 
            if (size() != unsigned(other.GetNbinsX())) { 
                FastHisto fh(other);
                swap(fh);
            } else CopyValues(other); 
            return *this;  
        }
        ~FastHisto() {  }
        void swap(FastHisto &other) {
            std::swap(size_, other.size_);
            std::swap(values_, other.values_);
            std::swap(binWidths_, other.binWidths_);
            std::swap(binEdges_, other.binEdges_);
        }
        T GetAt(const T &x) const ;
        int FindBin(const T &x) const ;
        const T & GetBinContent(int bin) const { return values_[bin]; }
        T IntegralWidth() const ;
        void Normalize() {
            T sum = IntegralWidth();
            if (sum > 0) Scale(1.0f/sum);
        }

        void Dump() const ;
    private:
        AT binEdges_;
        AT binWidths_;
    
};
class FastHisto2D : public FastTemplate {
    public:
        FastHisto2D() : FastTemplate(), binX_(0), binY_(0), binEdgesX_(), binEdgesY_(), binWidths_() {}
        FastHisto2D(const TH2 &hist, bool normXonly=false) ;
        FastHisto2D(const FastHisto2D &other) ;
        FastHisto2D & operator=(const FastHisto2D &other) { 
            if (binX() != other.binX() || binY() != other.binY()) {
                size_      = other.size_;
                values_    = other.values_;
                binWidths_ = other.binWidths_;
                binEdgesX_ = other.binEdgesX_;
                binEdgesY_ = other.binEdgesY_;
                binX_      = other.binX_;
                binY_      = other.binY_;
            } else CopyValues(other); 
            return *this; 
        }
        ~FastHisto2D() {  }
        void swap(FastHisto2D &other) {
            std::swap(size_, other.size_);
            std::swap(binX_, other.binX_);
            std::swap(binY_, other.binY_);
            std::swap(values_, other.values_);
            std::swap(binWidths_, other.binWidths_);
            std::swap(binEdgesX_, other.binEdgesX_);
            std::swap(binEdgesY_, other.binEdgesY_);
        }
        T GetAt(const T &x, const T &y) const ;
        T IntegralWidth() const ;
        unsigned int binX() const { return binX_; }
        unsigned int binY() const { return binY_; }
        void Normalize() {
            T sum = IntegralWidth();
            if (sum > 0) Scale(1.0f/sum);
        }
        /// For each X, normalize along Y
        void NormalizeXSlices() ;

        void Dump() const ;
    private:
        unsigned int binX_, binY_;
        AT binEdgesX_;
        AT binEdgesY_;
        AT binWidths_;
    
};

#endif
