#ifndef HiggsAnalysis_CombinedLimit_Accumulators_h
#define HiggsAnalysis_CombinedLimit_Accumulators_h

#include <vector>

class NaiveAccumulator {
    public:
        NaiveAccumulator(const double & value = 0) : sum_(value) {}
        NaiveAccumulator & operator+=(const double &inc) { sum_ += inc; return *this;  }
        NaiveAccumulator & operator-=(const double &inc) { sum_ -= inc; return *this; }
        double sum() const { return sum_; } 
    private:
        double sum_;
};

class KahanAccumulator {
    public:
        KahanAccumulator(const double & value = 0) : sum_(value), compensation_(0) {}
        void operator+=(const double &inc) { 
            double y = inc - compensation_;
            double sumnew = sum_ + y;
            double sumerr = ( sumnew - sum_ );
            compensation_ = sumerr - y;
            sum_ = sumnew; 
        }
        void operator-=(const double &inc) { this->operator+=(-inc); }
        double sum() const { return sum_; } 
    private:
        double sum_, compensation_;
};

template<typename A>
inline double sumWith(const std::vector<double> & vals) {
    A ret = 0;
    for (const double & v : vals) ret += v;
    return ret.sum();
}

typedef KahanAccumulator PreciseAccumulator;
typedef NaiveAccumulator FastAccumulator;
typedef PreciseAccumulator DefaultAccumulator;

inline double sumPrecise(const std::vector<double> & vals) {
    return sumWith<PreciseAccumulator>(vals);
}

inline double sumFast(const std::vector<double> & vals) {
    return sumWith<FastAccumulator>(vals);
}

inline double sumDefault(const std::vector<double> & vals) {
    return sumWith<DefaultAccumulator>(vals);
}

#endif
