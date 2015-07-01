#include "../../interface/Accumulators.h"
#include <cstdio>
#include <TRandom3.h>

void testTwo(int n) {
    n /= 2;
    std::vector<double> noise;
    std::vector<double> signal;
    std::vector<double> terms;
    for (unsigned int i = 0; i < n; ++i) {
        noise.push_back(gRandom->Gaus(0,1));
    }
    for (unsigned int i = 0; i < 2*n; ++i) {
        signal.push_back(gRandom->Gaus(0.001,0.001));
    }
    for (unsigned int i = 0; i < n; ++i) {
        terms.push_back(signal[i] + noise[i]);
    }
    for (unsigned int i = 0; i < n; ++i) {
        terms.push_back(signal[i+n] - noise[i]);
    }
    double best = sumPrecise(signal);
    printf("Naive   sum: %.7g \n", sumFast(terms)-best);
    printf("Precise sum: %.7g \n", sumPrecise(terms)-best);
}

void testOne() {
    double one = 1.0, eps = 1.0;
    do {
        double two = one + eps;
        if (two == one) break;
        eps /= 2;
    } while(true);
    printf("epsilon is %.7g\n",eps);
    std::vector<double> terms;
    std::vector<double> ok;
    terms.push_back(one);
    for (unsigned int i = 0; i < 37; ++i) { 
        terms.push_back(eps); 
        ok.push_back(eps); 
        if (i % 10 == 7) terms.push_back(one);
    }
    terms.push_back(-one);
    for (unsigned int i = 0; i < 37; ++i) { 
        if (i % 10 == 7) terms.push_back(-one);
    }
    printf("Naive   sum/eps: %.7g \n", sumFast(terms)/eps);
    printf("Precise sum/eps: %.7g \n", sumPrecise(terms)/eps);
    printf("True    sum/eps: %.7g \n", sumFast(ok)/eps);

}

int main(int argc, char **argv) {
    //testOne();
    testTwo(2000);
}
