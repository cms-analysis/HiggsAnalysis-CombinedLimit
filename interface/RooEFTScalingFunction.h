#ifndef ROO_EFTSCALINGFUNCTION
#define ROO_EFTSCALINGFUNCTION
#include <RooAbsReal.h>
#include <RooListProxy.h>
#include <TString.h>
#include <TObjString.h>

#include <vector>
#include <string>
#include <map>


class RooEFTScalingFunction : public RooAbsReal {
    public:
        RooEFTScalingFunction() {}
        RooEFTScalingFunction(const char *name, const char *title, const std::map<std::string,double> &coeffs, const RooArgList &terms);
        RooEFTScalingFunction(const RooEFTScalingFunction& other, const char* name=0);
        ~RooEFTScalingFunction() override {}
        TObject *clone(const char *newname) const override { return new RooEFTScalingFunction(*this,newname); } 
        const std::map<std::string,double> & coeffs() const { return coeffs_; }
        const RooArgList & terms() const { return terms_; }
    protected:
        std::map<std::string,double> coeffs_;
        RooListProxy terms_;
        std::map< std::vector<RooAbsReal *>, double> vcomponents_;
        double offset_;
        Double_t evaluate() const override ;
    private:
        ClassDefOverride(RooEFTScalingFunction,1)
};

#endif
