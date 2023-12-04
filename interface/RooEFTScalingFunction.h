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
        virtual ~RooEFTScalingFunction() {}
        virtual TObject *clone(const char *newname) const { return new RooEFTScalingFunction(*this,newname); } 
        const std::map<std::string,double> & coeffs() const { return coeffs_; }
        const RooArgList & terms() const { return terms_; }
    protected:
        std::map<std::string,double> coeffs_;
        RooListProxy terms_;
        std::map< std::vector<RooAbsReal *>, double> vcomponents_;
        double offset_;
        virtual Double_t evaluate() const ;
    private:
        ClassDef(RooEFTScalingFunction,1)
};

#endif
