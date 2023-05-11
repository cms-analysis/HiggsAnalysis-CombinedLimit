#ifndef ROO_EFTSCALINGFUNCTION_V2
#define ROO_EFTSCALINGFUNCTION_V2
#include <RooAbsReal.h>
#include <RooListProxy.h>
#include <TString.h>
#include <TObjString.h>

#include <vector>
#include <string>
#include <map>


class RooEFTScalingFunction_v2 : public RooAbsReal {
    public:
        RooEFTScalingFunction_v2() {}
        RooEFTScalingFunction_v2(const char *name, const char *title, const std::map<std::string,double> &coeffs, const RooArgList &terms);
        RooEFTScalingFunction_v2(const RooEFTScalingFunction_v2& other, const char* name=0);
        virtual ~RooEFTScalingFunction_v2() {}
        virtual TObject *clone(const char *newname) const { return new RooEFTScalingFunction_v2(*this,newname); } 
        const std::map<std::string,double> & coeffs() const { return coeffs_; }
        const RooArgList & terms() const { return terms_; }
    protected:
        std::map<std::string,double> coeffs_;
        RooListProxy terms_;
        std::map< std::vector<RooAbsReal *>, double> vcomponents_;
        double offset_;
        virtual Double_t evaluate() const ;
    private:
        ClassDef(RooEFTScalingFunction_v2,1)
};

#endif
