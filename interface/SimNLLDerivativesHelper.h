/*****************************************************************************
 * Authors:                                                                  *
 *  Andrea Carlo Marini,   CERN (CH),        andrea.carlo.marini@cern.ch        *
 *  Code based on roofit code.                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef SimNLLDerivHelper_H
#define SimNLLDerivHelper_H
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <exception>
#include "RooAddition.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProcessNormalization.h"

/* This class takes in input the CachingSimNLL and returns the derivatives of what it knows
 *
 */

// TODO: Make a common DerivativeAbstract class. Marke DerivativeLogNormal Inherits from this class
class DerivativeAbstract: public RooAbsReal
{
    protected: // Full access to derived classes
        const RooDataSet * data_{nullptr}; // not owned.
        // RooRealProxy?
        cacheutils::CachingAddNLL * pdf_{nullptr}; // CachingAddNLL
        RooRealProxy pdfproxy_;
    public:

        DerivativeAbstract(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data): 
            RooAbsReal(name, title),
            data_(data),
            pdf_(pdf),
            pdfproxy_( ( std::string("proxy_")+name).c_str() ,"",pdf) {};
        ~DerivativeAbstract(){};
        virtual Bool_t isDerived() const { return kTRUE; }
        const RooDataSet *data() const {return data_;}
        cacheutils::CachingAddNLL *pdf() { return pdf_; }
        bool verbose{true};

        virtual Double_t evaluate() const =0 ; // cacheutils::ReminderSum

};

// TODO: Inherit from DerivativeAbstract
class DerivativeLogNormal: public RooAbsReal
{
    const RooDataSet * data_{nullptr}; // not owned.
    // RooRealProxy?
    cacheutils::CachingAddNLL * pdf_{nullptr}; // CachingAddNLL
    RooRealProxy pdfproxy_;
    std::string thetaname_{""};

    RooRealProxy theta_;// useful for numerical derivatives -> terms

    public:
        DerivativeLogNormal(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data,const std::string& thetaname,int&found);
        ~DerivativeLogNormal(){};
        virtual Bool_t isDerived() const { return kTRUE; }
        virtual DerivativeLogNormal *clone(const char *name = 0) const ;
        const RooDataSet *data() const {return data_;}
        cacheutils::CachingAddNLL *pdf() { return pdf_; }

        virtual Double_t evaluate() const override ; // cacheutils::ReminderSum

        // name-> index association per process
        std::vector<int> kappa_pos_;

        bool verbose{true};

        //bool maskConstraints{false};

        // for debug purposes. Compute numerically the derivative corresponding to evaluate
        //Double_t numericalDerivative() const;
};

// Implementation of the derivative for rate params and rate POIs
class DerivativeRateParam : public DerivativeAbstract
{
    std::string ratename_{""};
    RooRealProxy rate_;// useful for numerical derivatives -> terms
    public:
        DerivativeRateParam(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data,const std::string& thetaname,int&found);
        ~DerivativeRateParam(){};
        virtual DerivativeRateParam *clone(const char *name = 0) const ;

        virtual Double_t evaluate() const override ; // cacheutils::ReminderSum
        // name-> index association per process in otherFactorList
        std::vector<int> rate_pos_;

};

//This class takes the nll and constructs the derivatives terms
class SimNLLDerivativesHelper
{
    private:
        //RooRealProxy?
        cacheutils::CachingSimNLL* nll_{nullptr}; //not owned, keep pointer
        RooRealVar unmask_;//("mask_derivative_constraint_","",1.);

        std::set<std::string> getServersVars(const RooAbsArg *node);
    public:
        SimNLLDerivativesHelper( cacheutils::CachingSimNLL * nll ): unmask_("mask_derivative_constraint_","",1.) {nll_=nll;}
        ~SimNLLDerivativesHelper(){
            for (auto x : derivatives_) x.second->Delete();
            derivatives_.clear();
        };
        void init();
        bool verbose{true}; // debug
        std::map<std::string,RooAbsReal*> derivatives_;

        void setMaskConstraint(int val=1){ 
            if (val ==0 or val ==1){
            unmask_.setVal( double(1-val) );
            }
            else throw std::invalid_argument("[SimNLLDerivativesHelper]::setMask called with value != 0,1");
        }
};
        
#endif

