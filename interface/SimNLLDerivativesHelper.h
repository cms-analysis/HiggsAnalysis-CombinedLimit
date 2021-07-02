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
#include "RooAddition.h"
#include "RooDataSet.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProcessNormalization.h"

/* This class takes in input the CachingSimNLL and returns the derivatives of what it knows
 *
 */

class DerivativeLogNormal: public RooAbsReal
{
    const RooDataSet * data_{nullptr}; // not owned.
    const cacheutils::CachingAddNLL * pdf_{nullptr}; // CachingAddNLL
    std::string kappaname_{""};

    public:
        DerivativeLogNormal(const char *name, const char *title, const cacheutils::CachingAddNLL *pdf, const RooDataSet *data,const std::string& kappaname,int&found);
        ~DerivativeLogNormal(){};
        virtual Bool_t isDerived() const { return kTRUE; }
        virtual DerivativeLogNormal *clone(const char *name = 0) const ;
        const RooDataSet *data() const {return data_;}
        const cacheutils::CachingAddNLL *pdf() const { return pdf_; }

        virtual Double_t evaluate() const ; // cacheutils::ReminderSum

        // name-> index association per process
        std::vector<int> kappa_pos_;

        bool verbose{true};
};


class SimNLLDerivativesHelper
{
    public:
        SimNLLDerivativesHelper( cacheutils::CachingSimNLL * nll ){nll_=nll;}
        ~SimNLLDerivativesHelper(){};
        void init();
        bool verbose{true}; // debug
        std::map<std::string,RooAbsReal*> derivatives_;
        //std::map<std::string,RooAbsReal*>& derivatives(){return derivatives;}
    private:
        cacheutils::CachingSimNLL* nll_{nullptr}; //not owned, keep pointer

        std::set<std::string> getServersVars(RooAbsArg *node);
};
        
#endif

