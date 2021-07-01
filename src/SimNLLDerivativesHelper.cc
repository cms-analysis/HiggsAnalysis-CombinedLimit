#include "HiggsAnalysis/CombinedLimit/interface/SimNLLDerivativesHelper.h"

DerivativeLogNormal::DerivativeLogNormal(const char *name, const char *title, const cacheutils::CachingAddNLL *pdf, const RooDataSet *data, const std::string& kappaname) :
    RooAbsReal(name, title),
    data_(data),
    pdf_(pdf),
    kappaname_(kappaname)
{
    for (unsigned int i=0;i<pdf_->pdfs_.size();++i)
    {
        int val=-1;
        ProcessNormalization *p = dynamic_cast<ProcessNormalization*> (pdf_->coeffs_[i]);
        for (unsigned int j=0 ; j<p->logKappa_.size();++j)
        {
            if ( p->thetaList_.at(j)->GetName() == kappaname_  ) {val=int(j);}
        }
        if (val==-1) { 
            std::cout<<"[ERROR]: Unable to find kappa corresponding to "<<kappaname_<<" in pdf "<<pdf_->GetName() <<std::endl;
            throw 1; 
        }// something
        kappa_pos_.push_back(val);
    }
}

DerivativeLogNormal* DerivativeLogNormal::clone(const char *name) const {
    return new DerivativeLogNormal( 
            (name)?name:GetName(), GetTitle(), pdf(), data(), kappaname_) ;
}

Double_t DerivativeLogNormal::evaluate() const {
    /* The derivative of a kappa in a channel is
     *   sum_bin lambdat_b - data_b lambdat_b/lambda_b
     *   lambda_b -> sum of expectations in the bin
     *   lambdat_b -> sum of expecations times logK
     */

    double sum=0.0;
    // caching: TODO
    for (int ib=0;ib<data_->numEntries();++ib)
    {
        data_->get(ib) ;
        double db = data_->weight();
        double lambda=pdf_->getVal();
        double lambdat=0.0;
        for (unsigned int i=0;i<pdf_->pdfs_.size();++i) // loop over processes
        { 
            ProcessNormalization *c = dynamic_cast<ProcessNormalization*>(pdf_->coeffs_[i]);
            double logK = c -> logKappa_ [kappa_pos_[i]];
            lambdat+= c->getVal() * pdf_->pdfs_.at(i).pdf()->getVal() * logK;
        } 
        sum += lambdat - db*lambdat/lambda;
    }
    return sum;

}

// look inside and check the parameters associated to kappas 
// compute the sums
// construct RooAbsReal 
// keep track of what I can't do like that.

void SimNLLDerivativesHelper::init(){
    if (verbose) std::cout<<"Init SimNLLDerivativesHelper"<<std::endl;
    // clean
    for (auto p : derivatives_) p.second->Delete(); 
    derivatives_.clear();
    //---

    std::set<std::string> logNormal;
    // Find LogNormal candidates
    for (auto f : nll_->constrainPdfsFast_)
    {
        if (verbose) std::cout<<"inserting candidate: "<< f->getX().GetName()<<std::endl;
        logNormal.insert( f->getX().GetName() );
    }

    // loop over sim components
    unsigned idx=0;
    for (std::vector<cacheutils::CachingAddNLL*>::const_iterator it = nll_->pdfs_.begin(), ed = nll_->pdfs_.end(); it != ed; ++it, ++idx) {
        cacheutils::CachingAddNLL * pdf = *it;
        if (pdf==nullptr) continue ; // ?!? needed?
        const RooAbsData* data= pdf->data();
        bool isWeighted = data->isWeighted(); // binned vs unbinned?!? TBC

        for(unsigned proc =0 ; proc < pdf ->pdfs_.size();++proc)
        {
            RooAbsReal * coeff = pdf -> coeffs_[proc] ;
            ProcessNormalization * pn = dynamic_cast<ProcessNormalization*> (coeff);
            if (pn == nullptr or not isWeighted) {  // remove all the lognormal candidates. Don't know how to deal with them
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing from the list of known derivatives because"
                        << ((pn==nullptr)? std::string(" coeff is not processNormalization but")+ std::string( coeff->ClassName() ): "")
                        << ((pn==nullptr and not isWeighted) ? " and":"" )
                        << ((not isWeighted)?" dataset is not weighted (unbinned?)":"")
                        << "the following variables"
                        <<std::endl;
                }
                RooArgList l = RooArgList(*coeff->getVariables());
                for (int idx=0; idx<l.getSize();++idx) { 
                    RooAbsArg*v = l.at(idx);
                    if (logNormal.find(v->GetName()) != logNormal.end()) {
                        logNormal.erase(v->GetName());
                        if (verbose) std::cout<<"|"<<v->GetName();
                    }
                }
                //--
                l = RooArgList(*pdf->getVariables());
                for (int idx=0; idx<l.getSize();++idx) { 
                    RooAbsArg*v = l.at(idx);
                    if (logNormal.find(v->GetName()) != logNormal.end()) {
                        logNormal.erase(v->GetName());
                        if (verbose) std::cout<<"|"<<v->GetName();
                    }
                }
                continue;
            } //  not a process normalization coefficient

            // remove all asymmThetaList
            //for (auto v : pn->asymmThetaList_)
            for(int idx=0; idx< pn->asymmThetaList_.getSize() ;++idx){
                RooAbsArg*v=pn->asymmThetaList_.at(idx);
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()<<" from the list of known derivatives because appears in an asymmThetaList "
                        <<std::endl;
                }
                if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
            }
            // remove all otherFactorList
            //for (auto v : pn->otherFactorList_)
            for(int idx=0; idx< pn->otherFactorList_.getSize() ;++idx){
                RooAbsArg*v=pn->otherFactorList_.at(idx); 
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()<<" from the list of known derivatives because appears in an otherFactors "
                        <<std::endl;
                }
                if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
            }

        }

        // cross check derivatives with logNormal
        //for (auto d : derivatives){
        //    if (logNormal.find(d.first) == logNormal.end()){
        //        //remove derivative
        //        derivatives.erase(d.first());
        //    }
        //}

    }

    for (auto name : logNormal){
        RooArgList list;
        // loop over sim components
        unsigned idx=0;
        for (std::vector<cacheutils::CachingAddNLL*>::const_iterator it = nll_->pdfs_.begin(), ed = nll_->pdfs_.end(); it != ed; ++it, ++idx) {
            // construct the derivative as sum
            cacheutils::CachingAddNLL * pdf = *it;
            if (pdf==nullptr) continue ; // ?!? needed?
            const RooAbsData* data= pdf->data();
            //(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data,std::string kappaname)
            DerivativeLogNormal der("","",pdf,dynamic_cast<const RooDataSet*>(data),name);
            list.add(der); // or addOwned?
        }
        // add derivative of constraint
        // nll_->constrainPdfsFast_
        for(unsigned i=0;i<nll_->constrainPdfsFast_.size();++i)
        {
            if (nll_->constrainPdfsFast_[i]->getX().GetName() != name) continue;
            RooFormulaVar constr("constr","(@0-@1)/(@2*@2)",RooArgList(nll_->constrainPdfsFast_[i]->getX(),nll_->constrainPdfsFast_[i]->getMean(), nll_->constrainPdfsFast_[i]->getSigma()));// (x-mean)/sigma^2
            list.add(constr) ; // or addOwned
        }
        derivatives_[name] = new RooAddition(("derivative_"+name).c_str(),"",list);
    } 


}

