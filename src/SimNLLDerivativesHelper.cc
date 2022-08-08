#include "HiggsAnalysis/CombinedLimit/interface/SimNLLDerivativesHelper.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"

#include "HiggsAnalysis/CombinedLimit/interface/CMSHistSum.h"
//#define DERIVATIVE_RATEPARAM_DEBUG 1
//#define DERIVATIVE_LOGNORMAL_DEBUG 1
#define DEBUG_CMSHISTSUM 1

DerivativeLogNormal::DerivativeLogNormal(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data, const std::string& thetaname,int&found) :
    RooAbsReal(name, title),
    data_(data),
    pdf_(pdf),
    pdfproxy_( ( std::string("proxy_")+name).c_str() ,"",pdf),
    thetaname_(thetaname)
{
    //setDirtyInhibit(true); // TODO: understand cache and proxies
    setOperMode(RooAbsArg::ADirty,false); // Recursive, bool

    found=0; // keep track if at least one depends on the given kappa
    for (unsigned int i=0;i<pdf_->pdfs_.size();++i) // process loop
    {
        if(verbose) std::cout<<"[DerivativeLogNormal]: loop "<< i <<"/"<<pdf_->pdfs_.size()<<" considering pdf "<<pdf_->GetName() <<std::endl;
        int val=-1;
        ProcessNormalization *p = dynamic_cast<ProcessNormalization*> (pdf_->coeffs_[i]);
        if (verbose) std::cout<<"[DerivativeLogNormal]: coeff is of class "<<pdf_->coeffs_[i]->ClassName()<<std::endl;
        // RooProduct
        RooProduct *pp = dynamic_cast<RooProduct*>(pdf_->coeffs_[i]);
        if (p == nullptr and pp != nullptr){
            for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
            {
                if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) p = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
            }
        }
        //
        if (p != nullptr) { // unable to add the sum for this process

            for (unsigned int j=0 ; j<p->logKappa_.size();++j)
            {
                if ( p->thetaList_.at(j)->GetName() == thetaname_  ) {val=int(j);}
            }
            if (val==-1) { 
                if(verbose) std::cout<<"[DerivativeLogNormal]: Unable to find kappa corresponding to "<<thetaname_<<" in pdf "<<pdf_->GetName() <<std::endl;
            }// something
            else{
                found=1;
            }
        }
        else{
            if(verbose) std::cout<<"[DerivativeLogNormal]: Unable to find kappa corresponding to "<<thetaname_<<" in pdf "<<pdf_->GetName() << "because not ProcessNormalization" <<std::endl;
        }

        kappa_pos_.push_back(val);
    }
    if(verbose) std::cout<<"[DerivativeLogNormal]: constructed derivative for "<<thetaname_<<" of pdf "<<pdf_->GetName() <<std::endl;
}

DerivativeLogNormal* DerivativeLogNormal::clone(const char *name) const {
    if(verbose) std::cout<<"[DerivativeLogNormal]: clone"<<std::endl;
    int found;
    return new DerivativeLogNormal( 
            (name)?name:GetName(), GetTitle(), pdf_, data(), thetaname_, found) ;
}

Double_t DerivativeLogNormal::evaluate() const {
#ifdef DERIVATIVE_LOGNORMAL_DEBUG
    if(verbose) std::cout<<"[DerivativeLogNormal]: evaluate"<<std::endl;
#endif
    /* The derivative of a kappa in a channel is
     *   sum_bin lambdat_b - data_b lambdat_b/lambda_b
     *   lambda_b -> sum of expectations in the bin
     *   lambdat_b -> sum of ( expecations times logK)
     */

    double sum=0.0;
    // caching: TODO
#ifdef DERIVATIVE_LOGNORMAL_DEBUG
    if(verbose) {
        std::cout<<"[DerivativeLogNormal]: data"<<std::endl;
        data_->Print("V");
        std::cout<<"[DerivativeLogNormal]: --"<<std::endl;
    }
#endif
    
    //const RooDataHist* dh = dynamic_cast<const RooDataHist*> (data_);
    //if (dh==nullptr) std::cout<<"[DerivativeLogNormal]:"<<" NO ROODATAHIST"<<std::endl;

    // FIX for ADDNLL_ROOREALSUM_KEEPZEROS 
    std::vector<double> reminder(pdf_->pdfs_.size(),0.0); // the sums must do 1, the ib==ndata does that with data=0
    int ndata=data_->numEntries();

    for (int ib=0;ib<=ndata;++ib) 
    //for (unsigned int ib=0;ib<pdf_->weights_.size();++ib)
    {
        auto x= (ib<ndata) ? data_->get(ib): data_->get(0) ;
        double db = (ib<ndata) ? data_->weight() : 0.; 
        //double bw = (dh != nullptr) ? dh->binVolume(): 1.0;
        double bw = 1; //(pdf_->binWidths_.size() > 1) ? pdf_->binWidths_[ib]: (pdf_->binWidths_.size()==1) ?pdf_->binWidths_[0]: 1.0;
        //        const RooArgSet *obs = data_->get();
        //                RooRealVar *xvar = dynamic_cast<RooRealVar *>(obs->first());
        if (x->getSize()==1) {
            RooRealVar *xvar=dynamic_cast<RooRealVar*>(x->first());
            const RooAbsBinning &bins = xvar->getBinning();
            bw=bins.binWidth(0); //  only costant supported. Need to figure out for zero bins
            //if (verbose){
            //    std::cout <<"BinWidth: ";
            //    for(int ii =0 ;ii<xvar->numBins();++ii){
            //        double bc2 = bins.binCenter(ii);
            //        std::cout<< bins.binWidth(ii);
            //    }
            //    std::cout <<std::endl;
            //}
        }

        //double db = pdf_->weights_[ib];
        //double bw = pdf_->binWidths_[ib];

        //double lambda=pdf_->getVal()*bw; // sum coeff*pdf
        double lambdat=0.0;
        double lambda=0.;// debug: recompute lambda with loop. not sure that the term before include everything I want

        //bool recursive=true;
        //double running_prod=1.0;
        for (unsigned int i=0;i<pdf_->pdfs_.size();++i) // loop over processes
        { 
            // find kappa
            ProcessNormalization *c = dynamic_cast<ProcessNormalization*>(pdf_->coeffs_[i]);
            
            // fix RooProduct
            RooProduct *pp = dynamic_cast<RooProduct*>(pdf_->coeffs_[i]);
            if (c == nullptr and pp != nullptr){
                for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
                {
                    if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) c = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                }
            }
            //double integral = (pdf_->isRooRealSum_ && pdf_->basicIntegrals_ < 2) ? pdf_->integrals_[i]->getVal(x): 1.0;
            //---
            double logK = (kappa_pos_[i]>=0)? c -> logKappa_ [kappa_pos_[i]] : 0.0;
            if (pdf_->isRooRealSum_ and verbose)std::cout<<"[DerivativeLogNormal]:DEBUG "<<" pdf is a RooRealSumPdf"<<std::endl;

            //pdf_->pdfs_.at(i).pdf()->setValueDirty();

            // TH1
            if (ib<ndata) reminder[i] += pdf_->pdfs_.at(i).pdf()->getVal(x)*bw;
            double expectedEvents= (ib<ndata)  ? pdf_->coeffs_[i]->getVal(x) * pdf_->pdfs_.at(i).pdf()->getVal(x)*bw : 
                        pdf_->coeffs_[i]->getVal(x) * (1.-reminder[i]);

            // sum coeff or sum coeff*integral
            //double expectedEvents= dynamic_cast<const RooAbsPdf*>(pdf_->pdfs_.at(i).pdf())->expectedEvents(x);
            //bool expEventsNoNorm=false;
            //double expectedEvents = (pdf_->isRooRealSum_ && !expEventsNoNorm ? dynamic_cast<const RooAbsPdf*>(pdf_->pdfs_.at(i).pdf())->getNorm(data_->get()) : pdf_->coeffs_[i]->getVal());
            lambdat+= expectedEvents * logK;
            lambda += expectedEvents;

            if (pdf_->basicIntegrals_ and verbose)std::cout<<"[DerivativeLogNormal]:DEBUG "<<" PDF has basic integrals:"<<pdf_->basicIntegrals_<<std::endl;

#ifdef DERIVATIVE_LOGNORMAL_DEBUG
            if(verbose /*and DERIVATIVE_LOGNORMAL_DEBUG > 1*/) std::cout<<"[DerivativeLogNormal]: *** ibin="<< ib<< " data="<<db <<" iproc="<< i<<" expEvents="<<expectedEvents<< " logK="<<logK<<" coeff="<<pdf_->coeffs_[i]->getVal(x) <<" pdf="<<pdf_->pdfs_.at(i).pdf()->getVal(x) <<" pdf2=<<"<<dynamic_cast<const RooAbsPdf*>(pdf_->pdfs_.at(i).pdf())->getValV(x)  << "pdfU="<<dynamic_cast<const RooAbsPdf*>(pdf_->pdfs_.at(i).pdf())->getValV(nullptr) <<"bw="<<bw  <<" norm="<<dynamic_cast<const RooAbsPdf*>(pdf_->pdfs_.at(i).pdf())->getNorm(x) <<" reminder(i)="<<reminder[i] << std::endl;
            //x->Print("V");
#endif
        } 
        sum += lambdat * (1. - db/lambda);

#ifdef DERIVATIVE_LOGNORMAL_DEBUG
        if(verbose) std::cout<<"[DerivativeLogNormal]: ibin="<< ib<< " data="<<db << " lambda="<<lambda << "( from pdf_" << pdf_->getVal()*bw <<")" <<" lambdat="<<lambdat<<" bw="<<bw<< " running sum="<<sum<<std::endl;
#endif
    }

#ifdef DERIVATIVE_LOGNORMAL_DEBUG
    if(verbose) std::cout<<"[DerivativeLogNormal]: thetaname="<<thetaname_<<" partial="<<sum<<std::endl;
#endif
    return sum;

}

// look inside and check the parameters associated to kappas 
// compute the sums
// construct RooAbsReal 
// keep track of what I can't do like that.

void SimNLLDerivativesHelper::init(){ // TODO. Will need to decompose this into list of logical operations. Split it in two: first find the parameters we know how to deal with, and second construct the derivatives
    if (verbose) std::cout<<"Init SimNLLDerivativesHelper"<<std::endl;
    // clean
    for (auto p : derivatives_) p.second->Delete(); 
    derivatives_.clear();
    //---

    std::set<std::string> logNormal;
    std::set<std::string> rateParams; // notice that logNormal cannot be in rateParams and vice-versa. Mixing could in principle be possible, but needs to be handled in logNormal to understand constraint part.

    // Find LogNormal candidates
    for (auto f : nll_->constrainPdfsFast_)
    {
        if (verbose) std::cout<<"inserting candidate: "<< f->getX().GetName()<<std::endl;
        logNormal.insert( f->getX().GetName() );
    }

    if (verbose) {
        std::cout<<"[SimNLLDerivativesHelper][init] List of constraint pdfs parameters:"<<std::endl;
        std::cout<<"     "; for( const auto& v : logNormal) std::cout<<v<<",";
        std::cout<<std::endl;
    }

    // find rateParams candidates
    RooArgList pl(nll_->params_);
    for (int ir=0;ir<pl.getSize();++ir)
    {
        auto& f = pl[ir];
        std::string name = f.GetName();
        if (logNormal.find(name) != logNormal.end()) continue; // it is not a rateParams candidate
        if (f.isConstant()) {  // Not elegant. Try to remove the constraint _In as well
            if (verbose) std::cout<<"rateparam candidate is constant: "<< f.GetName()<<". Not considering it"<<std::endl;
            continue;
        }
        if (verbose) std::cout<<"inserting candidate: "<< f.GetName()<<std::endl;
        rateParams.insert(name);
    }

    if (verbose) {
        std::cout<<"[SimNLLDerivativesHelper][init] List of unconstraint parameters candidates:"<<std::endl;
        std::cout<<"     "; for( const auto& v : rateParams) std::cout<<v<<",";
        std::cout<<std::endl;
    }


    // loop over sim components
    if (verbose) { std::cout << "[SimNLLDerivativesHelper][init] loop over sim components 1"<<std::endl;}
    unsigned idx=0;
    for (std::vector<cacheutils::CachingAddNLL*>::const_iterator it = nll_->pdfs_.begin(), ed = nll_->pdfs_.end(); it != ed; ++it, ++idx) {
        cacheutils::CachingAddNLL * pdf = *it;
        if (pdf==nullptr) continue ; // ?!? needed?
        const RooAbsData* data= pdf->data();
        bool isWeighted = data->isWeighted() ;//and (dynamic_cast<const RooDataHist*>(data) != nullptr); // binned vs unbinned?!? TBC

        if (verbose) { std::cout << "[SimNLLDerivativesHelper][init] > considering pdf "<<pdf->GetName()<< ((isWeighted)?" isWeighted":" notWeighted") <<std::endl;}

        //std::set<std::string> toRemove;

        for(unsigned proc =0 ; proc < pdf ->pdfs_.size();++proc)
        {
            RooAbsReal * coeff = pdf -> coeffs_[proc] ; // what to do for RooProduct?
            const RooAbsReal* pdfi= pdf->pdfs_[proc].pdf();

            ProcessNormalization * pn = dynamic_cast<ProcessNormalization*> (coeff);
            RooProduct *pp = dynamic_cast<RooProduct*>(coeff);
            RooRealVar *rv = dynamic_cast<RooRealVar*>(coeff);
            const CMSHistSum*hs = dynamic_cast<const CMSHistSum*>(pdfi);
            //if (rv !=nullptr) std::cout<<"[SimNLLDerivativesHelper]::[FIXME]"<<"derivative for "<rv->GetName()<<" is known, but discarded in the loop :( "<<std::endl;
#ifdef DEBUG_CMSHISTSUM
            std::cout<<"[SimNLLDerivativeHelper][init] DEBUG CMSHistSum: Coeff ClassName"<<coeff->ClassName() <<std::endl;
            std::cout<< "coeff name: "<<coeff->GetName()<<" val="<<coeff->getVal()<<std::endl;
            std::cout<<"[SimNLLDerivativeHelper][init] DEBUG CMSHistSum: PDF is "<<pdfi->ClassName() << " is CMS histSum" << ((hs==nullptr)?"no":"yes") <<std::endl;
#endif

            if (pn == nullptr and pp != nullptr){
                for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
                {
                    if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) pn = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                    else {
                            if (verbose) { std::cout <<"components:";}
                            std::set<std::string> servers= getServersVars(&pp->components()[ip]);
                            for(auto v : servers) {
                                    if (logNormal.find(v) != logNormal.end()) {
                                    logNormal.erase(v);
                                    if (verbose) std::cout<<"|(LN)"<<v;
                                    }
                                    if (rateParams.find(v) != rateParams.end() and rv == nullptr) {
                                        if (verbose) std::cout<< "|(R)"<<v;
                                        rateParams.erase(v);
                                    } // what if v is a RooRealVar in the RooProduct?
                            }
                            if (verbose) { std::cout <<std::endl;}

                    }
                } 
                if (verbose) {
                    if(pn) std::cout<<"[SimNLLDerivativesHelper][init]: RooProduct has inside a ProcessNormalization. using it"<<std::endl;
                    else std::cout<<"[SimNLLDerivativesHelper][init]: RooProduct has no product normalization inside"<<std::endl;
                }
            }

#ifdef DEBUG_CMSHISTSUM
            // THIS Guy below I need to do it only if it is not a hist sum. For them I need to do something else
#endif
            if (hs != nullptr) { // this is to establish what I know how to make a derivative
                // this is an CMSHistSum pdf
                // loop over the parameters and remove asymm theta variables and
                // things different from rate parameters
                std::cout<<"[SimNLLDerivativesHelper][init][HS]"<< ">> Removing all asymmThetaList from histsum "<<std::endl;
                // remove all asymmThetaList
                for (int i=0;i<hs->coefList().getSize() ;++i)
                {
                    std::cout<<"[SimNLLDerivativesHelper][init][HS]"<< "  Considering coeff n."<<i<<" with name"<< hs->coefList()[i].GetName()<<" of class "<< hs->coefList()[i].ClassName()<<std::endl;
                    ProcessNormalization *pn= dynamic_cast<ProcessNormalization*> ( &hs->coefList()[i]);
                    RooProduct *pp = dynamic_cast<RooProduct*>( &hs->coefList()[i]);
                    RooRealVar *rv = dynamic_cast<RooRealVar*>( &hs->coefList()[i]);

                    if (pn == nullptr and pp != nullptr){// removing strange components from RooProduct
                        for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
                        {
                            std::cout<<"[SimNLLDerivativesHelper][DEBUG][HS]: RooProduct component "<<ip<<"/"<<pp->components().getSize()<<" is a "<<pp->components()[ip].ClassName() <<std::endl;
                            if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) pn = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                            else {
                                if (verbose) { std::cout <<"components:";}
                                std::set<std::string> servers= getServersVars(&pp->components()[ip]);
                                for(auto v : servers) {
                                    if (logNormal.find(v) != logNormal.end()) {
                                        logNormal.erase(v);
                                        if (verbose) std::cout<<"|(LN)"<<v;
                                    }
                                    if (rateParams.find(v) != rateParams.end() and rv == nullptr) {
                                        if (verbose) std::cout<< "|(R)"<<v;
                                        rateParams.erase(v);
                                    } // what if v is a RooRealVar in the RooProduct?
                                }
                                if (verbose) { std::cout <<std::endl;}

                            }
                        } 
                        if (verbose) {
                            if(pn) std::cout<<"[SimNLLDerivativesHelper][init][HS]: RooProduct has inside a ProcessNormalization. using it"<<std::endl;
                            else std::cout<<"[SimNLLDerivativesHelper][init][HS]: RooProduct has no product normalization inside"<<std::endl;
                        }
                    }

                    if (pn != nullptr){
                        for(int idx=0; idx< pn->asymmThetaList_.getSize() ;++idx){
                            RooAbsArg*v=pn->asymmThetaList_.at(idx);
                            if (verbose) {
                                std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                                    <<"Removing "<<v->GetName()<<" from the list of known derivatives because appears in an asymmThetaList "
                                    <<std::endl;
                            }
                            if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
                            if (rateParams.find(v->GetName()) != rateParams.end()) {rateParams.erase(v->GetName());}  // probably not necessary, since this are constrained
                        }

                        for(int idx=0; idx< pn->thetaList_.getSize() ;++idx){
                            RooAbsArg*v=pn->thetaList_.at(idx);
                            if (verbose) {
                                std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                                    <<"Removing "<<v->GetName()
                                    <<" from the list of known derivatives because appears in an thetaList "
                                    <<std::endl;
                            }
                            if (rateParams.find(v->GetName()) != rateParams.end()) {rateParams.erase(v->GetName());}  // probably not necessary, since this are constrained
                            //if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName()); //FIXME: remove this when implemented. This are the logNormal
                        }
                    }// if there is a process normalization do something
                    // rv
                    if (pn == nullptr and rv == nullptr) { // remove everything. I don't know how to deal with it
                        if (verbose) { std::cout <<"coeff variables:";}
                        std::set<std::string> servers= getServersVars(coeff);
                        for(auto v : servers) {
                            if (logNormal.find(v) != logNormal.end()) {
                                logNormal.erase(v);
                                if (verbose) std::cout<<"|(LN)"<<v;
                            }
                            if (rateParams.find(v) != rateParams.end()) {
                                if (verbose) std::cout<<"|(R)"<<v;
                                rateParams.erase(v);
                            } 
                        }
                    }

                }//coeff loop

                if (verbose) std::cout<<"[SimNLLDerivativeHelper][init] Removing bin-by-bin parameters: ";
                for (const auto bp : hs->binpars_){ // this are constrained parameters for which I don't know the derivatives
                    if (logNormal.find( bp->GetName() ) != logNormal.end()) {
                            logNormal.erase(bp->GetName());
                            if (verbose) std::cout<<bp->GetName()<<", ";
                    }
                if (verbose) std::cout<<std::endl;

                }//binpars loop
            }// histsum

            if (hs == nullptr){
                // look inside pdfi if name is there -> likely shape 
                if (verbose) { std::cout <<"removing pdf variables (shape):";}
                std::set<std::string> servers= getServersVars(pdfi);
                for(auto v : servers) {
                        if (logNormal.find(v) != logNormal.end()) {
                        logNormal.erase(v);
                        if (verbose) std::cout<<"|(LN)"<<v;
                        }
                        if (rateParams.find(v) != rateParams.end()) {
                            if (verbose) std::cout<< "|(R)"<<v;
                            rateParams.erase(v);
                        } // 
                    }
                if (verbose) { std::cout<<std::endl;}
            }

            if ( (pn == nullptr and hs==nullptr) or not isWeighted ) {  // remove all the lognormal candidates. Don't know how to deal with them

                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][init]:"
                        <<"Removing from the list of known derivatives because"
                        << ((pn==nullptr)? std::string(" coeff is not processNormalization but ")+ std::string(coeff->ClassName()): "")
                        << ((pn==nullptr and not isWeighted) ? " and":"" )
                        << ((not isWeighted)?" dataset is not weighted (unbinned?)":"")
                        <<std::endl;
                }
#ifdef DERIVATIVE_LOGNORMAL_DEBUG
                if (verbose) { 
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " -- COEFF -- "<<std::endl;
                    coeff->Print("V");
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " -- PDF -- "<<std::endl;
                    pdfi->Print("V");
                    std::cout<<"-----------------"<<std::endl;
                }
#endif
                
                if (verbose) { std::cout <<"coeff variables:";}
                std::set<std::string> servers= getServersVars(coeff);
                for(auto v : servers) {
                        if (logNormal.find(v) != logNormal.end()) {
                        logNormal.erase(v);
                        if (verbose) std::cout<<"|(LN)"<<v;
                        }
                        if (rateParams.find(v) != rateParams.end()) {
                            if (verbose) std::cout<<"|(R)"<<v;
                            rateParams.erase(v);
                        } 
                }
                //--
                if (verbose) { std::cout<<std::endl <<"pdf variables:";}

                servers.clear();
                /*std::set<std::string>*/ servers= getServersVars(pdfi); // do I still need this block?
                for(auto v : servers) {
                        if (logNormal.find(v) != logNormal.end()) {
                        logNormal.erase(v);
                        if (verbose) std::cout<<"|(LN)"<<v;
                        }
                        if (rateParams.find(v) != rateParams.end()) {
                            rateParams.erase(v);
                            if (verbose) std::cout<<"|(R)"<<v;
                        } 
                }
                if (verbose) { std::cout<<std::endl;}

                continue; //  if pn==nullptr or not weighted   (A)
            } //  not a process normalization coefficient

            if (hs != nullptr) continue; // match the one above (A), since in this case all the things should have been already dealt wit

            std::cout<<"[SimNLLDerivativesHelper][init]"<< ">> Removing all asymmThetaList "<<std::endl;
            // remove all asymmThetaList
            for(int idx=0; idx< pn->asymmThetaList_.getSize() ;++idx){
                RooAbsArg*v=pn->asymmThetaList_.at(idx);
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()<<" from the list of known derivatives because appears in an asymmThetaList "
                        <<std::endl;
                }
                if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
                if (rateParams.find(v->GetName()) != rateParams.end()) {rateParams.erase(v->GetName());}  // probably not necessary, since this are constrained
            }

            // remove all otherFactorList
            std::cout<<"[SimNLLDerivativesHelper][init]"<< ">> Removing all otherfactorsList "<<std::endl;
            for(int idx=0; idx< pn->otherFactorList_.getSize() ;++idx){
                RooAbsArg*v=pn->otherFactorList_.at(idx); 
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()<<" from the list of known derivatives (LN) because appears in an otherFactors "
                        <<std::endl;
                }
                if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
                if (dynamic_cast<RooRealVar*>(v)==nullptr and rateParams.find(v->GetName()) != rateParams.end()) {
                    if (verbose) {
                        std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                            <<"Removing "<<v->GetName()<< " : "<<v->ClassName()<<" from the list of known derivatives (R) because appears in an otherFactors "
                            <<std::endl;
                    }
                    rateParams.erase(v->GetName());
                }  
            }

            // remove all symmThetaList
            for(int idx=0; idx< pn->thetaList_.getSize() ;++idx){
                RooAbsArg*v=pn->thetaList_.at(idx);
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()
                       // << ","<<v->GetX().GetName()
                       // << ","<<v->GetMean().GetName()
                       // << ","<<v->GetSigma().GetName()
                        <<" from the list of known derivatives because appears in an thetaList "
                        <<std::endl;
                }
                if (rateParams.find(v->GetName()) != rateParams.end()) {rateParams.erase(v->GetName());}  // probably not necessary, since this are constrained
                //if (rateParams.find(v->GetX().GetName()) != rateParams.end()) {rateParams.erase(v->GetX().GetName());}  // probably not necessary, since this are constrained
                //if (rateParams.find(v->GetMean().GetName()) != rateParams.end()) {rateParams.erase(v->GetMean().GetName());}  // probably not necessary, since this are constrained
                //if (rateParams.find(v->GetSigma().GetName()) != rateParams.end()) {rateParams.erase(v->GetSigma().GetName());}  // probably not necessary, since this are constrained
            }
            std::cout<<"[SimNLLDerivativesHelper][init]"<< ">> Continue loop "<<std::endl;

        }

        // cross check derivatives with logNormal
        //for (auto d : derivatives){
        //    if (logNormal.find(d.first) == logNormal.end()){
        //        //remove derivative
        //        derivatives.erase(d.first());
        //    }
        //}
        std::cout<<"[SimNLLDerivativesHelper][init]"<< "> Continue loop "<<std::endl;

    } // loop over sim components

    // Construct derivatives
    std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructing the derivatives sum for logNormal"<<std::endl;
    std::cout<<"[SimNLLDerivativesHelper][init]"<<"logNormals derivatives : "; for (auto name :logNormal) std::cout<<name<<","; std::cout<<std::endl;
    for (auto name : logNormal){
        RooArgList list;
        // loop over sim components
        unsigned idx=0;
        for (std::vector<cacheutils::CachingAddNLL*>::const_iterator it = nll_->pdfs_.begin(), ed = nll_->pdfs_.end(); it != ed; ++it, ++idx) {
            // construct the derivative as sum
            cacheutils::CachingAddNLL * pdf = *it;
            if (pdf==nullptr) continue ; // ?!? needed?
            const RooAbsData* data= pdf->data();
            //(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data,std::string thetaname)
            int found;
            DerivativeLogNormal *der= new DerivativeLogNormal((std::string("dln_")+name+Form("_%d",idx)).c_str(),"",pdf,dynamic_cast<const RooDataSet*>(data),name,found);
            if(found) list.add(*der); // or addOwned or add?
            else der->Delete();

            if (not found){
                DerivativeLogNormalCMSHistSum *der= new DerivativeLogNormalCMSHistSum((std::string("dln_")+name+Form("_%d",idx)).c_str(),"",pdf,dynamic_cast<const RooDataSet*>(data),name,found);
                if(found) list.add(*der); // or addOwned or add?
                else der->Delete();
            }
        }

        if (list.getSize() == 0){ // old fix for shapes. Leave it, probably useful for debug
            std::cout<<"[SimNLLDerivativesHelper][init][ERROR]"<< "Something went wrong: Unable to match  "<< name <<" with kappas"<<std::endl;
            continue;
        }
        // add derivative of constraint 
        // nll_->constrainPdfsFast_
        std::cout<<"[SimNLLDerivativesHelper][init]"<< "Adding constraint term for "<< name <<std::endl;
        for(unsigned i=0;i<nll_->constrainPdfsFast_.size();++i)
        {
            if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< " Considering "<< nll_->constrainPdfsFast_[i]->getX().GetName() <<std::endl;
            if (nll_->constrainPdfsFast_[i]->getX().GetName() != name) continue;

            //nll_->constrainPdfsFast_[i]->getMean().Print("V");
            if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Loop. found constraint "<< nll_->constrainPdfsFast_[i]->getX().GetName() <<" == "<<name<<". Building RooFormulaVar.  X="<<nll_->constrainPdfsFast_[i]->getX().GetName()<< ", "<<nll_->constrainPdfsFast_[i]->getMean().GetName()<<", "<<nll_->constrainPdfsFast_[i]->getSigma().GetName() <<std::endl;

            RooFormulaVar *constr;
            if ( std::string(nll_->constrainPdfsFast_[i]->getSigma().GetName()) == std::string("1")) // new RooFormulaVar get crazy if they parse 1
            {
                constr = new RooFormulaVar( (std::string("constr_")+name).c_str(),"@2*(@0-@1)",RooArgList(nll_->constrainPdfsFast_[i]->getX(),nll_->constrainPdfsFast_[i]->getMean(),unmask_));// (x-mean)/sigma^2
            }
            else{ // likely untested
                constr = new RooFormulaVar( (std::string("constr_")+name).c_str(),"@3*(@0-@1)/(@2*@2)",RooArgList(nll_->constrainPdfsFast_[i]->getX(),nll_->constrainPdfsFast_[i]->getMean(), nll_->constrainPdfsFast_[i]->getSigma(),unmask_));// (x-mean)/sigma^2
            }
            list.add(*constr) ; // or addOwned or add?
        }
        if (verbose) {std::cout <<"[SimNLLDerivativesHelper][init] --- LIST OF addition ---"<<std::endl;
            list.Print("V");
        }
        if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructing RooAddition for  "<< name <<std::endl;
        derivatives_[name] = new RooAddition(("derivative_"+name).c_str(),"",list);
        if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructed sum for "<< name <<std::endl;
    } 

    std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructing the derivatives sum for rateParams"<<std::endl;
    std::cout<<"[SimNLLDerivativesHelper][init]"<<"rate Params: "; for (auto name :rateParams) std::cout<<name<<","; std::cout<<std::endl;
    for (auto name : rateParams){
        RooArgList list;
        // loop over sim components
        unsigned idx=0;
        for (std::vector<cacheutils::CachingAddNLL*>::const_iterator it = nll_->pdfs_.begin(), ed = nll_->pdfs_.end(); it != ed; ++it, ++idx) {
            // construct the derivative as sum
            cacheutils::CachingAddNLL * pdf = *it;
            if (pdf==nullptr) continue ; // ?!? needed?
            const RooAbsData* data= pdf->data();
            //(const char *name, const char *title, RooAbsPdf *pdf, RooAbsData *data,std::string thetaname)
            const RooAbsReal* pdfi= (not pdf->pdfs_.empty()) ?  pdf->pdfs_[0].pdf() :nullptr;
            const CMSHistSum*hs = dynamic_cast<const CMSHistSum*>(pdfi);

            int found;
            if (hs == nullptr){
                DerivativeRateParam *der= new DerivativeRateParam((std::string("dln_")+name+Form("_%d",idx)).c_str(),"",pdf,dynamic_cast<const RooDataSet*>(data),name,found);
                if(found) {
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " adding DerivativeRateParam for parameter "<<name<<std::endl;
                    list.add(*der); // or addOwned or add?
                }
                else der->Delete();
            }

            if (hs) { // try HistSum
                DerivativeRateParamCMSHistSum *der= new DerivativeRateParamCMSHistSum((std::string("dln_")+name+Form("_%d",idx)).c_str(),"",pdf,dynamic_cast<const RooDataSet*>(data),name,found);
                if(found) {
                    list.add(*der); // or addOwned or add?
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " adding DerivativeRateParamiCMSHistSum for parameter "<<name<<std::endl;
                }
                else der->Delete();
            }
        }

        if (list.getSize() == 0){ // old fix for shapes. Leave it, probably useful for debug
            std::cout<<"[SimNLLDerivativesHelper][init][ERROR]"<< "Something went wrong: Unable to match  "<< name <<" with kappas"<<std::endl;
            continue;
        }
        // add derivative of constraint 
        // nll_->constrainPdfsFast_
        if (verbose) {std::cout <<"[SimNLLDerivativesHelper][init] --- LIST OF addition ---"<<std::endl;
            list.Print("V");
        }
        if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructing RooAddition for  "<< name <<std::endl;
        derivatives_[name] = new RooAddition(("derivative_"+name).c_str(),"",list);
        if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructed sum for "<< name <<std::endl;
    } 

    if (verbose){
        std::cout<<"[SimNLLDerivativesHelper][init]"<< "DONE INIT" <<std::endl;
        //std::cout<<"[SimNLLDerivativesHelper][init] DEBUG evaluate test"<<std::endl;
        //std::cout<< derivatives_.begin()->first <<" : "<<derivatives_.begin()->second->getVal(); // DEBUG
    }
}

//def getServers(node):
//    servers = []
//    iter = node.serverIterator()
//    while True:
//        server = iter.Next()
//        if server == None:
//            break
//        servers.append(server)
//    return servers

std::set<std::string>  SimNLLDerivativesHelper::getServersVars(const RooAbsArg *node){
    std::set<std::string> R;
    auto iter = node->serverIterator();
    while (true)
    {
        auto server = iter->Next(); 
        if (server==nullptr)break;

        if (dynamic_cast<RooRealVar*> (server) != nullptr)
        {
            R.insert(std::string(server->GetName()));
        }
        else{
            std::set<std::string> r1 = getServersVars((RooAbsArg*)server);
            R.merge(r1);
        }
    }
    return R;
}



/// ---
DerivativeRateParam::DerivativeRateParam(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data, const std::string& ratename,int&found) :
    DerivativeAbstract(name,title,pdf,data),
    ratename_(ratename)
{
    //setDirtyInhibit(true); // TODO: understand cache and proxies
    setOperMode(RooAbsArg::ADirty,false); // Recursive, bool

    found=0; // keep track if at least one depends on the given kappa
    for (unsigned int i=0;i<pdf_->pdfs_.size();++i) // process loop
    {
        if(verbose) std::cout<<"[DerivativeRateParam]: loop "<< i <<"/"<<pdf_->pdfs_.size()<<" considering pdf "<<pdf_->GetName() <<std::endl;
        int val=-1;
        ProcessNormalization *p = dynamic_cast<ProcessNormalization*> (pdf_->coeffs_[i]);
        if (verbose) std::cout<<"[DerivativeRateParam]: coeff is of class "<<pdf_->coeffs_[i]->ClassName()<<std::endl;
        // RooProduct
        RooProduct *pp = dynamic_cast<RooProduct*>(pdf_->coeffs_[i]);
        if (p == nullptr and pp != nullptr){
            for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
            {
                if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) p = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
            }
        }
        //
        RooRealVar *rv = dynamic_cast<RooRealVar*> (pdf_->coeffs_[i]);
        if (rv !=nullptr){
            val=-2; // the whole coeff is a RooRealVar
            found=1;
            if(verbose) std::cout<<"[DerivativeRateParam]: All Coeff corresponding to "<<ratename_<<" in pdf "<<pdf_->GetName() <<" is a RooRealVar"<<std::endl;
        }
        else if (p != nullptr) { // unable to add the sum for this process

            for (int j=0 ; j<p->otherFactorList_.getSize();++j)
            {
                if ( p->otherFactorList_.at(j)->GetName() == ratename_  ) {val=int(j);}
            }
            if (val==-1) { 
                if(verbose) std::cout<<"[DerivativeRateParam]: Unable to find rate corresponding to "<<ratename_<<" in pdf "<<pdf_->GetName() <<std::endl;
            }// something
            else{
                found=1;
            }
        }
        else{
            if(verbose) std::cout<<"[DerivativeRateParam]: Unable to find rate corresponding to "<<ratename_<<" in pdf "<<pdf_->GetName() << "because not ProcessNormalization" <<std::endl;
        }
        rate_pos_.push_back(val);
    } //for -- process loop
    if(verbose) std::cout<<"[DerivativeRateParam]: constructed derivative for "<<ratename_<<" of pdf "<<pdf_->GetName() <<std::endl;
}// constructor

Double_t DerivativeRateParam::evaluate() const {
#ifdef DERIVATIVE_RATEPARAM_DEBUG
    if(verbose) std::cout<<"[DerivativeRateParam]: evaluate"<<std::endl;
#endif
    /* The derivative of a kappa in a channel is
     *   sum_bin lambdat_b - data_b lambdat_b/lambda_b
     *   lambda_b -> sum of expectations in the bin
     *   lambdat_b -> sum of ( expectations with delta(rate param))
     */

    double sum=0.0;
    // caching: TODO
#ifdef DERIVATIVE_RATEPARAM_DEBUG
    if(verbose) {
        std::cout<<"[DerivativeRateParam]: data"<<std::endl;
        data_->Print("V");
        std::cout<<"[DerivativeRateParam]: --"<<std::endl;
    }
#endif
    
    //const RooDataHist* dh = dynamic_cast<const RooDataHist*> (data_);
    //if (dh==nullptr) std::cout<<"[DerivativeRateParam]:"<<" NO ROODATAHIST"<<std::endl;

    // FIX for ADDNLL_ROOREALSUM_KEEPZEROS 
    std::vector<double> reminder(pdf_->pdfs_.size(),0.0); // the sums must do 1, the ib==ndata does that with data=0
    int ndata=data_->numEntries();

    for (int ib=0;ib<=ndata;++ib) 
    //for (unsigned int ib=0;ib<pdf_->weights_.size();++ib)
    {
        auto x= (ib<ndata) ? data_->get(ib): data_->get(0) ;
        double db = (ib<ndata) ? data_->weight() : 0.; 
        double bw = 1; 
        if (x->getSize()==1) {
            RooRealVar *xvar=dynamic_cast<RooRealVar*>(x->first());
            const RooAbsBinning &bins = xvar->getBinning();
            bw=bins.binWidth(0); //  only costant supported. Need to figure out for zero bins
        }

        //double lambda=pdf_->getVal()*bw; // sum coeff*pdf
        double lambdat=0.0;
        double lambda=0.;// debug: recompute lambda with loop. not sure that the term before include everything I want

        for (unsigned int i=0;i<pdf_->pdfs_.size();++i) // loop over processes
        { 
            // find rate param. TODO make it more efficient using rate_pos
            ProcessNormalization *c = dynamic_cast<ProcessNormalization*>(pdf_->coeffs_[i]);
            RooRealVar *rv = dynamic_cast<RooRealVar*>(pdf_->coeffs_[i]);
            
            // fix RooProduct
            RooProduct *pp = dynamic_cast<RooProduct*>(pdf_->coeffs_[i]);
            if (c == nullptr and pp != nullptr){
                for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
                {
                    if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) c = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                }
            }
            RooRealVar *rate = dynamic_cast<RooRealVar*>( (rate_pos_[i]==-2)? rv: (rate_pos_[i]>=0) ? &c->otherFactorList_[rate_pos_[i]] : nullptr );

            //---
            int diracDelta = (rate_pos_[i]>=0 or rate_pos_[i]==-2)?1 : 0;
            if (pdf_->isRooRealSum_ and verbose)std::cout<<"[DerivativeRateParam]:DEBUG "<<" pdf is a RooRealSumPdf"<<std::endl;

            if (diracDelta and rate==nullptr){
                std::cout<<"[DerivativeRateParam]:[ERROR] RateParam"<<ratename_<<" is not a rooRealVar. rate_pos="<<rate_pos_[i]<<std::endl;
                // it will crash
            }

            // TH1
            if (ib<ndata) reminder[i] += pdf_->pdfs_.at(i).pdf()->getVal(x)*bw;
            double expectedEvents= (ib<ndata)  ? pdf_->coeffs_[i]->getVal(x) * pdf_->pdfs_.at(i).pdf()->getVal(x)*bw : 
                        pdf_->coeffs_[i]->getVal(x) * (1.-reminder[i]);

            if (diracDelta){ // only if process is affected by rate
                // for rate==0, we should probably think to set rate=1 and to get the expected events again. This will cancel the 0/0 in the right way.
                //lambdat+= expectedEvents * diracDelta/rate->getVal(); // over rate. Possible problem with rate=0. OK.

                double current=rate->getVal();
                rate->setVal(1.0);
                double expectedEventsTilde= (ib<ndata)  ? pdf_->coeffs_[i]->getVal(x) * pdf_->pdfs_.at(i).pdf()->getVal(x)*bw : pdf_->coeffs_[i]->getVal(x) * (1.-reminder[i]);
                rate->setVal(current);

                lambdat+= expectedEventsTilde; 
            }
            lambda += expectedEvents;

            if (pdf_->basicIntegrals_ and verbose)std::cout<<"[DerivativeRateParam]:DEBUG "<<" PDF has basic integrals:"<<pdf_->basicIntegrals_<<std::endl;

#ifdef DERIVATIVE_RATEPARAM_DEBUG
#endif
        } 
        sum += lambdat * (1. - db/lambda);

#ifdef DERIVATIVE_RATEPARAM_DEBUG
        if(verbose) std::cout<<"[DerivativeRateParam]: ibin="<< ib<< " data="<<db << " lambda="<<lambda << "( from pdf_" << pdf_->getVal()*bw <<")" <<" lambdat="<<lambdat<<" bw="<<bw<< " running sum="<<sum<<std::endl;
#endif
    }

#ifdef DERIVATIVE_RATEPARAM_DEBUG
    if(verbose) std::cout<<"[DerivativeRateParam]: ratename="<<ratename_<<" partial="<<sum<<std::endl;
#endif
    return sum;

}

DerivativeRateParam* DerivativeRateParam::clone(const char *name) const {
    if(verbose) std::cout<<"[DerivativeRateParam]: clone"<<std::endl;
    int found;
    return new DerivativeRateParam( 
            (name)?name:GetName(), GetTitle(), pdf_, data(), ratename_, found) ;
}


// ------------------------------------------------ CMS HIST SUM  part 

DerivativeRateParamCMSHistSum* DerivativeRateParamCMSHistSum::clone(const char *name) const {
    if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: clone"<<std::endl;
    int found;
    return new DerivativeRateParamCMSHistSum( 
            (name)?name:GetName(), GetTitle(), pdf_, data(), ratename_, found) ;
}

DerivativeRateParamCMSHistSum::DerivativeRateParamCMSHistSum(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data, const std::string& ratename,int&found) :
    DerivativeAbstract(name,title,pdf,data),
    ratename_(ratename)
{
    if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: constructor"<<std::endl;
    //setDirtyInhibit(true); // TODO: understand cache and proxies
    setOperMode(RooAbsArg::ADirty,false); // Recursive, bool
    if (this->pdf_->pdfs_.size() !=1 ) {
        std::cout<<"[DerivativeRateParamCMSHistSum]: ERROR expected only one component in CachingAddNLL. I got: " << this->pdf()->pdfs_.size() <<std::endl;
        // [TODO: raise exception here]
     }
    const CMSHistSum * histsum = dynamic_cast<const CMSHistSum*> (this->pdf()->pdfs_[0].pdf());
    if (histsum==nullptr)
    {
        std::cout<<"[DerivativeRateParamCMSHistSum]: ERROR. Pdf is not an HistSum but "<<this->pdf()->pdfs_[0].pdf()->ClassName() <<std::endl;
        // [TODO: raise exception here]
    }
    histsum->updateCache(); // we need to have an up-to-date  ?!?
    found=0; // keep track if at least one depends on the given kappa
    for (unsigned int i=0; i< histsum->coefList().size(); ++i){ // process loop
        if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: loop "<< i <<"/"<< histsum->coefList().size() <<std::endl;
        int val=-1;
        ProcessNormalization *p = dynamic_cast<ProcessNormalization*> (&histsum->coefList()[i]);
        if (verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: coeff is of class "<<histsum->coefList()[i].ClassName()<<std::endl;
        // RooProduct
        RooProduct *pp = dynamic_cast<RooProduct*>(&histsum->coefList()[i]);
        if (p == nullptr and pp != nullptr){
            std::cout<<"[DerivativeRateParamCMSHistSum]: looking for a ProcessNormalization in a RooProduct"<<std::endl;
            for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
            {
                if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) {
                    if (p!=nullptr) std::cout<<"[DerivativeRateParamCMSHistSum] ERROR. there is more than 1 ProcessNormalization in the RooProduct"<<std::endl;
                    p = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                }
            }
        }
        //
        RooRealVar *rv = dynamic_cast<RooRealVar*> (&histsum->coefList()[i]);
        if (rv !=nullptr){
            val=-2; // the whole coeff is a RooRealVar
            found=1;
            rate_.setArg(*rv);
            if(verbose) std::cout<<"[DerivativeRateParamCMSHistsum]: All Coeff corresponding to "<<ratename_<<" in pdf n."<<i <<" is a RooRealVar"<<std::endl;
        }
        else if (p != nullptr) { // unable to add the sum for this process
            if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: Looking inside otherFactors in process normalization " <<std::endl;
            for (int j=0 ; j<p->otherFactorList_.getSize();++j)
            {
                if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: otherFactor "<<j<<"/"<<p->otherFactorList_.getSize() <<std::endl;
                if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: the onthefactor pointer is "<< p->otherFactorList_.at(j)->ClassName()<<std::endl;
                if ( p->otherFactorList_.at(j)->GetName() == ratename_  ) {
                    val=int(j);
                    //rate_.setArg(*p->otherFactorListVec_.at(j));
                    rv_= dynamic_cast<RooAbsReal*>(p->otherFactorList_.at(j));
                    if (rv_ == nullptr) std::cout<<"[DerivativeRateParamCMSHistSum] ERROR: Rate is not a RooAbsReal?!?"<<std::endl;
                }
            }
            if (val==-1) { 
                if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: Unable to find rate corresponding to "<<ratename_<<" in pdf n."<<i <<std::endl;
            }// something
            else{
                found=1;
            }
        }
        else{
            if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: Unable to find rate corresponding to "<<ratename_<<" in pdf n."<<i << "because not ProcessNormalization" <<std::endl;
        }
        rate_pos_.push_back(val);
    }// process loop
    if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: -- construction done --"<<std::endl;
}// constructor

Double_t DerivativeRateParamCMSHistSum::evaluate() const {
#ifdef DEBUG_CMSHISTSUM
    if(verbose) std::cout<<"[DerivativeRateParamCMSHistSum]: evaluate"<<std::endl;
#endif
    /* The derivative of a kappa in a channel is
     *   sum_bin lambdat_b - data_b lambdat_b/lambda_b
     *   lambda_b -> sum of expectations in the bin
     *   lambdat_b -> sum of ( expectations with delta(rate param))
     */
    double sum=0.0;
    const CMSHistSum * histsum = dynamic_cast<const CMSHistSum*> ( pdf_->pdfs_[0].pdf()); // TODO save it somewhere
    histsum->updateCache(); // we need to have an up-to-date compcache_
    for (unsigned ib=0;ib<histsum->data_.size(); ++ib)
    {
        double db = histsum->data_[ib];
        double lambdat=0.0;
        double lambda=0.0; // expected events
        for (unsigned i = 0 ; i < histsum -> vcoeffpars_.size() ; ++i)
        {
            double coeffval = histsum->vcoeffpars_[i]->getVal();
            double frac = histsum->compcache_[i][ib];

            int diracDelta = (rate_pos_[i]>=0 or rate_pos_[i]==-2)?1 : 0;
            
#ifdef DEBUG_CMSHISTSUM
            /*
            std::cout<<"[DerivativeRateParamCMSHistSum][evaluate] Considering:"<<std::endl;
            std::cout<<"                               * data bin "<<ib<<") data="<<db<<std::endl;
            std::cout<<"                               * process "<<i<<") expected= " <<coeffval*frac<<std::endl;
            std::cout<<"                               * delta_p "<<diracDelta<<std::endl;
            */
#endif
            // I need to fetch the rate parameter and get its value

            lambda += coeffval * frac;
            if (diracDelta){
                double r = rv_->getVal();
#ifdef DEBUG_CMSHISTSUM
                std::cout<<"[DerivativeRateParamCMSHistSum] rate != 0 ? "<< ((abs(r)>1E-15)?"yes":"no")<<std::endl;
#endif
                // close to 0 this may have a problem; I could set rate to one and recompute the cache each time :( as done for the nonCMSHist version
                //  rate -> setVal = 1; expectedEventsTilde = coeff*frac; rate->setVal current
                lambdat += (abs(r)>1E-15) ? (coeffval * frac / r ) : 0.; 
            }
        }
#ifdef DEBUG_CMSHISTSUM
        std::cout<<"[DerivativeRateParamCMSHistSum]: lambda="<<lambda<<std::endl;
        std::cout<<"[DerivativeRateParamCMSHistSum]: lambdat="<<lambdat<<std::endl;
#endif
        sum += lambdat * (1. - db/lambda);
    }//bin loop
    return sum;
}

// -- CMS HIST SUM  part  logNormal
DerivativeLogNormalCMSHistSum* DerivativeLogNormalCMSHistSum::clone(const char *name) const {
    if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: clone"<<std::endl;
    int found;
    return new DerivativeLogNormalCMSHistSum( 
            (name)?name:GetName(), GetTitle(), pdf_, data(), thetaname_, found) ;
}

// constructor
DerivativeLogNormalCMSHistSum::DerivativeLogNormalCMSHistSum(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data,const std::string& thetaname,int&found) :
    DerivativeAbstract(name,title,pdf,data),
    thetaname_(thetaname)
{
    if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: constructor"<<std::endl;
    setOperMode(RooAbsArg::ADirty,false); // Recursive, bool

    if (this->pdf_->pdfs_.size() !=1 ) {
        std::cout<<"[DerivativeLogNormalCMSHistSum]: ERROR expected only one component in CachingAddNLL. I got: " << this->pdf()->pdfs_.size() <<std::endl;
        // [TODO: raise exception here]
     }
    const CMSHistSum * histsum = dynamic_cast<const CMSHistSum*> (this->pdf()->pdfs_[0].pdf());
    if (histsum==nullptr)
    {
        std::cout<<"[DerivativeLogNormalCMSHistSum]: ERROR. Pdf is not an HistSum but "<<this->pdf()->pdfs_[0].pdf()->ClassName() <<std::endl;
        // [TODO: raise exception here]
    }

    histsum->updateCache(); // we need to have an up-to-date  ?!?

    found=0; // keep track if at least one depends on the given kappa

    for (unsigned int i=0; i< histsum->coefList().size(); ++i){ // process loop
        if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: loop "<< i <<"/"<< histsum->coefList().size() <<std::endl;
        int val=-1;
        if (verbose) std::cout<<"[DerivativeLogNormal]: coeff is of class "<<histsum->coefList()[i].ClassName()<<std::endl;

        ProcessNormalization *p = dynamic_cast<ProcessNormalization*> (&histsum->coefList()[i]);
        RooProduct *pp = dynamic_cast<RooProduct*>(&histsum->coefList()[i]);

        if (p == nullptr and pp != nullptr){
            for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
            {
                    if (p!=nullptr) std::cout<<"[DerivativeLogNormalCMSHistSum] ERROR. there is more than 1 ProcessNormalization in the RooProduct"<<std::endl;
                if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) p = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
            }
        }

        //
        if (p != nullptr) { // unable to add the sum for this process

            for (unsigned int j=0 ; j<p->logKappa_.size();++j)
            {
                if ( p->thetaList_.at(j)->GetName() == thetaname_  ) {val=int(j);}
            }
            if (val==-1) { 
                if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: Unable to find kappa corresponding to "<<thetaname_<<" in pdfi "<< i <<std::endl;
            }// something
            else{
                found=1;
            }
        }
        else{
            if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: Unable to find kappa corresponding to "<<thetaname_<<" in pdfi "<<i << " because not ProcessNormalization" <<std::endl;
        }

        kappa_pos_.push_back(val);
    } // process loop
    if(verbose) std::cout<<"[DerivativeLogNormal]: constructed derivative for "<<thetaname_<<" of pdf "<<pdf_->GetName() <<std::endl;
}

Double_t DerivativeLogNormalCMSHistSum::evaluate() const {
#ifdef DEBUG_CMSHISTSUM
    if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: evaluate"<<std::endl;
#endif
    /* The derivative of a kappa in a channel is
     *   sum_bin lambdat_b - data_b lambdat_b/lambda_b
     *   lambda_b -> sum of expectations in the bin
     *   lambdat_b -> sum of ( expecations times logK)
     */

    double sum=0.0;
    const CMSHistSum * histsum = dynamic_cast<const CMSHistSum*> ( pdf_->pdfs_[0].pdf()); // TODO save it somewhere
    histsum->updateCache(); // we need to have an up-to-date compcache_

    for (unsigned ib=0;ib<histsum->data_.size(); ++ib)
    {

        double db = histsum->data_[ib];
        double lambdat=0.0;
        double lambda=0.0; // expected events
        for (unsigned i = 0 ; i < histsum -> vcoeffpars_.size() ; ++i)
        {
            double coeffval = histsum->vcoeffpars_[i]->getVal();
            double frac = histsum->compcache_[i][ib];
            ProcessNormalization *c = dynamic_cast<ProcessNormalization*> (&histsum->coefList()[i]);
            // find kappa
            RooProduct *pp = dynamic_cast<RooProduct*>(&histsum->coefList()[i]);
            if (c == nullptr and pp != nullptr){
                for ( int ip =0 ;ip<pp->components().getSize(); ++ip)
                {
                    if ( dynamic_cast<ProcessNormalization*>( &pp->components()[ip]) != nullptr) c = dynamic_cast<ProcessNormalization*>( &pp->components()[ip]);
                }
            }

            //---
            double logK = (kappa_pos_[i]>=0)? c -> logKappa_ [kappa_pos_[i]] : 0.0;

            double expectedEvents= coeffval * frac;
            
#ifdef DEBUG_CMSHISTSUM
            std::cout<<"[DerivativeLogNormalCMSHistSum][evaluate] Considering:"<<std::endl;
            std::cout<<"                               * data bin "<<ib<<") data="<<db<<std::endl;
            std::cout<<"                               * process "<<i<<") expected= " <<expectedEvents<<std::endl;
            std::cout<<"                               * logK "<<logK<<std::endl;
#endif
            lambdat+= expectedEvents * logK;
            lambda += expectedEvents;

        } // process loop
#ifdef DEBUG_CMSHISTSUM
            std::cout<<"[DerivativeLogNormalCMSHistSum][evaluate] lambda = "<<lambda<<", lambdat = "<<lambdat<<std::endl;
#endif
        sum += lambdat * (1. - db/lambda);

    } // data bin loop

#ifdef DEBUG_CMSHISTSUM
    if(verbose) std::cout<<"[DerivativeLogNormalCMSHistSum]: thetaname="<<thetaname_<<" partial="<<sum<<std::endl;
#endif
    return sum;
}

// -----------------------------------------------------------
