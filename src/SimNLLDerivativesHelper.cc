#include "HiggsAnalysis/CombinedLimit/interface/SimNLLDerivativesHelper.h"
#include "RooDataHist.h"
//#define DERIVATIVE_LOGNORMAL_DEBUG 1

DerivativeLogNormal::DerivativeLogNormal(const char *name, const char *title, cacheutils::CachingAddNLL *pdf, const RooDataSet *data, const std::string& thetaname,int&found) :
    RooAbsReal(name, title),
    data_(data),
    pdf_(pdf),
    pdfproxy_( ( std::string("proxy_")+name).c_str() ,"",pdf),
    thetaname_(thetaname)
{
    setDirtyInhibit(true); // TODO: understand cache and proxies

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
    
    const RooDataHist* dh = dynamic_cast<const RooDataHist*> (data_);

    // FIX for ADDNLL_ROOREALSUM_KEEPZEROS 
    std::vector<double> reminder(pdf_->pdfs_.size(),0.0); // the sums must do 1, the ib==ndata does that with data=0
    int ndata=data_->numEntries();

    for (int ib=0;ib<=ndata;++ib) 
    //for (unsigned int ib=0;ib<pdf_->weights_.size();++ib)
    {
        auto x= (ib<ndata) ? data_->get(ib): data_->get(0) ;
        double db = (ib<ndata) ? data_->weight() : 0.; 
        double bw = (dh != nullptr) ? dh->binVolume(): 1.0;

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
            ///--- from CachingAddNLL
            //static bool expEventsNoNorm = runtimedef::get("ADDNLL_ROOREALSUM_NONORM");
            //---
            double logK = (kappa_pos_[i]>=0)? c -> logKappa_ [kappa_pos_[i]] : 0.0;
            if (pdf_->isRooRealSum_ and verbose)std::cout<<"[DerivativeLogNormal]:DEBUG "<<" pdf is a RooRealSumPdf"<<std::endl;

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
            if(verbose and DERIVATIVE_LOGNORMAL_DEBUG > 1) std::cout<<"[DerivativeLogNormal]: *** ibin="<< ib<< " data="<<db <<" iproc="<< i<<" expEvents="<<expectedEvents<< " logK="<<logK<<" coeff="<<pdf_->coeffs_[i]->getVal() <<" pdf="<<pdf_->pdfs_.at(i).pdf()->getVal() <<" bw="<<bw <<std::endl;
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

    if (verbose) {
        std::cout<<"[SimNLLDerivativesHelper][init] List of constraint pdfs parameters:"<<std::endl;
        std::cout<<"     "; for( const auto& v : logNormal) std::cout<<v<<",";
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
                                    if (verbose) std::cout<<"|"<<v;
                                    }
                            }
                            if (verbose) { std::cout <<std::endl;}

                    }
                } 
                if (verbose) {
                    if(pn) std::cout<<"[SimNLLDerivativesHelper][init]: RooProduct has inside a ProcessNormalization. using it"<<std::endl;
                    else std::cout<<"[SimNLLDerivativesHelper][init]: RooProduct has no product normalization inside"<<std::endl;
                }
            }

            if (pn == nullptr or not isWeighted) {  // remove all the lognormal candidates. Don't know how to deal with them

                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][init]:"
                        <<"Removing from the list of known derivatives because"
                        << ((pn==nullptr)? std::string(" coeff is not processNormalization but ")+ std::string(coeff->ClassName()): "")
                        << ((pn==nullptr and not isWeighted) ? " and":"" )
                        << ((not isWeighted)?" dataset is not weighted (unbinned?)":"")
                        <<std::endl;
                }
                if (verbose) { 
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " -- COEFF -- "<<std::endl;
                    coeff->Print("V");
                    std::cout<<"[SimNLLDerivativesHelper][init]"<< " -- PDF -- "<<std::endl;
                    pdfi->Print("V");
                    std::cout<<"-----------------"<<std::endl;
                }
                
                if (verbose) { std::cout <<"coeff variables:";}
                std::set<std::string> servers= getServersVars(coeff);
                for(auto v : servers) {
                        if (logNormal.find(v) != logNormal.end()) {
                        logNormal.erase(v);
                        if (verbose) std::cout<<"|"<<v;
                        }
                }
                //--
                if (verbose) { std::cout<<std::endl <<"pdf variables:";}

                servers.clear();
                /*std::set<std::string>*/ servers= getServersVars(pdfi);
                for(auto v : servers) {
                        if (logNormal.find(v) != logNormal.end()) {
                        logNormal.erase(v);
                        if (verbose) std::cout<<"|"<<v;
                        }
                }
                if (verbose) { std::cout<<std::endl;}

                continue; //  if pn==nullptr or not weighted 
            } //  not a process normalization coefficient

            // TODO: remove all shapes

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
            }
            // remove all otherFactorList
            std::cout<<"[SimNLLDerivativesHelper][init]"<< ">> Removing all otherfactorsList "<<std::endl;
            for(int idx=0; idx< pn->otherFactorList_.getSize() ;++idx){
                RooAbsArg*v=pn->otherFactorList_.at(idx); 
                if (verbose) {
                    std::cout<<"[SimNLLDerivativesHelper][INFO]:"
                        <<"Removing "<<v->GetName()<<" from the list of known derivatives because appears in an otherFactors "
                        <<std::endl;
                }
                if (logNormal.find(v->GetName()) != logNormal.end()) logNormal.erase(v->GetName());
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

    }

    std::cout<<"[SimNLLDerivativesHelper][init]"<< "Constructing the derivatives sum "<<std::endl;
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
        }

        if (list.getSize() == 0){ // FIX for shapes TODO
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

            if(verbose) std::cout<<"[SimNLLDerivativesHelper][init]"<< "Loop. found constraint "<< nll_->constrainPdfsFast_[i]->getX().GetName() <<" == "<<name<<". Building RooFormulaVar" <<std::endl;
            RooFormulaVar *constr = new RooFormulaVar( (std::string("constr_")+name).c_str(),"@3*(@0-@1)/(@2*@2)",RooArgList(nll_->constrainPdfsFast_[i]->getX(),nll_->constrainPdfsFast_[i]->getMean(), nll_->constrainPdfsFast_[i]->getSigma(),unmask_));// (x-mean)/sigma^2
            list.add(*constr) ; // or addOwned or add?
        }
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
