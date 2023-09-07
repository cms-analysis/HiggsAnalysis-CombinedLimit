#include "HiggsAnalysis/CombinedLimit/interface/RooEFTScalingFunction.h"

ClassImp(RooEFTScalingFunction)

RooEFTScalingFunction::RooEFTScalingFunction(const char *name, const char *title, const std::map<std::string,double> &coeffs, const RooArgList &terms) :
    RooAbsReal(name,title),
    coeffs_(coeffs),
    terms_("!terms","RooArgList of Wilson coefficients",this),
    offset_(1.0)
{

    // Add Wilson coefficients to terms_ container
    RooFIter iter = terms.fwdIterator();
    for( RooAbsArg *a = iter.next(); a != 0; a = iter.next()) {
        RooAbsReal *rar = dynamic_cast<RooAbsReal *>(a);
        if (!rar) {
            throw std::invalid_argument(std::string("Term ")+a->GetName()+" of RooEFTScalingFunction is a "+a->ClassName());
        }
        terms_.add(*rar);
    }

    // Loop over elements in mapping: add components to vector depending on string
    for( auto const& x : coeffs ) {
        TString term_name = x.first;
        double term_prefactor = x.second;

        if( term_name.Contains("_") ) {
            TString first_term = dynamic_cast<TObjString *>( term_name.Tokenize("_")->At(0))->GetString();
            TString second_term = dynamic_cast<TObjString *>( term_name.Tokenize("_")->At(1))->GetString();

            // Squared-quadratic components
            if( second_term == "2" ){
                if( terms.find(first_term) ){
                    std::vector<RooAbsReal *> vterms;
                    RooAbsReal *rar1 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
                    RooAbsReal *rar2 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
                    vterms.push_back(rar1);
                    vterms.push_back(rar2);
                    vcomponents_.emplace(vterms,term_prefactor);
                }
            } else{
                // Cross-quadratic components
                if( terms.find(first_term) && terms.find(second_term) ){
                    std::vector<RooAbsReal *> vterms;
                    RooAbsReal *rar1 = dynamic_cast<RooAbsReal *>( terms.find(first_term) );
                    RooAbsReal *rar2 = dynamic_cast<RooAbsReal *>( terms.find(second_term) );
                    vterms.push_back(rar1);
                    vterms.push_back(rar2);
                    vcomponents_.emplace(vterms,term_prefactor);
                }
            }
        } else {

          if( terms.find(term_name) ) {
              std::vector<RooAbsReal *> vterms;
              RooAbsReal *rar = dynamic_cast<RooAbsReal *>( terms.find(term_name) );
              vterms.push_back(rar);
              vcomponents_.emplace(vterms,term_prefactor);
          } 
       }
    }
}

RooEFTScalingFunction::RooEFTScalingFunction(const RooEFTScalingFunction& other, const char* name) :
    RooAbsReal(other, name),
    coeffs_(other.coeffs_),
    terms_("!terms",this,other.terms_),
    vcomponents_(other.vcomponents_),
    offset_(other.offset_)
{
}

Double_t RooEFTScalingFunction::evaluate() const 
{
    if (vcomponents_.empty()) {
        return offset_;
    }
    double ret = offset_;
    for (auto const &x : vcomponents_) {
        double res = x.second;
        for( auto const &y : x.first ){
            res *= y->getVal();
        }
        ret += res;
    }
    return ret;
}
