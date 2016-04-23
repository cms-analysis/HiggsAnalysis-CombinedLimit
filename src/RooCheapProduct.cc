#include "HiggsAnalysis/CombinedLimit/interface/RooCheapProduct.h"
#include <RooConstVar.h>

RooCheapProduct::RooCheapProduct(const char *name, const char *title, const RooArgList &terms, bool pruneConstants) :
    RooAbsReal(name,title),
    terms_("!terms","Set of real product components",this),
    offset_(1.0)
{
    RooFIter iter = terms.fwdIterator();
    for (RooAbsArg *a = iter.next(); a != 0; a = iter.next()) {
        RooAbsReal *rar = dynamic_cast<RooAbsReal *>(a);
        if (!rar) { 
            throw std::invalid_argument(std::string("Component ")+a->GetName()+" of RooCheapProduct is a "+a->ClassName()); 
        }
        if (typeid(*rar) == typeid(RooConstVar) || (pruneConstants && rar->isConstant())) {
            offset_ *= rar->getVal();
            continue; 
        }
        terms_.add(*rar);
        vterms_.push_back(rar);
    }
}

RooCheapProduct::RooCheapProduct(const RooCheapProduct& other, const char* name) :
    RooAbsReal(other, name),
    terms_("!terms",this,other.terms_),
    vterms_(other.vterms_),
    offset_(other.offset_)
{
}

Double_t RooCheapProduct::evaluate() const 
{
    if (vterms_.empty()) {
        RooFIter iter = terms_.fwdIterator();
        std::vector<RooAbsReal *> & vterms = const_cast<std::vector<RooAbsReal *>&>(vterms_);
        vterms.reserve(terms_.getSize());
        for (RooAbsArg *a = iter.next(); a != 0; a = iter.next()) {
            vterms.push_back(dynamic_cast<RooAbsReal *>(a));
        }
    }
    double ret = offset_;
    for (const RooAbsReal *rar : vterms_) {
        ret *= rar->getVal();
    }
    return ret;
}
