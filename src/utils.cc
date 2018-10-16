#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <typeinfo>
#include <stdexcept>

#include <TIterator.h>
#include <TString.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooStats/ModelConfig.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <regex>

#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logger.h"

using namespace std;

// This is needed to be able to factorize products in factorizeFunc, since RooProduct::components() strips duplicates
namespace {
    class RooProductWithAccessors : public RooProduct {
        public:
            RooProductWithAccessors(const RooProduct &other) :
                RooProduct(other) {}
            RooArgList realTerms() const { 
                RooArgList ret;
                RooFIter compRIter = _compRSet.fwdIterator() ;
                RooAbsReal* rcomp;
                while((rcomp=(RooAbsReal*)compRIter.next())) {
                   ret.add(*rcomp); 
                }
                return ret;
            }
    };
}

void utils::printRDH(RooAbsData *data) {
  std::vector<std::string> varnames, catnames;
  const RooArgSet *b0 = data->get();
  TIterator *iter = b0->createIterator();
  for (RooAbsArg *a = 0; (a = (RooAbsArg *)iter->Next()) != 0; ) {
    if (a->InheritsFrom("RooRealVar")) {
      varnames.push_back(a->GetName());
    } else if (a->InheritsFrom("RooCategory")) {
      catnames.push_back(a->GetName());
    }
  }
  delete iter;
  size_t nv = varnames.size(), nc = catnames.size();
  printf(" bin  ");
  for (size_t j = 0; j < nv; ++j) { printf("%16.16s  ", varnames[j].c_str()); }
  for (size_t j = 0; j < nc; ++j) { printf("%16.16s  ", catnames[j].c_str()); }
  printf("  weight\n");
  for (int i = 0, nb = data->numEntries(); i < nb; ++i) {
    const RooArgSet *bin = data->get(i);
    printf("%4d  ",i);
    for (size_t j = 0; j < nv; ++j) { printf("%16g  ",    bin->getRealValue(varnames[j].c_str())); }
    for (size_t j = 0; j < nc; ++j) { printf("%16.16s  ", bin->getCatLabel(catnames[j].c_str())); }
    printf("%12.7f\n", data->weight());
  }
}

void utils::printRAD(const RooAbsData *d) {
  if (d->InheritsFrom("RooDataHist") || d->numEntries() != 1) printRDH(const_cast<RooAbsData*>(d));
  else d->get(0)->Print("V");
}

void utils::printPdf(const RooAbsReal *pdf) {
  std::cout << "Pdf " << pdf->GetName() << " parameters." << std::endl;
  std::auto_ptr<RooArgSet> params(pdf->getVariables());
  params->Print("V");
}


void utils::printPdf(RooStats::ModelConfig &mc) {
  std::cout << "ModelConfig " << mc.GetName() << " (" << mc.GetTitle() << "): pdf parameters." << std::endl;
  std::auto_ptr<RooArgSet> params(mc.GetPdf()->getVariables());
  params->Print("V");
}

void utils::printPdf(RooWorkspace *w, const char *pdfName) {
  std::cout << "PDF " << pdfName << " parameters." << std::endl;
  std::auto_ptr<RooArgSet> params(w->pdf(pdfName)->getVariables());
  params->Print("V");
}

RooAbsPdf *utils::factorizePdf(const RooArgSet &observables, RooAbsPdf &pdf, RooArgList &constraints) {
    //assert(&pdf);
    const std::type_info & id = typeid(pdf);
    if (id == typeid(RooProdPdf)) {
        //std::cout << " pdf is product pdf " << pdf.GetName() << std::endl;
        RooProdPdf *prod = dynamic_cast<RooProdPdf *>(&pdf);
        RooArgList newFactors; RooArgSet newOwned;
        RooArgList list(prod->pdfList());
        bool needNew = false;
        for (int i = 0, n = list.getSize(); i < n; ++i) {
            RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
            RooAbsPdf *newpdf = factorizePdf(observables, *pdfi, constraints);
            //std::cout << "    for " << pdfi->GetName() << "   newpdf  " << (newpdf == 0 ? "null" : (newpdf == pdfi ? "old" : "new"))  << std::endl;
            if (newpdf == 0) { needNew = true; continue; }
            if (newpdf != pdfi) { needNew = true; newOwned.add(*newpdf); }
            newFactors.add(*newpdf);
        }
        if (!needNew && newFactors.getSize() > 1) { copyAttributes(pdf, *prod); return prod; }
        else if (newFactors.getSize() == 0) return 0;
        else if (newFactors.getSize() == 1) {
            RooAbsPdf *ret = (RooAbsPdf *) newFactors.first()->Clone(TString::Format("%s_obsOnly", pdf.GetName()));
            copyAttributes(pdf, *ret);
            return ret;
        }
        RooProdPdf *ret = new RooProdPdf(TString::Format("%s_obsOnly", pdf.GetName()), "", newFactors);
        ret->addOwnedComponents(newOwned);
        copyAttributes(pdf, *ret);
        return ret;
    } else if (id == typeid(RooSimultaneous) || id == typeid(RooSimultaneousOpt)) {
        RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(&pdf);
        RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
        int nbins = cat->numBins((const char *)0);
        TObjArray factorizedPdfs(nbins); RooArgSet newOwned;
        bool needNew = false;
        for (int ic = 0, nc = nbins; ic < nc; ++ic) {
            cat->setBin(ic);
            RooAbsPdf *pdfi = sim->getPdf(cat->getLabel());
            RooAbsPdf *newpdf = factorizePdf(observables, *pdfi, constraints);
            factorizedPdfs[ic] = newpdf;
            if (newpdf == 0) { throw std::runtime_error(std::string("ERROR: channel ") + cat->getLabel() + " factorized to zero."); }
            if (newpdf != pdfi) { needNew = true; newOwned.add(*newpdf); }
        }
        if (id == typeid(RooSimultaneousOpt)) {
            RooSimultaneousOpt &o = dynamic_cast<RooSimultaneousOpt &>(pdf);
            RooLinkedListIter iter = o.extraConstraints().iterator();
            if (o.extraConstraints().getSize() > 0) needNew = true;
            for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next()) {
                if (!constraints.contains(*a) && (!a->getAttribute("ignoreConstraint")) ) constraints.add(*a);
            }
        }
        RooSimultaneous *ret = sim;
        if (needNew) {
            if (id == typeid(RooSimultaneousOpt)) {
                ret = new RooSimultaneousOpt(TString::Format("%s_obsOnly", pdf.GetName()), "", (RooAbsCategoryLValue&) sim->indexCat());
            } else {
                ret = new RooSimultaneous(TString::Format("%s_obsOnly", pdf.GetName()), "", (RooAbsCategoryLValue&) sim->indexCat());
            }
            for (int ic = 0, nc = nbins; ic < nc; ++ic) {
                cat->setBin(ic);
                RooAbsPdf *newpdf = (RooAbsPdf *) factorizedPdfs[ic];
                if (newpdf) ret->addPdf(*newpdf, cat->getLabel());
            }
            ret->addOwnedComponents(newOwned);
        }
        delete cat;
        copyAttributes(pdf, *ret);
        return ret;
    } else if (pdf.dependsOn(observables)) {
        return &pdf;
    } else {
        if (!constraints.contains(pdf) && (!pdf.getAttribute("ignoreConstraint"))) constraints.add(pdf);
        return 0;
    }

}

void utils::factorizePdf(RooStats::ModelConfig &model, RooAbsPdf &pdf, RooArgList &obsTerms, RooArgList &constraints, bool debug) {
    return factorizePdf(*model.GetObservables(), pdf, obsTerms, constraints, debug);
}
void utils::factorizePdf(const RooArgSet &observables, RooAbsPdf &pdf, RooArgList &obsTerms, RooArgList &constraints, bool debug) {
    //assert(&pdf); should be safe
    const std::type_info & id = typeid(pdf);
    if (id == typeid(RooProdPdf)) {
        RooProdPdf *prod = dynamic_cast<RooProdPdf *>(&pdf);
        RooArgList list(prod->pdfList());
        for (int i = 0, n = list.getSize(); i < n; ++i) {
            RooAbsPdf *pdfi = (RooAbsPdf *) list.at(i);
            factorizePdf(observables, *pdfi, obsTerms, constraints);
        }
    } else if (id == typeid(RooSimultaneous) || id == typeid(RooSimultaneousOpt)) {
        if (id == typeid(RooSimultaneousOpt)) {
            RooSimultaneousOpt &o = dynamic_cast<RooSimultaneousOpt &>(pdf);
            RooLinkedListIter iter = o.extraConstraints().iterator();
            for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next()) {
                if (!constraints.contains(*a) && (!a->getAttribute("ignoreConstraint"))) constraints.add(*a);
            }
        }
        RooSimultaneous *sim  = dynamic_cast<RooSimultaneous *>(&pdf);
        RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
        for (int ic = 0, nc = cat->numBins((const char *)0); ic < nc; ++ic) {
            cat->setBin(ic);
            RooAbsPdf *pdfi = sim->getPdf(cat->getLabel());
            if (pdfi != 0) factorizePdf(observables, *pdfi, obsTerms, constraints);
        }
        delete cat;
    } else if (pdf.dependsOn(observables)) {
        if (!obsTerms.contains(pdf)) obsTerms.add(pdf);
    } else {
        if (!constraints.contains(pdf) && (!pdf.getAttribute("ignoreConstraint")) ) constraints.add(pdf);
    }
}


RooArgList utils::factors(const RooProduct &prod) {
    return ::RooProductWithAccessors(prod).realTerms();
}
void utils::factorizeFunc(const RooArgSet &observables, RooAbsReal &func, RooArgList &obsTerms, RooArgList &constraints, bool keepDuplicate, bool debug) {
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(&func);
    if (pdf != 0) { 
        factorizePdf(observables, *pdf, obsTerms, constraints, debug); 
        return; 
    }
    const std::type_info & id = typeid(func);
    if (id == typeid(RooProduct)) {
        RooProduct *prod = dynamic_cast<RooProduct *>(&func);
        RooArgList components(utils::factors(*prod));
        //std::cout << "Function " << func.GetName() << " is a RooProduct with " << components.getSize() << " components." << std::endl;
        std::auto_ptr<TIterator> iter(components.createIterator());
        for (RooAbsReal *funci = (RooAbsReal *) iter->Next(); funci != 0; funci = (RooAbsReal *) iter->Next()) {
            //std::cout << "  component " << funci->GetName() << " of type " << funci->ClassName() << "(dep obs? " << funci->dependsOn(observables) << ")" << std::endl;
            factorizeFunc(observables, *funci, obsTerms, constraints, true);
        }
    } else if (func.dependsOn(observables)) {
        if (!obsTerms.contains(func) || keepDuplicate) obsTerms.add(func);
    } else {
        if (( !constraints.contains(func) && (!func.getAttribute("ignoreConstraint")) ) || keepDuplicate) constraints.add(func);
    }
}

RooAbsPdf *utils::makeNuisancePdf(RooStats::ModelConfig &model, const char *name) { 
    return utils::makeNuisancePdf(*model.GetPdf(), *model.GetObservables(), name);
}

RooAbsPdf *utils::makeNuisancePdf(RooAbsPdf &pdf, const RooArgSet &observables, const char *name) { 
    //assert(&pdf);
    RooArgList obsTerms, constraints;
    factorizePdf(observables, pdf, obsTerms, constraints);
    if (constraints.getSize() == 0) return 0;
    return new RooProdPdf(name,"", constraints);
}

RooAbsPdf *utils::makeObsOnlyPdf(RooStats::ModelConfig &model, const char *name) { 
    RooArgList obsTerms, constraints;
    factorizePdf(model, *model.GetPdf(), obsTerms, constraints);
    return new RooProdPdf(name,"", obsTerms);
}

RooAbsPdf *utils::fullClonePdf(const RooAbsPdf *pdf, RooArgSet &holder, bool cloneLeafNodes) {
  // Clone all FUNC compents by copying all branch nodes
  RooArgSet tmp("RealBranchNodeList") ;
  pdf->branchNodeServerList(&tmp);
  tmp.snapshot(holder, cloneLeafNodes); 
  // Find the top level FUNC in the snapshot list
  return (RooAbsPdf*) holder.find(pdf->GetName());
}
RooAbsReal *utils::fullCloneFunc(const RooAbsReal *pdf, RooArgSet &holder, bool cloneLeafNodes) {
  // Clone all FUNC compents by copying all branch nodes
  RooArgSet tmp("RealBranchNodeList") ;
  pdf->branchNodeServerList(&tmp);
  tmp.snapshot(holder, cloneLeafNodes); 
  // Find the top level FUNC in the snapshot list
  return (RooAbsReal*) holder.find(pdf->GetName());
}

RooAbsReal *utils::fullCloneFunc(const RooAbsReal *pdf, const RooArgSet &obs, RooArgSet &holder, bool cloneLeafNodes) {
  // Clone all FUNC compents by copying all branch nodes
  RooArgSet tmp("RealBranchNodeList"), toClone;
  pdf->branchNodeServerList(&tmp);
  unsigned int nitems = tmp.getSize();
  RooFIter iter = tmp.fwdIterator();
  for (RooAbsArg *a = iter.next(); a != 0; a = iter.next()) {
      if (a == pdf) toClone.add(*a);
      else if (a->dependsOn(obs)) toClone.add(*a);
  }
  unsigned int nobsitems = toClone.getSize();
  toClone.snapshot(holder, cloneLeafNodes); 
  if (runtimedef::get("fullCloneFunc_VERBOSE")) std::cout << "For PDF " << pdf->GetName() << ", cloned " << nobsitems << "/" << nitems << " items" << std::endl;
  // Find the top level FUNC in the snapshot list
  return (RooAbsReal*) holder.find(pdf->GetName());
}




void utils::getClients(const RooAbsCollection &values, const RooAbsCollection &allObjects, RooAbsCollection &clients) {
    std::auto_ptr<TIterator> iterAll(allObjects.createIterator());
    std::auto_ptr<TIterator> iterVal(values.createIterator());
    for (RooAbsArg *v = (RooAbsArg *) iterVal->Next(); v != 0; v = (RooAbsArg *) iterVal->Next()) {
        if (typeid(*v) != typeid(RooRealVar)) continue;
        std::auto_ptr<TIterator> clientIter(v->clientIterator());
        for (RooAbsArg *a = (RooAbsArg *) clientIter->Next(); a != 0; a = (RooAbsArg *) clientIter->Next()) {
            if (allObjects.containsInstance(*a) && !clients.containsInstance(*a)) clients.add(*a);
        }
    }
}

bool utils::setAllConstant(const RooAbsCollection &coll, bool constant) {
    bool changed = false;
    std::auto_ptr<TIterator> iter(coll.createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        RooRealVar *v = dynamic_cast<RooRealVar *>(a);
        RooCategory *cv = dynamic_cast<RooCategory *>(a);
        if (v && (v->isConstant() != constant)) {
            changed = true;
            v->setConstant(constant);
        }
        else if (cv && (cv->isConstant() != constant)) {
            changed = true;
            cv->setConstant(constant);
        }        
    }
    return changed;
}

TString utils::printRooArgAsString(RooAbsArg *a){
    RooRealVar *v = dynamic_cast<RooRealVar *>(a);
    if (v){
      return TString::Format("%s = %g %s",v->GetName(), v->getVal(), v->isConstant() ? "(constant)":"");
    } else {
      RooCategory *cv = dynamic_cast<RooCategory *>(a);
      if (cv){
        return TString::Format("%s = %d %s",cv->GetName(), cv->getIndex(), cv->isConstant() ? "(constant)":"") ;
      }
    }
    return "";
}

bool utils::checkModel(const RooStats::ModelConfig &model, bool throwOnFail) {
    bool ok = true; std::ostringstream errors; 
    std::auto_ptr<TIterator> iter;
    RooAbsPdf *pdf = model.GetPdf(); if (pdf == 0) throw std::invalid_argument("Model without Pdf");
    RooArgSet allowedToFloat; 
    if (model.GetObservables() == 0) { 
        ok = false; errors << "ERROR: model does not define observables.\n"; 
        std::cout << errors.str() << std::endl;
        if (throwOnFail) throw std::invalid_argument(errors.str()); else return false; 
    } else {
        allowedToFloat.add(*model.GetObservables());
    }
    if (model.GetParametersOfInterest() == 0) { 
        ok = false; errors << "ERROR: model does not define parameters of interest.\n";  
    } else {
        iter.reset(model.GetParametersOfInterest()->createIterator());
        for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
            RooRealVar *v = dynamic_cast<RooRealVar *>(a);
            if (!v) { ok = false; errors << "ERROR: parameter of interest " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; continue; }
            if (v->isConstant()) { ok = false; errors << "ERROR: parameter of interest " << a->GetName() << " is constant\n"; continue; }
            if (!pdf->dependsOn(*v)) { ok = false; errors << "ERROR: pdf does not depend on parameter of interest " << a->GetName() << "\n"; continue; }
            allowedToFloat.add(*v);
        }
    }
    if (model.GetNuisanceParameters() != 0) { 
        iter.reset(model.GetNuisanceParameters()->createIterator());
        for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
            RooRealVar *v = dynamic_cast<RooRealVar *>(a);
            if (!v) { ok = false; errors << "ERROR: nuisance parameter " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; continue; }
            if (v->isConstant()) { ok = false; errors << "ERROR: nuisance parameter " << a->GetName() << " is constant\n"; continue; }
            if (!pdf->dependsOn(*v)) { errors << "WARNING: pdf does not depend on nuisance parameter " << a->GetName() << "\n"; continue; }
            allowedToFloat.add(*v);
        }
    }
    if (model.GetGlobalObservables() != 0) { 
        iter.reset(model.GetGlobalObservables()->createIterator());
        for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
            RooRealVar *v = dynamic_cast<RooRealVar *>(a);
            if (!v) { ok = false; errors << "ERROR: global observable " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; continue; }
            if (!v->isConstant()) { ok = false; errors << "ERROR: global observable " << a->GetName() << " is not constant\n"; continue; }
            if (!pdf->dependsOn(*v)) { errors << "WARNING: pdf does not depend on global observable " << a->GetName() << "\n"; continue; }
        }
    }
    std::auto_ptr<RooArgSet> params(pdf->getParameters(*model.GetObservables()));
    iter.reset(params->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        if (a->getAttribute("flatParam") && a->isConstant()) {
            ok = false; errors << "ERROR: parameter " << a->GetName() << " is declared as flatParam but is constant.\n";
        }
        if (a->isConstant() || allowedToFloat.contains(*a)) continue;
        if (a->getAttribute("flatParam")) {
            RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
            if (rrv->getVal() > rrv->getMax() || rrv->getVal() < rrv->getMin()) {
                ok = false; errors << "ERROR: flatParam " << rrv->GetName() << " has a value " << rrv->getVal() << 
                                      " outside of the defined range [" << rrv->getMin() << ", " << rrv->getMax() << "]\n";
            }
        } else {
            errors << "WARNING: pdf parameter " << a->GetName() << " (type " << a->ClassName() << ") is not allowed to float (it's not nuisance, poi, observable or global observable\n"; 
        }
        RooRealVar *v = dynamic_cast<RooRealVar *>(a);
        if (v != 0 && !v->isConstant() && (v->getVal() < v->getMin() || v->getVal() > v->getMax())) {
            ok = false; errors << "ERROR: parameter " << a->GetName() << " has value " << v->getVal() << " outside range [ " << v->getMin() << " , " << v->getMax() << " ]\n"; 
        }
    }
    iter.reset();

    if (model.GetNuisanceParameters() != 0) {
        RooArgList constraints;
        std::auto_ptr<RooAbsPdf> factorizedPdf(utils::factorizePdf(*model.GetObservables(), *model.GetPdf(), constraints));
        if (factorizedPdf.get() == 0 || factorizedPdf.get() == model.GetPdf()) {
            factorizedPdf.release();
            ok = false; errors << "ERROR: have nuisance parameters, but can't factorize the pdf\n";
        }
        std::auto_ptr<RooArgSet> obsParams(factorizedPdf->getParameters(*model.GetObservables()));
        iter.reset(model.GetNuisanceParameters()->createIterator());
        for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
            if (!obsParams->contains(*a))  {
                errors << "WARNING: model pdf does not depend on nuisace parameter " << a->GetName() << "\n";
            }
        }
    }
    std::cout << errors.str() << std::endl;
    if (!ok && throwOnFail) throw std::invalid_argument(errors.str()); 
    return ok;
}

RooSimultaneous * utils::rebuildSimPdf(const RooArgSet &observables, RooSimultaneous *sim) {
    RooArgList constraints;
    RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
    int nbins = cat->numBins((const char *)0);
    TObjArray factorizedPdfs(nbins); 
    RooArgSet newOwned;
    for (int ic = 0, nc = nbins; ic < nc; ++ic) {
        cat->setBin(ic);
        RooAbsPdf *pdfi = sim->getPdf(cat->getLabel());
        if (pdfi == 0) { factorizedPdfs[ic] = 0; continue; }
        RooAbsPdf *newpdf = factorizePdf(observables, *pdfi, constraints);
        factorizedPdfs[ic] = newpdf;
        if (newpdf == 0) { continue; }
        if (newpdf != pdfi) { newOwned.add(*newpdf);  }
    }
    RooSimultaneous *ret = new RooSimultaneous(TString::Format("%s_reloaded", sim->GetName()), "", (RooAbsCategoryLValue&) sim->indexCat());
    for (int ic = 0, nc = nbins; ic < nc; ++ic) {
        cat->setBin(ic);
        RooAbsPdf *newpdf = (RooAbsPdf *) factorizedPdfs[ic];
        if (newpdf) {
            if (constraints.getSize() > 0) {
                RooArgList allfactors(constraints); allfactors.add(*newpdf);
                RooProdPdf *newerpdf = new RooProdPdf(TString::Format("%s_plus_constr", newpdf->GetName()), "", allfactors);
                ret->addPdf(*newerpdf, cat->getLabel());
                copyAttributes(*newpdf, *newerpdf);
                newOwned.add(*newerpdf);
            } else {
                ret->addPdf(*newpdf, cat->getLabel());
            }
        }
    }
    ret->addOwnedComponents(newOwned);
    copyAttributes(*sim, *ret);
    delete cat;
    return ret;
}

void utils::copyAttributes(const RooAbsArg &from, RooAbsArg &to) {
    if (&from == &to) return;
    const std::set<std::string> attribs = from.attributes();
    if (!attribs.empty()) {
        for (std::set<std::string>::const_iterator it = attribs.begin(), ed = attribs.end(); it != ed; ++it) to.setAttribute(it->c_str());
    }
    const std::map
<std::string, std::string> strattribs = from.stringAttributes();
    if (!strattribs.empty()) {
        for (std::map
        <std::string,std::string>::const_iterator it = strattribs.begin(), ed = strattribs.end(); it != ed; ++it) to.setStringAttribute(it->first.c_str(), it->second.c_str());
    }
}

void utils::guessChannelMode(RooSimultaneous &simPdf, RooAbsData &simData, bool verbose) 
{
    RooAbsCategoryLValue &cat = const_cast<RooAbsCategoryLValue &>(simPdf.indexCat());
    TList *split = simData.split(cat, kTRUE);
    for (int i = 0, n = cat.numBins((const char *)0); i < n; ++i) {
        cat.setBin(i);
        RooAbsPdf *pdf = simPdf.getPdf(cat.getLabel());
        if (pdf->getAttribute("forceGenBinned") || pdf->getAttribute("forceGenUnbinned")) {
            if (verbose) std::cout << " - " << cat.getLabel() << " has generation mode already set" << std::endl;
            continue;
        }
        RooAbsData *spl = (RooAbsData *) split->FindObject(cat.getLabel());
        if (spl == 0) { 
            if (verbose) std::cout << " - " << cat.getLabel() << " has no dataset, cannot guess" << std::endl; 
            continue;
        }
        if (spl->numEntries() != spl->sumEntries()) {
            if (verbose) std::cout << " - " << cat.getLabel() << " has " << spl->numEntries() << " num entries of sum " << spl->sumEntries() << ", mark as binned" << std::endl;
            pdf->setAttribute("forceGenBinned");
        } else {
            if (verbose) std::cout << " - " << cat.getLabel() << " has " << spl->numEntries() << " num entries of sum " << spl->sumEntries() << ", mark as unbinned" << std::endl;
            pdf->setAttribute("forceGenUnbinned");
        }
    }
}

void
utils::setChannelGenModes(RooSimultaneous &simPdf, const std::string &binned, const std::string &unbinned, int verbose) {
    std::regex r_binned(binned, std::regex::ECMAScript);
    std::regex r_unbinned(unbinned, std::regex::ECMAScript);
    std::smatch match;

    RooAbsCategoryLValue &cat = const_cast<RooAbsCategoryLValue &>(simPdf.indexCat());
    for (int i = 0, n = cat.numBins((const char *)0); i < n; ++i) {
        cat.setBin(i);
        const std::string & label = cat.getLabel();
        RooAbsPdf *pdf = simPdf.getPdf(label.c_str());
        if (!binned.empty() && std::regex_match(label, match, r_binned)) {
            if (pdf->getAttribute("forceGenUnbinned")) {
                if (verbose) std::cout << "Overriding generation mode for " << pdf->GetName() << " in " << label << " from unbinned to binned." << std::endl;
                pdf->setAttribute("forceGenUnbinned", false);
            } else {
                if (verbose) std::cout << "Setting generation mode for " << pdf->GetName() << " in " << label << " to binned." << std::endl;
            }
            pdf->setAttribute("forceGenBinned");
        }
        if (!unbinned.empty() && std::regex_match(label, match, r_unbinned)) {
            if (pdf->getAttribute("forceGenBinned")) {
                if (verbose) std::cout << "Overriding generation mode for " << pdf->GetName() << " in " << label << " from binned to unbinned." << std::endl;
                pdf->setAttribute("forceGenBinned", false);
            } else {
                if (verbose) std::cout << "Setting generation mode for " << pdf->GetName() << " in " << label << " to unbinned." << std::endl;
            }
            pdf->setAttribute("forceGenUnbinned");
        }
        if (!pdf->getAttribute("forceGenUnbinned") && !pdf->getAttribute("forceGenBinned")) {
            if (verbose > -1) std::cout << "Warning: pdf generation mode for " << pdf->GetName() << " in " << label << " is not set" << std::endl;
        }
    }
}
TGraphAsymmErrors * utils::makeDataGraph(TH1 * dataHist, bool asDensity)
{
    // Properly normalise 
    TGraphAsymmErrors * dataGraph = new TGraphAsymmErrors(dataHist->GetNbinsX());

    dataHist->Sumw2(false);
    dataHist->SetBinErrorOption(TH1::kPoisson);

    for (int b=1;b <= dataHist->GetNbinsX();b++){
	double yv = dataHist->GetBinContent(b);
	double bw = dataHist->GetBinWidth(b);

	//double up; 
	//double dn;

	//RooHistError::instance().getPoissonInterval(yv,dn,up,1);

	//double errlow  = (yv-dn);
	//double errhigh = (up-yv);
    //
    double errlow = dataHist->GetBinErrorLow(b);
    double errhigh = dataHist->GetBinErrorUp(b);


	if (asDensity) {
		yv/=bw;
		errlow/=bw;
		errhigh/=bw;
	}

	dataGraph->SetPoint(b-1,dataHist->GetBinCenter(b),yv);
	dataGraph->SetPointError(b-1,bw/2,bw/2,errlow,errhigh);
    }

    return dataGraph;
}

std::vector<RooPlot *>
utils::makePlots(const RooAbsPdf &pdf, const RooAbsData &data, const char *signalSel, const char *backgroundSel, float rebinFactor, RooFitResult *fitRes) {
    std::vector<RooPlot *> ret;
    RooArgList constraints;
    RooAbsPdf *facpdf = factorizePdf(*data.get(0), const_cast<RooAbsPdf &>(pdf), constraints);
    
    const std::type_info & id = typeid(*facpdf);
    if (id == typeid(RooSimultaneous) || id == typeid(RooSimultaneousOpt)) {
        const RooSimultaneous *sim  = dynamic_cast<const RooSimultaneous *>(&pdf);
        const RooAbsCategoryLValue &cat = (RooAbsCategoryLValue &) sim->indexCat();
        TList *datasets = data.split(cat, true);
        TIter next(datasets);
        for (RooAbsData *ds = (RooAbsData *) next(); ds != 0; ds = (RooAbsData *) next()) {
            RooAbsPdf *pdfi  = sim->getPdf(ds->GetName());
            std::auto_ptr<RooArgSet> obs(pdfi->getObservables(ds));
            if (obs->getSize() == 0) break;
	    TIterator *obs_iter = obs->createIterator();
	    std::cout << " PDF CHECKING " << std::endl; 
	    pdfi->Print("v");
	    ds->Print("v");
	    std::cout << " ------------ " << std::endl; 

	    //for (int iobs=0;iobs<obs->getSize();iobs++){
  	    for (RooAbsArg *a = 0; (a = (RooAbsArg *)obs_iter->Next()) != 0; ) {
	      RooRealVar *x = dynamic_cast<RooRealVar *>(a);
	      if (x == 0) continue;
	      int nbins = x->numBins(); if (nbins == 0) nbins = 100;
	      if (nbins/rebinFactor > 6) nbins = ceil(nbins/rebinFactor);
	      ret.push_back(x->frame(RooFit::Title(ds->GetName()), RooFit::Bins(nbins)));
	      ret.back()->SetName(Form("%s_%s",ds->GetName(),x->GetName()));
	      ds->plotOn(ret.back(), RooFit::DataError(RooAbsData::Poisson));
	      if (fitRes)pdfi->plotOn(ret.back(),RooFit::Normalization(pdfi->expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent),RooFit::VisualizeError(*fitRes,1) ,RooFit::FillColor(kOrange));
	      if (signalSel && strlen(signalSel))         pdfi->plotOn(ret.back(), RooFit::LineColor(209), RooFit::Components(signalSel),RooFit::Normalization(pdfi->expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent));
	      if (signalSel && strlen(signalSel))         pdfi->plotOn(ret.back(), RooFit::LineColor(209), RooFit::Components(signalSel));
	      if (backgroundSel && strlen(backgroundSel)) pdfi->plotOn(ret.back(), RooFit::LineColor(206), RooFit::Components(backgroundSel),RooFit::Normalization(pdfi->expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent));
	      std::cout << "[utils::makePlots] Number of events for pdf in " << ret.back()->GetName() << ", pdf " << pdfi->GetName() << " = " << pdfi->expectedEvents(RooArgSet(*x)) << std::endl;  
	      pdfi->plotOn(ret.back(),RooFit::Normalization(pdfi->expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent));
	      ds->plotOn(ret.back(), RooFit::DataError(RooAbsData::Poisson));
	    }
            delete ds;
        }
        delete datasets;
    } else if (pdf.canBeExtended()) {
        std::auto_ptr<RooArgSet> obs(pdf.getObservables(&data));

	//for (int iobs=0;iobs<obs->getSize();iobs++){
	TIterator *obs_iter = obs->createIterator();
  	for (RooAbsArg *a = 0; (a = (RooAbsArg *)obs_iter->Next()) != 0; ) {
          RooRealVar *x = dynamic_cast<RooRealVar *>(a);
	  if (x != 0) {
	      ret.push_back(x->frame());
	      ret.back()->SetName(Form("data_%s",x->GetName()));
	      data.plotOn(ret.back(), RooFit::DataError(RooAbsData::Poisson));
	      if (fitRes) pdf.plotOn(ret.back(),RooFit::Normalization(pdf.expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent),RooFit::VisualizeError(*fitRes,1) ,RooFit::FillColor(kOrange));
	      if (signalSel && strlen(signalSel))         pdf.plotOn(ret.back(), RooFit::LineColor(209), RooFit::Components(signalSel),RooFit::Normalization(pdf.expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent));
	      if (backgroundSel && strlen(backgroundSel)) pdf.plotOn(ret.back(), RooFit::LineColor(206), RooFit::Components(backgroundSel),RooFit::Normalization(pdf.expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent));
	      std::cout << "[utils::makePlots] Number of events for pdf in " << ret.back()->GetName() << ", pdf " << pdf.GetName() << " = " << pdf.expectedEvents(RooArgSet(*x)) << std::endl;  
	      pdf.plotOn(ret.back(),RooFit::Normalization(pdf.expectedEvents(RooArgSet(*x)),RooAbsReal::NumEvent),RooFit::VisualizeError(*fitRes,1));
	      data.plotOn(ret.back(), RooFit::DataError(RooAbsData::Poisson));
	  }
	}
    }
    if (facpdf != &pdf) { delete facpdf; }
    return ret;

}

void utils::CheapValueSnapshot::readFrom(const RooAbsCollection &src) {
    if (&src != src_) {
        src_ = &src;
        values_.resize(src.getSize());
    }
    RooLinkedListIter iter = src.iterator(); int i = 0;
    for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next(), ++i) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv == 0) {
          RooCategory *rc = dynamic_cast<RooCategory *>(a);
	  if (rc == 0){
		throw std::invalid_argument("Collection to read from contains a non-RooRealVar/RooCategory");
	  } 
	  values_[i] = (double)rc->getIndex();
	} else {
          values_[i] = rrv->getVal();
	}
    }
}

void utils::CheapValueSnapshot::writeTo(const RooAbsCollection &src) const {
    if (&src == src_) {
        RooLinkedListIter iter = src.iterator();  int i = 0;
        for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next(), ++i) {
            RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
	    if (rrv!=0) rrv->setVal(values_[i]);
	    else {
		RooCategory *rc = dynamic_cast<RooCategory *>(a);
		rc->setIndex((int)values_[i]);
	    }
        }
    } else {
        RooLinkedListIter iter = src_->iterator();  int i = 0;
        for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next(), ++i) {
            RooAbsArg *a2 = src.find(a->GetName()); if (a2 == 0) continue;
            RooRealVar *rrv = dynamic_cast<RooRealVar *>(a2);
            if (rrv!=0) rrv->setVal(values_[i]);
	    else {
		RooCategory *rc = dynamic_cast<RooCategory *>(a2);
		rc->setIndex((int)values_[i]);
	    }
        }
    }
}

void utils::CheapValueSnapshot::Print(const char *fmt) const {
    if (src_ == 0) { printf("<NIL>\n"); return; }
    if (fmt[0] == 'V') {
        RooLinkedListIter iter = src_->iterator(); int i = 0;
        for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next(), ++i) {
            printf(" %3d) %-30s = %9.6g\n", i, a->GetName(), values_[i]);
        }
        printf("\n");
    } else {
        src_->Print(fmt);
    }
}

void utils::setModelParameters( const std::string & setPhysicsModelParameterExpression, const RooArgSet & params) {

  vector<string> SetParameterExpressionList;  
  boost::split(SetParameterExpressionList, setPhysicsModelParameterExpression, boost::is_any_of(","));
  for (UInt_t p = 0; p < SetParameterExpressionList.size(); ++p) {
    vector<string> SetParameterExpression;
    boost::split(SetParameterExpression, SetParameterExpressionList[p], boost::is_any_of("="));
      
    if (SetParameterExpression.size() != 2) {
      std::cout << "Error parsing physics model parameter expression : " << SetParameterExpressionList[p] << endl;
    } 
    // check for regex syntax: rgx{regex}                                                                                                                                     
    else if (boost::starts_with(SetParameterExpression[0], "rgx{") && boost::ends_with(SetParameterExpression[0], "}")) {

        std::string reg_esp = SetParameterExpression[0].substr(4, SetParameterExpression[0].size()-5);
        std::cout<<"interpreting "<<reg_esp<<" as regex "<<std::endl;
        std::regex rgx( reg_esp, std::regex::ECMAScript);

        std::auto_ptr<TIterator> iter(params.createIterator());
        for (RooAbsArg *tmp = (RooAbsArg*) iter->Next(); tmp != 0; tmp = (RooAbsArg*) iter->Next()) {

            bool isrvar = tmp->IsA()->InheritsFrom(RooRealVar::Class());  // check its type    

            if (isrvar) {
                RooRealVar *tmpParameter = dynamic_cast<RooRealVar *>(tmp);
                const std::string &target = tmpParameter->GetName();
                std::smatch match;
                if (std::regex_match(target, match, rgx)) {
                    double PhysicsParameterValue = atof(SetParameterExpression[1].c_str());
                    cout << "Set Default Value of Parameter " << target
                     << " To : " << PhysicsParameterValue << "\n";
                    tmpParameter->setVal(PhysicsParameterValue);
                }
            } else {
                RooCategory *tmpCategory  = dynamic_cast<RooCategory*>(tmp);
                const std::string &target = tmpCategory->GetName();
                std::smatch match;
                if (std::regex_match(target, match, rgx)) {
                    int PhysicsParameterValue = atoi(SetParameterExpression[1].c_str());
                    cout << "Set Default Index of Parameter " << target
                         << " To : " << PhysicsParameterValue << "\n";
                    tmpCategory->setIndex(PhysicsParameterValue);
                }
            }
        }

    }
    else {
      
      RooAbsArg  *tmp = (RooAbsArg*)params.find(SetParameterExpression[0].c_str());     
      if (tmp){

        bool isrvar = tmp->IsA()->InheritsFrom(RooRealVar::Class());  // check its type

        if (isrvar) {
          RooRealVar *tmpParameter = dynamic_cast<RooRealVar*>(tmp);
          double PhysicsParameterValue = atof(SetParameterExpression[1].c_str());
          cout << "Set Default Value of Parameter " << SetParameterExpression[0] 
               << " To : " << PhysicsParameterValue << "\n";
          tmpParameter->setVal(PhysicsParameterValue);
        } else {
          RooCategory *tmpCategory  = dynamic_cast<RooCategory*>(tmp);
          int PhysicsParameterValue = atoi(SetParameterExpression[1].c_str());
          cout << "Set Default Index of Parameter " << SetParameterExpression[0] 
               << " To : " << PhysicsParameterValue << "\n";
          tmpCategory->setIndex(PhysicsParameterValue);
	}
      }
        else {
        std::cout << "Warning: Did not find a parameter with name " << SetParameterExpression[0] << endl;
      }
    }
  }

}

void utils::setModelParameterRanges( const std::string & setPhysicsModelParameterRangeExpression, const RooArgSet & params) {

  vector<string> SetParameterRangeExpressionList;  
  boost::split(SetParameterRangeExpressionList, setPhysicsModelParameterRangeExpression, boost::is_any_of(":"));
  for (UInt_t p = 0; p < SetParameterRangeExpressionList.size(); ++p) {
    vector<string> SetParameterRangeExpression;
    boost::split(SetParameterRangeExpression, SetParameterRangeExpressionList[p], boost::is_any_of("=,"));
      
    if (SetParameterRangeExpression.size() != 3) {
      std::cout << "Error parsing physics model parameter expression : " << SetParameterRangeExpressionList[p] << endl;
    } else if (SetParameterRangeExpression[1] == "-inf" && (
                    SetParameterRangeExpression[2] == "inf" || 
                    SetParameterRangeExpression[2] == "+inf") ) {
      RooRealVar *tmpParameter = (RooRealVar*)params.find(SetParameterRangeExpression[0].c_str());            
      if (tmpParameter) {
        cout << "Leave Parameter " << SetParameterRangeExpression[0] 
             << " freely floating, with no range\n";
        tmpParameter->removeRange();
      } else {
        std::cout << "Warning: Did not find a parameter with name " << SetParameterRangeExpression[0] << endl;
      }
 
    } else {
      double PhysicsParameterRangeLow = atof(SetParameterRangeExpression[1].c_str());
      double PhysicsParameterRangeHigh = atof(SetParameterRangeExpression[2].c_str());

      // check for regex syntax: rgx{regex}
      if (boost::starts_with(SetParameterRangeExpression[0], "rgx{") && boost::ends_with(SetParameterRangeExpression[0], "}")) {
          
          std::string reg_esp = SetParameterRangeExpression[0].substr(4, SetParameterRangeExpression[0].size()-5);
          std::cout<<"interpreting "<<reg_esp<<" as regex "<<std::endl;
          std::regex rgx( reg_esp, std::regex::ECMAScript);

          std::auto_ptr<TIterator> iter(params.createIterator());
          for (RooAbsArg *a = (RooAbsArg*) iter->Next(); a != 0; a = (RooAbsArg*) iter->Next()) {
              RooRealVar *tmpParameter = dynamic_cast<RooRealVar *>(a);
              const std::string &target = tmpParameter->GetName();
              std::smatch match;
              if (std::regex_match(target, match, rgx)) {
                  if (tmpParameter->isConstant()) continue;
                  cout << "Set Range of Parameter " << target
                       << " To : (" << PhysicsParameterRangeLow << "," << PhysicsParameterRangeHigh << ")\n";
                  tmpParameter->setRange(PhysicsParameterRangeLow,PhysicsParameterRangeHigh);
              }
          }

      }         
      else {
          RooRealVar *tmpParameter = (RooRealVar*)params.find(SetParameterRangeExpression[0].c_str());            
          if (tmpParameter) {
              cout << "Set Range of Parameter " << SetParameterRangeExpression[0] 
                   << " To : (" << PhysicsParameterRangeLow << "," << PhysicsParameterRangeHigh << ")\n";
              tmpParameter->setRange(PhysicsParameterRangeLow,PhysicsParameterRangeHigh);
          } else {
              std::cout << "Warning: Did not find a parameter with name " << SetParameterRangeExpression[0] << endl;
          }
      }
    }

  }
}

void utils::createSnapshotFromString( const std::string expression, const RooArgSet &allvars, RooArgSet &output, const char *context) {
    if (expression.find("=") == std::string::npos) {
        if (allvars.getSize() != 1) throw std::invalid_argument(std::string("Error: the argument to ")+context+" is a single value, but there are multiple variables to choose from");
        allvars.snapshot(output);
        errno = 0; // check for errors in str->float conversion
        ((RooRealVar*)output.first())->setVal(strtod(expression.c_str(),NULL));
        if (errno != 0) std::invalid_argument(std::string("Error: the argument to ")+context+" is not a valid number.");
    } else {
        std::string::size_type eqidx = 0, colidx = 0, colidx2;
        do {
            eqidx   = expression.find("=", colidx);
            colidx2 = expression.find(",", colidx+1);
            if (eqidx == std::string::npos || (colidx2 != std::string::npos && colidx2 < eqidx)) {
                throw std::invalid_argument(std::string("Error: the argument to ")+context+" is not in the forms 'value' or 'name1=value1,name2=value2,...'\n");
            }
            std::string poiName = expression.substr(colidx, eqidx-colidx);
            std::string poiVal  = expression.substr(eqidx+1, (colidx2 == std::string::npos ? std::string::npos : colidx2 - eqidx - 1));
            RooAbsArg *poi = allvars.find(poiName.c_str());
            if (poi == 0) throw std::invalid_argument(std::string("Error: unknown parameter '")+poiName+"' passed to "+context+".");
            output.addClone(*poi);
            errno = 0;
            output.setRealValue(poi->GetName(), strtod(poiVal.c_str(),NULL));
            if (errno != 0) throw std::invalid_argument(std::string("Error: invalid value '")+poiVal+"' for parameter '"+poiName+"' passed to "+context+".");
            colidx = colidx2+1;
        } while (colidx2 != std::string::npos);
    }
}

void utils::reorderCombinations(std::vector<std::vector<int> > &vec, const std::vector<int> &max, const std::vector<int> &pos){
	
    // Assuming everyone starts from 0, add the start position modulo number of positions
    std::vector<std::vector<int> >::iterator vit = vec.begin();
    for (;vit!=vec.end();vit++){
      std::vector<int>::iterator pit=(*vit).begin();
      int cp=0;
      for (;pit!=(*vit).end();pit++){
  	int newpos = (*pit+pos[cp])%max[cp];
	*pit=newpos;
	cp++;
      }
    }
}

std::vector<std::vector<int> > utils::generateOrthogonalCombinations(const std::vector<int> &vec) {
    int n = vec.size();
    std::vector< std::vector<int> > result;
    std::vector<int> idx(n,0);
    result.push_back(idx);
    for (int i = 0; i < n; ++i) {    
      std::vector<int> idx(n,0);
      bool exit_loop = false;
      while(exit_loop == false){
	if (idx[i]+1==vec[i]){
		exit_loop=true;
	} else {
	   ++(idx[i]);
	   result.push_back(idx);
	}
      }                  
    }                        
  return result;                                    
}


std::vector<std::vector<int> > utils::generateCombinations(const std::vector<int> &vec) {
    // Generate all combinations choosing from N sets each of which have Mn values
    int n = vec.size();
    std::vector<int> idx(n,0);
    std::vector< std::vector<int> > result;
    result.push_back(idx);
    bool exit_loop = false;
    while(exit_loop == false){
      // Find the first index we can increment (if possible)
      for (int i = 0; i < n; ++i) {
        if ((idx[i]+1) == vec[i]) {
          if (i != n-1) {
            idx[i] = 0;
          } else {
            // we're done
            exit_loop = true;
            break;
          }
        } else {
          ++(idx[i]);
          result.push_back(idx);
          break;
        }
      }
    }

  return result;
}


bool utils::isParameterAtBoundary( const RooRealVar &param ){

    double vMin = param.getMin();
    double vMax = param.getMax();
    double val = param.getVal();
    double errLo = -1.0 * param.getErrorLo(); 
    double errHi = param.getErrorHi(); 

    double pullMin = (val-vMin) / (errLo);
    double pullMax = (vMax-val) / (errHi);

    float nSigma=1.0;

    if(pullMin < nSigma || pullMax < nSigma){
        return true;
    }
    
    return false;
}


bool utils::anyParameterAtBoundaries( const RooArgSet &params, int verbosity ){

    static std::unordered_map<std::string, unsigned char> timesFoundAtBoundary;
    bool isAnyBad = false;

    RooLinkedListIter iter = params.iterator(); int i = 0;
    for (RooRealVar *a = (RooRealVar *) iter.Next(); a != 0; a = (RooRealVar *) iter.Next(), ++i) {

        bool isBad = isParameterAtBoundary(*a);

        if(isBad){
            std::string varName((*a).GetName());

            if( verbosity >= 9 || (timesFoundAtBoundary[varName] < 3 && verbosity > -1) ){
                fprintf(CloseCoutSentry::trueStdOutGlobal(),"  [WARNING] Found [%s] at boundary. \n", (*a).GetName());
                std::cout << "       "; (*a).Print();
            }
	    
            if( verbosity > 0 ){
            	Logger::instance().log(std::string(Form("utils.cc: %d -- Found parameter %s at boundary (within ~1sigma): %g+/-%g",__LINE__,(*a).GetName(),(*a).getVal(),(*a).getError())),Logger::kLogLevelInfo,__func__);
	    }

            timesFoundAtBoundary[varName]++;
        }

        isAnyBad |= isBad;
    }

    // for( std::unordered_map<std::string, unsigned char>::value_type e : timesFoundAtBoundary ){
    //     printf("e %s -> %i\n", e.first.c_str(), e.second);
    // }
    
    return isAnyBad;
}

int utils::countFloating(const RooArgSet &params){
	int count=0;
        RooLinkedListIter iter = params.iterator(); int i = 0;
        for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next(), ++i) {
		if (!a->isConstant()) count++;
        }
	return count;
}

RooArgSet utils::returnAllVars(RooWorkspace *w){
	// Helper function to retun *all* vars, including RooCategories in workspace
	RooArgSet args(w->allVars());
	args.add(w->allCats());
	return args;
}
