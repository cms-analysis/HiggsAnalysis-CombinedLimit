#include "../interface/SimpleProdPdf.h"
#include "../interface/RooCheapProduct.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooHistPdf.h"
#include "RooNameReg.h"

#include <ROOT/RConfig.hxx>  // for ROOT_VERSION

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TString.h"

#include <memory>
#include <stdexcept>
#include <vector>

SimpleProdPdf::SimpleProdPdf() : _cacheMgr(this, 10) {}

SimpleProdPdf::SimpleProdPdf(const char* name, const char* title, RooArgList const& pdfList)
    : RooAbsPdf(name, title), _pdfList("!pdfs", "List of PDFs", this), _cacheMgr(this, 10) {
  for (RooAbsArg *arg : pdfList) {
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(arg);
    if (!pdf) {
      coutW(InputArguments) << "SimpleProdPdf::SimpleProdPdf(" << GetName() << ") list arg "
                            << arg->GetName() << " is not a PDF, ignored" << std::endl;
      continue;
    }
    _pdfList.add(*pdf);
    _pdfSettings.push_back(0);
  }
}

SimpleProdPdf::SimpleProdPdf(const char* name,
                             const char* title,
                             RooArgList& pdfList,
                             std::vector<int> const& pdfSettings)
    : RooAbsPdf(name, title), _pdfList("!pdfs", "List of PDFs", this), _cacheMgr(this, 10) {
  int iSetting = -1;
  for (RooAbsArg *arg : pdfList) {
    ++iSetting;
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(arg);
    if (!pdf) {
      coutW(InputArguments) << "SimpleProdPdf::SimpleProdPdf(" << GetName() << ") list arg "
                            << arg->GetName() << " is not a PDF, ignored" << std::endl;
      continue;
    }
    _pdfList.add(*pdf);
    _pdfSettings.push_back(iSetting);
  }
}

SimpleProdPdf::SimpleProdPdf(const char* name, const char* title, RooProdPdf const& prodPdf)
    : SimpleProdPdf(name, title, prodPdf.pdfList()) {
  for (int i = 0; i < _pdfList.getSize(); ++i) {
    RooArgSet* nset = prodPdf.findPdfNSet(*static_cast<RooAbsPdf*>(_pdfList.at(i)));
    if (nset && TString(nset->GetName()) == "nset" && nset->getSize() > 0) {
      _pdfSettings.at(i) = 1;
    }
  }
}

SimpleProdPdf::SimpleProdPdf(const SimpleProdPdf& other, const char* name)
    : RooAbsPdf(other, name),
      _pdfList("!pdfs", this, other._pdfList),
      _pdfSettings(other._pdfSettings),
      _cacheMgr(other._cacheMgr, this) {}

Double_t SimpleProdPdf::evaluate() const {
    Double_t value = 1.0;
    for (int i = 0; i < _pdfList.getSize(); ++i) {
        value *= static_cast<RooAbsPdf*>(_pdfList.at(i))->getVal();
    }
    return value;
}

bool SimpleProdPdf::forceAnalyticalInt(const RooAbsArg& /*dep*/) const { return true; }

Int_t SimpleProdPdf::getAnalyticalIntegralWN(RooArgSet& allVars,
                                             RooArgSet& analVars,
                                             const RooArgSet* normSet,
                                             const char* rangeName) const {
  if (_forceNumInt)
    return 0;
  // All observables are integrated as a product of the integrals of the components
  analVars.add(allVars);
  return getIntegralCode(allVars, normSet, rangeName) + 1;
}

double SimpleProdPdf::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const {
  if (code == 0)
    return getVal(normSet);

  auto* cache = dynamic_cast<CacheElem*>(_cacheMgr.getObjByIndex(code - 1));
  if (cache == nullptr) {
    // The cache was sterilized (e.g. after a server redirection), so rebuild
    // the integral for this code from the sets of variables stored by name.
    std::unique_ptr<RooArgSet> vars{getParameters(RooArgSet())};
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 26, 0)
    RooArgSet nset = _cacheMgr.selectFromSet1(*vars, code - 1);
    RooArgSet iset = _cacheMgr.selectFromSet2(*vars, code - 1);
#else
    RooArgSet nset, iset;
    nset.add(*std::unique_ptr<RooArgSet>{_cacheMgr.nameSet1ByIndex(code - 1)->select(*vars)});
    iset.add(*std::unique_ptr<RooArgSet>{_cacheMgr.nameSet2ByIndex(code - 1)->select(*vars)});
#endif
    Int_t code2 = getIntegralCode(iset, nset.empty() ? nullptr : &nset, rangeName);
    cache = dynamic_cast<CacheElem*>(_cacheMgr.getObjByIndex(code2));
  }
  return cache->integral->getVal();
}

Int_t SimpleProdPdf::getIntegralCode(const RooArgSet& iset, const RooArgSet* nset, const char* rangeName) const {
  Int_t sterileIdx = -1;
  if (_cacheMgr.getObj(nset, &iset, &sterileIdx, rangeName)) {
    return _cacheMgr.lastIndex();
  }
  auto cache = std::make_unique<CacheElem>();
  cache->integral = makeIntegral(iset, nset, rangeName);
  return _cacheMgr.setObj(nset, &iset, cache.release(), RooNameReg::ptr(rangeName));
}

// Build the integral over iset as the product of the integrals of the components
std::unique_ptr<RooAbsReal> SimpleProdPdf::makeIntegral(const RooArgSet& iset,
                                                        const RooArgSet* nset,
                                                        const char* rangeName) const {
  RooArgList terms;

  std::vector<RooAbsPdf*> pdfVec;
  std::vector<std::unique_ptr<RooArgSet>> iVarsVec;
  for (int i = 0; i < _pdfList.getSize(); ++i) {
    auto* iPdf = dynamic_cast<RooAbsPdf*>(_pdfList.at(i));
    std::unique_ptr<RooArgSet> iVars(iPdf->getVariables());
    pdfVec.push_back(iPdf);
    // The cast is only redundant since ROOT 6.34, before that selectCommon() returned a RooAbsCollection*
    iVarsVec.emplace_back(static_cast<RooArgSet*>(iset.selectCommon(*iVars)));  // NOLINT(readability-redundant-casting)
  }
  for (unsigned i = 0; i < pdfVec.size(); ++i) {
    if (_pdfSettings.at(i) == 1) {
      for (unsigned j = 0; j < pdfVec.size(); ++j) {
        if (j == i)
          continue;
        iVarsVec.at(i)->remove(*iVarsVec.at(j));
      }
    }
  }
  for (unsigned i = 0; i < pdfVec.size(); ++i) {
    TString intName = pdfVec.at(i)->GetName();
    intName.Append(this->integralNameSuffix(*iVarsVec.at(i), nset, rangeName));
    RooAbsArg* prexistingInt = nullptr;
    if (pdfVec.at(i)->ownedComponents() != nullptr) {
      prexistingInt = pdfVec.at(i)->ownedComponents()->find(intName);
      if (prexistingInt) {
        // std::cout << "Found pre-existing integral " << prexistingInt->GetName() << ", recycling...\n";
        terms.add(*prexistingInt);
      }
    }
    if (prexistingInt == nullptr) {
      RooAbsReal* iInt = pdfVec.at(i)->createIntegral(*iVarsVec.at(i), nset, nullptr, rangeName);
      iInt->SetName(intName);
      if (dynamic_cast<RooHistPdf*>(pdfVec.at(i))) {
        pdfVec.at(i)->addOwnedComponents(*iInt);
      }
      terms.add(*iInt);
    }
  }
  TString name(this->GetName());
  name.Append(this->integralNameSuffix(iset, nset, rangeName));
  return std::make_unique<RooCheapProduct>(name, /*title=*/"", terms, /*pruneConstants=*/true);
}

void SimpleProdPdf::printMetaArgs(std::ostream& os) const

{
  Bool_t first(kTRUE);

  for (int i = 0; i < _pdfList.getSize(); ++i) {
    if (!first) {
      os << " * ";
    } else {
      first = kFALSE;
    }
    os << _pdfList.at(i)->GetName();
  }
  os << " ";
}

std::list<Double_t>* SimpleProdPdf::binBoundaries(RooAbsRealLValue& obs, Double_t xlo,
                                                  Double_t xhi) const {
  for (RooAbsArg *pdf : _pdfList) {
    std::list<Double_t>* hint = static_cast<RooAbsPdf*>(pdf)->binBoundaries(obs, xlo, xhi);
    if (hint) {
      return hint;
    }
  }

  return 0;
}

ClassImp(SimpleProdPdf)
