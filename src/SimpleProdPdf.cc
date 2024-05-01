#include "../interface/SimpleProdPdf.h"
#include "../interface/RooCheapProduct.h"

#include <string>
#include <memory>
#include <stdexcept>
#include "RooHistPdf.h"

SimpleProdPdf::SimpleProdPdf() : RooAbsPdf() {}

SimpleProdPdf::SimpleProdPdf(const char* name, const char* title, RooArgList const& pdfList)
    : RooAbsPdf(name, title), _pdfList("!pdfs", "List of PDFs", this) {
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

SimpleProdPdf::SimpleProdPdf(const char* name, const char* title, RooArgList& pdfList,
                             std::vector<int> const& pdfSettings)
    : RooAbsPdf(name, title), _pdfList("!pdfs", "List of PDFs", this) {
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
      _pdfSettings(other._pdfSettings) {}

Double_t SimpleProdPdf::evaluate() const {
    Double_t value = 1.0;
    for (int i = 0; i < _pdfList.getSize(); ++i) {
        value *= static_cast<RooAbsPdf*>(_pdfList.at(i))->getVal();
    }
    return value;
}

RooAbsReal* SimpleProdPdf::createIntegral(const RooArgSet& iset, const RooArgSet* nset,
                           const RooNumIntConfig* cfg, const char* rangeName) const {
    RooArgList terms;

    std::vector<RooAbsPdf*> pdfVec;
    std::vector<std::unique_ptr<RooArgSet>> iVarsVec;
    for (int i = 0; i < _pdfList.getSize(); ++i) {
      RooAbsPdf *iPdf = static_cast<RooAbsPdf*>(_pdfList.at(i));
      std::unique_ptr<RooArgSet> iVars(iPdf->getVariables());
      pdfVec.push_back(iPdf);
      iVarsVec.emplace_back(static_cast<RooArgSet*>(iset.selectCommon(*iVars)));
    }
    for (unsigned i = 0; i < pdfVec.size(); ++i) {
      if (_pdfSettings.at(i) == 1) {
        for (unsigned j = 0; j < pdfVec.size(); ++j) {
          if (j == i) continue;
          iVarsVec.at(i)->remove(*iVarsVec.at(j));
        }
      }
    }
    for (unsigned i = 0; i < pdfVec.size(); ++i) {
        TString intName = pdfVec.at(i)->GetName();
        intName.Append(this->integralNameSuffix(*iVarsVec.at(i),nset,rangeName));
        RooAbsArg *prexistingInt = nullptr;
        if (pdfVec.at(i)->ownedComponents() != nullptr) {
          prexistingInt = pdfVec.at(i)->ownedComponents()->find(intName);
          if (prexistingInt) {
            // std::cout << "Found pre-existing integral " << prexistingInt->GetName() << ", recycling...\n";
            terms.add(*prexistingInt);
          }
        }
        if (prexistingInt == nullptr) {
          RooAbsReal *iInt = pdfVec.at(i)->createIntegral(*iVarsVec.at(i), nset, cfg, rangeName);
          iInt->SetName(intName);
          if (dynamic_cast<RooHistPdf*>(pdfVec.at(i))) {
            pdfVec.at(i)->addOwnedComponents(*iInt);
          }
          terms.add(*iInt);
        }
    }
    TString name(this->GetName()) ;
    name.Append(this->integralNameSuffix(iset,nset,rangeName));
    RooCheapProduct *prod = new RooCheapProduct(name, "", terms, true);
    return prod;
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
