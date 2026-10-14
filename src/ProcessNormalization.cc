#include "../interface/ProcessNormalization.h"

#include "../interface/CombineMathFuncs.h"

#include <cmath>
#include <cassert>
#include <cstdio>

ProcessNormalization::ProcessNormalization(const char *name, const char *title, double nominal) :
        RooAbsReal(name,title),
        nominalValue_(nominal),
        thetaList_("thetaList","List of nuisances for symmetric kappas", this), 
        asymmThetaList_("asymmThetaList","List of nuisances for asymmetric kappas", this), 
        otherFactorList_("otherFactorList","Other multiplicative terms", this)
{ 
}

ProcessNormalization::ProcessNormalization(const char *name, const char *title, RooAbsReal &nominal)
  : ProcessNormalization{name, title, 1.0}
{
   otherFactorList_.add(nominal);
}

ProcessNormalization::ProcessNormalization(const ProcessNormalization &other, const char *newname) :
        RooAbsReal(other, newname ? newname : other.GetName()),
        nominalValue_(other.nominalValue_),
        logKappa_(other.logKappa_),
        thetaList_("thetaList", this, other.thetaList_), 
        logAsymmKappa_(other.logAsymmKappa_),
        asymmThetaList_("asymmThetaList", this, other.asymmThetaList_), 
        otherFactorList_("otherFactorList", this, other.otherFactorList_)
{
}

void ProcessNormalization::addLogNormal(double kappa, RooAbsReal &theta) {
    if (kappa != 0.0 && kappa != 1.0) {
        logKappa_.push_back(std::log(kappa));
        thetaList_.add(theta);
    }
}

void ProcessNormalization::addAsymmLogNormal(double kappaLo, double kappaHi, RooAbsReal &theta) {
    if (fabs(kappaLo*kappaHi - 1) < 1e-5) {
        addLogNormal(kappaHi, theta);
    } else {
        logAsymmKappa_.push_back(std::make_pair(std::log(kappaLo), std::log(kappaHi)));
        asymmThetaList_.add(theta);
    }
}

void ProcessNormalization::addOtherFactor(RooAbsReal &factor) {
    otherFactorList_.add(factor);
}

void ProcessNormalization::fillAsymmKappaVecs() const
{
    if (logAsymmKappaLow_.size() != logAsymmKappa_.size()) {
       logAsymmKappaLow_.reserve(logAsymmKappa_.size());
       logAsymmKappaHigh_.reserve(logAsymmKappa_.size());
       for (auto [lo, hi] : logAsymmKappa_) {
          logAsymmKappaLow_.push_back(lo);
          logAsymmKappaHigh_.push_back(hi);
       }
    }
}

Double_t ProcessNormalization::evaluate() const
{
    thetaListVec_.resize(thetaList_.size());
    asymmThetaListVec_.resize(asymmThetaList_.size());
    otherFactorListVec_.resize(otherFactorList_.size());
    for (std::size_t i = 0; i < thetaList_.size(); ++i) {
        thetaListVec_[i] = static_cast<RooAbsReal const&>(thetaList_[i]).getVal();
    }
    for (std::size_t i = 0; i < asymmThetaList_.size(); ++i) {
        asymmThetaListVec_[i] = static_cast<RooAbsReal const&>(asymmThetaList_[i]).getVal();
    }
    for (std::size_t i = 0; i < otherFactorList_.size(); ++i) {
        otherFactorListVec_[i] = static_cast<RooAbsReal const&>(otherFactorList_[i]).getVal();
    }

    fillAsymmKappaVecs();
    return RooFit::Detail::MathFuncs::processNormalization(nominalValue_,
            thetaList_.size(), asymmThetaList_.size(), otherFactorList_.size(),
            thetaListVec_.data(), logKappa_.data(), asymmThetaListVec_.data(),
            logAsymmKappaLow_.data(), logAsymmKappaHigh_.data(),
            otherFactorListVec_.data());
}

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
// Override RooAbsReal::doEval to support RooFits vectorized likelihood
// evaluation. It's not strictly needed, as it can also fallback to
// RooAbsReal::evaluate(), but this comes at a performance price.
void ProcessNormalization::doEval(RooFit::EvalContext& ctx) const {
  thetaListVec_.resize(thetaList_.size());
  asymmThetaListVec_.resize(asymmThetaList_.size());
  otherFactorListVec_.resize(otherFactorList_.size());
  for (std::size_t i = 0; i < thetaList_.size(); ++i) {
    thetaListVec_[i] = ctx.at(&thetaList_[i])[0];
  }
  for (std::size_t i = 0; i < asymmThetaList_.size(); ++i) {
    asymmThetaListVec_[i] = ctx.at(&asymmThetaList_[i])[0];
  }
  for (std::size_t i = 0; i < otherFactorList_.size(); ++i) {
    otherFactorListVec_[i] = ctx.at(&otherFactorList_[i])[0];
  }

  fillAsymmKappaVecs();

  using RooFit::Detail::MathFuncs::processNormalization;

  double out = processNormalization(nominalValue_,
                                    thetaList_.size(),
                                    asymmThetaList_.size(),
                                    otherFactorList_.size(),
                                    thetaListVec_.data(),
                                    logKappa_.data(),
                                    asymmThetaListVec_.data(),
                                    logAsymmKappaLow_.data(),
                                    logAsymmKappaHigh_.data(),
                                    otherFactorListVec_.data());

  ctx.output()[0] = out;
}
#endif

void ProcessNormalization::dump() const {
    std::cout << "Dumping ProcessNormalization " << GetName() << " @ " << (void*)this << std::endl;
    std::cout << "\tnominal value: " << nominalValue_ << std::endl;
    std::cout << "\tlog-normals (" << logKappa_.size() << "):"  << std::endl;
    for (unsigned int i = 0; i < logKappa_.size(); ++i) {
        std::cout << "\t\t kappa = " << exp(logKappa_[i]) << ", logKappa = " << logKappa_[i] << 
                     ", theta = " << thetaList_.at(i)->GetName() << " = " << ((RooAbsReal*)thetaList_.at(i))->getVal() << std::endl;
    }
    std::cout << "\tasymm log-normals (" << logAsymmKappa_.size() << "):"  << std::endl;
    for (unsigned int i = 0; i < logAsymmKappa_.size(); ++i) {
        std::cout << "\t\t kappaLo = " << exp(logAsymmKappa_[i].first) << ", logKappaLo = " << logAsymmKappa_[i].first << 
                     ", kappaHi = " << exp(logAsymmKappa_[i].second) << ", logKappaHi = " << logAsymmKappa_[i].second << 
                     ", theta = " << asymmThetaList_.at(i)->GetName() << " = " << ((RooAbsReal*)asymmThetaList_.at(i))->getVal() << std::endl;
    }
    std::cout << "\tother terms (" << otherFactorList_.getSize() << "):"  << std::endl;
    for (int i = 0; i < otherFactorList_.getSize(); ++i) {  
        std::cout << "\t\t term " << otherFactorList_.at(i)->GetName() <<
                     " (class " << otherFactorList_.at(i)->ClassName() << 
                     "), value = " << ((RooAbsReal*)otherFactorList_.at(i))->getVal() << std::endl;
    }
    std::cout << std::endl;
}

std::vector<std::tuple<std::string,double,double>> ProcessNormalization::sigmaVariationsAll() const
{
    std::vector<std::tuple<std::string,double,double>> out;
    out.reserve(thetaList_.size() + asymmThetaList_.size());

    // Symmetric list: thetaList_[i] corresponds to logKappa_[i]
    for (std::size_t i = 0; i < logKappa_.size(); ++i) {
        const double lk = logKappa_[i];
        const double up = std::exp(lk);
        const double down = std::exp(-lk);
        const char *tname = thetaList_.at(i)->GetName();
        out.emplace_back(std::string(tname), down, up);
    }

    // Asymmetric list: asymmThetaList_[i] corresponds to logAsymmKappa_[i]
    for (std::size_t i = 0; i < logAsymmKappa_.size(); ++i) {
        const auto &p = logAsymmKappa_[i];
        const double down = std::exp(p.first);
        const double up   = std::exp(p.second);
        const char *tname = asymmThetaList_.at(i)->GetName();
        out.emplace_back(std::string(tname), down, up);
    }

    return out;
}






ClassImp(ProcessNormalization)


#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFitHS3/JSONIO.h>
#include <RooFit/Detail/JSONInterface.h>

#include "RooArgList.h"
#include "RooRealVar.h"

#include "static_execute.h"

namespace {

using RooFit::Detail::JSONNode;

class ProcessNormalizationFactory : public RooFit::JSONIO::Importer {
    public:
        bool importArg(RooJSONFactoryWSTool *tool, const JSONNode &p) const override
        {
            std::string name(RooJSONFactoryWSTool::name(p));

            if (!p.has_child("nominalValue")) {
                RooJSONFactoryWSTool::error("no nominalValue given in '" + name + "'");
            }
            double nominal_value(p["nominalValue"].val_double());
            tool->wsEmplace<ProcessNormalization>(name, nominal_value);
            RooArgList arg_list = tool->requestArgList<RooAbsReal>(p, "otherFactorList");
            RooWorkspace &ws = *tool->workspace();
            auto& process_normalization = static_cast<ProcessNormalization&>(*ws.function(name));
            for (auto* val : static_range_cast<RooAbsReal*>(arg_list)) {
                process_normalization.addOtherFactor(*val);
            } 
            return true;
        }
};

class ProcessNormalizationStreamer : public RooFit::JSONIO::Exporter {
    public:
        std::string const &key() const override {
	  static std::string key_ = "CMS::ProcessNormalization";
	  return key_;
	}
    
        bool exportObject(RooJSONFactoryWSTool *, const RooAbsArg *func, JSONNode &elem) const override
        {
            const ProcessNormalization *norm = static_cast<const ProcessNormalization *>(func);
            elem["type"] << key();

            elem["expression"] << norm->GetName();

            elem["nominalValue"] << norm->nominalValue();

            RooJSONFactoryWSTool::exportArray(norm->logKappa().size(), norm->logKappa().data(), elem["logKappa"]);

            RooJSONFactoryWSTool::fillSeq(elem["thetaList"], norm->thetaList());

            // Convert std::vector<std::pair<double, double>> to RooArgList
            // try elem["logKappa"].appendChild() << val; instead
            // Should be avoided with ROOT 6.32
            RooArgList logAsymmKappaList;
            for (const auto& pair : norm->logAsymmKappa()) {
                // Create RooRealVars for the pair
                RooRealVar* firstVal = new RooRealVar("first", "First value", pair.first);
                RooRealVar* secondVal = new RooRealVar("second", "Second value", pair.second);

                // Create a sublist for the pair
                RooArgList* pairList = new RooArgList();
                pairList->add(*firstVal);
                pairList->add(*secondVal);

                // Add the sublist to the main list
                logAsymmKappaList.add(*pairList);
            }
            RooJSONFactoryWSTool::fillSeq(elem["logAsymmKappa"], logAsymmKappaList);

            RooJSONFactoryWSTool::fillSeq(elem["asymmThetaList"], norm->asymmThetaList());

            RooJSONFactoryWSTool::fillSeq(elem["otherFactorList"], norm->otherFactorList());

            return true;
        }
};
}

STATIC_EXECUTE([]() {
    using namespace RooFit::JSONIO;

    // Register the importer and exporter
    registerImporter<ProcessNormalizationFactory>("CMS::ProcessNormalization", false);
    registerExporter<ProcessNormalizationStreamer>("ProcessNormalization", false);
});
