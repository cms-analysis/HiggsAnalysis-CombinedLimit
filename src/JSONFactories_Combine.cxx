#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFitHS3/JSONIO.h>
#include <RooFit/Detail/JSONInterface.h>
#include "../interface/static_execute.h"

#include "RooArgList.h"
#include "RooFit.h"
#include "RooRealVar.h"

#include "../interface/ProcessNormalization.h"

using RooFit::Detail::JSONNode;

class ProcessNormalizationFactory : public RooFit::JSONIO::Importer {
    public:
        bool importArg(RooJSONFactoryWSTool *tool, const JSONNode &p) const override
        {
            std::cout << "Importing ProcessNormalization" << std::endl;
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
        std::string const &key() const override;
        bool exportObject(RooJSONFactoryWSTool *, const RooAbsArg *func, JSONNode &elem) const override
        {
            std::cout << "Exporting ProcessNormalization" << std::endl;
            const ProcessNormalization *norm = static_cast<const ProcessNormalization *>(func);
            elem["type"] << key();

            elem["expression"] << norm->GetName();

            elem["nominalValue"] << norm->nominalValue();

            // Convert std::vector<double> to RooArgList
            // try elem["logKappa"].appendChild() << val; instead
            // Should be avoided with ROOT 6.32
            RooArgList logKappaList;
            for (double val : norm->logKappa()) {
                logKappaList.add(*new RooRealVar("", "", val));  // Create RooRealVar for each value
            }
            RooJSONFactoryWSTool::fillSeq(elem["logKappa"], logKappaList);

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

#define DEFINE_EXPORTER_KEY(class_name, name)    \
   std::string const &class_name::key() const    \
   {                                             \
      const static std::string keystring = name; \
      return keystring;                          \
   }
DEFINE_EXPORTER_KEY(ProcessNormalizationStreamer, "ProcessNormalization");

STATIC_EXECUTE([]() {
    using namespace RooFit::JSONIO;

    // Register the importer and exporter
    registerImporter<ProcessNormalizationFactory>("ProcessNormalization", false);
    registerExporter<ProcessNormalizationStreamer>(ProcessNormalization::Class(), false);
});
