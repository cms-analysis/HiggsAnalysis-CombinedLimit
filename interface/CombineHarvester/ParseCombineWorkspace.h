#ifndef CombineTools_ParseCombineWorkspace_h
#define CombineTools_ParseCombineWorkspace_h
#include <string>
#include "RooWorkspace.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/CombineHarvester.h"

class RooAddPdf;
class RooAbsPdf;

namespace ch {
void ParseCombineWorkspace(CombineHarvester& cb, RooWorkspace& ws,
                           std::string const& modelcfg,
                           std::string const& data, bool verbose = false);

void ParseCombineWorkspacePy(CombineHarvester& cb, RooWorkspace const& ws,
                           std::string const& modelcfg,
                           std::string const& data, bool verbose = false);

RooAbsReal* FindAddPdf(RooAbsReal* input);
}
#endif
