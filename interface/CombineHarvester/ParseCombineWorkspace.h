#ifndef CombineTools_CopyTools_h
#define CombineTools_CopyTools_h
#include <string>
#include "RooWorkspace.h"
#include "CombineTools/interface/CombineHarvester.h"

class RooAddPdf;
class RooAbsPdf;

namespace ch {
void ParseCombineWorkspace(CombineHarvester& cb, RooWorkspace& ws,
                           std::string const& modelcfg,
                           std::string const& data, bool verbose = false);

RooAddPdf* FindAddPdf(RooAbsPdf* input);
}
#endif
