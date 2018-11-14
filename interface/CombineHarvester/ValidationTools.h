#ifndef CombineTools_ValidationTools_h
#define CombineTools_ValidationTools_h
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "TGraph.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"
#include "CombineHarvester/CombineTools/interface/json.hpp"

namespace ch { class CombineHarvester; }

namespace ch {
using json = nlohmann::json;
void PrintSystematic(ch::Systematic *syst);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb, json &jsobj);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb);
void CheckEmptyShapes(CombineHarvester& cb, json &jsobj);
void CheckEmptyShapes(CombineHarvester& cb);
void CheckNormEff(CombineHarvester& cb, double maxNormEff, json &jsobj);
void CheckNormEff(CombineHarvester& cb, double maxNormEff);
void CheckSizeOfShapeEffect(CombineHarvester& cb, json& jsobj);
void CheckSizeOfShapeEffect(CombineHarvester& cb);
void ValidateCards(CombineHarvester& cb, std::string const& filename, double maxNormEff);

}

#endif
