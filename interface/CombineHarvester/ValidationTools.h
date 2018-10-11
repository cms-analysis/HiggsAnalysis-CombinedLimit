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

namespace ch { class CombineHarvester; }

namespace ch {
void PrintSystematic(ch::Systematic *syst);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb);
void CheckEmptyShapes(CombineHarvester& cb);

}

#endif
