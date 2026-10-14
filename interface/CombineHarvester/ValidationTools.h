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
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/Parameter.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/Systematic.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/Process.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/Logging.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/json.hpp"

namespace ch { class CombineHarvester; }

namespace ch {
using json = nlohmann::json;
void PrintSystematic(ch::Systematic *syst);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb, json &jsobj);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb);
void CheckEmptyShapes(CombineHarvester& cb, json &jsobj);
void CheckEmptyShapes(CombineHarvester& cb);
void CheckEmptyBins(CombineHarvester& cb, json &jsobj);
void CheckEmptyBins(CombineHarvester& cb);
bool HistErrorsAreSqrtN(const TH1* hist, double rate);
void CheckBinErrors(CombineHarvester& cb, double maxRelBinErr, json& jsobj);
void CheckBinErrors(CombineHarvester& cb, double maxRelBinErr);
void CheckNormEff(CombineHarvester& cb, double maxNormEff, json &jsobj);
void CheckNormEff(CombineHarvester& cb, double maxNormEff);
void CheckSizeOfShapeEffect(CombineHarvester& cb, json& jsobj);
void CheckSizeOfShapeEffect(CombineHarvester& cb);
void CheckSmallSignals(CombineHarvester& cb, double minSigFrac, json& jsobj);
void CheckSmallSignals(CombineHarvester& cb, double minSigFrac);
void ValidateShapeTemplates(CombineHarvester& cb, json &jsobj);
void ValidateShapeTemplates(CombineHarvester& cb);
void ValidateCards(CombineHarvester& cb, std::string const& filename, double maxNormEff, double minSigFrac, double maxRelBinErr);

}

#endif
