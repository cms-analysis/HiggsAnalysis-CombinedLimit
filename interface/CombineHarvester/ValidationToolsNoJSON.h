#ifndef CombineTools_ValidationTools_h
#define CombineTools_ValidationTools_h
/* Alternate version of ValidationTools.h
 * with no forward declaration of JSON-enabled methods
 * The problem is that the type declaration of
 * nlohmann::json clashes with another inside ROOT
 * source code "include/ROOT/REveElement.hxx" and so
 * cppyy spits lots of errors because both are in global
 * namespace :(
 */
#include <string>
#include "CombineHarvester/CombineTools/interface/Systematic.h"

namespace ch { class CombineHarvester; }

namespace ch {
void PrintSystematic(ch::Systematic *syst);
void ValidateShapeUncertaintyDirection(CombineHarvester& cb);
void CheckEmptyShapes(CombineHarvester& cb);
void CheckEmptyBins(CombineHarvester& cb);
void CheckNormEff(CombineHarvester& cb, double maxNormEff);
void CheckSizeOfShapeEffect(CombineHarvester& cb);
void CheckSmallSignals(CombineHarvester& cb, double minSigFrac);
void ValidateShapeTemplates(CombineHarvester& cb);
void ValidateCards(CombineHarvester& cb, std::string const& filename, double maxNormEff, double minSigFrac);

}

#endif
