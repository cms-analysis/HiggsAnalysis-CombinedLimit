#include "CombineHarvester/CombineTools/interface/CopyTools.h"
#include <string>
#include <vector>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"


namespace ch {
void SplitSyst(ch::CombineHarvester& cb, std::string const& syst_in,
               std::string const& split1, std::string const& split2,
               double val1, double val2) {
  CloneSysts(cb.cp().syst_name({syst_in}), cb, [&](Systematic *s) {
    s->set_name(split1);
    s->set_value_u(val1);
  });
  CloneSysts(cb.cp().syst_name({syst_in}), cb, [&](Systematic *s) {
    s->set_name(split2);
    s->set_value_u(val2);
  });
  cb.syst_name({syst_in}, false);
}
}
