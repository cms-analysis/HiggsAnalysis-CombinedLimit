#ifndef CombineTools_CopyTools_h
#define CombineTools_CopyTools_h
#include <string>
#include <map>
#include <vector>
#include <set>
#include "CombineTools/interface/CombineHarvester.h"
#include "CombineTools/interface/Object.h"

namespace ch {

template <typename Function>
void CloneObs(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachObs([&](ch::Observation *obs) {
    ch::Observation cpy = *obs;
    func(&cpy);
    dest.InsertObservation(cpy);
  });
}

template <typename Function>
void CloneProcs(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachProc([&](ch::Process *proc) {
    ch::Process cpy = *proc;
    func(&cpy);
    dest.InsertProcess(cpy);
  });
}

template <typename Function>
void CloneSysts(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachSyst([&](ch::Systematic *sys) {
    ch::Systematic cpy = *sys;
    func(&cpy);
    dest.InsertSystematic(cpy);
  });
}

template <typename Function>
void CloneProcsAndSysts(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  CloneProcs(src, dest, func);
  CloneSysts(src, dest, func);
}

void SplitSyst(ch::CombineHarvester& cb, std::string const& syst_in,
               std::string const& split1, std::string const& split2,
               double val1, double val2);
}

#endif
