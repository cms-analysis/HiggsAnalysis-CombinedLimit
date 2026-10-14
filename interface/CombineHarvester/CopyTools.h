#ifndef CombineTools_CopyTools_h
#define CombineTools_CopyTools_h
#include <string>
#include <map>
#include <vector>
#include <set>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Object.h"

namespace ch {

/**
 * Duplicate Observation objects in one instance, which are then modified by a
 * user-supplied function, before being inserted into another instance
 *
 * All entries in **src** will be copied into **dest** (which could well be
 * the same instance). The third argument is a function or callable that
 * accepts a single `ch::Observation *` argument, and that will be applied for
 * each duplicated object. 
 */
template <typename Function>
void CloneObs(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachObs([&](ch::Observation *obs) {
    ch::Observation cpy = *obs;
    func(&cpy);
    dest.InsertObservation(cpy);
  });
}

/**
 * Duplicate Process objects in one instance, which are then modified by a
 * user-supplied function, before being inserted into another instance
 *
 * All entries in **src** will be copied into **dest** (which could well be
 * the same instance). The third argument is a function or callable that
 * accepts a single `ch::Process *` argument, and that will be applied for
 * each duplicated object.
 */
template <typename Function>
void CloneProcs(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachProc([&](ch::Process *proc) {
    ch::Process cpy = *proc;
    func(&cpy);
    dest.InsertProcess(cpy);
  });
}

/**
 * Duplicate Systematic objects in one instance, which are then modified by a
 * user-supplied function, before being inserted into another instance
 *
 * All entries in **src** will be copied into **dest** (which could well be
 * the same instance). The third argument is a function or callable that
 * accepts a single `ch::Systematic *` argument, and that will be applied for
 * each duplicated object.
 */
template <typename Function>
void CloneSysts(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  src.ForEachSyst([&](ch::Systematic *sys) {
    ch::Systematic cpy = *sys;
    func(&cpy);
    dest.InsertSystematic(cpy);
  });
}

/**
 * Duplicate the Process and Systematic objects in one instance, which are
 * then modified by a user-supplied function, before being inserted into
 * another instance
 *
 * All entries in **src** will be copied into **dest** (which could well be
 * the same instance). The third argument is a function or callable that
 * accepts a single `ch::Object *` argument, and that will be applied for each
 * duplicated object.
 *
 * This function is useful for making a full copy of an existing process. As
 * an example, in the Run I SM HTT analysis a search for a second Higgs boson
 * was performed. This meant taking the mH=125 GeV signal process and
 * converting it into a background for all signal mass hypotheses. This can be
 * achieved with
 * 
 *     ch::CloneProcsAndSysts(cb.cp().signals().mass({"125"}), cb,
 *                       [](ch::Object* p) {
 *      p->set_process(p->process() + "_SM125");
 *      p->set_signal(false);
 *      p->set_mass("*");
 *     });
 */
template <typename Function>
void CloneProcsAndSysts(ch::CombineHarvester& src, ch::CombineHarvester& dest,
                Function func) {
  CloneProcs(src, dest, func);
  CloneSysts(src, dest, func);
}

/**
 * Replace all instances of an existing systematic with two copies having
 * modified values
 *
 * The function copies every Systematic object with a given name twice,
 * modifies the name and value of each copy, before removing all instances of
 * the original.
 *
 * @param cb CombineHarvester instance to modify
 * @param syst_in name of systematic that should be removed and split into two
 * @param split1 the new name for the first duplicate
 * @param split2 the new name for the second duplicate
 * @param val1 the new uncertainty value for the first duplicate
 * @param val2 the new uncertainty value for the second duplicate
 */
void SplitSyst(ch::CombineHarvester& cb, std::string const& syst_in,
               std::string const& split1, std::string const& split2,
               double val1, double val2);
}

#endif
