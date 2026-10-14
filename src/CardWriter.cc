#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include "boost/format.hpp"
#include "CombineHarvester/CombineTools/interface/Logging.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"

namespace ch {

CardWriter::CardWriter(std::string const& text_pattern,
                       std::string const& root_pattern)
    : text_pattern_(text_pattern),
      root_pattern_(root_pattern),
      wildcard_masses_({"*"}),
      v_(0),
      create_dirs_(true) {}

CardWriter& CardWriter::SetVerbosity(unsigned v) {
  v_ = v;
  return *this;
}

CardWriter& CardWriter::CreateDirectories(bool flag) {
  create_dirs_ = flag;
  return *this;
}


CardWriter& CardWriter::SetWildcardMasses(
    std::vector<std::string> const& masses) {
  wildcard_masses_ = masses;
  return *this;
}

auto CardWriter::BuildMap(std::string const& pattern,
                          ch::CombineHarvester& cmb) const -> PatternMap {
  PatternMap f_map;
  // We first filter Objects having a mass value in the wildcard list
  cmb.cp().mass(wildcard_masses_, false)
    .ForEachObj([&](ch::Object const* obj) {
      // Build the fully-compiled key
      std::string key = Compile(pattern, obj);
      if (f_map.count(key)) return;
      std::set<std::string> mappings;
      // The fully-compiled pattern always goes in
      mappings.insert(key);
      // Compile again, but this time skipping the $MASS substitution
      std::string maps_pattern = Compile(pattern, obj, true);
      // Create a set of patterns by substituting each mass wildcard
      for (auto m : wildcard_masses_) {
        std::string tmp = maps_pattern;
        boost::replace_all(tmp, "$MASS", m);
        mappings.insert(tmp);
      }
      f_map[key] = mappings;
    });
  auto masses = cmb.cp().mass(wildcard_masses_, false).mass_set();
  cmb.cp().mass(wildcard_masses_, true)
    .ForEachObj([&](ch::Object const* obj) {
      // Build the fully-compiled key
      std::string key = Compile(pattern, obj, true);
      for (auto m : masses) {
        std::string full_key = key;
        boost::replace_all(full_key, "$MASS", m);
        if (f_map.count(full_key)) continue;
        std::set<std::string> mappings;
        // The fully-compiled pattern always goes in
        std::string wildcard_key = key;
        boost::replace_all(wildcard_key, "$MASS", obj->mass());
        mappings.insert(wildcard_key);
        f_map[full_key] = mappings;
      }
    });
  /*for (auto const& it : f_map) {
    std::cout << it.first << ": ";
    for (auto const& it2 : it.second) std::cout << it2 << " ";
    std::cout << "\n";
  }*/
  return f_map;
}

void CardWriter::MakeDirs(PatternMap const& map) const {
  std::set<boost::filesystem::path> mk_dirs;
  for (auto const& f : map) {
    mk_dirs.insert(boost::filesystem::path(f.first).parent_path());
  }
  for (auto & dir : mk_dirs) {
    FNLOGC(std::cout, v_ > 0) << "Creating dir " << dir.string() << "\n";

    boost::filesystem::create_directories(dir);
  }
}

std::map<std::string, CombineHarvester> CardWriter::WriteCards(std::string const& tag,
                       ch::CombineHarvester& cmb) const {
  #ifdef TIME_FUNCTIONS
    LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  #endif

  // Update the current "$TAG" value
  tag_ = tag;

  // Build the set of files we need to create (the map keys), and associate a
  // set of compiled Object patterns that are associated to each file (the map
  // values). One entry will typically look like:
  //
  // ["htt_mt_125.root"] = {"htt_mt_125.root", "htt_mt_*.root"}
  //
  // i.e. in addition to getting a pattern identical to the filename, we also
  // get a pattern where each $MASS placeholder has been replaced by a wildcard
  // expression.
  auto f_map = BuildMap(root_pattern_, cmb);

  if (create_dirs_) MakeDirs(f_map);

  // It turns out that the Compile method takes a small but appreciable amount
  // of CPU time. To evaluate it every time it's needed could be O(10^7) calls,
  // equivalent to tens of seconds for a complex model. To avoid this we just
  // calculate once for each ch::Object and store the result in a map, which we
  // use as a look-up later.
  std::map<Object const*, std::string> root_map;
  cmb.ForEachObj([&](ch::Object const* obj) {
      root_map[obj] = Compile(root_pattern_, obj);
    });
  std::map<Object const*, std::string> text_map;
  cmb.ForEachObj([&](ch::Object const* obj) {
      text_map[obj] = Compile(text_pattern_, obj);
    });

  std::map<std::string, CombineHarvester> datacards;
  for (auto const& f : f_map) {
    // Create each ROOT file (overwrite pre-existing)
    FNLOGC(std::cout, v_ > 0) << "Creating file " << f.first << "\n";
    TFile file(f.first.c_str(), "RECREATE");

    // Filter CH instance to leave only the objects that will be written into
    // this file
    CombineHarvester f_cmb = cmb.cp().FilterAll([&](ch::Object const* obj) {
      return !ch::contains(f.second, root_map.at(obj));
    });

    // Call BuildMap again - this time to figure out which text datacards to
    // create
    auto d_map = BuildMap(text_pattern_, f_cmb);

    // Create dirs if we're allowed to
    if (create_dirs_) MakeDirs(d_map);

    // Loop through each datacard
    for (auto const& d : d_map) {
      // Filter CH instance to leave only the objects that will be written into
      // this text datacard
      CombineHarvester d_cmb = f_cmb.cp().FilterAll([&](ch::Object const* obj) {
        return !ch::contains(d.second, text_map.at(obj));
      });
      FNLOGC(std::cout, v_ > 0) << "Creating datacard " << d.first << "\n";
      d_cmb.WriteDatacard(d.first, file);
      datacards[d.first] = d_cmb;
    }
  };
  return datacards;
}

std::string CardWriter::Compile(std::string pattern, ch::Object const* obj,
                                bool skip_mass) const {
  #ifdef TIME_FUNCTIONS
    LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  #endif
  boost::replace_all(pattern, "$TAG", tag_);
  boost::replace_all(pattern, "$BINID",
                     boost::lexical_cast<std::string>(obj->bin_id()));
  boost::replace_all(pattern, "$BIN", obj->bin());
  boost::replace_all(pattern, "$PROCESS", obj->process());
  if (!skip_mass) boost::replace_all(pattern, "$MASS", obj->mass());
  boost::replace_all(pattern, "$ERA", obj->era());
  boost::replace_all(pattern, "$CHANNEL", obj->channel());
  boost::replace_all(pattern, "$ANALYSIS", obj->analysis());
  return pattern;
}
}
