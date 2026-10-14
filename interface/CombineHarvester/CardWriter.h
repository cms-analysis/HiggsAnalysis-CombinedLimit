#ifndef CombineTools_CardWriter_h
#define CombineTools_CardWriter_h
#include <string>
#include <map>
#include <vector>
#include <set>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Object.h"

namespace ch {

/**
 * Automates the writing of datacards into directory structures
 *
 * A CardWriter is constructed with two strings which give patterns for where
 * the datacard text and ROOT files should be located and how they should be
 * named. These patterns can contain placeholders for the object metadata,
 * e.g.:
 *
 *     ch::CardWriter writer(
 *        "$TAG/$MASS/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt",
 *        "$TAG/common/$ANALYSIS_$CHANNEL.input_$ERA.root");
 *
 * All of the standard metadata placeholders are supported. Then, assuming a
 * CombineHarvester instance `cb` has already been configured, the cards are
 * written with:
 *
 *     writer.WriteCards("LIMITS/cmb", cb);
 *
 * The CardWriter will determine the sets of root files and datacards to
 * create given the metadata of the objects in `cb`. The first argument of
 * `WriteCards` is substituted into the `$TAG` placeholder. In this way
 * multiple calls to \ref WriteCards can produce output in different folders.
 *
 * By default CardWriter will create any directories that are needed, but this
 * behaviour can be switched off with the `CreateDirectories` method. Also
 * note that any existing files with the same name will be over-written.
 *
 * \note The $MASS variable is treated differently than the others when
 * constructing the list of output files. It is anticipated that a mass value
 * of `*` should not give rise to actual files containing the `*` character,
 * but rather implies an Object that is valid for any mass-point. It is
 * commonly used to denote background processes that do not depend on the
 * signal mass hypothesis. Therefore this class ignores objects with a mass
 * value of `*` when determining which files to create, and when creating the
 * files treats this object as matching any mass value. It is possible to
 * alter or remove this behaviour by supplying a new list of wildcard values
 * with the \ref SetWildcardMasses method.
 */
class CardWriter {
 public:
  /// Must be constructed with text and ROOT file patterns
  CardWriter(std::string const& text_pattern, std::string const& root_pattern);
  /// Write datacards according to patterns, substituting `$TAG` for `tag`
  std::map<std::string, CombineHarvester> WriteCards(std::string const& tag, ch::CombineHarvester& cmb) const;
  /// Set >= 1 for verbose output, otherwise silent
  CardWriter& SetVerbosity(unsigned v);
  /// Control whether directories can be created if missing
  CardWriter& CreateDirectories(bool flag);
  /// Redefine the mass values that should be treated as wildcards
  CardWriter& SetWildcardMasses(std::vector<std::string> const& masses);

 private:
  typedef std::map<std::string, std::set<std::string>> PatternMap;
  std::string text_pattern_;
  std::string root_pattern_;
  mutable std::string tag_;
  std::vector<std::string> wildcard_masses_;
  unsigned v_;
  bool create_dirs_;

  std::string Compile(std::string pattern, ch::Object const* obj,
                      bool skip_mass = false) const;
  PatternMap BuildMap(std::string const& pattern,
                      ch::CombineHarvester& cmb) const;
  void MakeDirs(PatternMap const& map) const;
};
}

#endif
