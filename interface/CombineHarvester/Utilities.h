#ifndef CombineTools_Utilities_h
#define CombineTools_Utilities_h
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "TGraph.h"
#include "RooFitResult.h"
#include "CombineTools/interface/Parameter.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Logging.h"

namespace ch { class CombineHarvester; }

namespace ch {
// ---------------------------------------------------------------------------
// Paramter extraction/manipulation
// ---------------------------------------------------------------------------
std::vector<ch::Parameter> ExtractFitParameters(RooFitResult const& res);
std::vector<ch::Parameter> ExtractSampledFitParameters(RooFitResult const& res);


// ---------------------------------------------------------------------------
// Property matching & editing
// ---------------------------------------------------------------------------
void SetStandardBinNames(CombineHarvester & cb);

template<class T>
void SetStandardBinName(T *input) {
  std::string newname;
  newname = input->analysis() + "_"
    + input->channel() + "_"
    + boost::lexical_cast<std::string>(input->bin_id()) + "_"
    + input->era();
  input->set_bin(newname);
}

template<class T, class U>
bool MatchingProcess(T const& first, U const& second) {
  if (first.bin()         == second.bin()         &&
      first.process()     == second.process()     &&
      first.signal()      == second.signal()      &&
      first.analysis()    == second.analysis()    &&
      first.era()         == second.era()         &&
      first.channel()     == second.channel()     &&
      first.bin_id()      == second.bin_id()      &&
      first.mass()        == second.mass()) {
    return true;
  } else {
    return false;
  }
}

template<class T, class U>
void SetProperties(T * first, U const* second) {
  first->set_bin(second->bin());
  first->set_process(second->process());
  first->set_signal(second->signal());
  first->set_analysis(second->analysis());
  first->set_era(second->era());
  first->set_channel(second->channel());
  first->set_bin_id(second->bin_id());
  first->set_mass(second->mass());
}


template<class T>
void SetFromBinName(T *input, std::string parse_rules) {
  boost::replace_all(parse_rules, "$ANALYSIS",  "(?<ANALYSIS>\\w+)");
  boost::replace_all(parse_rules, "$ERA",       "(?<ERA>\\w+)");
  boost::replace_all(parse_rules, "$CHANNEL",   "(?<CHANNEL>\\w+)");
  boost::replace_all(parse_rules, "$BINID",     "(?<BINID>\\w+)");
  boost::replace_all(parse_rules, "$MASS",      "(?<MASS>\\w+)");
  boost::regex rgx(parse_rules);
  boost::smatch matches;
  boost::regex_search(input->bin(), matches, rgx);
  if (matches.str("ANALYSIS").length())
    input->set_analysis(matches.str("ANALYSIS"));
  if (matches.str("ERA").length())
    input->set_era(matches.str("ERA"));
  if (matches.str("CHANNEL").length())
    input->set_channel(matches.str("CHANNEL"));
  if (matches.str("BINID").length())
    input->set_bin_id(boost::lexical_cast<int>(matches.str("BINID")));
  if (matches.str("MASS").length())
    input->set_mass(matches.str("MASS"));
}

// ---------------------------------------------------------------------------
// Rate scaling
// ---------------------------------------------------------------------------

void ParseTable(std::map<std::string, TGraph>* graphs, std::string const& file,
                std::vector<std::string> const& fields);

void ScaleProcessRate(ch::Process* p,
                      std::map<std::string, TGraph> const* graphs,
                      std::string const& prod, std::string const& decay);



// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------
template<class T>
std::vector<T> Join(std::vector<std::vector<T>> const& in) {
  unsigned size = 0;
  for (auto const& x : in) size += x.size();
  std::vector<T> res(size);
  unsigned i = 0;
  for (auto const& x : in) {
    for (auto const& y : x) {
      res[i] = y;
      ++i;
    }
  }
  return res;
}

std::vector<std::string> JoinStr(
    std::vector<std::vector<std::string>> const& in);

template<class T>
std::vector<T> Set2Vec(std::set<T> const& in) {
  return std::vector<T>(in.begin(), in.end());
}

RooDataHist TH1F2Data(TH1F const& hist, RooRealVar const& x,
                      std::string const& name);

TH1F RebinHist(TH1F const& hist);

std::vector<std::vector<unsigned>> GenerateCombinations(
    std::vector<unsigned> vec);

std::vector<std::string> ParseFileLines(std::string const& file_name);

std::vector<std::string> MassesFromRange(std::string const& input,
                                         std::string const& fmt = "%.0f");
}

#endif
