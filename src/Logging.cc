#include "CombineTools/interface/Logging.h"
#include <string>
#include <iostream>
#include "boost/lexical_cast.hpp"

namespace ch {

std::string GetQualififedName(std::string const& str) {
  std::size_t bpos = str.find_first_of('(');
  if (bpos == std::string::npos) return str;
  std::size_t apos = bpos;
  bool break_cond = false;
  bool in_clause = false;
  while (!break_cond && apos > 0) {
    --apos;
    if (str[apos] == '>') in_clause = true;
    if (str[apos] == '<') in_clause = false;
    if (str[apos] == ' ' && !in_clause) {
      break_cond = true;
      ++apos;
    }
  }
  std::string cpy = str.substr(apos, bpos - apos);
  return cpy;
}

void LogLine(std::ostream& stream, std::string const& func,
             std::string const& message) {
  stream << "[" << func << "] " << message << "\n";
}

std::string FnError(std::string const& message, std::string const& file,
                    unsigned line, std::string const& fn) {
  std::string res;
  std::string banner(79, '*');
  res += "\n" + banner;
  res += "\nContext: Function " + GetQualififedName(fn) + " at ";
  res += "\n  " + file + ":" + boost::lexical_cast<std::string>(line);
  res += "\nProblem: " + message;
  res += "\n" + banner;
  res +=
      "\nPlease report issues at\n  "
      "https://github.com/cms-analysis/HiggsAnalysis-HiggsToTauTau/issues/new";
  res += "\n" + banner;
  return res;
}
}
