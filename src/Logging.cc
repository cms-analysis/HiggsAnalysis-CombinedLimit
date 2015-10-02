#include "CombineHarvester/CombineTools/interface/Logging.h"
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
      "https://github.com/cms-analysis/CombineHarvester/issues";
  res += "\n" + banner;
  return res;
}

// Implementation of FnTimer ("Function Timer") class
FnTimer::FnTimer(std::string name) : name_(name), calls_(0), elapsed_(0.) {}
FnTimer::~FnTimer() {
  printf(
      "[Timer] %-40s Calls: %-20u Total [s]: %-20.5g Per-call [s]: %-20.5g\n",
      name_.c_str(), calls_, elapsed_, elapsed_ / double(calls_));
}
FnTimer::Token FnTimer::Inc() {
  ++calls_;
  return Token(this);
}
void FnTimer::StartTimer() { start_ = std::chrono::system_clock::now(); }
void FnTimer::StopTimer() {
  end_ = std::chrono::system_clock::now();
  elapsed_ += std::chrono::duration<double>(end_ - start_).count();
}
FnTimer::Token::Token(FnTimer* src) : src_(src) { src_->StartTimer(); }
FnTimer::Token::~Token() { src_->StopTimer(); }
}
