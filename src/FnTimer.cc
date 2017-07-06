#include "HiggsAnalysis/CombinedLimit/interface/FnTimer.h"
#include <string>
#include <iostream>
#include "boost/lexical_cast.hpp"


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


// Implementation of FnTimer ("Function Timer") class
FnTimer::FnTimer(std::string name) : name_(name), calls_(0), elapsed_(0.), elapsed_overhead_(0.) {}
FnTimer::~FnTimer() {
  printf(
      "[Timer] %-40s Calls: %-20u Total [s]: %-20.5g Per-call [s]: %-20.5g\n",
      name_.c_str(), calls_, elapsed_ - elapsed_overhead_, (elapsed_ - elapsed_overhead_) / double(calls_));

}
FnTimer::Token FnTimer::Inc() {
  ++calls_;
  return Token(this);
}
FnTimer::Token::Token(FnTimer* src) : src_(src) { src_->StartTimer(); }
FnTimer::Token::~Token() { src_->StopTimer(); }

