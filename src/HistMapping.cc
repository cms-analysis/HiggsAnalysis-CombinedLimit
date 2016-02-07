#include "CombineHarvester/CombineTools/interface/HistMapping.h"
#include <string>
#include "boost/format.hpp"

namespace ch {

HistMapping::HistMapping() : is_fake(false) {}

HistMapping::HistMapping(std::string const& p, std::string const& c,
                         std::string const& pat, std::string const& s_pat)
    : process(p), category(c), pattern(pat), syst_pattern(s_pat), is_fake(false) {}

bool HistMapping::IsHist() const {
  if (pattern.size() > 0 && pattern.find(':') == pattern.npos) {
    return true;
  } else {
    return false;
  }
}

bool HistMapping::IsPdf() const {
  if (pattern.size() > 0 && pattern.find(':') != pattern.npos) {
    return true;
  } else {
    return false;
  }
}

bool HistMapping::IsData() const {
  if (pattern.size() > 0 && pattern.find(':') != pattern.npos) {
    return true;
  } else {
    return false;
  }
}

std::string HistMapping::WorkspaceName() const {
  std::size_t colon = pattern.find_last_of(':');
  if (colon != pattern.npos) {
    return pattern.substr(0, colon);
  } else {
    return std::string();
  }
}
std::string HistMapping::WorkspaceObj() const {
  std::size_t colon = pattern.find_last_of(':');
  if (colon != pattern.npos) {
    return pattern.substr(colon + 1);
  } else {
    return std::string();
  }
}
std::string HistMapping::SystWorkspaceName() const {
  std::size_t colon = syst_pattern.find_last_of(':');
  if (colon != syst_pattern.npos) {
    return syst_pattern.substr(0, colon);
  } else {
    return std::string();
  }
}
std::string HistMapping::SystWorkspaceObj() const {
  std::size_t colon = syst_pattern.find_last_of(':');
  if (colon != syst_pattern.npos) {
    return syst_pattern.substr(colon + 1);
  } else {
    return std::string();
  }
}

std::ostream& operator<<(std::ostream& out, HistMapping const& val) {
  out << boost::format("%-7s %-7s %-20s %-20s %-20s") % val.process %
             val.category % (val.file.get() ? val.file->GetName() : "0") %
             val.pattern % val.syst_pattern;
  return out;
}
}
