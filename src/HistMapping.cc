#include "CombineTools/interface/HistMapping.h"
#include <string>
#include "boost/format.hpp"

namespace ch {

  bool HistMapping::IsHist() const {
    if (pattern.find(':') == pattern.npos) {
      return true;
    } else {
      return false;
    }
  }

  bool HistMapping::IsPdf() const {
    if (pattern.find(':') != pattern.npos) {
      return true;
    } else {
      return false;
    }
  }

  bool HistMapping::IsData() const {
    if (pattern.find(':') != pattern.npos) {
      return true;
    } else {
      return false;
    }
  }

  std::ostream& operator<<(std::ostream& out, HistMapping const& val) {
    out << boost::format("%-7s %-7s %-20s %-20s %-20s")
    % val.process
    % val.category
    % (val.file.get() ? val.file->GetName() : "0")
    % val.pattern
    % val.syst_pattern;
    return out;
  }
}
