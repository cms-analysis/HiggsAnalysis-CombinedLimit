#ifndef CombineTools_HistMapping_h
#define CombineTools_HistMapping_h
#include <string>
#include <memory>
#include "TFile.h"

namespace ch {

  struct HistMapping {
    std::string process;
    std::string category;
    std::shared_ptr<TFile> file;
    std::string pattern;
    std::string syst_pattern;
    bool IsHist() const;
    bool IsPdf() const;
    bool IsData() const;

    friend std::ostream& operator<< (std::ostream &out, HistMapping const& val);
  };
}

#endif
