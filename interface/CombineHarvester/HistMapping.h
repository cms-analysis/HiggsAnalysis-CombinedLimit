#ifndef CombineTools_HistMapping_h
#define CombineTools_HistMapping_h
#include <string>
#include <memory>
#include "TFile.h"
#include "RooWorkspace.h"

namespace ch {

  struct HistMapping {
    std::string process;
    std::string category;
    std::string pattern;
    std::string syst_pattern;
    std::shared_ptr<TFile> file;
    std::shared_ptr<RooWorkspace> ws;
    std::shared_ptr<RooWorkspace> sys_ws;
    bool is_fake;
    bool IsHist() const;
    bool IsPdf() const;
    bool IsData() const;
    std::string WorkspaceName() const;
    std::string WorkspaceObj() const;
    std::string SystWorkspaceName() const;
    std::string SystWorkspaceObj() const;

    HistMapping();
    HistMapping(std::string const& p, std::string const& c,
                std::string const& pat, std::string const& s_pat);

    friend std::ostream& operator<< (std::ostream &out, HistMapping const& val);
  };
}

#endif
