#ifndef CombineTools_HistMapping_h
#define CombineTools_HistMapping_h
#include <map>
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
    // The path of the file as given in the shapes line, joined with the path
    // of the datacard. It can contain $CHANNEL, $BIN, $PROCESS and $MASS, in
    // which case file is only opened when the shapes of an object are loaded
    std::string file_name;
    // Files opened for the mappings of one datacard, by path
    std::shared_ptr<std::map<std::string, std::shared_ptr<TFile>>> file_store;
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
    HistMapping(std::string const& p, std::string const& c, std::string const& pat, std::string const& s_pat);

    friend std::ostream& operator<<(std::ostream& out, HistMapping const& val);
  };
}  // namespace ch

#endif
