#include "CombineHarvester/CombineTools/interface/TFileIO.h"
#include <memory>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace ch {

std::unique_ptr<TH1> GetClonedTH1(TFile* file, std::string const& path) {
  if (!file) {
    throw std::runtime_error(FNERROR("Supplied ROOT file pointer is null"));
  }
  TDirectory* backup_dir = gDirectory;
  file->cd();
  bool cur_status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TObject *obj = gDirectory->Get(path.c_str());
  if (!obj) {
    gDirectory = backup_dir;
    throw std::runtime_error(
        FNERROR("TH1 " + path + " not found in " + file->GetName()));
  }
  std::unique_ptr<TH1> res(
      dynamic_cast<TH1*>(obj));
  TH1::AddDirectory(cur_status);
  if (!res) {
    gDirectory = backup_dir;
    throw std::runtime_error(FNERROR("Object " + path + " in " +
                                     file->GetName() + " is not of type TH1"));
  }
  // res->SetDirectory(0);
  gDirectory = backup_dir;
  return res;
}
}
