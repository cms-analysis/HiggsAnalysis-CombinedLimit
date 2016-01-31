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
  if (!gDirectory->Get(path.c_str())) {
    gDirectory = backup_dir;
    throw std::runtime_error(
        FNERROR("TH1 " + path + " not found in " + file->GetName()));
  }
  std::unique_ptr<TH1> res(
      dynamic_cast<TH1*>(gDirectory->Get(path.c_str())->Clone()));
  if (!res) {
    gDirectory = backup_dir;
    throw std::runtime_error(FNERROR("Object " + path + " in " +
                                     file->GetName() + " is not of type TH1"));
  }
  res->SetDirectory(0);
  gDirectory = backup_dir;
  return res;
}


void WriteHistToTFileAndSetName(TH1F * ptr, TFile* file, std::string const& path) {
  #ifdef TIME_FUNCTIONS
    LAUNCH_FUNCTION_TIMER(__timer__, __token__)
  #endif
  file->cd();
  std::vector<std::string> as_vec;
  boost::split(as_vec, path, boost::is_any_of("/"));
  if (as_vec.size() >= 1) {
    for (unsigned i = 0; i < as_vec.size() - 1; ++i) {
      if (!gDirectory->GetDirectory(as_vec[i].c_str())) {
        gDirectory->mkdir(as_vec[i].c_str());
      }
      gDirectory->cd(as_vec[i].c_str());
    }
    if (!gDirectory->FindKey(as_vec.back().c_str())) {
      ptr->SetName((as_vec.back()).c_str());
      gDirectory->WriteTObject(ptr, as_vec.back().c_str());
    }
    gDirectory->cd("/");
  }
}

}
