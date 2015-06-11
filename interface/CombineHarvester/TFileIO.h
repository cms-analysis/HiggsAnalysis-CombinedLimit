#ifndef CombineTools_TFileIO_h
#define CombineTools_TFileIO_h
#include <memory>
#include <string>
#include <vector>
#include "boost/algorithm/string.hpp"
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"
#include "CombineTools/interface/Logging.h"

namespace ch {

std::unique_ptr<TH1> GetClonedTH1(TFile* file, std::string const& path);

template <class T>
void WriteToTFile(T const* ptr, TFile* file, std::string const& path);

// Extracts objects from the form:
// "path/to/a/file.root:path/to/object"
template <class T>
T OpenFromTFile(std::string const& fullpath);

template <class T>
T OpenFromTFile(TFile* file, std::string const& path);
}

// Template function implementation
// ------------------------------------------------------------------
template <class T>
void ch::WriteToTFile(T const* ptr, TFile* file, std::string const& path) {
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
      gDirectory->WriteTObject(ptr, as_vec.back().c_str());
    }
    gDirectory->cd("/");
  }
}

template <class T>
T ch::OpenFromTFile(std::string const& fullpath) {
  std::size_t pos = fullpath.find(':');
  std::string filepath = "";
  std::string objectpath = "";
  if (pos == std::string::npos || pos == 0 || pos == (fullpath.size() - 1)) {
    throw std::runtime_error(
        FNERROR("Input path must of the format file.root:object"));
  } else {
    filepath = fullpath.substr(0, pos);
    objectpath = fullpath.substr(pos + 1);
  }
  TFile file(filepath.c_str());
  if (!file.IsOpen()) {
    throw std::runtime_error(FNERROR("File is invalid"));
  }
  file.cd();
  T* obj_ptr = dynamic_cast<T*>(gDirectory->Get(objectpath.c_str()));
  if (!obj_ptr) {
    throw std::runtime_error(
        FNERROR("Object " + objectpath + " is missing or of wrong type"));
  }
  return *obj_ptr;
}

template <class T>
T ch::OpenFromTFile(TFile* file, std::string const& path) {
  if (!file || !file->IsOpen()) {
    throw std::runtime_error(FNERROR("File is null or invalid"));
  }
  file->cd();
  T* obj_ptr = dynamic_cast<T*>(gDirectory->Get(path.c_str()));
  if (!obj_ptr) {
    throw std::runtime_error(
        FNERROR("Object " + path + " is missing or of wrong type"));
  }
  return *obj_ptr;
}

#endif
