#include "CombineHarvester/CombineTools/interface/JsonTools.h"

#include <fstream>
#include <string>
#include "boost/algorithm/string.hpp"

namespace ch {
Json::Value ExtractJsonFromFile(std::string const& file) {
  Json::Value js;
  Json::Reader json_reader;
  std::fstream input;
  input.open(file);
  json_reader.parse(input, js);
  return js;
}

Json::Value ExtractJsonFromString(std::string const& str) {
  Json::Value js;
  Json::Reader reader(Json::Features::all());
  reader.parse(str, js);
  return js;
}

void UpdateJson(Json::Value& a, Json::Value const& b) {
  if (!a.isObject() || !b.isObject()) return;
  for (auto const& key : b.getMemberNames()) {
    if (a[key].isObject()) {
      UpdateJson(a[key], b[key]);
    } else {
      a[key] = b[key];
    }
  }
}

Json::Value MergedJson(int argc, char* argv[]) {
  Json::Value res(Json::objectValue);
  if (argc >= 1) {
    for (int i = 1; i < argc; ++i) {
      if (boost::algorithm::iends_with(argv[i], ".json")) {
        Json::Value val = ExtractJsonFromFile(argv[i]);
        UpdateJson(res, val);
      } else {
        Json::Value val = ExtractJsonFromString(argv[i]);
        UpdateJson(res, val);
      }
    }
  }
  return res;
}

Json::Value MergedJson(std::vector<std::string> const& vec) {
  Json::Value res(Json::objectValue);
  for (auto const& str : vec) {
    if (boost::algorithm::iends_with(str, ".json")) {
      Json::Value val = ExtractJsonFromFile(str);
      UpdateJson(res, val);
    } else {
      Json::Value val = ExtractJsonFromString(str);
      UpdateJson(res, val);
    }
  }
  return res;
}
}
