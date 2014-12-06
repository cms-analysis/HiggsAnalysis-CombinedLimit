#include "CombineTools/interface/Utilities.h"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <map>
#include "boost/format.hpp"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "CombineTools/interface/CombineHarvester.h"

namespace ch {
// ---------------------------------------------------------------------------
// Paramter extraction/manipulation
// ---------------------------------------------------------------------------
std::vector<ch::Parameter> ExtractFitParameters(RooFitResult const& res) {
  std::vector<ch::Parameter> params;
  params.resize(res.floatParsFinal().getSize());
  for (int i = 0; i < res.floatParsFinal().getSize(); ++i) {
    RooRealVar const* var =
        dynamic_cast<RooRealVar const*>(res.floatParsFinal().at(i));
    params[i].set_name(std::string(var->GetName()));
    params[i].set_val(var->getVal());
    params[i].set_err_d(var->getErrorLo());
    params[i].set_err_u(var->getErrorHi());
  }
  return params;
}

std::vector<ch::Parameter> ExtractSampledFitParameters(
    RooFitResult const& res) {
  std::vector<ch::Parameter> params;
  params.resize(res.floatParsFinal().getSize());
  RooArgList const& rands = res.randomizePars();
  for (int i = 0; i < res.floatParsFinal().getSize(); ++i) {
    RooRealVar const* var = dynamic_cast<RooRealVar const*>(rands.at(i));
    params[i].set_name(std::string(var->GetName()));
    params[i].set_val(var->getVal());
    params[i].set_err_d(var->getErrorLo());
    params[i].set_err_u(var->getErrorHi());
  }
  return params;
}

// ---------------------------------------------------------------------------
// Property matching & editing
// ---------------------------------------------------------------------------
void SetStandardBinNames(CombineHarvester & cb) {
  cb.ForEachObs(ch::SetStandardBinName<ch::Observation>);
  cb.ForEachProc(ch::SetStandardBinName<ch::Process>);
  cb.ForEachSyst(ch::SetStandardBinName<ch::Systematic>);
}

// ---------------------------------------------------------------------------
// Rate scaling
// ---------------------------------------------------------------------------
void ParseTable(std::map<std::string, TGraph>* graphs, std::string const& file,
                std::vector<std::string> const& fields) {
  auto lines = ch::ParseFileLines(file);
  for (auto const& f : fields) {
    (*graphs)[f] = TGraph(lines.size());
  }
  for (unsigned i = 0; i < lines.size(); ++i) {
    std::vector<std::string> words;
    boost::split(words, lines[i], boost::is_any_of("\t "),
                 boost::token_compress_on);
    double m = boost::lexical_cast<double>(words.at(0));
    for (unsigned f = 0; f < fields.size(); ++f) {
      (*graphs)[fields[f]]
          .SetPoint(i, m, boost::lexical_cast<double>(words.at(f + 1)));
    }
  }
}

void ScaleProcessRate(ch::Process* p,
                      std::map<std::string, TGraph> const* graphs,
                      std::string const& prod, std::string const& decay,
                      std::string const& force_mass) {
  double mass =
      boost::lexical_cast<double>(force_mass.size() ? force_mass : p->mass());
  double scale = 1.0;
  if (prod != "") {
    if (!graphs->count(prod))
      throw std::runtime_error(
          FNERROR("Requested TGraph " + prod + " not found in map"));
    scale *= graphs->find(prod)->second.Eval(mass);
  }
  if (decay != "") {
    if (!graphs->count(decay))
      throw std::runtime_error(
          FNERROR("Requested TGraph " + decay + " not found in map"));
    scale *= graphs->find(decay)->second.Eval(mass);
  }
  p->set_rate(p->rate() * scale);
}

// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------
std::vector<std::string> JoinStr(
    std::vector<std::vector<std::string>> const& in) {
  return Join<std::string>(in);
}

RooDataHist TH1F2Data(TH1F const& hist, RooRealVar const& x,
                      std::string const& name) {
  TH1F shape("tmp", "tmp", hist.GetNbinsX(), 0.,
             static_cast<float>(hist.GetNbinsX()));
  for (int i = 1; i <= hist.GetNbinsX(); ++i) {
    shape.SetBinContent(i, hist.GetBinContent(i));
  }
  RooDataHist dh(name.c_str(), name.c_str(),
                 RooArgList(x), RooFit::Import(shape, false));
  return dh;
}

TH1F RebinHist(TH1F const& hist) {
  TH1::AddDirectory(0);
  TH1F shape("tmp", "tmp", hist.GetNbinsX(), 0.,
             static_cast<float>(hist.GetNbinsX()));
  for (int i = 1; i <= hist.GetNbinsX(); ++i) {
    shape.SetBinContent(i, hist.GetBinContent(i));
  }
  return shape;
}

TH1F RestoreBinning(TH1F const& src, TH1F const& ref) {
  TH1F res = ref;
  res.Reset();
  for (int x = 1; x <= res.GetNbinsX(); ++x) {
    res.SetBinContent(x, src.GetBinContent(x));
    res.SetBinError(x, src.GetBinError(x));
  }
  return res;
}

std::vector<std::vector<unsigned>> GenerateCombinations(
    std::vector<unsigned> vec) {
  unsigned n = vec.size();
  std::vector<unsigned> idx(n, 0);
  std::vector<std::vector<unsigned>> result;
  // if any one of the elements is zero there are no
  // combinations to build
  if (std::find(vec.begin(), vec.end(), 0) != vec.end()) {
    return result;
  }
  result.push_back(idx);
  bool exit_loop = false;
  while (exit_loop == false) {
    // Find the first index we can increment (if possible)
    for (unsigned i = 0; i < n; ++i) {
      if ((idx[i] + 1) == vec[i]) {
        if (i != n - 1) {
          idx[i] = 0;
        } else {
          // we're done
          exit_loop = true;
          break;
        }
      } else {
        ++(idx[i]);
        result.push_back(idx);
        break;
      }
    }
  }
  return result;
}

std::vector<std::string> ParseFileLines(std::string const& file_name) {
  // Build a vector of input files
  std::vector<std::string> files;
  std::ifstream file;
  file.open(file_name.c_str());
  if (!file.is_open()) {
    throw std::runtime_error(
        FNERROR("File " + file_name + " could not be opened"));
  }
  std::string line = "";
  while (std::getline(file, line)) {  // while loop through lines
    files.push_back(line);
  }
  file.close();
  return files;
}

std::vector<std::string> MassesFromRange(std::string const& input,
                                         std::string const& fmt) {
  std::set<double> mass_set;
  std::vector<std::string> tokens;
  boost::split(tokens, input, boost::is_any_of(","));
  for (auto const& t : tokens) {
    std::vector<std::string> sub_tokens;
    boost::split(sub_tokens, t, boost::is_any_of("-:"));
    if (sub_tokens.size() == 1) {
      double mass_val = boost::lexical_cast<double>(sub_tokens[0]);
      mass_set.insert(mass_val);
    } else if (sub_tokens.size() == 3) {
      double lo = boost::lexical_cast<double>(sub_tokens[0]);
      double hi = boost::lexical_cast<double>(sub_tokens[1]);
      double step = boost::lexical_cast<double>(sub_tokens[2]);
      if (hi <= lo)
        throw std::runtime_error(
            "[MassesFromRange] High mass is smaller than low mass!");
      double start = lo;
      while (start < hi + 0.001) {
        mass_set.insert(start);
        start += step;
      }
    }
  }
  std::vector<std::string> result;
  for (auto const& m : mass_set) {
    result.push_back((boost::format(fmt) % m).str());
  }
  return result;
}

boost::filesystem::path make_relative(boost::filesystem::path p_from,
                                      boost::filesystem::path p_to) {
  p_from = boost::filesystem::absolute(p_from);
  p_to = boost::filesystem::absolute(p_to);
  boost::filesystem::path ret;
  boost::filesystem::path::const_iterator itrFrom(p_from.begin()),
      itrTo(p_to.begin());
  // Find common base
  for (boost::filesystem::path::const_iterator toEnd(p_to.end()),
       fromEnd(p_from.end());
       itrFrom != fromEnd && itrTo != toEnd && *itrFrom == *itrTo;
       ++itrFrom, ++itrTo);
  // Navigate backwards in directory to reach previously found base
  for (boost::filesystem::path::const_iterator fromEnd(p_from.end());
       itrFrom != fromEnd; ++itrFrom) {
    if ((*itrFrom) != ".") ret /= "..";
  }
  // Now navigate down the directory branch
  // ret.append(itrTo, p_to.end());
  for (boost::filesystem::path::const_iterator toEnd(p_to.end());
       itrTo != toEnd; ++itrTo) {
    ret /= *itrTo;
  }
  return ret;
}
}
