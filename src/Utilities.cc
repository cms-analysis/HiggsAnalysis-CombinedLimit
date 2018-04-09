#include "CombineHarvester/CombineTools/interface/Utilities.h"
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
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {

RooArgSet ParametersByName(RooAbsReal const* pdf, RooArgSet const* dat_vars) {
  // Get all pdf parameters first
  // We are expected to manage the memory of the RooArgSet pointer we're given,
  // so let's use a unique_ptr to ensure it gets cleaned up
  std::unique_ptr<RooArgSet> all_vars(pdf->getParameters(RooArgSet()));
  // Get the data variables and fill a set with all the names
  std::set<std::string> names;
  RooFIter dat_it = dat_vars->fwdIterator();
  RooAbsArg *dat_arg = nullptr;
  while((dat_arg = dat_it.next())) {
    names.insert(dat_arg->GetName());
  }

  // Build a new RooArgSet from all_vars, excluding any in names
  RooArgSet result_set;
  RooFIter vars_it = all_vars->fwdIterator();
  RooAbsArg *vars_arg = nullptr;
  while((vars_arg = vars_it.next())) {
    if (!names.count(vars_arg->GetName())) {
      result_set.add(*vars_arg);
    }
  }
  return result_set;
}

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
void SetStandardBinNames(CombineHarvester& cb, std::string const& pattern) {
  cb.ForEachObj([&](ch::Object* obj) {
    ch::SetStandardBinName(obj, pattern);
  });
}

void SetStandardBinName(ch::Object* obj, std::string pattern) {
  boost::replace_all(pattern, "$BINID",
                     boost::lexical_cast<std::string>(obj->bin_id()));
  boost::replace_all(pattern, "$BIN", obj->bin());
  boost::replace_all(pattern, "$PROCESS", obj->process());
  boost::replace_all(pattern, "$MASS", obj->mass());
  boost::replace_all(pattern, "$ERA", obj->era());
  boost::replace_all(pattern, "$CHANNEL", obj->channel());
  boost::replace_all(pattern, "$ANALYSIS", obj->analysis());
  obj->set_bin(pattern);
}

void SetFromBinName(ch::Object *input, std::string parse_rules) {
  boost::replace_all(parse_rules, "$ANALYSIS",  "(?<ANALYSIS>\\w+)");
  boost::replace_all(parse_rules, "$ERA",       "(?<ERA>\\w+)");
  boost::replace_all(parse_rules, "$CHANNEL",   "(?<CHANNEL>\\w+)");
  boost::replace_all(parse_rules, "$BINID",     "(?<BINID>\\w+)");
  boost::replace_all(parse_rules, "$MASS",      "(?<MASS>\\w+)");
  boost::regex rgx(parse_rules);
  boost::smatch matches;
  boost::regex_search(input->bin(), matches, rgx);
  if (matches.str("ANALYSIS").length())
    input->set_analysis(matches.str("ANALYSIS"));
  if (matches.str("ERA").length())
    input->set_era(matches.str("ERA"));
  if (matches.str("CHANNEL").length())
    input->set_channel(matches.str("CHANNEL"));
  if (matches.str("BINID").length())
    input->set_bin_id(boost::lexical_cast<int>(matches.str("BINID")));
  if (matches.str("MASS").length())
    input->set_mass(matches.str("MASS"));
}


// ---------------------------------------------------------------------------
// Rate scaling
// ---------------------------------------------------------------------------
TGraph TGraphFromTable(std::string filename, std::string const& x_column, std::string const& y_column) {
  auto lines = ch::ParseFileLines(filename);
  if (lines.size() == 0) throw std::runtime_error(FNERROR("Table is empty"));
  std::vector<std::string> fields;
  boost::split(fields, lines[0], boost::is_any_of("\t "),
               boost::token_compress_on);
  auto x_it = std::find(fields.begin(), fields.end(), x_column);
  if (x_it == fields.end())
    throw std::runtime_error(FNERROR("No column with label " + x_column));
  auto y_it = std::find(fields.begin(), fields.end(), y_column);
  if (y_it == fields.end())
    throw std::runtime_error(FNERROR("No column with label " + y_column));
  unsigned x_col(x_it - fields.begin());
  unsigned y_col(y_it - fields.begin());
  TGraph res(lines.size() - 1);
  for (unsigned i = 1; i < lines.size(); ++i) {
    std::vector<std::string> words;
    boost::split(words, lines[i], boost::is_any_of("\t "),
                 boost::token_compress_on);
    if (words.size() != fields.size()) {
      FNLOG(std::cout) << "Skipped this line:\n" << lines[i] << "\n";
      continue;
    }
    res.SetPoint(i, boost::lexical_cast<double>(words.at(x_col)),
                    boost::lexical_cast<double>(words.at(y_col)));
  }
  return res;
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
    shape.SetBinError(i, hist.GetBinError(i));
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

bool is_float(std::string const& str) {
  std::istringstream iss(str);
  float f;
  iss >> std::noskipws >> f;  // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail();
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

std::vector<std::string> ValsFromRange(std::string const& input,
                                       std::string const& fmt) {
  std::set<double> mass_set;
  std::vector<std::string> tokens;
  boost::split(tokens, input, boost::is_any_of(","));
  for (auto const& t : tokens) {
    std::vector<std::string> sub_tokens;
    boost::split(sub_tokens, t, boost::is_any_of(":|"));
    if (sub_tokens.size() == 1) {
      double mass_val = boost::lexical_cast<double>(sub_tokens[0]);
      mass_set.insert(mass_val);
    } else if (sub_tokens.size() == 3) {
      double lo = boost::lexical_cast<double>(sub_tokens[0]);
      double hi = boost::lexical_cast<double>(sub_tokens[1]);
      double step = boost::lexical_cast<double>(sub_tokens[2]);
      if (hi <= lo)
        throw std::runtime_error(
            "[ValsFromRange] High mass is smaller than low mass!");
      double start = lo;
      while (start < hi + 1E-4) {
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

bool HasNegativeBins(TH1 const* h) {
  bool has_negative = false;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) < 0.) {
      has_negative = true;
    }
  }
  return has_negative;
}

void ZeroNegativeBins(TH1 *h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinContent(i) < 0.) {
      h->SetBinContent(i, 0.);
    }
  }
}
}
