#include <string>
#include <iostream>
#include <vector>
#include "boost/format.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include "TTree.h"
#include "TFile.h"
/*
Warning: this program is still in the prototype stage! There is currently no way
to configure its behaviour without modifying the code directly.
*/

using namespace std;
using boost::format;
namespace po = boost::program_options;

// Struct holds info on each column in the table
struct Vals {
  std::string parameter;
  double par_best;
  double par_lo;
  double par_hi;
  double par_best_post;
  double par_lo_post;
  double par_hi_post;
  double par_best_pre;
  double par_lo_pre;
  double par_hi_pre;
  double impact;
  double impact_post;
  double impact_pre;
  int rank_impact;
  int rank_impact_post;
  int rank_impact_pre;
};


int main(int argc, char* argv[]) {
  string input_file = "";
  unsigned max = 99999;
  bool do_latex = false;
  bool do_text = false;
  std::vector<string> filters;

  po::variables_map vm;
  po::notify(vm);
  po::options_description config("configuration");
  config.add_options()
    ("input,i", po::value<string>(&input_file)->required())
    ("max,m", po::value<unsigned>(&max)->default_value(max))
    ("latex,l", po::value<bool>(&do_latex)->implicit_value(true))
    ("text,t", po::value<bool>(&do_text)->implicit_value(true))
    ("filter", po::value<vector<string>>(&filters)->multitoken());
  po::store(po::command_line_parser(argc, argv)
    .options(config).allow_unregistered().run(), vm);
  po::notify(vm);

  std::vector<boost::regex> rgx;
  for (auto const& ele : filters) rgx.emplace_back(ele);


  TFile f(input_file.c_str());
  TTree *t = (TTree*)f.Get("impact");

  Vals reader;
  TString *param_str = nullptr;
  t->SetBranchAddress("parameter",        &param_str);
  t->SetBranchAddress("par_best",         &reader.par_best);
  t->SetBranchAddress("par_lo",           &reader.par_lo);
  t->SetBranchAddress("par_hi",           &reader.par_hi);
  t->SetBranchAddress("par_best_post",    &reader.par_best_post);
  t->SetBranchAddress("par_lo_post",      &reader.par_lo_post);
  t->SetBranchAddress("par_hi_post",      &reader.par_hi_post);
  t->SetBranchAddress("par_best_pre",     &reader.par_best_pre);
  t->SetBranchAddress("par_lo_pre",       &reader.par_lo_pre);
  t->SetBranchAddress("par_hi_pre",       &reader.par_hi_pre);
  t->SetBranchAddress("impact",           &reader.impact);
  t->SetBranchAddress("impact_post",      &reader.impact_post);
  t->SetBranchAddress("impact_pre",       &reader.impact_pre);
  t->SetBranchAddress("rank_impact",      &reader.rank_impact);
  t->SetBranchAddress("rank_impact_post", &reader.rank_impact_post);
  t->SetBranchAddress("rank_impact_pre",  &reader.rank_impact_pre);
  std::vector<Vals> val_vec;
  for (unsigned i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    reader.parameter = *param_str;
    std::cout << reader.parameter << "\n";
    bool matches = false;
    for (unsigned j = 0; j < rgx.size(); ++j) {
      if (boost::regex_match((*param_str).Data(), rgx[j])) {
        matches = true;
        break;
      }
    }
    if (rgx.size() && !matches) continue;
    val_vec.push_back(reader);
  }

  double tot_obs = 0;
  double tot_post = 0;
  double tot_pre = 0;
  for (unsigned i = 0; i < max && i < val_vec.size(); ++i) {
    Vals &val = val_vec[i];
    tot_obs += (val.impact * val.impact);
    tot_post += (val.impact_post * val.impact_post);
    tot_pre += (val.impact_pre * val.impact_pre);
  }
  tot_obs = std::sqrt(tot_obs);
  tot_post = std::sqrt(tot_post);
  tot_pre = std::sqrt(tot_pre);

  string fmt = "& %-5.2f & %-5.2f";

  if (do_text) {
    cout << format("%-70s %-15s | %-15s | %-15s | %-13s | "
          "%-13s | %-13s\n")
      % "Parameter"
      % "Best (obs)"
      % "Best (post)"
      % "Best (pre)"
      % "Impact (obs)"
      % "Impact (post)"
      % "Impact (pre)";

    for (unsigned i = 0; i < max && i < val_vec.size(); ++i) {
      Vals &val = val_vec[i];

      cout << format(
          "%-70s %7.2f %7.2f | %7.2f %7.2f | %7.2f %7.2f | %7.3f %5i | "
          "%7.3f %5i | %7.3f %5i\n")
        % val.parameter
        % val.par_best
        % ((val.par_hi - val.par_lo)/2.)
        % val.par_best_post
        % ((val.par_hi_post - val.par_lo_post)/2.)
        % val.par_best_pre
        % ((val.par_hi_pre - val.par_lo_pre)/2.)
        % val.impact % val.rank_impact
        % val.impact_post % val.rank_impact_post
        % val.impact_pre % val.rank_impact_pre;
    }
  }

  // Here we'll make use of C++11 "string literals", to avoid having to escape
  // all the latex '\' characters. More info on the format here:
  // https://solarianprogrammer.com/2011/10/16/cpp-11-raw-strings-literals-tutorial/
  if (do_latex) {
    cout << R"(\small)";
    cout << R"(\begin{tabular}{l)";
    cout << R"(r@{$ \,\,\pm\,\, $}l)";
    cout << R"(r@{$ \,\,\pm\,\, $}l)";
    cout << R"(rrrr)";
    cout << R"(})" << "\n";


    cout << format("%-60s") % "Nuisance Parameter"
         << R"(& \multicolumn{2}{c}{$\hat{\theta}$ (obs)})"
         << R"(& \multicolumn{2}{c}{$\hat{\theta}$ (asimov)})"
         << R"(& $\Delta\hat{\mu}$ (obs) & Rank)"
         << R"(& $\Delta\hat{\mu}$ (asimov) & Rank)"
         << R"(\\)" << "\n" << R"(\hline)" << "\n";


    for (unsigned i = 0; i < max && i < val_vec.size(); ++i) {
      Vals & val = val_vec[i];

      boost::replace_all(val.parameter, "_", "\\_");

      cout << format("%-60s") % val.parameter
           << format(fmt) % val.par_best % ((val.par_hi - val.par_lo)/2.)
           << format(fmt) % val.par_best_post % ((val.par_hi_post - val.par_lo_post)/2.)
           << format("& %-5.3f & %-4i") % val.impact % val.rank_impact
           << format("& %-5.3f & %-4i") % val.impact_post % val.rank_impact_post
           << R"(\\)" << "\n";
    }

    cout << R"(\hline)" << "\n";
    cout << R"(\end{tabular})" << "\n";
  }
  std::cout << "Total impact obs:  " << tot_obs << "\n";
  std::cout << "Total impact post: " << tot_post << "\n";
  std::cout << "Total impact pre:  " << tot_pre << "\n";
  return 0;
}
