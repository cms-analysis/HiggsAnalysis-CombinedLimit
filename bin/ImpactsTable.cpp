#include <string>
#include <iostream>
#include <vector>
#include "boost/format.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "TTree.h"
#include "TFile.h"
/*
Warning: this program is still in the prototype stage! There is currently no way
to configure its behaviour without modifying the code directly.
*/

using namespace std;
using boost::format;

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
  TFile f("impact.root");
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
    val_vec.push_back(reader);
    std::cout << val_vec.back().parameter << "\n";
  }

  unsigned max = 30;
  string fmt = "& %-5.2f & %-5.2f";

  // Here we'll make use of C++11 "string literals", to avoid having to escape
  // all the latex '\' characters. More info on the format here:
  // https://solarianprogrammer.com/2011/10/16/cpp-11-raw-strings-literals-tutorial/
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
  return 0;
}
