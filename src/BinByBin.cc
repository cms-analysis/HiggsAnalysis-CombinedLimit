#include "CombineTools/interface/BinByBin.h"
#include <iostream>
#include <string>
#include <vector>
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

namespace ch {

BinByBinFactory::BinByBinFactory()
    : pattern_("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bin_$#"),
      v_(0),
      bbb_threshold_(0.),
      merge_threshold_(0.),
      fix_norm_(true) {}


void BinByBinFactory::MergeBinErrors(CombineHarvester &cb) {
  // Reduce merge_threshold very slightly to avoid numerical issues
  // E.g. two backgrounds each with bin error 1.0. merge_threshold of
  // 0.5 should not result in merging - but can do depending on
  // machine and compiler
  double merge_threshold = merge_threshold_ - 1E-9 * merge_threshold_;
  auto bins = cb.bin_set();
  for (auto const& bin : bins) {
    unsigned bbb_added = 0;
    unsigned bbb_removed = 0;
    CombineHarvester tmp = std::move(cb.cp().bin({bin}).histograms());
    std::vector<Process *> procs;
    tmp.ForEachProc([&](Process *p) {
      if (p->shape()->GetSumw2N() == 0) {
        std::cout << "Process " << p->process()
                  << " does not continue the weights information needed for "
                     "valid errors, skipping\n";
      } else {
        procs.push_back(p);
      }
    });
    if (procs.size() == 0) continue;

    std::vector<std::unique_ptr<TH1>> h_copies(procs.size());
    for (unsigned i = 0; i < h_copies.size(); ++i) {
      h_copies[i] = procs[i]->ClonedScaledShape();
    }

    for (int i = 1; i <= h_copies[0]->GetNbinsX(); ++i) {
      double tot_bbb_added = 0.0;
      std::vector<std::pair<double, TH1 *>> result;
      for (unsigned j = 0; j < h_copies.size(); ++j) {
        double val = h_copies[j]->GetBinContent(i);
        double err = h_copies[j]->GetBinError(i);
        if (val == 0.0 &&  err == 0.0) continue;
        if (val == 0 || (err/val) > bbb_threshold_) {
          bbb_added += 1;
          tot_bbb_added += (err * err);
          result.push_back(std::make_pair(err*err, h_copies[j].get()));
        }
      }
      if (tot_bbb_added == 0.0) continue;
      std::sort(result.begin(), result.end());
      double removed = 0.0;
      for (unsigned r = 0; r < result.size(); ++r) {
        if ((result[r].first + removed) < (merge_threshold * tot_bbb_added) &&
            r < (result.size() - 1)) {
          bbb_removed += 1;
          removed += result[r].first;
          result[r].second->SetBinError(i, 0.0);
        }
      }
      double expand = std::sqrt(1. / (1. - (removed / tot_bbb_added)));
      for (unsigned r = 0; r < result.size(); ++r) {
        result[r]
            .second->SetBinError(i, result[r].second->GetBinError(i) * expand);
      }
    }
    for (unsigned i = 0; i < h_copies.size(); ++i) {
      procs[i]->set_shape(std::move(h_copies[i]), false);
    }
    if (v_ > 0) {
      std::cout << "BIN: " << bin << std::endl;
      std::cout << "Total bbb added:    " << bbb_added << "\n";
      std::cout << "Total bbb removed:  " << bbb_removed << "\n";
      std::cout << "Total bbb =======>: " << bbb_added-bbb_removed << "\n";
    }
  }
}

void BinByBinFactory::AddBinByBin(CombineHarvester &src, CombineHarvester &dest) {
  unsigned bbb_added = 0;
  std::vector<Process *> procs;
  src.ForEachProc([&](Process *p) { 
    procs.push_back(p);
  });
  for (unsigned i = 0; i < procs.size(); ++i) {
    if (!procs[i]->shape()) continue;
    TH1 const* h = procs[i]->shape();
    if (h->GetSumw2N() == 0) {
      std::cout << "Process " << procs[i]->process()
                << " does not continue the weights information needed for "
                   "valid errors, skipping\n";
      continue;
    }
    unsigned n_pop_bins = 0;
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      if (h->GetBinContent(j) > 0.0) ++n_pop_bins;
    }
    if (n_pop_bins <= 1 && fix_norm_) {
      if (v_ >= 1) {
        std::cout << "Requested fixed_norm but template has <= 1 populated "
                     "bins, skipping\n";
        std::cout << Process::PrintHeader << *(procs[i]) << "\n";
      }
      continue;
    }
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      bool do_bbb = false;
      double val = h->GetBinContent(j);
      double err = h->GetBinError(j);
      if (val == 0. && err > 0.) do_bbb = true;
      if (val > 0. && (err / val) > bbb_threshold_) do_bbb = true;
      // if (h->GetBinContent(j) <= 0.0) {
      //   if (h->GetBinError(j) > 0.0) {
      //     std::cout << *(procs_[i]) << "\n";
      //     std::cout << "Bin with content <= 0 and error > 0 found, skipping\n";
      //   }
      //   continue;
      // }
      if (do_bbb) {
        ++bbb_added;
        ch::Systematic sys;
        ch::SetProperties(&sys, procs[i]);
        sys.set_type("shape");
        std::string name = pattern_;
        boost::replace_all(name, "$ANALYSIS", sys.analysis());
        boost::replace_all(name, "$CHANNEL", sys.channel());
        boost::replace_all(name, "$BIN", sys.bin());
        boost::replace_all(name, "$BINID", boost::lexical_cast<std::string>(sys.bin_id()));
        boost::replace_all(name, "$ERA", sys.era());
        boost::replace_all(name, "$PROCESS", sys.process());
        boost::replace_all(name, "$MASS", sys.mass());
        boost::replace_all(name, "$#", boost::lexical_cast<std::string>(j));
        sys.set_name(name);
        sys.set_asymm(true);
        std::unique_ptr<TH1> h_d(static_cast<TH1 *>(h->Clone()));
        std::unique_ptr<TH1> h_u(static_cast<TH1 *>(h->Clone()));
        h_d->SetBinContent(j, val - err);
        if (h_d->GetBinContent(j) < 0.) h_d->SetBinContent(j, 0.);
        h_u->SetBinContent(j, val + err);
        if (fix_norm_) {
          sys.set_value_d(1.0);
          sys.set_value_u(1.0);
        } else {
          sys.set_value_d(h_d->Integral()/h->Integral());
          sys.set_value_u(h_u->Integral()/h->Integral());
        }
        sys.set_shapes(std::move(h_u), std::move(h_d), nullptr);
        dest.CreateParameterIfEmpty(sys.name());
        dest.InsertSystematic(sys);
      }
    }
  }
  // std::cout << "bbb added: " << bbb_added << std::endl;
}

void BinByBinFactory::MergeAndAdd(CombineHarvester &src, CombineHarvester &dest) {
  MergeBinErrors(src);
  AddBinByBin(src, dest);
}
}
