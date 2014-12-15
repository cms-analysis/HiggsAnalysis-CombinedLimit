#include "CombinePdfs/interface/MorphFunctions.h"
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"
#include "boost/multi_array.hpp"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooProduct.h"
#include "RooConstVar.h"
#include "CombineTools/interface/Logging.h"

namespace ch {


void BuildRooMorphing(RooWorkspace& ws, CombineHarvester& cb,
                      std::string const& bin, std::string const& process,
                      RooAbsReal& mass_var, std::string const& mass_min,
                      std::string const& mass_max, bool allow_morph, bool verbose) {
  using std::set;
  using std::vector;
  using std::string;
  using boost::lexical_cast;
  using boost::multi_array;
  using boost::extents;

  RooFit::MsgLevel backup_msg_level =
      RooMsgService::instance().globalKillBelow();
  if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if (verbose)
    std::cout << ">> Bin: " << bin << "  Process: " << process << "\n";
  TString key = bin + "_" + process;

  CombineHarvester cb_bp =
      std::move(cb.cp().bin({bin}).process({process}));

  vector<string> masses_all = Set2Vec(cb_bp.GenerateSetFromProcs<string>(
      std::mem_fn(&ch::Process::mass)));

  std::sort(masses_all.begin(), masses_all.end(),
    [](string const& s1, string const& s2) {
      return lexical_cast<double>(s1) < lexical_cast<double>(s2);
    });

  auto m_it_lo = std::find(masses_all.begin(), masses_all.end(), mass_min);
  auto m_it_hi = std::find(masses_all.begin(), masses_all.end(), mass_max);

  if (m_it_lo == masses_all.end() || m_it_hi == masses_all.end()) {
    throw std::runtime_error(
        FNERROR("Bin " + bin + ", process " + process +
                " does not have entries for the min/max mass points"));
  }

  // Define mass points = n
  //  -- masses_str[n]
  //  -- masses[n]
  vector<string> m_str_vec(m_it_lo, ++m_it_hi);
  vector<double> m_vec;
  for (auto const& s : m_str_vec) {
    // if (verbose) std::cout << ">>>> Mass point: " << s << "\n";
    m_vec.push_back(lexical_cast<double>(s));
  }
  unsigned m = m_vec.size();

  // ss = "shape systematic"
  vector<string> ss_vec =
      Set2Vec(cb_bp.cp().syst_type({"shape"}).syst_name_set());
  unsigned ss = ss_vec.size();

  // ls = "lnN systematic"
  vector<string> ls_vec =
      Set2Vec(cb_bp.cp().syst_type({"lnN"}).syst_name_set());
  unsigned ls = ls_vec.size();

  multi_array<ch::Process *,  1> pr_arr(extents[m]);
  multi_array<ch::Systematic *, 2> ss_arr(extents[ss][m]);
  multi_array<double, 1> ss_scale_arr(extents[ss]);
  multi_array<bool, 1> ss_must_scale_arr(extents[ss]);
  multi_array<ch::Systematic *, 2> ls_arr(extents[ls][m]);

  multi_array<std::shared_ptr<RooRealVar>, 1> ss_scale_var_arr(extents[ss]);
  multi_array<std::shared_ptr<RooConstVar>, 1> ss_scale_fac_arr(extents[ss]);
  multi_array<std::shared_ptr<RooProduct>, 1> ss_scale_prod_arr(extents[ss]);


  for (unsigned mi = 0; mi < m; ++mi) {
    cb_bp.cp().mass({m_str_vec[mi]}).ForEachProc([&](ch::Process *p) {
      pr_arr[mi] = p;
    });
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      cb_bp.cp().mass({m_str_vec[mi]}).syst_name({ss_vec[ssi]})
        .ForEachSyst([&](ch::Systematic *n) {
            ss_arr[ssi][mi] = n;
        });
    }
    for (unsigned lsi = 0; lsi < ls; ++lsi) {
      cb_bp.cp().mass({m_str_vec[mi]}).syst_name({ls_vec[lsi]})
        .ForEachSyst([&](ch::Systematic *n) {
            ls_arr[lsi][mi] = n;
        });
    }
  }

  RooArgList ss_list;
  for (unsigned ssi = 0; ssi < ss; ++ssi) {
    ss_scale_var_arr[ssi] =
        std::make_shared<RooRealVar>(ss_vec[ssi].c_str(), "", 0);
    set<double> scales;
    for (unsigned mi = 0; mi < m; ++mi) {
      scales.insert(ss_arr[ssi][mi]->scale());
    }
    if (scales.size() > 1) {
      std::runtime_error(FNERROR(
          "Shape morphing parameters that vary with mass are not allowed"));
    } else {
      ss_scale_arr[ssi] = *(scales.begin());
      if (std::fabs(ss_scale_arr[ssi] - 1.0) > 1E-6) {
        ss_must_scale_arr[ssi] = true;
        ss_scale_fac_arr[ssi] = std::make_shared<RooConstVar>(
            TString::Format("%g", ss_scale_arr[ssi]), "",
            ss_scale_arr[ssi]);
        ss_scale_prod_arr[ssi] = std::make_shared<RooProduct>(
            ss_vec[ssi] + "_scaled_" + key, "",
            RooArgList(*(ss_scale_var_arr[ssi]), *(ss_scale_fac_arr[ssi])));
        ss_list.add(*(ss_scale_prod_arr[ssi]));
      } else {
        ss_list.add(*(ss_scale_var_arr[ssi]));
      }
    }
  }
  // ss_list.Print();

  if (verbose) {
    std::cout << ">> Shape systematics: " << ss << "\n";
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      std::cout << boost::format("%-50s %-5i %-8.3g\n")
        % ss_vec[ssi] % ss_must_scale_arr[ssi] % ss_scale_arr[ssi];
    }
  }


  // lms = "lnN morphing systematic"
  vector<string> lms_vec;
  set<string> lms_set;
  vector<unsigned > lms_vec_idx;
  for (unsigned lsi = 0; lsi < ls; ++lsi) {
    set<double> k_hi;
    set<double> k_lo;
    for (unsigned mi = 0; mi < m; ++mi) {
      Systematic *n = ls_arr[lsi][mi];
      k_hi.insert(n->value_u());
      if (n->asymm()) {
        k_lo.insert(n->value_d());
      }
    }
    if (k_hi.size() > 1 || k_lo.size() > 1) {
      lms_vec.push_back(ls_vec[lsi]);
      lms_set.insert(ls_vec[lsi]);
      lms_vec_idx.push_back(lsi);
    }
  }
  unsigned lms = lms_vec.size();
  multi_array<ch::Systematic *, 2> lms_arr(extents[lms][m]);
  multi_array<std::shared_ptr<RooRealVar>, 1> lms_var_arr(extents[lms]);

  for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
    lms_var_arr[lmsi] =
        std::make_shared<RooRealVar>(lms_vec[lmsi].c_str(), "", 0);
  }

  for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
    for (unsigned mi = 0; mi < m; ++mi) {
      lms_arr[lmsi][mi] = ls_arr[lms_vec_idx[lmsi]][mi];
    }
  }
  if (verbose) {
    std::cout << ">> lnN morphing systematics: " << lms << "\n";
    for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
      std::cout << ">>>> " << lms_vec[lmsi] << "\n";
    }
  }

  multi_array<std::shared_ptr<TH1F>, 2> hist_arr(extents[m][1+ss*2]);
  multi_array<double, 1> rate_arr(extents[m]);
  multi_array<std::shared_ptr<TList>, 1> list_arr(extents[m]);
  multi_array<double, 2> ss_k_hi_arr(extents[ss][m]);
  multi_array<double, 2> ss_k_lo_arr(extents[ss][m]);
  multi_array<std::shared_ptr<RooSpline1D>, 1> ss_spl_hi_arr(extents[ss]);
  multi_array<std::shared_ptr<RooSpline1D>, 1> ss_spl_lo_arr(extents[ss]);
  multi_array<std::shared_ptr<AsymPow>, 1> ss_asy_arr(extents[ss]);

  multi_array<double, 2> lms_k_hi_arr(extents[ss][m]);
  multi_array<double, 2> lms_k_lo_arr(extents[ss][m]);
  multi_array<std::shared_ptr<RooSpline1D>, 1> lms_spl_hi_arr(extents[ss]);
  multi_array<std::shared_ptr<RooSpline1D>, 1> lms_spl_lo_arr(extents[ss]);
  multi_array<std::shared_ptr<AsymPow>, 1> lms_asy_arr(extents[ss]);

  for (unsigned mi = 0; mi < m; ++mi) {
    hist_arr[mi][0] =
        std::make_shared<TH1F>(RebinHist(AsTH1F(pr_arr[mi]->shape())));
    rate_arr[mi] = pr_arr[mi]->rate();
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      hist_arr[mi][1 + 2 * ssi] =
          std::make_shared<TH1F>(RebinHist(AsTH1F(ss_arr[ssi][mi]->shape_u())));
      hist_arr[mi][2 + 2 * ssi] =
          std::make_shared<TH1F>(RebinHist(AsTH1F(ss_arr[ssi][mi]->shape_d())));
      ss_k_hi_arr[ssi][mi] = ss_arr[ssi][mi]->value_u();
      ss_k_lo_arr[ssi][mi] = ss_arr[ssi][mi]->value_d();
    }
    for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
      lms_k_hi_arr[lmsi][mi] = lms_arr[lmsi][mi]->value_u();
      if (lms_arr[lmsi][mi]->asymm()) {
        lms_k_lo_arr[lmsi][mi] = lms_arr[lmsi][mi]->value_d();
      } else {
        lms_k_lo_arr[lmsi][mi] = 1. / lms_arr[lmsi][mi]->value_u();
      }
    }
  }

  for (unsigned mi = 0; mi < m; ++mi) {
    list_arr[mi] = std::make_shared<TList>();
    for (unsigned xi = 0; xi < (1 + ss * 2); ++xi) {
      list_arr[mi]->Add(hist_arr[mi][xi].get());
    }
  }

  if (verbose) {
    for (unsigned mi = 0; mi < m; ++mi) {
      std::cout << boost::format("%-10s") % m_str_vec[mi];
    }
    std::cout << "\n";
    for (unsigned mi = 0; mi < m; ++mi) {
      std::cout << boost::format("%-10.5g") % rate_arr[mi];
    }
    std::cout << "\n";
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      for (unsigned mi = 0; mi < m; ++mi) {
        std::cout << boost::format("%-10.5g") % ss_k_hi_arr[ssi][mi];
      }
      std::cout << "\n";
      for (unsigned mi = 0; mi < m; ++mi) {
        std::cout << boost::format("%-10.5g") % ss_k_lo_arr[ssi][mi];
      }
      std::cout << "\n";
    }
  }

  TString interp = "LINEAR";
  RooSpline1D rate_spline("interp_rate_"+key, "", mass_var, m, m_vec.data(),
                          rate_arr.data(), interp);
  RooArgList rate_prod(rate_spline);
  for (unsigned ssi = 0; ssi < ss; ++ssi) {
    ss_spl_hi_arr[ssi] = std::make_shared<RooSpline1D>("spline_hi_" +
        key + "_" + ss_vec[ssi], "", mass_var, m, m_vec.data(),
        ss_k_hi_arr[ssi].origin(), interp);
    ss_spl_lo_arr[ssi] = std::make_shared<RooSpline1D>("spline_lo_" +
        key + "_" + ss_vec[ssi], "", mass_var, m, m_vec.data(),
        ss_k_lo_arr[ssi].origin(), interp);
    ss_asy_arr[ssi] = std::make_shared<AsymPow>("systeff_" +
        key + "_" + ss_vec[ssi], "",
        *(ss_spl_lo_arr[ssi]), *(ss_spl_hi_arr[ssi]),
         *(reinterpret_cast<RooAbsReal*>(ss_list.at(ssi))));
    rate_prod.add(*(ss_asy_arr[ssi]));
  }
  for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
    lms_spl_hi_arr[lmsi] = std::make_shared<RooSpline1D>("spline_hi_" +
        key + "_" + lms_vec[lmsi], "", mass_var, m, m_vec.data(),
        lms_k_hi_arr[lmsi].origin(), interp);
    lms_spl_lo_arr[lmsi] = std::make_shared<RooSpline1D>("spline_lo_" +
        key + "_" + lms_vec[lmsi], "", mass_var, m, m_vec.data(),
        lms_k_lo_arr[lmsi].origin(), interp);
    lms_asy_arr[lmsi] = std::make_shared<AsymPow>("systeff_" +
        key + "_" + lms_vec[lmsi], "", *(lms_spl_lo_arr[lmsi]),
        *(lms_spl_hi_arr[lmsi]), *(lms_var_arr[lmsi]));
    rate_prod.add(*(lms_asy_arr[lmsi]));
  }

  TH1F data_hist = cb_bp.GetObservedShape();
  TH1F proc_hist = cb_bp.GetShape();
  RooRealVar xvar("CMS_th1x", "CMS_th1x", 0,
                 static_cast<float>(data_hist.GetNbinsX()));
  xvar.setBins(data_hist.GetNbinsX());

  RooRealVar morph_xvar(("CMS_th1x_"+bin).c_str(), "", 0,
                 static_cast<float>(proc_hist.GetNbinsX()));
  morph_xvar.setConstant();
  morph_xvar.setBins(proc_hist.GetNbinsX());
  if (verbose) {
    xvar.Print();
    morph_xvar.Print();
  }
  multi_array<std::shared_ptr<FastVerticalInterpHistPdf2>, 1> vpdf_arr(
      extents[m]);
  RooArgList vpdf_list;

  TString vert_name = key + "_";

  for (unsigned mi = 0; mi < m; ++mi) {
    vpdf_arr[mi] = std::make_shared<FastVerticalInterpHistPdf2>(
        vert_name + m_str_vec[mi] + "_vmorph", "", morph_xvar, *(list_arr[mi]),
        ss_list, 1, 0);
    vpdf_list.add(*(vpdf_arr[mi]));
  }
  TString morph_name = key + "_morph";
  RooMorphingPdf morph_pdf(morph_name, "", xvar, mass_var, vpdf_list,
                           m_vec, allow_morph, *(data_hist.GetXaxis()),
                           *(proc_hist.GetXaxis()));
  RooProduct morph_rate(morph_name + "_norm", "", rate_prod);

  ws.import(morph_pdf, RooFit::RecycleConflictNodes());
  ws.import(morph_rate, RooFit::RecycleConflictNodes());

  if (!verbose) RooMsgService::instance().setGlobalKillBelow(backup_msg_level);

  // Now we can cleanup the CB instance a bit
  cb.FilterProcs([&](ch::Process const* p) {
    return p->bin() == bin && p->process() == process && p->mass() != mass_min;
  });
  cb.FilterSysts([&](ch::Systematic const* n) {
    return (n->bin() == bin && n->process() == process) &&
           ((n->mass() != mass_min) || (n->type() == "shape") ||
            (lms_set.count(n->name())));
  });
  cb.ForEachProc([&](ch::Process * p) {
    if (p->bin() == bin && p->process() == process) {
      p->set_mass("*");
      p->set_shape(nullptr, false);
      p->set_rate(1.0);
    }
  });
  cb.ForEachSyst([&](ch::Systematic * n) {
    if (n->bin() == bin && n->process() == process) {
      n->set_mass("*");
    }
  });
}

void BuildRooMorphing(RooWorkspace& ws, CombineHarvester& cb, RooAbsReal& mh,
                      bool verbose, std::string norm_postfix) {

  // First step is to figure out what we can build with the CombineHarvester
  // instance we've been given by the user. The target will be to produce a
  // RooMorphingPdf for each process in each bin, using whichever mass points
  // are available in each case. The only requirements are that the `min` and
  // `max` points are available for each process.

  // Generate a list of bin/process combinations

  // Iterate through this list
    // Extract the set of mass points
    // check that `min` and `max` are present, and that we can covert
    // to a sorted list of doubles ok
    //
    // Get the list of systematics for this bin/process combination
    // Treatment is as follows, depending on type
    // "shape": must be present for all mass points. At each mass point
    // a list of Vars will be created, multiplied by "scale" if scale != 1. Will
    // "lnN": will check if it varies with mh. if it does: then we'll build another
    // Spline
  using std::set;
  using std::string;
  using std::vector;
  using boost::lexical_cast;

  // Get a vector of mass points and sort the string numerically
  vector<string> mass_str_vec = Set2Vec(cb.mass_set());
  std::sort(mass_str_vec.begin(), mass_str_vec.end(),
    [](string const& s1, string const& s2) {
      return lexical_cast<double>(s1) < lexical_cast<double>(s2);
    });

  // Fill a vector of mass points as doubles, aligned with the string version
  vector<double> mass_vec;
  for (auto const& s : mass_str_vec) {
    if (verbose) std::cout << ">> Mass point: " << s << "\n";
    mass_vec.push_back(lexical_cast<double>(s));
  }
  unsigned n = mass_vec.size();

  double mass_min = mass_vec.front();
  double mass_max = mass_vec.back();

  if (verbose) {
    std::cout << ">> Found " << mass_vec.size() << " mass points with min "
              << mass_min << " and max " << mass_max << "\n";
    mh.Print();
  }

  auto bins = cb.bin_set();

  for (auto const& b : bins) {
    auto sigs = cb.process_set();
    for (auto s : sigs) {
      if (verbose) std::cout << ">> bin: " << b << " process: " << s << "\n";
      ch::CombineHarvester tmp =
          std::move(cb.cp().bin({b}).process({s}).syst_type({"shape"}));
      TH1F data_hist = tmp.GetObservedShape();
      // tmp2.PrintAll();
      RooRealVar mtt("CMS_th1x", "CMS_th1x", 0,
                     static_cast<float>(data_hist.GetNbinsX()));
      mtt.setBins(data_hist.GetNbinsX());
      RooRealVar morph_mtt(("CMS_th1x_"+b).c_str(), "", 0,
                     static_cast<float>(tmp.GetShape().GetNbinsX()));
      morph_mtt.setConstant();
      morph_mtt.setBins(tmp.GetShape().GetNbinsX());
      if (verbose) {
        mtt.Print();
        morph_mtt.Print();
      }
      auto systs = ch::Set2Vec(tmp.syst_name_set());
      if (verbose) std::cout << ">> Found shape systematics:\n";
      RooArgList syst_list;
      for (auto const& syst: systs) {
        if (verbose) std::cout << ">>>> " << syst << "\n";
        syst_list.addClone(RooRealVar(syst.c_str(), syst.c_str(), 0));
      }
      vector<vector<TH1F>> h_vec(n);
      vector<TList> list_vec(n);
      vector<double> yield_vec(n);
      vector<FastVerticalInterpHistPdf2> v_pdfs;

      vector<vector<double>> k_vals_hi(systs.size(), vector<double>(n));
      vector<vector<double>> k_vals_lo(systs.size(), vector<double>(n));
      vector<RooSpline1D> k_splines_hi;
      vector<RooSpline1D> k_splines_lo;
      vector<AsymPow> k_asym;
      RooArgList k_list;

      std::string key = b + "_" + s;
      for (unsigned m = 0; m < n; ++m) {
        ch::CombineHarvester tmp3 =
            std::move(tmp.cp().mass({mass_str_vec[m]}));
        yield_vec[m] = tmp3.GetRate();
        h_vec[m].push_back(RebinHist(tmp3.GetShape()));
        for (unsigned s = 0; s < systs.size(); ++s) {
          tmp3.cp().syst_name({systs[s]}).ForEachSyst([&](Systematic const* n) {
            h_vec[m].push_back(RebinHist(*(TH1F*)(n->shape_u())));
            h_vec[m].push_back(RebinHist(*(TH1F*)(n->shape_d())));
            k_vals_hi[s][m] = n->value_u();
            k_vals_lo[s][m] = n->value_d();
            // if (std::fabs(n->scale() - 1.) >= 1E-6) {
            //   string scale_str = (boost::format("%g") % n->scale()).str();
            //   RooConstVar scale(scale_str.c_str(), "", n->scale());
            //   syst_list.addClone(RooProduct( (systs[s] + "_scaled_" + key).c_str(), "", RooArgList  ))
            // }
          });
        }
      }

      RooSpline1D yield((key + "_yield").c_str(), "", mh,
                        n, &(mass_vec[0]), &(yield_vec[0]), "LINEAR");
      k_list.add(yield);
      for (unsigned s = 0; s < systs.size(); ++s) {
        k_splines_hi.push_back(RooSpline1D((key + "_" + systs[s] + "_hi").c_str(), "", mh, n,
                        &(mass_vec[0]), &(k_vals_hi[s][0]), "LINEAR"));
        k_splines_lo.push_back(RooSpline1D((key + "_" + systs[s] + "_lo").c_str(), "", mh, n,
                        &(mass_vec[0]), &(k_vals_lo[s][0]), "LINEAR"));
      }
      for (unsigned s = 0; s < systs.size(); ++s) {
        k_asym.push_back(
            AsymPow((key + "_" + systs[s] + "_lnN").c_str(), "",
                    k_splines_lo[s], k_splines_hi[s],
                    *(reinterpret_cast<RooRealVar*>(syst_list.at(s)))));
      }
      for (unsigned s = 0; s < systs.size(); ++s) k_list.add(k_asym[s]);

      for (unsigned m = 0; m < n; ++m) {
        for (unsigned h = 0; h < h_vec[m].size(); ++h) {
          list_vec[m].Add(&h_vec[m][h]);
        }
      }

      for (unsigned m = 0; m < n; ++m) {
        v_pdfs.push_back(FastVerticalInterpHistPdf2(
            (b + "_" + s + mass_str_vec[m] + "_vpdf").c_str(), "", morph_mtt,
            list_vec[m], syst_list, 1, 0));
        // v_pdfs.back().Print();
      }
      RooArgList v_pdf_list;
      for (unsigned m = 0; m < n; ++m) v_pdf_list.add(v_pdfs[m]);

      RooMorphingPdf morph((b + "_" + s + "_mpdf").c_str(), "", mtt, mh,
                           v_pdf_list, mass_vec, true, *(data_hist.GetXaxis()),
                           *(tmp.GetShape().GetXaxis()));

      RooProduct full_yield((b + "_" + s + "_mpdf" + norm_postfix).c_str(), "",
                            k_list);

      ws.import(morph, RooFit::RecycleConflictNodes());
      ws.import(full_yield, RooFit::RecycleConflictNodes());
    }
  }
  // ws.Print();
}
}
