#include "CombineHarvester/CombinePdfs/interface/CMSHistFuncFactory.h"
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"
#include "boost/multi_array.hpp"
#include "TVector.h"
#include "TGraphErrors.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooProduct.h"
#include "RooConstVar.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace ch {

/*
TODO:
2) Morphing options
3) Remove morphing proc errors
*/

CMSHistFuncFactory::CMSHistFuncFactory() : v_(1), hist_mode_(0), rebin_(false) {}

void CMSHistFuncFactory::Run(ch::CombineHarvester &cb, RooWorkspace &ws, std::map<std::string, std::string> process_vs_norm_postfix_map) {
  for (auto const& bin : cb.bin_set()) {
    for (auto const& proc : cb.cp().bin({bin}).process_set()) {
      if (v_) {
        std::cout << ">> Processing " << bin << "," << proc << "\n";
      }
      if (process_vs_norm_postfix_map.find(proc) == process_vs_norm_postfix_map.end())
      {
        norm_postfix_ = "norm";
      }
      else
      {
        norm_postfix_ = process_vs_norm_postfix_map[proc];
      }
      RunSingleProc(cb, ws, bin, proc);
    }
    TH1F data_hist = cb.cp().bin({bin}).GetObservedShape();
    if (rebin_) data_hist = RebinHist(data_hist);
    // data_hist.Print("range");

    RooRealVar *xvar = ws.var(TString::Format("CMS_x_%s", bin.c_str()));
    RooDataHist rdh_dat(TString::Format("%s_data_obs", bin.c_str()), "",
                        RooArgList(*xvar), &data_hist);

    ws.import(rdh_dat);

    cb.cp().bin({bin}).ForEachObs([&](ch::Observation * p) {
      p->set_shape(nullptr, false);
      p->set_rate(1.0);
    });
  }
}

void CMSHistFuncFactory::Run(ch::CombineHarvester &cb, RooWorkspace &ws) {
  CMSHistFuncFactory::Run(cb, ws,  {});
}

void CMSHistFuncFactory::RunSingleProc(CombineHarvester& cb, RooWorkspace& ws,
                                       std::string bin, std::string process) {
  using std::vector;
  using std::set;
  using std::string;
  using boost::lexical_cast;
  using boost::multi_array;
  using boost::extents;

  TString key = bin + "_" + process;


  CombineHarvester cbp = cb.cp().bin({bin}).process({process});
  vector<string> m_str_vec = Set2Vec(cbp.SetFromProcs(
      std::mem_fn(&ch::Process::mass)));
  unsigned m = m_str_vec.size();
  vector<double> m_vec;
  if (m_str_vec.size() == 1) {
  } else {
    std::sort(m_str_vec.begin(), m_str_vec.end(),
              [](string const& s1, string const& s2) {
                return lexical_cast<double>(s1) < lexical_cast<double>(s2);
              });
    for (auto const& s : m_str_vec) {
      if (v_) std::cout << ">>>> Mass point: " << s << "\n";
      m_vec.push_back(lexical_cast<double>(s));
    }
  }

  // ss = "shape systematic"
  // Make a list of the names of shape systematics affecting this process
  vector<string> ss_vec =
      Set2Vec(cbp.cp().syst_type({"shape", "shapeU"}).syst_name_set());
  // Now check if all shape systematics are present for all mass points
  for (auto const& s : m_str_vec) {
    if (cbp.cp().syst_type({"shape", "shapeU"}).mass({s}).syst_name_set().size() !=
        ss_vec.size()) {
      throw std::runtime_error(FNERROR(
          "Some mass points do not have the full set of shape systematics, "
          "this is currently unsupported"));
    }
  }
  unsigned ss = ss_vec.size();  // number of shape systematics

  vector<string> ls_vec =
      Set2Vec(cbp.cp().syst_type({"lnN"}).syst_name_set());
  unsigned ls = ls_vec.size();  // number of lnN systematics

  // Store pointers to each ch::Process (one per mass point) in the CH instance
  multi_array<ch::Process *,  1> pr_arr(extents[m]);
  // Store pointers pointers to each ch::Systematic for each mass point (hence
  // an ss * m 2D array)
  multi_array<ch::Systematic *, 2> ss_arr(extents[ss][m]);
  // With shape systematics we have to support cases where the value in the
  // datacard is != 1.0, i.e. we are scaling the parameter that goes into the
  // Gaussian constraint PDF - we have to pass this factor on when we build the
  // normal vertical-interpolation PDF for each mass point. We will store these
  // scale factors in this array, one per shape systematic.
  multi_array<double, 1> ss_scale_arr(extents[ss]);
  // Really just for book-keeping, we'll set this flag to true when the shape
  // systematic scale factor != 1
  multi_array<bool, 1> ss_must_scale_arr(extents[ss]);
  // Similar to the ss_arr above, store the lnN ch::Systematic objects. Note
  // implicit in all this is the assumption that the processes at each mass
  // point have exactly the same list of systematic uncertainties.
  multi_array<ch::Systematic *, 2> ls_arr(extents[ls][m]);

  // This array holds pointers to the RooRealVar objects that will become our
  // shape nuisance parameters, e.g. "CMS_scale_t_mutau_8TeV"
  multi_array<std::shared_ptr<RooRealVar>, 1> ss_scale_var_arr(extents[ss]);
  // And this array holds the constant scale factors that we'll build from
  // ss_scale_arr above
  multi_array<std::shared_ptr<RooConstVar>, 1> ss_scale_fac_arr(extents[ss]);
  // Finally this array will contain the scale_var * scale_fac product where we
  // need to provide a scaled value of the nuisance parameter instead of the
  // parameter itself in building the vertical-interp. PDF
  multi_array<std::shared_ptr<RooProduct>, 1> ss_scale_prod_arr(extents[ss]);

    // Now let's fill some of these arrays...
    for (unsigned mi = 0; mi < m; ++mi) {
      // The ch::Process pointers
      cbp.cp().mass({m_str_vec[mi]}).ForEachProc([&](ch::Process *p) {
        pr_arr[mi] = p;
      });
      for (unsigned ssi = 0; ssi < ss; ++ssi) {
        // The ch::Systematic pointers for shape systematics
        cbp.cp().mass({m_str_vec[mi]}).syst_name({ss_vec[ssi]})
          .ForEachSyst([&](ch::Systematic *n) {
              ss_arr[ssi][mi] = n;
          });
      }
      for (unsigned lsi = 0; lsi < ls; ++lsi) {
        // The ch::Systematic pointers for lnN systematics
        cbp.cp().mass({m_str_vec[mi]}).syst_name({ls_vec[lsi]})
          .ForEachSyst([&](ch::Systematic *n) {
              ls_arr[lsi][mi] = n;
          });
      }
    }

    //! [part3]
    // We need to build a RooArgList of the vertical morphing parameters for the
    // vertical-interpolation pdf - this will be the same for each mass point so
    // we only build it once
    RooArgList ss_list;
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      // Make the nuisance parameter var. We use shared_ptrs here as a convenience
      // because they will take care of automatically cleaning up the memory at
      // the end
      ss_scale_var_arr[ssi] =
          std::make_shared<RooRealVar>(ss_vec[ssi].c_str(), "", 0);
      // We'll make a quick check that the scale factor for this systematic is the
      // same for all mass points. We could do a separate scaling at each mass
      // point but this would create a lot of complications
      set<double> scales;
      // Insert the scale from each mass point into the set, if it has a size
      // larger than one at the end then we have a problem!
      for (unsigned mi = 0; mi < m; ++mi) {
        scales.insert(ss_arr[ssi][mi]->scale());
      }
      if (scales.size() > 1) {
        // Don't let the user proceed, we can't build the model they want
        std::runtime_error(FNERROR(
            "Shape morphing parameters that vary with mass are not allowed"));
      } else {
        // Everything ok, set the scale value in its array
        ss_scale_arr[ssi] = *(scales.begin());
        // Handle the case where the scale factor is != 1
        if (std::fabs(ss_scale_arr[ssi] - 1.0) > 1E-6) {
          ss_must_scale_arr[ssi] = true;
          // Build the RooConstVar with the value of the scale factor
          ss_scale_fac_arr[ssi] = std::make_shared<RooConstVar>(
              TString::Format("%g", ss_scale_arr[ssi]), "",
              ss_scale_arr[ssi]);
          // Create the product of the floating nuisance parameter and the
          // constant scale factor
          ss_scale_prod_arr[ssi] = std::make_shared<RooProduct>(
              ss_vec[ssi] + "_scaled_" + key, "",
              RooArgList(*(ss_scale_var_arr[ssi]), *(ss_scale_fac_arr[ssi])));
          // Add this to the list
          ss_list.add(*(ss_scale_prod_arr[ssi]));
        } else {
          // If the scale factor is 1.0 then we just add the nuisance parameter
          // directly to our list
          ss_list.add(*(ss_scale_var_arr[ssi]));
        }
      }
    }
    //! [part3]
  // Summarise the info on the shape systematics and scale factors
  if (v_) {
    std::cout << ">> Shape systematics: " << ss << "\n";
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      std::cout << boost::format("%-50s %-5i %-8.3g\n")
        % ss_vec[ssi] % ss_must_scale_arr[ssi] % ss_scale_arr[ssi];
    }
  }

    // lms = "lnN morphing systematic"
    // Now we have some work to do with the lnN systematics. We can consider two cases:
    //  a) The uncertainty is the same for each mass point => we can leave it in
    //     the datacard as is and let text2workspace do its normal thing
    //  b) The uncertainty varies between mass points => we can't capture this
    //     information in the text datacard in the usual way, so we'll build a RooFit
    //     object that effectively makes the lnN uncertainty a function of the mass
    //     variable
    // We'll use "lms" to refer to case b), which we'll try to figure out now...
    vector<string> lms_vec;  // vec of systematic names
    set<string> lms_set;     // set of systematic names
    // index positions in our full ls_arr array for the lms systematics
    vector<unsigned > lms_vec_idx;
    for (unsigned lsi = 0; lsi < ls; ++lsi) {
      // Extra complication is that the user might have been evil and mixed
      // symmetric and asymmetric lnN values, we'll try and detect changes in
      // either
      set<double> k_hi;
      set<double> k_lo;
      // Go through each mass point for this systematic and add the uncertainty
      // values (so-called "kappa" values)
      for (unsigned mi = 0; mi < m; ++mi) {
        Systematic *n = ls_arr[lsi][mi];
        k_hi.insert(n->value_u());
        if (n->asymm()) {
          k_lo.insert(n->value_d());
        }
      }
      // If either of these sets has more than one entry then this is a lms case
      if (k_hi.size() > 1 || k_lo.size() > 1) {
        lms_vec.push_back(ls_vec[lsi]);
        lms_set.insert(ls_vec[lsi]);
        lms_vec_idx.push_back(lsi);
      }
    }
    unsigned lms = lms_vec.size();
    // New array for the pointers to the lms Systematic objects
    multi_array<ch::Systematic *, 2> lms_arr(extents[lms][m]);
    for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
      for (unsigned mi = 0; mi < m; ++mi) {
        lms_arr[lmsi][mi] = ls_arr[lms_vec_idx[lmsi]][mi];
      }
    }

    // We will need to create the nuisance parameters for these now
    multi_array<std::shared_ptr<RooRealVar>, 1> lms_var_arr(extents[lms]);
    for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
      lms_var_arr[lmsi] =
          std::make_shared<RooRealVar>(lms_vec[lmsi].c_str(), "", 0);
    }

    // Give a summary of the lms systematics to the user
    if (v_) {
      std::cout << ">> lnN morphing systematics: " << lms << "\n";
      for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
        std::cout << ">>>> " << lms_vec[lmsi] << "\n";
      }
    }

    // 2D array of all input histograms, size is (mass points * (nominal +
    // 2*shape-systs)). The factor of 2 needed for the Up and Down shapes
    multi_array<std::shared_ptr<TH1F>, 2> hist_arr(extents[m][1+ss*2]);
    // We also need the array of process yields vs mass, because this will have to
    // be interpolated too
    multi_array<double, 1> rate_arr(extents[m]);
    // Also store the yield uncertainty - we don't actually need this for the signal
    // model, but will include it in the debug TGraph
    multi_array<double, 1> rate_unc_arr(extents[m]);
    // The vertical-interpolation PDF needs the TH1 inputs in the format of a TList
    multi_array<std::shared_ptr<TList>, 1> list_arr(extents[m]);
    // Combine always treats the normalisation part of shape systematics as
    // distinct from the actual shape morphing. Essentially the norm part is
    // treated as an asymmetric lnN. We have to make the kappa_hi and kappa_lo a
    // function of the mass parameter too, so we need two more arrays in (ss * m)
    multi_array<double, 2> ss_k_hi_arr(extents[ss][m]);
    multi_array<double, 2> ss_k_lo_arr(extents[ss][m]);
    // For each shape systematic we will build a RooSpline1D, configured to
    // interpolate linearly between the kappa values
    multi_array<std::shared_ptr<RooAbsReal>, 1> ss_spl_hi_arr(extents[ss]);
    multi_array<std::shared_ptr<RooAbsReal>, 1> ss_spl_lo_arr(extents[ss]);
    // To define the actually process scaling as a function of these kappa values,
    // we need an AsymPow object per mass point
    multi_array<std::shared_ptr<AsymPow>, 1> ss_asy_arr(extents[ss]);

    // Similar set of objects needed for the lms normalisation systematics
    multi_array<double, 2> lms_k_hi_arr(extents[lms][m]);
    multi_array<double, 2> lms_k_lo_arr(extents[lms][m]);
    multi_array<std::shared_ptr<RooSpline1D>, 1> lms_spl_hi_arr(extents[lms]);
    multi_array<std::shared_ptr<RooSpline1D>, 1> lms_spl_lo_arr(extents[lms]);
    multi_array<std::shared_ptr<AsymPow>, 1> lms_asy_arr(extents[lms]);

    for (unsigned mi = 0; mi < m; ++mi) {
      // Grab the nominal process histograms. We also have to convert every
      // histogram to a uniform integer binning, because this is what
      // text2workspace will do for all the non-morphed processes in our datacard,
      // and we need the binning of these to be in sync.
      hist_arr[mi][0] =
          std::make_shared<TH1F>(AsTH1F(pr_arr[mi]->ClonedScaledShape().get()));
      if (rebin_) *hist_arr[mi][0] = RebinHist(*hist_arr[mi][0]);
      if (m > 1) {
        for (int b = 1; b < hist_arr[mi][0]->GetNbinsX() + 1; ++b) {
          hist_arr[mi][0]->SetBinError(b, 0.);
        }
      }
      // If the user supplied a TFile pointer we'll dump a bunch of info into it
      // for debugging
      // if (file) {
      //   file->WriteTObject(pr_arr[mi]->shape(), key + "_" + m_str_vec[mi]);
      // }
      // Store the process rate
      rate_arr[mi] = 1.;
      // auto proc_hist = pr_arr[mi]->ClonedScaledShape();
      // proc_hist->IntegralAndError(1, proc_hist->GetNbinsX(), rate_unc_arr[mi]);
      // Do the same for the Up and Down shapes
      for (unsigned ssi = 0; ssi < ss; ++ssi) {
        hist_arr[mi][1 + 2 * ssi] =
            std::make_shared<TH1F>(AsTH1F(ss_arr[ssi][mi]->shape_u()));
        hist_arr[mi][2 + 2 * ssi] =
            std::make_shared<TH1F>(AsTH1F(ss_arr[ssi][mi]->shape_d()));
        if (rebin_) *hist_arr[mi][1 + 2 * ssi] = RebinHist(*hist_arr[mi][1 + 2 * ssi]);
        if (rebin_) *hist_arr[mi][2 + 2 * ssi] = RebinHist(*hist_arr[mi][2 + 2 * ssi]);
        TH1F* h_hi = hist_arr[mi][1 + 2 * ssi].get();
        TH1F* h_lo = hist_arr[mi][2 + 2 * ssi].get();
        if (h_hi->Integral() > 0.) {
          h_hi->Scale(hist_arr[mi][0]->Integral() / h_hi->Integral());
        }
        if (h_lo->Integral() > 0.) {
          h_lo->Scale(hist_arr[mi][0]->Integral() / h_lo->Integral());
        }

        // TODO: Need to implement the logic of normalising shape variations here

        // if (file) {
        //   file->WriteTObject(ss_arr[ssi][mi]->shape_u(),
        //                      key + "_" + m_str_vec[mi] + "_" + ss_vec[ssi] + "Up");
        //   file->WriteTObject(ss_arr[ssi][mi]->shape_d(),
        //                      key + "_" + m_str_vec[mi] + "_" + ss_vec[ssi] + "Down");
        // }
        // Store the uncertainty ("kappa") values for the shape systematics
        ss_k_hi_arr[ssi][mi] = ss_arr[ssi][mi]->value_u();
        ss_k_lo_arr[ssi][mi] = ss_arr[ssi][mi]->value_d();
        // For the normalisation we scale the kappa instead of putting the scaling
        // parameter as the variable
        if (std::fabs(ss_scale_arr[ssi] - 1.0) > 1E-6) {
          ss_k_hi_arr[ssi][mi] = std::pow(ss_arr[ssi][mi]->value_u(), ss_scale_arr[ssi]);
          ss_k_lo_arr[ssi][mi] = std::pow(ss_arr[ssi][mi]->value_d(), ss_scale_arr[ssi]);
        }
      }
      // And now the uncertainty values for the lnN systematics that vary with mass
      // We'll force these to be asymmetric even if they're not
      for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
        lms_k_hi_arr[lmsi][mi] = lms_arr[lmsi][mi]->value_u();
        if (lms_arr[lmsi][mi]->asymm()) {
          lms_k_lo_arr[lmsi][mi] = lms_arr[lmsi][mi]->value_d();
        } else {
          lms_k_lo_arr[lmsi][mi] = 1. / lms_arr[lmsi][mi]->value_u();
        }
      }
    }
    // Now we'll fill these objects..


    // Now we've made all our histograms, we'll put them in the TList format
    // of [nominal, syst_1_Up, syst_1_Down, ... , syst_N_Up, syst_N_Down]
    // for (unsigned mi = 0; mi < m; ++mi) {
    //   list_arr[mi] = std::make_shared<TList>();
    //   for (unsigned xi = 0; xi < (1 + ss * 2); ++xi) {
    //     list_arr[mi]->Add(hist_arr[mi][xi].get());
    //   }
    // }

    // Print the values of the yields and kappa factors that will be inputs
    // to our spline interpolation objects
    if (v_) {
      for (unsigned mi = 0; mi < m; ++mi) {
        std::cout << boost::format("%-10s") % m_str_vec[mi];
      }
      std::cout << "\n";
      for (unsigned mi = 0; mi < m; ++mi) {
        std::cout << boost::format("%-10.5g") % rate_arr[mi];
      }
      std::cout << "\n";
      for (unsigned ssi = 0; ssi < ss; ++ssi) {
        std::cout << ss_vec[ssi] << " Up" << std::endl;
        for (unsigned mi = 0; mi < m; ++mi) {
          std::cout << boost::format("%-10.5g") % ss_k_hi_arr[ssi][mi];
        }
        std::cout << "\n";
        std::cout << ss_vec[ssi] << " Down" << std::endl;
        for (unsigned mi = 0; mi < m; ++mi) {
          std::cout << boost::format("%-10.5g") % ss_k_lo_arr[ssi][mi];
        }
        std::cout << "\n";
      }
      for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
        std::cout << lms_vec[lmsi] << " Up" << std::endl;
        for (unsigned mi = 0; mi < m; ++mi) {
          std::cout << boost::format("%-10.5g") % lms_k_hi_arr[lmsi][mi];
        }
        std::cout << "\n";
        std::cout << lms_vec[lmsi] << " Down" << std::endl;
        for (unsigned mi = 0; mi < m; ++mi) {
          std::cout << boost::format("%-10.5g") % lms_k_lo_arr[lmsi][mi];
        }
        std::cout << "\n";
      }
    }

    // Can do more sophistical spline interpolation if we want, but let's use
    // simple LINEAR interpolation for now
    TString interp = "LINEAR";

    // Here when the force_template_limit option is requested 
    // we have to add some extra terms to the vector of masses to ensure that 
    // the signal pdf goes to 0 outside of the MC template range.

    vector<double> new_m_vec(m_vec);
    // Insert an entry at either end of the vector for a mass just slightly
    // outside of the range
    if (m > 1) {
      new_m_vec.insert(new_m_vec.begin(),m_vec[0]-1E-6);
      new_m_vec.push_back(m_vec[m-1]+1E-6);
    }
    // Create a corresponding rate array with 0 entries for these new masses
    multi_array<double, 1> new_rate_arr(extents[m+2]);
    new_rate_arr[0] = 0.0;
    // for(unsigned i = 0; i < m; ++i) new_rate_arr[i+1] = rate_arr[i] ;

    // With CMSHistFunc the normalisation is handled internally
    for(unsigned i = 0; i < m; ++i) new_rate_arr[i+1] = 1.0;
    new_rate_arr[m+1] = 0.0;

    bool force_template_limit = false;

    // if (v_ && force_template_limit) {
    //    std::cout << ">>>> Forcing rate to 0 outside of template range:" << "\n";
    //    for(unsigned mi = 0; mi < m+2; ++mi) {
    //       std::cout << boost::format("%-10.5g") % new_m_vec[mi];
    //    }
    //    std::cout << "\n";
    //    for(unsigned mi = 0; mi < m+2; ++mi) {
    //       std::cout << boost::format("%-10.5g") % new_rate_arr[mi];
    //    }
    //    std::cout << "\n";
    // }
    // Create the 1D spline directly from the rate array
    //! [part4]
    std::shared_ptr<RooSpline1D> rate_spline;
    if (m > 1) rate_spline = std::make_shared<RooSpline1D>("interp_rate_"+key, "", *mass_var[process],
                          force_template_limit ? m+2 : m,
                          force_template_limit ? new_m_vec.data() : m_vec.data(),
                          force_template_limit ? new_rate_arr.data() : rate_arr.data(),
                          interp);
    // if (file) {
    //   TGraphErrors tmp(m, m_vec.data(), rate_arr.data(), nullptr, rate_unc_arr.data());
    //   file->WriteTObject(&tmp, "interp_rate_"+key);
    // }
    // Collect all terms that will go into the total normalisation:
    //   nominal * systeff_1 * systeff_2 * ... * systeff_N
    RooArgList rate_prod;
    if (m > 1) {
      rate_prod.add(*rate_spline);
    }
    // For each shape systematic build a 1D spline for kappa_hi and kappa_lo
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      if (m > 1) {
        ss_spl_hi_arr[ssi] = std::make_shared<RooSpline1D>("spline_hi_" +
            key + "_" + ss_vec[ssi], "", *mass_var[process], m, m_vec.data(),
            ss_k_hi_arr[ssi].origin(), interp);
        ss_spl_lo_arr[ssi] = std::make_shared<RooSpline1D>("spline_lo_" +
            key + "_" + ss_vec[ssi], "", *mass_var[process], m, m_vec.data(),
            ss_k_lo_arr[ssi].origin(), interp);
      } else {
        ss_spl_hi_arr[ssi] = std::make_shared<RooConstVar>(TString::Format("%g", ss_k_hi_arr[ssi][0]), "", ss_k_hi_arr[ssi][0]);
        ss_spl_lo_arr[ssi] = std::make_shared<RooConstVar>(TString::Format("%g", ss_k_lo_arr[ssi][0]), "", ss_k_lo_arr[ssi][0]);
      }
      ss_asy_arr[ssi] = std::make_shared<AsymPow>("systeff_" +
          key + "_" + ss_vec[ssi], "",
          *(ss_spl_lo_arr[ssi]), *(ss_spl_hi_arr[ssi]),
           *(ss_scale_var_arr[ssi]));
      // if (file) {
      //   TGraph tmp_hi(m, m_vec.data(), ss_k_hi_arr[ssi].origin());
      //   file->WriteTObject(&tmp_hi, "spline_hi_" + key + "_" + ss_vec[ssi]);
      //   TGraph tmp_lo(m, m_vec.data(), ss_k_lo_arr[ssi].origin());
      //   file->WriteTObject(&tmp_lo, "spline_lo_" + key + "_" + ss_vec[ssi]);
      // }
      // Then build the AsymPow object for each systematic as a function of the
      // kappas and the nuisance parameter
      rate_prod.add(*(ss_asy_arr[ssi]));
    }
    // Same procedure for the lms normalisation systematics: build the splines
    // then the AsymPows and add to the rate_prod list
    for (unsigned lmsi = 0; lmsi < lms; ++lmsi) {
      lms_spl_hi_arr[lmsi] = std::make_shared<RooSpline1D>("spline_hi_" +
          key + "_" + lms_vec[lmsi], "", *mass_var[process], m, m_vec.data(),
          lms_k_hi_arr[lmsi].origin(), interp);
      lms_spl_lo_arr[lmsi] = std::make_shared<RooSpline1D>("spline_lo_" +
          key + "_" + lms_vec[lmsi], "", *mass_var[process], m, m_vec.data(),
          lms_k_lo_arr[lmsi].origin(), interp);
      lms_asy_arr[lmsi] = std::make_shared<AsymPow>("systeff_" +
          key + "_" + lms_vec[lmsi], "", *(lms_spl_lo_arr[lmsi]),
          *(lms_spl_hi_arr[lmsi]), *(lms_var_arr[lmsi]));
      rate_prod.add(*(lms_asy_arr[lmsi]));
    }
    // We'll come back to this rate_prod a bit later.

    // Now some fun with binning. We anticipate that the process histograms we
    // have been supplied could have a finer binning than is actually wanted for
    // the analysis (and the fit), in order to avoid known problems with the RMS
    // of a peaking distribution not being morphed smoothly from mass point to
    // mass point if the binning is too wide. The RooMorphingPdf will handle
    // re-binning on the fly, but we have to tell it how to rebin. To do this we
    // assume the observed data histogram has the target binning.
    TH1F data_hist = cbp.GetObservedShape();
    TH1F proc_hist = cbp.GetShape();
    // The x-axis variable has to be called "CMS_th1x", as this is what
    // text2workspace will use for all the normal processes
    // data_hist.Print("range");
    TString xvar_name = TString::Format("CMS_x_%s", bin.c_str());
    RooRealVar xvar = VarFromHist(xvar_name, xvar_name, data_hist);
    // RooRealVar xvar("CMS_th1x", "CMS_th1x", 0,
    //                static_cast<float>(data_hist.GetNbinsX()));
    // xvar.setBins(data_hist.GetNbinsX());

    // Create a second x-axis variable, named specific to the bin, that will be
    // for the finer-binned input
    // RooRealVar morph_xvar(("CMS_th1x_"+bin).c_str(), "", 0,
    //                static_cast<float>(proc_hist.GetNbinsX()));
    // We're not going to need roofit to evaluate anything as a function of this
    // morphing x-axis variable, so we set it constant
    // morph_xvar.setConstant();
    // morph_xvar.setBins(proc_hist.GetNbinsX());
    if (v_) {
      xvar.Print();
      // morph_xvar.Print();
    }

    // Follow what ShapeTools.py does and set the smoothing region
    // to the minimum of all of the shape scales
    double qrange = 1.;
    for (unsigned ssi = 0; ssi < ss; ++ssi) {
      if (ss_scale_arr[ssi] < qrange) qrange = ss_scale_arr[ssi];
    }


    TString morph_name = key + "_morph";
    // At long last, we can build our pdf, giving it:
    //   xvar :       the fixed "CMS_th1x" x-axis variable with uniform integer binning
    //   mass_var:    the floating mass value for the interpolation
    //   vpdf_list:   the list of vertical-interpolation pdfs
    //   m_vec:       the corresponding list of mass points
    //   allow_morph: if false will just evaluate to the closest pdf in mass
    //   data_hist.GetXaxis(): The original (non-uniform) target binning
    //   proc_hist.GetXaxis(): The original (non-uniform) morphing binning
    //! [part5]
    CMSHistFunc morph_func(morph_name, "", xvar, *hist_arr[0][0]);
    morph_func.setVerticalSmoothRegion(qrange);
    morph_func.setHorizontalType(CMSHistFunc::HorizontalType::Integral);
    morph_func.setVerticalMorphs(ss_list);
    if (m > 1) {
      morph_func.addHorizontalMorph(*mass_var[process], TVectorD(m, m_vec.data()));
    }
    morph_func.prepareStorage();

    for (unsigned mi = 0; mi < m; ++mi) {
      morph_func.setShape(0, mi, 0, 0, *hist_arr[mi][0]);
      for (unsigned ssi = 0; ssi < ss; ++ssi) {
        morph_func.setShape(0, mi, ssi + 1, 0, *hist_arr[mi][2 + 2 * ssi]);
        morph_func.setShape(0, mi, ssi + 1, 1, *hist_arr[mi][1 + 2 * ssi]);
      }
    }
    // And we can make the final normalisation product
    // The value of norm_postfix is very important. text2workspace will only look
    // for for a term with the pdf name + "_norm". But it might be the user wants
    // to add even more terms to this total normalisation, so we give them the option
    // of using some other suffix.
    //! [part6]
    RooProduct morph_rate(morph_name + "_" + TString(norm_postfix_), "",
                          rate_prod);
    //! [part6]

    // Dump even more plots
    // if (file) MakeMorphDebugPlots(&morph_pdf, &mass_var, m_vec, file, &data_hist);

    // Load our pdf and norm objects into the workspace
    ws.import(morph_func, RooFit::RecycleConflictNodes());
    if (rate_prod.getSize() > 0) {
      ws.import(morph_rate, RooFit::RecycleConflictNodes());
    } else {
      std::cout << "No normalisation term to import!\n";
    }

    // if (!verbose) RooMsgService::instance().setGlobalKillBelow(backup_msg_level);
// 
    // Now we can clean up the CH instance a bit
    // We only need one entry for each Process or Systematic now, not one per mass point
    // We'll modify the first mass point to house our new pdf and norm, and drop
    // the rest.
    std::string mass_min = m_str_vec.at(0);

    // Dump Process entries for other mass points
    cb.FilterProcs([&](ch::Process const* p) {
      return p->bin() == bin && p->process() == process && p->mass() != mass_min;
    });
    // Dump Systematic entries for other mass points, but only if the type is shape
    // or a lnN that we found varied as a function of the mass. We've already built
    // these uncertainties into our normalisation term. Constant lnN systematics can
    // remain in the datacard and be added by text2workspace
    cb.FilterSysts([&](ch::Systematic const* n) {
      return (n->bin() == bin && n->process() == process) &&
             ((n->mass() != mass_min) || (n->type() == "shape" || n->type() == "shapeU") ||
              (lms_set.count(n->name())));
    });
    // With the remaining Process entry (should only be one if we did this right),
    // Make the mass generic ("*"), drop the TH1 and set the rate to 1.0, as this
    // will now be read from our norm object
    cb.ForEachProc([&](ch::Process * p) {
      if (p->bin() == bin && p->process() == process) {
        if (m > 1) {
          p->set_mass("*");
        }
        p->set_shape(nullptr, false);
        p->set_rate(1.0);
      }
    });
    // Just declare the mass to be generic for the remaining systematics
    if (m > 1) {
      cb.ForEachSyst([&](ch::Systematic * n) {
        if (n->bin() == bin && n->process() == process) {
          n->set_mass("*");
        }
      });
    }

    // And we're done, but note that we haven't populated the Process entry with
    // the PDF or norm objects we created, as we assume the user has more work to
    // do on their RooWorkspace before copying into the CH instance. Once they've
    // imported the workspace:
    //
    //    cb.AddWorkspace(ws);
    //
    // They can populate the Process entries, e.g. if all signals were replaced
    // with Morphing PDFs:
    // 
    //   cb.cp().signals().ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");
    //! [part4]
}

// }

// void MakeMorphDebugPlots(RooMorphingPdf* pdf, RooAbsReal* mass,
//                          std::vector<double> const& masses, TFile* f, TH1 *ref_bins) {
//   RooRealVar *rrv = dynamic_cast<RooRealVar*>(mass);
//   if (!rrv) return;
//   f->cd();
//   f->mkdir(pdf->GetName());
//   gDirectory->cd(pdf->GetName());
//   for (unsigned i = 0; i < masses.size(); ++i) {
//     rrv->setVal(masses[i]);
//     TH1 * h = pdf->createHistogram("CMS_th1x");
//     h->AddDirectory(false);
//     TH1 * h2 = nullptr;
//     if (ref_bins) { 
//       h2 = (TH1*)ref_bins->Clone();
//       h2->Reset();
//       for (int x = 1; x <= h2->GetNbinsX(); ++x) {
//         h2->SetBinContent(x, h->GetBinContent(x));
//       }
//     } else {
//       h2 = h;
//       h = nullptr;
//     }
//     h2->AddDirectory(false);
//     h2->SetName(TString::Format("actual_point_%g", masses[i]));
//     gDirectory->WriteTObject(h2);
//     if (h) delete h;
//     if (h2) delete h2;
//   }
//   double m = masses.front();
//   double step = 1;
//   if (((masses.back() - masses.front()) / step) > 100.) step = step * 10.;
//   if (((masses.back() - masses.front()) / step) > 100.) step = step * 10.;
//   if (((masses.back() - masses.front()) / step) < 10.) step = step/10.;
//   while (m <= masses.back()) {
//     rrv->setVal(m);
//     TH1* hm = pdf->createHistogram("CMS_th1x");
//     hm->AddDirectory(false);
//     TH1 * hm2 = nullptr;
//     if (ref_bins) { 
//       hm2 = (TH1*)ref_bins->Clone();
//       hm2->Reset();
//       for (int x = 1; x <= hm2->GetNbinsX(); ++x) {
//         hm2->SetBinContent(x, hm->GetBinContent(x));
//       }
//     } else {
//       hm2 = hm;
//       hm = nullptr;
//     }
//     hm2->AddDirectory(false);
//     hm2->SetName(TString::Format("morph_point_%g", m));
//     //It struggles with m=0, instead adds really small number close to 0
//     if(fabs(m) < 1E-5) hm2->SetName(TString::Format("morph_point_0"));
//     gDirectory->WriteTObject(hm2);
//     if (hm) delete hm;
//     if (hm2) delete hm2;
//     m += step;
//   }
//   f->cd();
// }
}
