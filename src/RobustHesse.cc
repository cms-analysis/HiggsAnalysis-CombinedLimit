#include "HiggsAnalysis/CombinedLimit/interface/RobustHesse.h"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <typeinfo>
#include <stdexcept>

#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooWorkspace.h"
#include "TDecompBK.h"
#include "TMatrixDSymEigen.h"


RobustHesse::RobustHesse(RooAbsReal &nll) : nll_(&nll) {
  targetNllForStencils_ = 0.1;
  minNllForStencils_ = 0.095;
  maxNllForStencils_ = 0.105;
  doRescale_ = true;
  initialize();
}

void RobustHesse::initialize() {
  nll0_ = nll_->getVal();
  std::unique_ptr<RooArgSet> allpars(nll_->getParameters(RooArgSet()));
  RooFIter iter = allpars->fwdIterator();
  RooAbsArg *item;
  while ((item = iter.next())) {
    RooRealVar *rrv = dynamic_cast<RooRealVar*>(item);
    if (rrv && !rrv->isConstant()) {
      v_.push_back(rrv);
      nominal_.push_back(rrv->getVal());
    }
  }
  std::cout << ">> Found " << v_.size() << " floating parameters\n";
  for (unsigned i = 0; i < v_.size(); ++i) {
    v_[i]->Print();
  }

  stencils_.resize(v_.size());
  d1coeffs_.resize(v_.size());
  d2coeffs_.resize(v_.size());
  rescales_.resize(v_.size());
}

double RobustHesse::deltaNLL() {
  return nll_->getVal() - nll0_;
}

double RobustHesse::deltaNLL(std::vector<unsigned> const& indices, std::vector<double> const& vals) {
  if (indices.size() != vals.size()) {
    std::cout << ">> Error in deltaNLL: number of indices != number of vals\n";
    return 0.;
  }
  if (indices.size() == 1) {
    auto it = nllcache_.find(std::make_pair(indices[0], vals[0]));
    if (it == nllcache_.end()) {
      for (unsigned i = 0; i < indices.size(); ++i) {
        v_[indices[i]]->setVal(vals[i]);
      }
      it = nllcache_.emplace(std::make_pair(indices[0], vals[0]), deltaNLL()).first;
      for (unsigned i = 0; i < indices.size(); ++i) {
        v_[indices[i]]->setVal(nominal_[indices[i]]);
      }
    }
    return it->second;
  }
  for (unsigned i = 0; i < indices.size(); ++i) {
    v_[indices[i]]->setVal(vals[i]);
  }
  double result = deltaNLL();

  for (unsigned i = 0; i < indices.size(); ++i) {
    v_[indices[i]]->setVal(nominal_[indices[i]]);
  }
  return result;
}


int RobustHesse::setParameterStencil(unsigned i) {
  double x = nominal_[i];
  RooRealVar *rrv = v_[i];

  double valLo = x - rrv->getError();
  double valHi = x + rrv->getError();

  // Am I near a boundary?
  double boundLo = rrv->getMin();
  double boundHi = rrv->getMax();

  // Do I have a boundary?
  // bool hasBoundLo = rrv->hasMin();
  // bool hasBoundHi = rrv->hasMax();

  bool closeToLo = valLo < boundLo;
  bool closeToHi = valHi > boundHi;

  if (closeToLo) {
    std::cout << ">> Parameter " << rrv->GetName() << " is too close to the lower bound: \n";
    rrv->Print();
    valLo = boundLo + 1E-7;
  }
  if (closeToHi) {
    std::cout << ">> Parameter " << rrv->GetName() << " is too close to the upper bound: \n";
    rrv->Print();
    valHi = boundHi + 1E-7;
  }

  double dNllLo = deltaNLL({i}, {valLo});
  double dNllHi = deltaNLL({i}, {valHi});

  // What sort of problems can we have here?
  //  a) dNll can be negative => we weren't quite at the minimum,
  //     but maybe we did just take a very small step
  //  b) dNLL is too large or too small. Should perform a search to get in the target range
  //  c) dNLL could be zero, or zero-ish. Might be a non-active parameter in the multipdf
  // Solution to cases (a) and (b) is the same: perform a search in the allowed range to see if
  // we can find a better val. If we don't on one side only, then we can use an asymmetric stencil
  // if we fail on both sides we have a problem.
  // Note that a very large dNll could just mean the parameter has a very strong anti-correlation
  // with another.

  std::cout << rrv->GetName() << "\t" << dNllLo << "\t" << dNllHi << "\n";

  bool notViableLo = false;
  bool notViableHi = false;

  if (dNllLo < minNllForStencils_) {
    if (closeToLo) {
      notViableLo = true;
    } else {
      auto bound_result = findBound(i, x, valLo - x, 2., 2., targetNllForStencils_, boundLo, 100);
      valLo = bound_result.second;
      dNllLo = deltaNLL({i}, {valLo});
      if (dNllLo < minNllForStencils_) {
        std::cout << ">> dNllLo is still too small after [findBound]\n";
        notViableLo = true;
      } else {
        std::cout << ">> dNllLo is now " << dNllLo << "\n";
      }
    }
    if (notViableLo) std::cout << ">> dNllLo is not viable\n";
  }

  if (dNllHi < minNllForStencils_) {
    if (closeToHi) {
      notViableHi = true;
    } else {
      auto bound_result = findBound(i, x, valHi - x, 2., 2., targetNllForStencils_, boundHi, 100);
      valHi = bound_result.second;
      dNllHi = deltaNLL({i}, {valHi});
      if (dNllHi < minNllForStencils_) {
        std::cout << ">> dNllHi is still too small after [findBound]\n";
        notViableHi = true;
      } else {
        std::cout << ">> dNllHi is now " << dNllHi << "\n";
      }
    }
    if (notViableHi) std::cout << ">> dNllHi is not viable\n";
  }

  if (dNllLo > maxNllForStencils_) {
    valLo = improveWithBisect(i, x, valLo, targetNllForStencils_, minNllForStencils_,  maxNllForStencils_, 20);
    dNllLo = deltaNLL({i}, {valLo});
  }
  if (dNllHi > maxNllForStencils_) {
    valHi = improveWithBisect(i, x, valHi, targetNllForStencils_, minNllForStencils_,  maxNllForStencils_, 20);
    dNllHi = deltaNLL({i}, {valHi});
  }


  std::cout << "[setParameterStencil] " << rrv->GetName() << "\t" << x << "\t" << valLo << " (" << dNllLo << ")\t" << valHi << " (" << dNllHi << ")\n";


  rescales_[i] = 1.0;
  if (!notViableLo && !notViableHi) {
    if (doRescale_) {
      rescales_[i] = (valHi - valLo) / 2.;
    }
    stencils_[i] = {(valLo - x) / rescales_[i], 0., (valHi - x) / rescales_[i]};
  } else if (notViableLo && !notViableHi) {
    if (doRescale_) {
      rescales_[i] = (valHi - x);
    }
    stencils_[i] = {0.,  (valHi - x) / (2. * rescales_[i]) ,(valHi - x) / rescales_[i]};
  } else if (!notViableLo && notViableHi) {
    if (doRescale_) {
      rescales_[i] = (x - valLo);
    }
    stencils_[i] = {(valLo - x) / rescales_[i], (valLo - x) / (2. * rescales_[i]), 0.};
  } else {
    std::cout << "[setParameterStencil] " << rrv->GetName() << "not viable at all!\n";
  }


  return 0;
}

std::pair<int, double> RobustHesse::findBound(unsigned i, double x, double initialDelta, double initialMult, double scaleMult, double threshold, double hardBound, unsigned maxIters) {
  std::cout << ">> [findBound] for parameter " << i << ", best-fit = " << x << ", initialDelta = " << initialDelta << "\n";
  double mult = initialMult;
  double trial = x + initialDelta;
  for (unsigned j = 0; j < maxIters; ++j) {
    bool trial_hit_bound = false;
    trial = x + initialDelta * mult;
    if (hardBound > x && trial > hardBound) {
      std::cout << "Step would go past bound, reducing\n";
      trial = hardBound - 1E-7;
      trial_hit_bound = true;
    }
    if (hardBound < x && trial < hardBound) {
      std::cout << "Step would go past bound, reducing\n";
      trial = hardBound + 1E-7;
      trial_hit_bound = true;
    }

    double dNll = deltaNLL({i}, {trial});
    std::cout << ">> dNLL at " << trial << " (x" << mult << " of original delta) = " << dNll << "\n";
    if ((mult > 1. && dNll > threshold) || (mult <= 1. && dNll < threshold)) {
      std::cout << ">> This is past the threshold, returning\n";
      return std::make_pair(0, trial);
    } else if (trial_hit_bound) {
      return std::make_pair(1, trial);
    }
    mult *= scaleMult;
  }
  return std::make_pair(2, trial);
}


double RobustHesse::improveWithBisect(unsigned i, double x, double max, double target, double targetLo, double targetHi, unsigned maxIters) {
  double cmin = x;
  double cmax = max;
  double trial = max;
  for (unsigned j = 0; j < maxIters; ++j) {
    trial = (cmin + cmax) / 2.;
    double dNll = deltaNLL({i}, {trial});
    std::cout << "[Bisect] " << j << " cmin = " << cmin << ", cmax = " << cmax << ", trial = " << trial << " dNll = " << dNll << "\n";
    if (dNll > targetLo && dNll < targetHi) {
      break;
    }
    if (dNll < target) {
      cmin = trial;
    } else {
      cmax = trial;
    }
  }
  return trial;
}

std::vector<double> RobustHesse::getFdCoeffs(unsigned n, std::vector<double> const& stencil) {
  std::vector<double> result;
  unsigned N = stencil.size();
  TMatrixD smatrix(N, N);
  TVectorD dvec(N);
  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < N; ++j) {
      smatrix(i, j) = std::pow(stencil[j], double(i));
    }
    dvec[i] = double((i == n) ? factorial(n) : 0);
  }
  TMatrixD smatrix2 = smatrix;
  smatrix2.Invert();
  TVectorD res = (smatrix2 * dvec);
  for (unsigned i = 0; i < N; ++i) {
    result.push_back(res[i]);
  }
  return result;
}


int RobustHesse::hesse() {
  // ---- Step 1: figure out stencils for FD
  // -- Step 1a: determine reasonable uncertainties
  for (unsigned i = 0; i < v_.size(); ++i) {
    setParameterStencil(i);
  }

  std::vector<RooRealVar *> v_new;
  std::vector<double> nominal_new;
  std::vector<std::vector<double>> stencils_new;
  std::vector<std::vector<double>> d1coeffs_new;
  std::vector<std::vector<double>> d2coeffs_new;
  std::vector<double> rescales_new;

  for (unsigned i = 0; i < v_.size(); ++i) {
    if (stencils_[i].size() == 3) {
      v_new.push_back(v_[i]);
      nominal_new.push_back(nominal_[i]);
      stencils_new.push_back(stencils_[i]);
      rescales_new.push_back(rescales_[i]);
    }
  }
  d1coeffs_.resize(v_new.size());
  d2coeffs_.resize(v_new.size());

  v_ = v_new;
  nominal_ = nominal_new;
  stencils_ = stencils_new;
  rescales_ = rescales_new;

  for (unsigned i = 0; i < v_.size(); ++i) {
    if (stencils_[i].size() == 3) {
      d1coeffs_[i] = getFdCoeffs(1, stencils_[i]);
      d2coeffs_[i] = getFdCoeffs(2, stencils_[i]);
      printf("%-80s %-10.3f %-10.3f %-10.3f | %-10.3f %-10.3f %-10.3f | %-10.3f %-10.3f %-10.3f\n",
        v_[i]->GetName(), stencils_[i][0], stencils_[i][1], stencils_[i][2],
        d1coeffs_[i][0], d1coeffs_[i][1], d1coeffs_[i][2],
        d2coeffs_[i][0], d2coeffs_[i][1], d2coeffs_[i][2]
        );
    }
  }



  // Step 2: Calculate and populate hessian
  TMatrixDSym hessian(v_.size());
  unsigned ntotal = (((v_.size() * v_.size()) - v_.size()) / 2) + v_.size();
  unsigned idx = 0;

  if (loadFile_ != "") {
    TFile fin(loadFile_.c_str());
    hessian = *((TMatrixDSym*)gDirectory->Get("hessian"));
  } else {
    for (unsigned i = 0; i < v_.size(); ++i) {
      for (unsigned j = i; j < v_.size(); ++j) {
        if (idx % 100 == 0) {
          std::cout << " - Done " << idx << "/" << ntotal << " terms\n";
        }
        double term = 0.;
        if (i == j) {
          for (unsigned k = 0; k < stencils_[i].size(); ++k) {
            if (stencils_[i][k] != 0.) {
              term += deltaNLL({i}, {nominal_[i] + rescales_[i] * stencils_[i][k]}) * d2coeffs_[i][k];
            }
          }
        } else {
          for (unsigned k = 0; k < stencils_[i].size(); ++k) {
            double c1 = d1coeffs_[i][k];
            double c2 = 0.;
            for (unsigned l = 0; l < stencils_[j].size();++l) {
              if (stencils_[i][k] == 0. && stencils_[j][l] == 0.) {
                continue;
              } else if (stencils_[i][k] == 0.) {
                c2 += deltaNLL({j}, {nominal_[j] + stencils_[j][l] * rescales_[j]}) * d1coeffs_[j][l];
              } else if (stencils_[j][l] == 0.) {
                c2 += deltaNLL({i}, {nominal_[i] + stencils_[i][k] * rescales_[i]}) * d1coeffs_[j][l];
              } else {
                c2 += deltaNLL({i, j}, {nominal_[i] + stencils_[i][k] * rescales_[i], nominal_[j] + stencils_[j][l] * rescales_[j]}) * d1coeffs_[j][l];
              }
            }
            term += (c2 * c1);
          }
        }
        hessian[i][j] = term;
        hessian[j][i] = term;
        ++idx;
      }
    }
  }
  if (saveFile_ != "") {
    TFile fout(saveFile_.c_str(), "RECREATE");
    gDirectory->WriteObject(&hessian, "hessian");
  }

  std::vector<std::string> remove_vars = {"CMS_fake_ele_wh3l", "QCDscale_tt", "CMS_fake_mu_stat_wh3l", "CMS_scale_j_RelativeJERHF_13TeV", "CMS_scale_j_RelativeJEREC2_13TeV", "CMS_ttHl_fakes", "CMS_ttH_ddQCD", "CMS_vhbb_bTagWeightJES_pt3_eta3"};

  std::vector<RooRealVar *> v_new2;
  std::vector<double> nominal_new2;
  std::vector<std::vector<double>> stencils_new2;
  // std::vector<std::vector<double>> d1coeffs_new2;
  // std::vector<std::vector<double>> d2coeffs_new2;
  std::vector<double> rescales_new2;

  TMatrixDSym hessian2(v_.size() - remove_vars.size());

  unsigned inew = 0;
  unsigned jnew = 0;
  std::vector<bool> ikeep(v_.size());

  for (unsigned i = 0; i < v_.size(); ++i) {
    ikeep[i] = (std::find(remove_vars.begin(), remove_vars.end(), v_[i]->GetName()) == remove_vars.end());
  }

  for (unsigned i = 0; i < v_.size(); ++i) {
    if (ikeep[i]) {
      v_new2.push_back(v_[i]);
      nominal_new2.push_back(nominal_[i]);
      stencils_new2.push_back(stencils_[i]);
      rescales_new2.push_back(rescales_[i]);
    }
    jnew = inew;
    if (ikeep[i]) {
      for (unsigned j = i; j < v_.size(); ++j) {
        if (ikeep[j]) {
          hessian2[inew][jnew] = hessian[i][j];
          hessian2[jnew][inew] = hessian[j][i];
          ++jnew;
        }
      }
      ++inew;
    }
  }

  v_ = v_new2;
  nominal_ = nominal_new2;
  stencils_ = stencils_new2;
  rescales_ = rescales_new2;
  hessian.ResizeTo(hessian2);
  hessian = hessian2;

  // hessian.Print();

  bool print_only_negative = true;
  unsigned max = 10;
  for (unsigned ai = 0; ai < max; ++ai) {
    bool have_negative_eigenvalues = false;
    TMatrixDSymEigen eigen(hessian);
    auto eigenvals = eigen.GetEigenValues();
    auto eigenvecs = eigen.GetEigenVectors();

    // for (int ei =0; ei < hessian.GetNrows(); ++ei) {
    //   double diag = hessian[ei][ei];
    //   double off_diag = 0.;
    //   for (int ej = 0; ej < hessian.GetNrows(); ++ej) {
    //     if (ej != ei) {
    //       off_diag += hessian[ei][ej];
    //     }
    //   }
    //   std::cout << "Diag = " << diag << "\tOff-diag = " << off_diag << " (" << v_[ei]->GetName() << ")\n";
    // }


    for (int ei = 0; ei < eigenvals.GetNrows(); ++ei) {
      if (eigenvals[ei] < 0.) have_negative_eigenvalues = true;


      std::vector<std::pair<unsigned, double>> eigenvec;
      double mag2 = 0.;
      for (int ej = 0; ej < eigenvecs.GetNrows(); ++ej) {
        mag2 += (eigenvecs[ej][ei]) * (eigenvecs[ej][ei]);
        eigenvec.push_back(std::make_pair(unsigned(ej), eigenvecs[ej][ei]));
      }
      if ((!print_only_negative) || (print_only_negative && eigenvals[ei] < 0.)) {
        std::cout << "Eigenvalue " << ei << " : " << eigenvals[ei] << "\n";
        std::sort(eigenvec.begin(), eigenvec.end(), [](std::pair<unsigned, double> const& p1, std::pair<unsigned, double> const& p2) {
          return std::fabs(p1.second) > std::fabs(p2.second);
        });
        unsigned n_show = 50;
        double trial = 0.002;
        for (unsigned es = 0; es < TMath::Min(n_show, unsigned(eigenvec.size())); ++es) {
          // unsigned int hi = eigenvec[es].first;
          // hessian[hi][hi] = hessian[hi][hi] * (1. + trial * std::pow(std::fabs(eigenvec[es].second), 2.));
          std::cout << " - " << eigenvec[es].second << "\t" << v_[eigenvec[es].first]->GetName() << "\n";
          for (unsigned es2 = es; es2 < TMath::Min(n_show, unsigned(eigenvec.size())); ++es2) {
            if (es != es2) {
              unsigned int hi = eigenvec[es].first;
              unsigned int hj = eigenvec[es2].first;
              double current = hessian[hi][hj];
              hessian[hi][hj] = current * (1. - trial);
              hessian[hj][hi] = current * (1. - trial);
            }
          }
        }
      }
    }
    if (!have_negative_eigenvalues) {
      break;
    } else {
      std::cout << ">> Attempt " << ai << " failed, increasing hessian diagonal elements...\n";
      // // Try and scale up the diagonal??
      // for (int ei =0; ei < hessian.GetNrows(); ++ei) {
      //   hessian[ei][ei] = hessian[ei][ei] * 1.001;
      // }
    }
  }

  // eigenvecs.Print();

  // Try to fix...
  // TMatrixD eigenvals_matrix(eigenvals.GetNrows(), eigenvals.GetNrows());
  // for (int ei = 0; ei < eigenvals.GetNrows(); ++ei) {
  //   if (eigenvals[ei] < 1E-6) {
  //     eigenvals_matrix[ei][ei] = 1E-6;
  //   } else {
  //     eigenvals_matrix[ei][ei] = eigenvals[ei];
  //   }
  // }
  // auto invert_eigenvecs = eigenvecs.Invert();

  // TMatrixD new_hessian(eigenvecs * (eigenvals_matrix * invert_eigenvecs));
  // for (int ei = 0; ei < eigenvals.GetNrows(); ++ei) {
  //   for (int ej = 0; ej < eigenvals.GetNrows(); ++ej) {
  //     hessian[ei][ej] = new_hessian[ei][ej];
  //   }
  // }

  TMatrixDSym covariance = hessian;
  std::cout << " - Inverting Hessian...\n";
  covariance.Invert();

  // TDecompBK decomp(hessian);
  // std::cout << " - Inverting Hessian...\n";
  // TMatrixDSym covariance = decomp.Invert();
  // covariance.Invert();


  for (unsigned i = 0; i < v_.size(); ++i) {
    std::cout << v_[i]->GetName() << "\t" << (std::sqrt(covariance[i][i]) * rescales_[i]) << "\n";
    for (unsigned j = i+1; j < v_.size(); ++j) {
      double corr = covariance[i][j] / (std::sqrt(covariance[i][i]) * std::sqrt(covariance[j][j]));
      if (std::abs(corr) > 0.8) {
        std::cout << " - correlation with " << v_[j]->GetName() << " = " << corr << "\n";
      }
    }
  }

  if (saveFile_ != "") {
    TFile fout2("covariance_output.root", "RECREATE");
    TH2F h_cov("covariance", "covariance", v_.size(), 0, v_.size(), v_.size(), 0, v_.size());
    RooWorkspace w("w", "w");
    for (unsigned i = 0; i < v_.size(); ++i) {
      w.import(*(v_[i]));
      h_cov.GetXaxis()->SetBinLabel(i + 1, v_[i]->GetName());
      h_cov.GetYaxis()->SetBinLabel(i + 1, v_[i]->GetName());
      for (unsigned j = i; j < v_.size(); ++j) {
        h_cov.SetBinContent(i + 1, j + 1, covariance[i][j] * rescales_[i] * rescales_[j]);
        h_cov.SetBinContent(j + 1, i + 1, covariance[i][j] * rescales_[i] * rescales_[j]);
      }
    }
    gDirectory->WriteObject(&hessian, "hessian");
    gDirectory->WriteObject(&covariance, "covariance");
    gDirectory->WriteObject(&h_cov, "h_covariance");
    w.Write(); 
    fout2.Close();
  }



  // Step 3: Test inversion

  return 0;

}

void RobustHesse::SaveHessianToFile(std::string const& filename) {
  saveFile_ = filename;
}
void RobustHesse::LoadHessianFromFile(std::string const& filename) {
  loadFile_ = filename;
}


