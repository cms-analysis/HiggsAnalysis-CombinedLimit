#include "HiggsAnalysis/CombinedLimit/interface/GoodnessOfFit.h"
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooRandom.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooSimultaneous.h>
#include <RooAddPdf.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>
#include <RooCategory.h>
#include <RooStats/ModelConfig.h>
#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"
#include "HiggsAnalysis/CombinedLimit/interface/Significance.h"
#include "HiggsAnalysis/CombinedLimit/interface/CascadeMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/CloseCoutSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"

#include <numeric>
#include <memory>


#include <Math/MinimizerOptions.h>

using namespace RooStats;

std::string GoodnessOfFit::algo_;
//std::string GoodnessOfFit::minimizerAlgo_ = "Minuit2";
//float       GoodnessOfFit::minimizerTolerance_ = 1e-4;
//int         GoodnessOfFit::minimizerStrategy_  = 1;
float       GoodnessOfFit::mu_ = 0.0;
bool        GoodnessOfFit::fixedMu_ = false;
bool        GoodnessOfFit::makePlots_ = false;
TDirectory* GoodnessOfFit::plotDir_ = nullptr;
std::vector<std::string>  GoodnessOfFit::binNames_;
std::vector<float>        GoodnessOfFit::qVals_;
std::string GoodnessOfFit::setParametersForFit_ = "";
std::string GoodnessOfFit::setParametersForEval_ = "";

GoodnessOfFit::GoodnessOfFit() :
    LimitAlgo("GoodnessOfFit specific options")
{
    options_.add_options()
        ("algorithm",          boost::program_options::value<std::string>(&algo_), "Goodness of fit algorithm. Supported algorithms are 'saturated', 'KS' and 'AD'.")
        ("setParametersForFit",   boost::program_options::value<std::string>(&setParametersForFit_)->default_value(""), "Set parameters values for the saturated model fitting step")
        ("setParametersForEval",   boost::program_options::value<std::string>(&setParametersForEval_)->default_value(""), "Set parameter values for the saturated model NLL eval step")
  //      ("minimizerAlgo",      boost::program_options::value<std::string>(&minimizerAlgo_)->default_value(minimizerAlgo_), "Choice of minimizer (Minuit vs Minuit2)")
  //      ("minimizerTolerance", boost::program_options::value<float>(&minimizerTolerance_)->default_value(minimizerTolerance_),  "Tolerance for minimizer")
  //      ("minimizerStrategy",  boost::program_options::value<int>(&minimizerStrategy_)->default_value(minimizerStrategy_),      "Stragegy for minimizer")
        ("fixedSignalStrength", boost::program_options::value<float>(&mu_)->default_value(mu_),  "Compute the goodness of fit for a fixed signal strength. If not specified, it's left floating")
        ("plots",  "Make plots containing information of the computation of the Anderson-Darling or Kolmogorov-Smirnov test statistic")
    ;
}

void GoodnessOfFit::applyOptions(const boost::program_options::variables_map &vm) 
{
    fixedMu_ = !vm["fixedSignalStrength"].defaulted();
    if (algo_ == "saturated") {
      std::cout << "Will use saturated models to compute goodness of fit for a binned likelihood" << std::endl;
    } else if (algo_ == "AD") {
      std::cout << "Will use the Anderson-Darling test to compute goodness of fit for a binned likelihood" << std::endl;
    } else if (algo_ == "KS") {
      std::cout << "Will use the Kolmogorov-Smirnov test to compute goodness of fit for a binned likelihood" << std::endl;
    } else {
      throw std::invalid_argument("GoodnessOfFit: algorithm "+algo_+" not supported");
    }
    makePlots_ = vm.count("plots");
}

bool GoodnessOfFit::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) { 
  //double minimizerTolerance_  = ROOT::Math::MinimizerOptions::DefaultTolerance();
  //std::string minimizerAlgo_       = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
  //Significance::MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);

  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
  if (fixedMu_) { r->setVal(mu_); r->setConstant(true); }
  if (algo_ == "saturated") return runSaturatedModel(w, mc_s, mc_b, data, limit, limitErr, hint);
  static bool is_init = false;
  if (algo_ == "AD" || algo_ == "KS") {
    if (!is_init) {
      initKSandAD(mc_s);
      if (makePlots_) {
        plotDir_ = outputFile ? outputFile->mkdir("GoodnessOfFit") : 0;
      }
      is_init = true;
    }
    if (algo_ == "AD") return runKSandAD(w, mc_s, mc_b, data, limit, limitErr, hint, 0);
    if (algo_ == "KS") return runKSandAD(w, mc_s, mc_b, data, limit, limitErr, hint, 1);
  }

  return false;  
}

void GoodnessOfFit::initKSandAD(RooStats::ModelConfig *mc_s) {
  RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(mc_s->GetPdf());
  if (sim) {
    std::auto_ptr<RooAbsCategoryLValue> cat(
        static_cast<RooAbsCategoryLValue *>(sim->indexCat().Clone()));
    int nbins = cat->numBins((const char *)0);
    binNames_.resize(nbins);
    qVals_.resize(nbins);
    for (int i = 0; i < nbins; ++i) {
      cat->setBin(i);
      binNames_[i] = cat->getLabel();
      qVals_[i] = 0.;
      Combine::addBranch(binNames_[i].c_str(), &qVals_[i], (binNames_[i]+"/F").c_str());
    }
  }
}

bool GoodnessOfFit::runSaturatedModel(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) { 
  RooAbsPdf *pdf_nominal = mc_s->GetPdf();
  // now I need to make the saturated pdf
  std::auto_ptr<RooAbsPdf> saturated;
  // factorize away constraints anyway
  RooArgList constraints;
  RooAbsPdf *obsOnlyPdf = utils::factorizePdf(*mc_s->GetObservables(), *pdf_nominal, constraints);
  // case 1:
  RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(obsOnlyPdf);
  if (sim) {
      RooAbsCategoryLValue *cat = (RooAbsCategoryLValue *) sim->indexCat().Clone();
      std::auto_ptr<TList> datasets(data.split(*cat, true));
      int nbins = cat->numBins((const char *)0);
      RooArgSet newPdfs;
      TString satname = TString::Format("%s_saturated", sim->GetName());
      RooSimultaneous *satsim = (typeid(*sim) == typeid(RooSimultaneousOpt)) ? new RooSimultaneousOpt(satname, "", *cat) : new RooSimultaneous(satname, "", *cat); 
      for (int ic = 0, nc = nbins; ic < nc; ++ic) {
          cat->setBin(ic);
          RooAbsPdf *pdfi = sim->getPdf(cat->getLabel());
          if (pdfi == 0) continue;
          RooAbsData *datai = (RooAbsData *) datasets->FindObject(cat->getLabel());
          if (datai == 0) throw std::runtime_error(std::string("Error: missing dataset for category label ")+cat->getLabel());
          std::unique_ptr<RooArgSet> data_observables(pdfi->getObservables(datai));
          std::unique_ptr<RooAbsData> data_reduced(datai->reduce(*data_observables));
          RooAbsPdf *saturatedPdfi = makeSaturatedPdf(*data_reduced);
          //delete datai;
          if (constraints.getSize() > 0) {
            RooArgList terms(constraints); terms.add(*saturatedPdfi);
            RooProdPdf *prodpdf = new RooProdPdf(TString::Format("%s_constr", saturatedPdfi->GetName()), "", terms);
            prodpdf->addOwnedComponents(RooArgSet(*saturatedPdfi));
            saturatedPdfi = prodpdf;
          }
          satsim->addPdf(*saturatedPdfi, cat->getLabel());
          satsim->addOwnedComponents(RooArgSet(*saturatedPdfi));
      }
      // Transfer the channel masks manually
      RooSimultaneousOpt* satsimopt = dynamic_cast<RooSimultaneousOpt*>(satsim);
      RooSimultaneousOpt* simopt = dynamic_cast<RooSimultaneousOpt*>(pdf_nominal);
      if (satsimopt && simopt) {
        satsimopt->addChannelMasks(simopt->channelMasks());
      }
      saturated.reset(satsim);
  } else {
      RooAbsPdf *saturatedPdfi = makeSaturatedPdf(data);
      if (constraints.getSize() > 0) {
          RooArgList terms(constraints); terms.add(*saturatedPdfi);
          RooProdPdf *prodpdf = new RooProdPdf(TString::Format("%s_constr", saturatedPdfi->GetName()), "", terms);
          prodpdf->addOwnedComponents(RooArgSet(*saturatedPdfi));
          saturatedPdfi = prodpdf;
      }
      saturated.reset(saturatedPdfi);
  }

  CloseCoutSentry sentry(verbose < 2);

  const RooCmdArg &constrainCmdArg = withSystematics  ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooCmdArg();
  std::auto_ptr<RooAbsReal> nominal_nll(pdf_nominal->createNLL(data, constrainCmdArg));
  std::auto_ptr<RooAbsReal> saturated_nll(saturated->createNLL(data, constrainCmdArg));

  if (setParametersForFit_ != "") {
    utils::setModelParameters(setParametersForFit_, w->allVars());
  }
  CascadeMinimizer minimn(*nominal_nll, CascadeMinimizer::Unconstrained);
 // minimn.setStrategy(minimizerStrategy_);
  minimn.minimize(verbose-2);
  // This test is a special case where we are comparing the likelihoods of two
  // different models and so we can't re-zero the NLL with respect to the
  // initial parameters.
  if (dynamic_cast<cacheutils::CachingSimNLL*>(nominal_nll.get())) {
    static_cast<cacheutils::CachingSimNLL*>(nominal_nll.get())->clearConstantZeroPoint();
  }

  if (setParametersForEval_ != "") {
    utils::setModelParameters(setParametersForEval_, w->allVars());
  }
  double nll_nominal = nominal_nll->getVal();

  if (setParametersForFit_ != "") {
    utils::setModelParameters(setParametersForFit_, w->allVars());
  }
  CascadeMinimizer minims(*saturated_nll, CascadeMinimizer::Unconstrained);
  //minims.setStrategy(minimizerStrategy_);
  minims.minimize(verbose-2);
  if (dynamic_cast<cacheutils::CachingSimNLL*>(saturated_nll.get())) {
    static_cast<cacheutils::CachingSimNLL*>(saturated_nll.get())->clearConstantZeroPoint();
  }

  if (setParametersForEval_ != "") {
    utils::setModelParameters(setParametersForEval_, w->allVars());
  }
  double nll_saturated = saturated_nll->getVal();

  sentry.clear();

  saturated.reset();
  for (int i = 0, n = tempData_.size(); i < n; ++i) delete tempData_[i]; 
  tempData_.clear();

  if (fabs(nll_nominal) > 1e10 || fabs(nll_saturated) > 1e10) return false;
  limit = 2*(nll_nominal-nll_saturated);

  std::cout << "\n --- GoodnessOfFit --- " << std::endl;
  std::cout << "Best fit test statistic: " << limit << std::endl;
  return true;
}

// Code for the Anderson-Darling test originates from https://gist.github.com/neggert/4586791
bool GoodnessOfFit::runKSandAD(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint, bool kolmo) { 
  RooAbsPdf *pdf = mc_s->GetPdf();

  // Don't want the constraints here
  RooArgList constraints;
  RooAbsPdf *obsOnlyPdf = utils::factorizePdf(*mc_s->GetObservables(), *pdf, constraints);

  //First, find the best fit values
  CloseCoutSentry sentry(verbose < 2);
  
  /*
  const RooCmdArg &minim = RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                             ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  int minimizerStrategy_  = ROOT::Math::MinimizerOptions::DefaultStrategy();

  std::auto_ptr<RooFitResult> result(pdf->fitTo(data, RooFit::Save(1), minim, RooFit::Strategy(minimizerStrategy_), RooFit::Hesse(0), RooFit::Constrain(*mc_s->GetNuisanceParameters())));
  sentry.clear();
  */

  const RooCmdArg &constrainCmdArg = withSystematics  ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooCmdArg();
  std::auto_ptr<RooAbsReal> nll(pdf->createNLL(data, constrainCmdArg));
  CascadeMinimizer minim(*nll, CascadeMinimizer::Unconstrained);
  //minims.setStrategy(minimizerStrategy_);
  minim.minimize(verbose-2);

  sentry.clear();

  RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(obsOnlyPdf);

  limit = 0;

  // now check if we have categories to deal with
  if (binNames_.size()) {
    RooAbsCategoryLValue const &cat = sim->indexCat();    
    // Split the data using the option "true" to include empty datasets
    std::auto_ptr<TList> datasets(data.split(cat, true));
    
    // Number of categories should always equal the number of datasets
    if (datasets->GetSize() != int(binNames_.size())) {
      throw std::runtime_error(
          "Numbers of categories and datasets do not match");
    }

    for (unsigned i = 0; i < binNames_.size(); i++) {
      RooAbsData *cat_data = dynamic_cast<RooAbsData *>(datasets->At(i));
      RooAbsPdf *cat_pdf = sim->getPdf(binNames_[i].c_str());
      std::auto_ptr<RooArgSet> observables(cat_pdf->getObservables(cat_data));
      if (observables->getSize() > 1) {
        std::cout << "Warning, KS and AD statistics are not well defined for "
                     "models with more than one observable\n";
      }
      RooRealVar *x = dynamic_cast<RooRealVar *>(observables->first());
      qVals_[i] = EvaluateADDistance(*cat_pdf, *cat_data, *x, kolmo);
    }
    limit = std::accumulate(qVals_.begin(), qVals_.end(), 0.);
  } else {
    RooRealVar *obs = dynamic_cast<RooRealVar *>(
        std::auto_ptr<RooArgSet>(pdf->getObservables(data))->first());
    limit = EvaluateADDistance(*pdf, data, *obs, kolmo);
  }


  std::cout << "\n --- GoodnessOfFit --- " << std::endl;
  if(kolmo){
      std::cout << "Kolmogorov-Smirnov test statistic: " << limit << std::endl;
  }else{
      std::cout << "Anderson-Darling test statistic: " << limit << std::endl;
  }

  makePlots_ = false;

  return true;
}

Double_t GoodnessOfFit::EvaluateADDistance(RooAbsPdf& pdf, RooAbsData& data, RooRealVar& observable, bool kolmo) {
    typedef std::pair<double, double> double_pair;
    std::vector<double_pair> data_points;
    Int_t n_data = data.numEntries();
    Double_t s_data = data.sumEntries();

    const RooArgSet* datavals;
    RooRealVar* observable_val;
    for (int i = 0; i < n_data; i++) {
        datavals = data.get(i);
        observable_val = (RooRealVar*)(datavals->find(observable.GetName()));
        data_points.push_back(std::make_pair(observable_val->getVal(), data.weight()));
    }

    std::stable_sort(data_points.begin(), data_points.end(),
         [](double_pair i, double_pair j) {
           return i.first < j.first;
         });

    double test_stat = 0.;
    double current_cdf_val = 0.;
    double last_cdf_val = 0.;
    double empirical_df = 0.;
    double observableval = 0.;
    double bin_prob = 0.;
    double distance = 0.;

    // CDF of the PDF
    // If RooFit needs to use the scanning technique then increase the number
    // of sampled bins from 1000 to 10000
    std::auto_ptr<RooAbsReal> cdf(pdf.createCdf(observable, RooFit::ScanAllCdf(), RooFit::ScanParameters(10000, 2)));
    TH1 * hCdf = nullptr;
    TH1 * hEdf = nullptr;
    TH1 * hDiff = nullptr;
    if (plotDir_ && makePlots_) {
      hCdf = observable.createHistogram("cdf");
      hEdf = observable.createHistogram("edf");
      hDiff = observable.createHistogram("diff");
      hCdf->SetName((std::string(data.GetName())+"_cdf").c_str());
      hEdf->SetName((std::string(data.GetName())+"_edf").c_str());
      hDiff->SetName((std::string(data.GetName())+"_diff").c_str());
    }

    int bin = 0;
    for (std::vector<double_pair>::const_iterator d = data_points.begin();
         d != data_points.end(); d++, ++bin) {
        // observableval = ((d+1)->first + d->first)/2.; // d->first is middle of bin, want upper edge.
      
        // This is a better way to get the upper bin edge in the case where we
        // have variable bin widths (I hope)
        observableval = observable.getBinning().binHigh(
            observable.getBinning().binNumber(d->first));
        observable.setVal(observableval);
        // observable.bin
        current_cdf_val = cdf->getVal();
        empirical_df += d->second/s_data;

        if (plotDir_ && makePlots_) {
          hCdf->SetBinContent(bin+1, current_cdf_val);
          hEdf->SetBinContent(bin+1, empirical_df);
        }

        if (kolmo){
            distance = std::abs(empirical_df-current_cdf_val);
            if (verbose >= 3) {
              std::cout << "Observable: " << observableval << "\tdata: " << d->second << "\tedf: " << empirical_df << "\tcdf: " << current_cdf_val << "\tdistance: " << distance << "\n";
            }
            if (distance > test_stat) test_stat = distance;
        }else{
            bin_prob = current_cdf_val-last_cdf_val;
            distance = s_data*pow((empirical_df-current_cdf_val), 2)/current_cdf_val/(1.-current_cdf_val)*bin_prob;
            if (current_cdf_val >= 1.0) {
              distance = 0.;
            }
            if (verbose >= 3) {
              std::cout << "Observable: " << observableval << "\tdata: " << d->second << "\tedf: " << empirical_df << "\tcdf: " << current_cdf_val << "\tdistance: " << distance << "\n";
            }
            // from L. Demortier, CDF/ANAL/JET/CDFR/3419
            test_stat += distance;
        }
        if (plotDir_ && makePlots_) {
          hDiff->SetBinContent(bin+1, distance);
        }

        last_cdf_val = current_cdf_val;
    }

    if (plotDir_ && makePlots_) {
      plotDir_->WriteTObject(hCdf);
      plotDir_->WriteTObject(hEdf);
      plotDir_->WriteTObject(hDiff);
      delete hCdf;
      delete hEdf;
      delete hDiff;
    }

    // if (kolmo){
    //     if (test_stat < oldlimit) test_stat = oldlimit; // The test statistic of the Kolmogorov-Smirnov test is the maximum of the test statistics of all individual PDFs.
    // }else{
    //     test_stat += oldlimit; // The test statistic of the Anderson-Darling test is the sum of the test statistics of all individual PDFs.
    // }
    return test_stat;
}

RooAbsPdf * GoodnessOfFit::makeSaturatedPdf(RooAbsData &data) {
  if (verbose > 1) std::cout << "Generating saturated model for " << data.GetName() << std::endl;
  RooDataHist *rdh = new RooDataHist(TString::Format("%s_binned", data.GetName()), "", *data.get(), data); tempData_.push_back(rdh);
  if (verbose > 1) utils::printRDH(rdh);
  RooHistPdf *hpdf = new RooHistPdf(TString::Format("%s_shape", data.GetName()), "", *rdh->get(), *rdh);
  RooConstVar *norm = new RooConstVar(TString::Format("%s_norm", data.GetName()), "", data.sumEntries());
  // we use RooAddPdf because this works with CachingNLL
  RooAddPdf *ret = new RooAddPdf(TString::Format("%s_saturated", data.GetName()), "", RooArgList(*hpdf), RooArgList(*norm));
  ret->addOwnedComponents(RooArgSet(*norm));
  ret->addOwnedComponents(RooArgSet(*hpdf));
  return ret;
}


