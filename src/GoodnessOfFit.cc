#include "../interface/GoodnessOfFit.h"
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
#include <TLegend.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>
#include <RooStats/ModelConfig.h>
#include "RooStats/RooStatsUtils.h"
#include "RooCategory.h"
#include "../interface/Combine.h"
#include "../interface/ProfileLikelihood.h"
#include "../interface/CloseCoutSentry.h"
#include "../interface/RooSimultaneousOpt.h"
#include "../interface/utils.h"
#include <fstream>


#include <Math/MinimizerOptions.h>

using namespace RooStats;

std::string GoodnessOfFit::algo_;
std::string GoodnessOfFit::minimizerAlgo_ = "Minuit2";
float       GoodnessOfFit::minimizerTolerance_ = 1e-4;
int         GoodnessOfFit::minimizerStrategy_  = 1;
float       GoodnessOfFit::mu_ = 0.0;
bool        GoodnessOfFit::fixedMu_ = false;
bool        GoodnessOfFit::nofit_ = false;
bool        GoodnessOfFit::makeplot_ = false;

GoodnessOfFit::GoodnessOfFit() :
    LimitAlgo("GoodnessOfFit specific options")
{
    options_.add_options()
        ("algorithm",          boost::program_options::value<std::string>(&algo_), "Goodness of fit algorithm. Supported algorithms are 'saturated', 'kolmogorovsmirnov' and 'andersondarling'.")
        ("minimizerAlgo",      boost::program_options::value<std::string>(&minimizerAlgo_)->default_value(minimizerAlgo_), "Choice of minimizer (Minuit vs Minuit2)")
        ("minimizerTolerance", boost::program_options::value<float>(&minimizerTolerance_)->default_value(minimizerTolerance_),  "Tolerance for minimizer")
        ("minimizerStrategy",  boost::program_options::value<int>(&minimizerStrategy_)->default_value(minimizerStrategy_),      "Stragegy for minimizer")
        ("fixedSignalStrength", boost::program_options::value<float>(&mu_)->default_value(mu_),  "Compute the goodness of fit for a fixed signal strength. If not specified, it's left floating")
        ("nofit", boost::program_options::value<bool>(&nofit_)->default_value(nofit_),  "Doesn't take the NLL values directly from the fits (for saturated model). Default: False")
        ("makeplot", boost::program_options::value<bool>(&makeplot_)->default_value(makeplot_),  "Make a plot containing information of the computation of the Anderson-Darling or Kolmogorov-Smirnov test statistic")
    ;
}

void GoodnessOfFit::applyOptions(const boost::program_options::variables_map &vm) 
{
    fixedMu_ = !vm["fixedSignalStrength"].defaulted();
    if (algo_ == "saturated") std::cout << "Will use saturated models to compute goodness of fit for a binned likelihood" << std::endl;
    else if (algo_ == "andersondarling") std::cout << "Will use the Anderson-Darling test to compute goodness of fit for a binned likelihood" << std::endl;
else if (algo_ == "kolmogorovsmirnov") std::cout << "Will use the Kolmogorov-Smirnov test to compute goodness of fit for a binned likelihood" << std::endl;
    else throw std::invalid_argument("GoodnessOfFit: algorithm "+algo_+" not supported");
}

bool GoodnessOfFit::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) { 
  ProfileLikelihood::MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);

  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
  if (fixedMu_) { r->setVal(mu_); r->setConstant(true); }
  if (algo_ == "saturated") return runSaturatedModel(w, mc_s, mc_b, data, limit, limitErr, hint);
  if (algo_ == "andersondarling") return runAndeDarlKolmSmir(w, mc_s, mc_b, data, limit, limitErr, hint, 0);
  if (algo_ == "kolmogorovsmirnov") return runAndeDarlKolmSmir(w, mc_s, mc_b, data, limit, limitErr, hint, 1);
  return false;  
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
          RooAbsPdf *saturatedPdfi = makeSaturatedPdf(*datai);
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

  double nll_nominal, nll_saturated;

  CloseCoutSentry sentry(verbose < 2);
  // let's assume fits converge, for a while
  const RooCmdArg &minim = RooFit::Minimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                             ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  std::auto_ptr<RooFitResult> result_nominal(pdf_nominal->fitTo(data, RooFit::Save(1), minim, RooFit::Strategy(minimizerStrategy_), RooFit::Hesse(0), RooFit::Constrain(*mc_s->GetNuisanceParameters())));
  std::auto_ptr<RooFitResult> result_saturated(saturated->fitTo(data, RooFit::Save(1), minim, RooFit::Strategy(minimizerStrategy_), RooFit::Hesse(0), RooFit::Constrain(*mc_s->GetNuisanceParameters())));
  sentry.clear();

  if(nofit_){
    std::auto_ptr<RooAbsReal> nominal_nll(pdf_nominal->createNLL(data, RooFit::Constrain(*mc_s->GetNuisanceParameters())));
    std::auto_ptr<RooAbsReal> saturated_nll(saturated->createNLL(data, RooFit::Constrain(*mc_s->GetNuisanceParameters())));
    //if (nominal_nll.get()   == 0 or saturated_nll.get() == 0) return false;
    nll_nominal = nominal_nll->getVal();
    nll_saturated = saturated_nll->getVal();
  }else{
    if (result_nominal.get()   == 0 or result_saturated.get() == 0) return false;
    nll_nominal   = result_nominal->minNll();
    nll_saturated = result_saturated->minNll();
  }

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
bool GoodnessOfFit::runAndeDarlKolmSmir(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint, bool kolmo) { 
   /*
  Find the Anderson-Darling difference between fPdf and the data. This is just
  the maximum over data points of the distance between the CDF and the emprical
  distribution function of the data.
  */

  const RooArgSet * POI = mc_s->GetParametersOfInterest();
  RooRealVar* poi = dynamic_cast<RooRealVar*>(POI->first());
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Set POI constant so the fit doesn't adjust them
  // RooRealVar* param;
  // TIterator *param_iter = params.createIterator();
  // while (1) {
  //     param = dynamic_cast<RooRealVar*>(param_iter->Next());
  //     if (!param) break;
  //     param->setConstant(kTRUE);
  // }

  //First, find the best fit values
  RooAbsPdf *fPdf = mc_s->GetPdf();
  RooArgSet* pdfParams = fPdf->getParameters(data);
  RooRealVar* pdfPOI= dynamic_cast<RooRealVar*>(pdfParams->find(poi->GetName()));
  pdfPOI->setVal(poi->getVal());
  pdfPOI->setConstant();
  RemoveConstantParameters(pdfParams);
  fPdf->fitTo(data, RooFit::Constrain(*pdfParams), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE));

  // Unset POI constant so the fit doesn't adjust them
  // param_iter = params.createIterator();
  // while (1) {
  //     param = dynamic_cast<RooRealVar*>(param_iter->Next());
  //     if (!param) break;
  //     param->setConstant(kFALSE);
  // }

  limit = 0;
  RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(fPdf);

  // now check if we have categories to deal with
  if (sim) {
      const RooCategory* cat = dynamic_cast<const RooCategory*>(&(sim->indexCat()));
      // split the data
      TList* datasets = data.split(*cat);
      Int_t num_cats = datasets->GetEntries();
      for(Int_t i= 0; i < num_cats; i++) {
          // set to current category
          const Text_t* cat_name = cat->lookupType(i)->GetName();
          RooAbsData *cat_data = dynamic_cast<RooAbsData*>(datasets->At(i));
          RooAbsPdf *cat_pdf = sim->getPdf(cat_name);
          RooArgSet* observables = cat_pdf->getObservables(cat_data);
          RooRealVar* x = dynamic_cast<RooRealVar*>(observables->first());
          limit = EvaluateADDistance(*cat_pdf, *cat_data, *x, limit, kolmo);
      }
  } else {
      RooRealVar* obs = dynamic_cast<RooRealVar*>(fPdf->getObservables(data)->first());
      limit = EvaluateADDistance(*fPdf, data, *obs, limit, kolmo);
  }

  std::cout << "\n --- GoodnessOfFit --- " << std::endl;
  if(kolmo){
      std::cout << "Kolmogorov-Smirlov test statistic: " << limit << std::endl;
  }else{
      std::cout << "Anderson-Darling test statistic: " << limit << std::endl;
  }

  RooMsgService::instance().setGlobalKillBelow(msglevel);
  return true;
}


typedef std::pair<Double_t, Double_t> double_pair;

bool sort_pairs (double_pair i, double_pair j) {return (i.first < j.first);}

Double_t GoodnessOfFit::EvaluateADDistance( RooAbsPdf& pdf, RooAbsData& data, RooRealVar& observable, double& oldlimit, bool kolmo) {

    std::vector<double_pair > data_points;
    Int_t n_data = data.numEntries();
    Double_t s_data = data.sumEntries();

    const RooArgSet* datavals;
    RooRealVar* observable_val;

    for (int i = 0; i < n_data; i++) {
        datavals = data.get(i);
        observable_val = (RooRealVar*)(datavals->find(observable.GetName()));
        data_points.push_back(double_pair(observable_val->getVal(), data.weight()));
    }

    sort(data_points.begin(), data_points.end(), sort_pairs);

    Double_t test_stat = 0.;
    Double_t current_cdf_val = 0.;
    Double_t last_cdf_val = 0.;
    Double_t empirical_df = 0.;
    Double_t observableval, bin_prob;
    Double_t distance = 0.;

    //CDF of the PDF
    RooAbsReal* cdf;

    // Declarations for the saveplot option
    Int_t count=-1;
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.14);
    TCanvas* c = new TCanvas("c","c",600,600);
    TH1F* hcdf = new TH1F("cdf","",n_data-1,0,n_data-1);
    TH1F* hedf = new TH1F("edf","",n_data-1,0,n_data-1);
    TH1F* hcdif = new TH1F("cdif","",n_data-1,0,n_data-1);
    TH1F* hpdif = new TH1F("pdif","",n_data-1,0,n_data-1);
    TH1F* hADadd = new TH1F("ADadd","",n_data-1,0,n_data-1);
    TH1F* hADtrend = new TH1F("ADtrend","",n_data-1,0,n_data-1);
    TLegend* leg1 = new TLegend(0.14,0.91,0.9,1.0);
    TLegend* leg2 = new TLegend(0.14,0.91,0.9,1.0);
    TPad* pad1 = new TPad("pad1","pad1",0,0.22,1,0.98);
    TPad* pad2 = new TPad("pad2","pad2",0,0.15,1,0.3);
    TPad* pad3 = new TPad("pad3","pad3",0,0,1,0.15);

    if (makeplot_){
        pad1->SetFillStyle(4000);
        pad1->Draw();
        pad2->SetFillStyle(4000);
        pad2->Draw();
        pad3->SetFillStyle(4000);
        pad3->Draw();
        pad1->cd();

        hcdf->SetMinimum(0);
        hedf->SetMinimum(0);
        leg1->AddEntry(hcdf,"Cumulative distribution function","l");
        leg1->AddEntry(hedf,"Empirical distribution function","l");
        leg1->SetTextSize(0.04);

        hcdif->SetOption("P0");
        hcdif->SetMarkerStyle(20);
        hcdif->SetMarkerSize(1);
        hcdif->SetYTitle("CDF-EDF");
        hcdif->GetYaxis()->SetTitleSize(0.20);
        hcdif->GetYaxis()->SetTitleOffset(0.25);
        hcdif->GetYaxis()->SetNdivisions(306);
        hcdif->GetYaxis()->SetLabelSize(0.15);
        hcdif->GetXaxis()->SetLabelSize(0.15);

        hpdif->SetOption("P0");
        hpdif->SetMarkerStyle(20);
        hpdif->SetMarkerSize(1);
        hpdif->SetYTitle("Data/MC");
        hpdif->GetYaxis()->SetTitleSize(0.20);
        hpdif->GetYaxis()->SetTitleOffset(0.25);
        hpdif->GetYaxis()->SetNdivisions(306);
        hpdif->GetYaxis()->SetLabelSize(0.15);
        hpdif->GetXaxis()->SetLabelSize(0.15);

        hADadd->SetMinimum(0);
        hADtrend->SetMinimum(0);
        leg2->AddEntry(hADadd,"Gets added to the test statistic","l");
        leg2->AddEntry(hADtrend,"Trend of the test statistic","l");
        leg2->SetTextSize(0.04);
    }

    for(std::vector<double_pair>::const_iterator d = data_points.begin(); d != data_points.end()-1; d++) {
        observableval = ((d+1)->first + d->first)/2.; // d->first is middle of bin, want upper edge.
        observable.setVal(observableval);
        cdf = pdf.createCdf(observable, RooFit::ScanAllCdf());
        current_cdf_val = cdf->getVal();
        empirical_df += d->second/s_data;

        count++;
        hcdif->Fill(count,current_cdf_val-empirical_df);
        hpdif->Fill(count,(d->second/s_data)/(current_cdf_val-last_cdf_val));

        if (kolmo){
            distance = std::abs(empirical_df-current_cdf_val);
            hcdf->Fill(count,current_cdf_val);
            hedf->Fill(count,empirical_df);
            if (distance > test_stat) test_stat = distance;
        }else{
            bin_prob = current_cdf_val-last_cdf_val;
            // from L. Demortier, CDF/ANAL/JET/CDFR/3419
            test_stat += s_data*pow((empirical_df-current_cdf_val), 2)/current_cdf_val/(1.-current_cdf_val)*bin_prob;
            hADadd->Fill(count,s_data*pow((empirical_df-current_cdf_val), 2)/current_cdf_val/(1.-current_cdf_val)*bin_prob);
            hADtrend->Fill(count,test_stat);
        }
        last_cdf_val = current_cdf_val;

        if (last_cdf_val >= 1.0)
            break;
    }

    if(makeplot_){
        if (hcdif->GetMaximum() > std::abs(hcdif->GetMinimum())){
            hcdif->SetMinimum(-1.0*(hcdif->GetMaximum()));//hcdif->GetYaxis()->SetRange(-1.05*(hcdif->GetMaximum()),1.05*(hcdif->GetMaximum()));//
        }else{
            hcdif->SetMaximum(-1.0*(hcdif->GetMinimum()));//hcdif->GetYaxis()->SetRange(1.05*(hcdif->GetMinimum()),-1.05*(hcdif->GetMinimum()));//
        }
        if ((hpdif->GetMaximum())-1 > 1-(hpdif->GetMinimum())){
            hpdif->SetMinimum(2-(hpdif->GetMaximum()));//hpdif->GetYaxis()->SetRange(-1.05*(hpdif->GetMaximum())+2,1.05*(hpdif->GetMaximum()));//
        }else{
            hpdif->SetMaximum(2-(hpdif->GetMinimum()));//hpdif->GetYaxis()->SetRange(1.05*(hpdif->GetMinimum()),-1.05*(hpdif->GetMinimum())+2);//
        }

        if (kolmo){
            hcdf->Draw();
            hedf->SetLineColor(kRed);
            hedf->Draw("sames");
            leg1->Draw();
        }else{
            hADtrend->Draw();
            hADadd->SetLineColor(kRed);
            hADadd->Draw("sames");
            leg2->Draw();
        }
        pad2->cd();
        gPad->SetGridy();
        hcdif->Draw();
        pad3->cd();
        gPad->SetGridy();
        hpdif->Draw();

        Int_t filecount=0;
        std::string filename;
        std::fstream ifile;
        do{
            if (ifile != 0) ifile.close();
            filecount++;
            filename = "Gof_plot_"+std::to_string(filecount)+".png";
            ifile.open(filename);
        }while(ifile != 0);
        ifile.close();
        TString tfilename=filename;
        c->SaveAs(tfilename);
    }

    delete c;
    delete hcdf;
    delete hedf;
    delete hcdif;
    delete hpdif;
    delete hADadd;
    delete hADtrend;

    if (kolmo){
        if (test_stat < oldlimit) test_stat = oldlimit; // The test statistic of the Kolmogorov-Smirnov test is the maximum of the test statistics of all individual PDFs.
    }else{
        test_stat += oldlimit; // The test statistic of the Anderson-Darling test is the sum of the test statistics of all individual PDFs.
    }
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


