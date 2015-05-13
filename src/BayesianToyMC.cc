#include <stdexcept>
#include <cmath>
#include "HiggsAnalysis/CombinedLimit/interface/BayesianToyMC.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include <Math/DistFuncMathCore.h>

#include "HiggsAnalysis/CombinedLimit/interface/Combine.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include "HiggsAnalysis/CombinedLimit/interface/utils.h"

using namespace RooStats;

int BayesianToyMC::numIters_ = 1000;
std::string BayesianToyMC::integrationType_ = "toymc";
unsigned int BayesianToyMC::tries_ = 1;
float BayesianToyMC::hintSafetyFactor_ = 5.;
std::vector<std::string> BayesianToyMC::twoPoints_;

BayesianToyMC::BayesianToyMC() :
    LimitAlgo("BayesianToyMC specific options")
{
    options_.add_options()
        ("integrationType", boost::program_options::value<std::string>(&integrationType_)->default_value(integrationType_), "Integration algorithm to use")
        ("tries",      boost::program_options::value<unsigned int>(&tries_)->default_value(tries_), "Number of times to run the ToyMC on the same data")
        ("numIters,i", boost::program_options::value<int>(&numIters_)->default_value(numIters_),    "Number of iterations or calls used within iteration (0=ROOT Default)")
        ("hintSafetyFactor",
                boost::program_options::value<float>(&hintSafetyFactor_)->default_value(hintSafetyFactor_),
                "set range of integration equal to this number of times the hinted limit")
        ("twoPoints",
                boost::program_options::value<std::vector<std::string> >(&twoPoints_)->multitoken(), "Compute BF comparing two points in parameter space");
        ;
}

void BayesianToyMC::applyOptions(const boost::program_options::variables_map &vm) {
    if (!withSystematics) {
        std::cout << "BayesianToyMC: when running with no systematics, BayesianToyMC is identical to BayesianSimple." << std::endl;
        tries_    = 1;
        numIters_ = 0;
        integrationType_ = "";
    }
    if (!twoPoints_.empty() && twoPoints_.size() != 2) throw std::logic_error("twoPoints option requires exactly two points\n");
    if (!twoPoints_.empty() && !doSignificance_) throw std::logic_error("twoPoints option works with --significance\n"); 
}
bool BayesianToyMC::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
  if (doSignificance_) return runBayesFactor(w,mc_s,mc_b,data,limit,limitErr,hint);

  RooArgSet  poi(*mc_s->GetParametersOfInterest());

  RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
  if ((hint != 0) && (*hint > r->getMin())) {
    r->setMax(hintSafetyFactor_*(*hint));
  }
  double rMax = r->getMax();
  
  std::auto_ptr<RooStats::ModelConfig> mc_noNuis(0);
  if (!withSystematics && mc_s->GetNuisanceParameters() != 0) {
    mc_noNuis.reset(new RooStats::ModelConfig(w));
    mc_noNuis->SetPdf(*mc_s->GetPdf());
    mc_noNuis->SetObservables(*mc_s->GetObservables());
    mc_noNuis->SetParametersOfInterest(*mc_s->GetParametersOfInterest());
    mc_noNuis->SetPriorPdf(*mc_s->GetPriorPdf());
    mc_s = mc_noNuis.get();
  }
  std::auto_ptr<RooAbsPdf> nuisancePdf;
  for (;;) {
    limit = 0; limitErr = 0; bool rerun = false;
    for (unsigned int i = 0; i < tries_; ++i) {
        BayesianCalculator bcalc(data, *mc_s);
        bcalc.SetLeftSideTailFraction(0);
        bcalc.SetConfidenceLevel(cl); 
        if (!integrationType_.empty()) bcalc.SetIntegrationType(integrationType_.c_str());
        if (numIters_) bcalc.SetNumIters(numIters_);
        if (integrationType_ == "toymc") {
            if (nuisancePdf.get() == 0)  nuisancePdf.reset(utils::makeNuisancePdf(*mc_s));
            bcalc.ForceNuisancePdf(*nuisancePdf);
            // turn off optimization of constraint terms in the NLL, otherwise
            // it does not divide properly by the nuisance pdf
            cacheutils::CachingSimNLL::forceUnoptimizedConstraints();
        }
        // get the interval
        std::auto_ptr<SimpleInterval> bcInterval(bcalc.GetInterval());
        if (bcInterval.get() == 0) return false;
        double lim = bcInterval->UpperLimit();
        // check against bound
        if (lim >= 0.5*r->getMax()) { 
            std::cout << "Limit " << r->GetName() << " < " << lim << "; " << r->GetName() << " max < " << r->getMax() << std::endl;
            if (r->getMax()/rMax > 20) return false;
            r->setMax(r->getMax()*2); 
            rerun = true; break;
        }
        // add to running sum(x) and sum(x2)
        limit    += lim;
        limitErr += lim*lim;
        if (tries_ > 1 && verbose > 1) std::cout << " - limit from try " << i << ": " << lim << std::endl;
    }
    if (rerun) continue;
    limit /= tries_; 
    limitErr = (tries_ > 1 ? std::sqrt((limitErr/tries_ - limit*limit)/((tries_-1)*tries_)) : 0);
    if (verbose > -1) {
        std::cout << "\n -- BayesianToyMC -- " << "\n";
        if (limitErr > 0) {
            std::cout << "Limit: " << r->GetName() << " < " << limit << " +/- " << limitErr << " @ " << cl * 100 << "% CL" << std::endl;
        } else {
            std::cout << "Limit: " << r->GetName() << " < " << limit << " @ " << cl * 100 << "% CL" << std::endl;
        }
    }
    break;
  }
  return true;
}

bool BayesianToyMC::runBayesFactor(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    std::pair<double,double> ppS, ppB;
    double offset = std::numeric_limits<double>::quiet_NaN();
    if (twoPoints_.empty()) {
        ppS = priorPredictiveDistribution(mc_s,data,0,&offset);
        ppB = priorPredictiveDistribution(mc_b,data,0,&offset);
    } else {
        RooArgSet POI(*mc_s->GetParametersOfInterest());
        RooArgSet points[2]; 
        for (int i=0; i<2; ++i) {
            utils::createSnapshotFromString(twoPoints_[i], POI, points[i], "--twoPoints");
            if (verbose > 1) { std::cout << "Point " << i+1 <<" : " <<std::endl; points[i].Print("V"); }
        }
        ppS = priorPredictiveDistribution(mc_s,data,&points[0],&offset);
        ppB = priorPredictiveDistribution(mc_s,data,&points[1],&offset);
    }
    limit = ppS.first/ppB.first;
    limitErr = std::hypot(ppS.second/ppB.first, ppS.first*ppB.second/(ppB.first*ppB.first));
    if (verbose > -1) {
        std::cout << "\n -- BayesianToyMC -- " << "\n";
        std::cout << "Bayes factor: " << limit << " +/- " << limitErr << std::endl;
        std::cout << "    (min(1/BF,BF) expressed as sigma's: " << ROOT::Math::normal_quantile_c(std::min(limit,1.0/limit),1.0) << ")" << std::endl;
    }
    return true;
}

std::pair<double,double> BayesianToyMC::priorPredictiveDistribution(RooStats::ModelConfig *mc, RooAbsData &data, const RooArgSet *point, double *offset) {
  // factorize away nuisance pdf
  RooAbsPdf *pdf = mc->GetPdf();
  std::auto_ptr<RooAbsPdf>  nuisancePdf, nonNuisancePdf; 
  if (withSystematics) {
    RooArgList constraints;
    nonNuisancePdf.reset(utils::factorizePdf(*data.get(), *pdf, constraints));
    if (constraints.getSize() > 0) {
        nuisancePdf.reset(new RooProdPdf("nuis","",constraints));
        pdf = nonNuisancePdf.get();
    } else nonNuisancePdf.release();
  }
  std::cout << "Factorized PDF, now creating NLL" << std::endl;
  // create NLL
  const RooCmdArg &constrain  = withSystematics  ? RooFit::Constrain(*mc->GetNuisanceParameters()) : RooCmdArg::none();
  std::auto_ptr<RooAbsReal> nll(pdf->createNLL(data, constrain, RooFit::Extended(pdf->canBeExtended())));
  std::auto_ptr<RooArgSet>  params(nll->getParameters(data));

  // Determine which POIs we have to generate, if any
  RooArgSet poiToGen; 
  if (mc->GetParametersOfInterest()) poiToGen.add(*mc->GetParametersOfInterest());
  if (point != 0) poiToGen.remove(*point, false, true);
  if (verbose && poiToGen.getSize()) { std::cout << "POI to generate:"; poiToGen.Print(); }

  // Determine what other floating parameters besides POI and nuisances we have to generate
  RooArgSet otherParams(*params);
  RooStats::RemoveConstantParameters(&otherParams);
  if (mc->GetParametersOfInterest()) otherParams.remove(*mc->GetParametersOfInterest()); 
  if (withSystematics && mc->GetNuisanceParameters()) otherParams.remove(*mc->GetNuisanceParameters());
  if (verbose && otherParams.getSize()) { std::cout << "Other unnamed parameters to generate:"; otherParams.Print(); }

  // Set the point we're running at
  if (point != 0) params->assignValueOnly(*point);

  // start running
  std::vector<double> results; double sum = 0;
  for (unsigned int t = 0; t < tries_; ++t) {
      std::auto_ptr<RooDataSet> nuisanceValues, poiValues;
      if (withSystematics) nuisanceValues.reset(nuisancePdf->generate(*mc->GetNuisanceParameters(), numIters_));
      if (poiToGen.getSize() > 0) {
        if (mc->GetPriorPdf() == 0) throw std::logic_error(std::string("Missing prior in model: ")+ mc->GetName());
        poiValues.reset(mc->GetPriorPdf()->generate(poiToGen, numIters_));
      }
      for (int i = 0; i < numIters_; ++i) {
        if (nuisanceValues.get() != 0) *params = *nuisanceValues->get(i);
        if (poiValues.get() != 0) *params = *poiValues->get(i);
        if (otherParams.getSize()) RooStats::RandomizeCollection(otherParams);
        if (verbose > 2) { std::cout << "\n\n==== POINT "<< t << ","<<i<<" ====" << std::endl; params->Print("V"); }
        double nllVal = nll->getVal();
        if (offset) { 
            if (isnan(*offset)) *offset = nllVal; 
            nllVal -= *offset; 
        }
        if (verbose > 1) std::cout << "nll[" << t << ","<<i<<"] = " << nllVal << ", p = " << std::exp(-nllVal) << std::endl;
        results.push_back(std::exp(-nllVal));
        sum += results.back();
      }
  }
  double n = results.size();
  sum /= n; 
  double sumd = 0, sumd2 = 0;
  for (int i = 0, ni = results.size(); i < ni; ++i) {
      sumd  += results[i] - sum;
      sumd2 += std::pow(results[i] - sum,2);
  }
  sum += sumd/numIters_; 
  double err = std::sqrt((sumd2/numIters_ - std::pow(sumd/numIters_,2))/numIters_);
  return std::make_pair(sum,err);
}

