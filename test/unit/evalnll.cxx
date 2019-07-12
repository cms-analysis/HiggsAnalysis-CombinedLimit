#include <cstdlib>
#include <unistd.h>
#include "TFile.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooStats/ModelConfig.h"
#include "Math/MinimizerOptions.h"
#include "Math/IOptions.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProfilingTools.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"

void init_rtd() {
  // CMSDAS Defaults (you can turn off with --X-rtd <name>=0
  runtimedef::set("OPTIMIZE_BOUNDS", 1);
  runtimedef::set("ADDNLL_RECURSIVE", 1);
  runtimedef::set("ADDNLL_GAUSSNLL", 1);
  runtimedef::set("ADDNLL_HISTNLL", 1);
  runtimedef::set("ADDNLL_CBNLL", 1);
  runtimedef::set("TMCSO_AdaptivePseudoAsimov", 1);
  // Optimization for bare RooFit likelihoods (--optimizeSimPdf=0)
  runtimedef::set("MINIMIZER_optimizeConst", 2); 
  runtimedef::set("MINIMIZER_rooFitOffset", 1); 
  // Optimization for ATLAS HistFactory likelihoods
  runtimedef::set("ADDNLL_ROOREALSUM_FACTOR",1);
  runtimedef::set("ADDNLL_ROOREALSUM_NONORM",1);
  runtimedef::set("ADDNLL_ROOREALSUM_BASICINT",1);
  runtimedef::set("ADDNLL_ROOREALSUM_KEEPZEROS",1);
  runtimedef::set("ADDNLL_PRODNLL",1);
  runtimedef::set("ADDNLL_HFNLL",1);
  runtimedef::set("ADDNLL_HISTFUNCNLL",1);
  runtimedef::set("ADDNLL_ROOREALSUM_CHEAPPROD",1);
  // these are not default in combine yet
  runtimedef::set("MINIMIZER_analytic",1);
  runtimedef::set("FAST_VERTICAL_MORPH",1);
}

void set_rtd(const char *rtd) {
  std::string rtds(rtd);
  std::string::size_type idx = rtds.find('=');
  if (idx == std::string::npos) {
      runtimedef::set(rtd, 1);
      std::cout << "Turning on runtime-define " << rtd << std::endl;
  } else {
      std::string name  = rtds.substr(0, idx);
      std::string svalue = rtds.substr(idx+1);
      int ivalue = atoi( svalue.c_str() );
      std::cout << "Setting runtime-define " << name << " to " << ivalue << std::endl;
      runtimedef::set(name, ivalue);
  }
}

int main(int argc, char **argv) {
    if (argc <= 1) { printf("Usage: %s file -w workspace(=w) -c modelConfig(=ModelConfig) -D dataset(=data_obs) -m mass(=125)  -S snapshot  -s strategy(=0) -N evals(=10000)  \n",argv[0]); return 1; }
    const char *workspace = "w"; // -w
    const char *dataset   = "data_obs"; // -D
    const char *modelConfig = "ModelConfig"; // -c
    const char *snapshot    = NULL; // -S
    int   strategy = 0; // -s
    float mass     = 125; // -m
    int   maxevals = 10000; // -N
    init_rtd(); 
    do {
        int opt = getopt(argc, argv, "w:D:c:S:s:m:R:N:");
        switch (opt) {
            case 'w': workspace = optarg; break;
            case 'D': dataset = optarg; break;
            case 'c': modelConfig = optarg; break;
            case 'S': snapshot = optarg; break;
            case 's': strategy = atoi(optarg); break;
            case 'm': mass = atof(optarg); break;
            case 'R': set_rtd(optarg); break;
            case 'N': maxevals = atoi(optarg); break;
            case '?': std::cerr << "Unsupported option. Please see the code. " << std::endl; return 1; break;
        }
        if (opt == -1) break;
    } while (true);
    TStopwatch timer;
    TFile *f = TFile::Open(argv[optind]);
    if (!f) { std::cerr << "ERROR: could not open " << argv[optind] << std::endl; return 2; }

    timer.Stop();
    printf("File opened in %.2f min (cpu), %.2f min (real)\n", timer.CpuTime()/60., timer.RealTime()/60.); fflush(stdout);
    timer.Start();
    
    RooWorkspace *w = (RooWorkspace *) f->Get(workspace);
    if (!w)  { std::cerr << "ERROR: could not find workspace '" << workspace << "' in " << argv[optind] << std::endl; f->ls(); return 2; }
    RooStats::ModelConfig *mc = (RooStats::ModelConfig *) w->genobj(modelConfig);
    if (!mc) { std::cerr << "ERROR: could not find ModelConfig '" << modelConfig << "' in worskpace" << std::endl; return 2; }
    RooAbsData *d = w->data(dataset);
    if (!d) { 
        std::cerr << "ERROR: could not find dataset '" << dataset << "' in workspace. Available datasets are: " << std::endl;
        std::list<RooAbsData*> datasets = w->allData();
        for (std::list<RooAbsData*>::const_iterator it = datasets.begin(), ed = datasets.end(); it != ed; ++it) {
            (*it)->Print("");
        }
        return 2;
    }
    if (snapshot) w->loadSnapshot(snapshot);
    if (mass) {
        if (w->var("MH")) w->var("MH")->setVal(mass);
    }
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
    RooAbsPdf *pdf = mc->GetPdf();
    const RooArgSet *nuisances = mc->GetNuisanceParameters();

    timer.Stop();
    printf("Workspace opened in %.2f min (cpu), %.2f min (real)\n", timer.CpuTime()/60., timer.RealTime()/60.);  fflush(stdout);
    timer.Start();

    RooAbsReal *nll = pdf->createNLL(*d, RooFit::Constrain(*nuisances), RooFit::Extended(pdf->canBeExtended()), RooFit::Offset(true));
    cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(nll);
    if (!simnll) {
        std::cerr << "ERROR: not a cacheutils::CachingSimNLL !" << std::endl;
        return 1;
    }
    if (runtimedef::get(std::string("MINIMIZER_analytic"))) {
        simnll->setAnalyticBarlowBeeston(true);
    }
    if (runtimedef::get("FAST_VERTICAL_MORPH")) {
        CMSHistFunc::EnableFastVertical();
    }

    if (simnll) {
        simnll->setZeroPoint();
    }
    timer.Stop(); 
    printf("NLL created in %.2f min (cpu), %.2f min (real)\n", timer.CpuTime()/60., timer.RealTime()/60.);  fflush(stdout);
    std::vector<RooRealVar *> var; std::vector<std::pair<float,float>> vrange;
    RooLinkedListIter iter = nuisances->iterator();
    for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next()) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv == 0 || rrv->isConstant()) continue;
        var.push_back(rrv);
        vrange.emplace_back(rrv->getMin(), rrv->getMax());
    }
    timer.Start();
    int ieval = 0;
    do {
        float nll0 = nll->getVal(); ++ieval;
        switch (strategy) {
            case -1: break;
            case 0:
            case 1:
                for (unsigned int i = 0, nvar = var.size(); i < nvar; ++i) {
                    RooRealVar *rrv = var[i];
                    float wval = rrv->getVal(), wnext = (wval == vrange[i].second ? vrange[i].first : vrange[i].second);
                    rrv->setVal(0.999 * wval + 0.001 * wnext);
                    float nll1 = nll->getVal(); ++ieval;
                    rrv->setVal(wval);
                    if (ieval >= maxevals) break;
                }
                break;
        }
        for (unsigned int i = 0, nvar = var.size(); i < nvar; ++i) {
            RooRealVar *rrv = var[i];
            float wval = rrv->getVal(), wnext = (wval == vrange[i].second ? vrange[i].first : vrange[i].second);
            rrv->setVal(0.95 * wval + 0.05 * wnext);
        }
    } while ( ieval < maxevals );
    printf("Done %d evals in %.2f min (cpu), %.2f min (real)\nAverage  eval time %.2f ms, eval time %.2f ms\n", ieval, timer.CpuTime()/60., timer.RealTime()/60., timer.CpuTime()*1000./ieval, timer.RealTime()*1000./ieval);
}
