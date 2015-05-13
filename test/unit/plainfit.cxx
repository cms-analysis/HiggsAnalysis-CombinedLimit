#include <cstdlib>
#include <unistd.h>
#include "TFile.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooStats/ModelConfig.h"
#include "Math/MinimizerOptions.h"

int main(int argc, char **argv) {
    if (argc <= 1) { printf("Usage: %s file -w workspace(=w) -c modelConfig(=ModelConfig) -D dataset(=data_obs)  -S snapshot  -s strategy(=0) -t tolerance(=1) -M param_to_run_minos_on  \n",argv[0]); return 1; }
    const char *workspace = "w"; // -w
    const char *dataset   = "data_obs"; // -D
    const char *modelConfig = "ModelConfig"; // -c
    const char *snapshot    = NULL; // -S
    int   strategy    = 0; // -s
    float tolerance   = 1; // -t
    const char *minos = NULL; // -M
    int   optimize    = 2; // -O
    bool  useFitTo    = false; // -f ( use pdf->fitTo instead of Minimizer )
    do {
        int opt = getopt(argc, argv, "w:D:c:S:s:t:M:O:f");
        switch (opt) {
            case 'w': workspace = optarg; break;
            case 'D': dataset = optarg; break;
            case 'c': modelConfig = optarg; break;
            case 'S': snapshot = optarg; break;
            case 's': strategy = atoi(optarg); break;
            case 't': tolerance = atof(optarg); break;
            case 'M': minos = optarg; break;
            case 'O': optimize = atoi(optarg); break;
            case 'f': useFitTo = true; break;
            case '?': std::cerr << "Unsupported option. Please see the code. " << std::endl; return 1; break;
        }
        if (opt == -1) break;
    } while (true);
    TFile *f = TFile::Open(argv[optind]);
    if (!f) { std::cerr << "ERROR: could not open " << argv[optind] << std::endl; return 2; }
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
    if (minos && !w->var(minos)) {
        std::cerr << "ERROR: POI '" << minos << "' not found. Available POIs are: " << std::endl;
        mc->GetParametersOfInterest()->Print("");
        return 2;
    }
    if (snapshot) w->loadSnapshot(snapshot);
    RooArgSet poi; if (minos) poi.add(*w->var(minos));
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(tolerance);
    TStopwatch timer;
    if (useFitTo) {
        const RooCmdArg & maybeMinos = (minos ? RooFit::Minos(poi) : RooCmdArg::none());
        mc->GetPdf()->fitTo(*d, 
                RooFit::Constrain(*mc->GetNuisanceParameters()),
                RooFit::Minimizer("Minuit2","minimize"),
                RooFit::Offset(true),
                RooFit::Optimize(optimize > 0),
                RooFit::Strategy(strategy),
                RooFit::Hesse(kFALSE),
                maybeMinos);
        if (minos) w->var(minos)->Print("");
    } else {
        RooAbsReal *nll = mc->GetPdf()->createNLL(*d, RooFit::Constrain(*mc->GetNuisanceParameters()));
        RooMinimizer minim(*nll);
        minim.setStrategy(strategy);
        minim.setEps(tolerance);
        minim.setOffsetting(1);
        minim.optimizeConst(optimize);
        minim.minimize("Minuit2","minimize");
        if (minos) { minim.minos(poi); w->var(minos)->Print(""); }
    }
    timer.Stop(); printf("Done in %.2f min (cpu), %.2f min (real)\n", timer.CpuTime()/60., timer.RealTime()/60.);
}
