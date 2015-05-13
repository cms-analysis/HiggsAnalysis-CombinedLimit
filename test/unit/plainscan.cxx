#include <cstdlib>
#include <unistd.h>
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooStats/ModelConfig.h"
#include "Math/MinimizerOptions.h"

int main(int argc, char **argv) {
    if (argc <= 1) { printf("Usage: %s file -w workspace(=w) -c modelConfig(=ModelConfig) -D dataset(=data_obs)  -S snapshot  -s strategy(=0) -t tolerance(=1) -o outfile(=scan.root) -P poi(=r) -m min -M max -n points(=10)  \n",argv[0]); return 1; }
    const char *workspace = "w"; // -w
    const char *dataset   = "data_obs"; // -D
    const char *modelConfig = "ModelConfig"; // -c
    const char *snapshot    = NULL; // -S
    int   strategy    = 0; // -s
    float tolerance   = 1; // -t
    float minval      = NAN;
    float maxval      = NAN;
    const char *poi = "r"; // -P
    const char *out = "out.root"; // -P
    int   optimize    = 2; // -O
    int   points      = 10;
    bool  allvars = false;
    do {
        int opt = getopt(argc, argv, "w:D:c:S:s:t:M:O:m:n:o:P:A");
        switch (opt) {
            case 'w': workspace = optarg; break;
            case 'D': dataset = optarg; break;
            case 'c': modelConfig = optarg; break;
            case 'S': snapshot = optarg; break;
            case 's': strategy = atoi(optarg); break;
            case 't': tolerance = atof(optarg); break;
            case 'n': points = atoi(optarg); break;
            case 'm': minval = atof(optarg); break;
            case 'M': maxval = atof(optarg); break;
            case 'P': poi = optarg; break;
            case 'o': out = optarg; break;
            case 'O': optimize = atoi(optarg); break;
            case 'A': allvars = true; break;
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
    if (!w->var(poi)) {
        std::cerr << "ERROR: POI '" << poi << "' not found. Available POIs are: " << std::endl;
        mc->GetParametersOfInterest()->Print("");
        return 2;
    }
    if (snapshot) w->loadSnapshot(snapshot);
    RooRealVar *r = w->var(poi);
    if (!std::isnan(minval)) r->setMin(minval);
    if (!std::isnan(maxval)) r->setMax(maxval);

    TFile *fOut = TFile::Open(out,"RECREATE");
    TTree *limit = new TTree("limit","limit");
    float xvar, deltaNll; 
    limit->Branch(poi, &xvar, (std::string(poi)+"/F").c_str());
    limit->Branch("deltaNLL", &deltaNll, "deltaNLL/F");
    std::vector<std::pair<float,RooRealVar*> > vars;
    if (allvars) {
        RooLinkedListIter iter = w->allVars().iterator();
        for (RooAbsArg *arg = (RooAbsArg*) iter.Next(); arg; arg = (RooAbsArg*) iter.Next()) {
            RooRealVar *rrv = dynamic_cast<RooRealVar *>(arg);
            if (rrv != 0 && !rrv->isConstant() && rrv != r) {
                vars.push_back(std::make_pair(r->getVal(),r));
            }
        }
        for (std::vector<std::pair<float,RooRealVar*> >::iterator it = vars.begin(), ed = vars.end(); it != ed; ++it) {
            limit->Branch(it->second->GetName(), & it->first, (std::string(it->second->GetName())+"/F").c_str());
        }
    }
    
    TStopwatch timer;
    RooAbsReal *nll = mc->GetPdf()->createNLL(*d, RooFit::Constrain(*mc->GetNuisanceParameters()));
    RooMinimizer minim(*nll);
    minim.setStrategy(strategy);
    minim.setEps(tolerance);
    minim.setOffsetting(1);
    minim.optimizeConst(optimize);
    minim.minimize("Minuit2","minimize");
    xvar = r->getVal(); deltaNll = 0; 
    for (std::vector<std::pair<float,RooRealVar*> >::iterator it = vars.begin(), ed = vars.end(); it != ed; ++it) {
        it->first = it->second->getVal();
    }
    limit->Fill();
    double nll0 = nll->getVal();
    for (int i = 0; i < points; ++i) {
        xvar = r->getMin() + ((i+0.5)/points)*(r->getMax()-r->getMin());
        r->setVal(xvar);
        minim.minimize("Minuit2","minimize");
        std::cout << "Point " << i << "/" << points << ": " << poi << " = " << xvar << std::endl;
        deltaNll = nll->getVal() - nll0; 
        for (std::vector<std::pair<float,RooRealVar*> >::iterator it = vars.begin(), ed = vars.end(); it != ed; ++it) {
            it->first = it->second->getVal();
        }
        limit->Fill();
    }
    timer.Stop(); printf("Done in %.2f min (cpu), %.2f min (real)\n", timer.CpuTime()/60., timer.RealTime()/60.);
    fOut->Close();

}
