#include "RooWorkspace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TEnv.h"
#include "RooAbsPdf.h"
#include <dlfcn.h>
#include <RooStats/ModelConfig.h>
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"

void (*dump_)(const char *);


int main(int argc, char *argv[]) {
	if (void *sym = dlsym(0, "igprof_dump_now")) {
	  dump_ = __extension__ (void(*)(const char *)) sym;
	} else {
	  dump_=0;
	  std::cout << "Heap profile requested but application is not"
	            << " currently being profiled with igprof" << std::endl;
	}

	gROOT->GetName();

	if (argc < 3) return 1;

	TFile f(argv[1]);
	if(dump_) {
		dump_("profdump_file.out.gz");
	}

	RooWorkspace *w = (RooWorkspace*)gDirectory->Get(argv[2]);
	// w->Print();
	auto allfuncs = w->allFunctions();
	auto it =allfuncs.fwdIterator();
	for (RooAbsArg *a = it.next(); a != 0; a = it.next()) {
		auto rrv = dynamic_cast<RooAbsReal*>(a);
		if (rrv) {
			rrv->getVal();
		}
	}
	if(dump_) {
		dump_("profdump_wsp.out.gz");
	}

        RooStats::ModelConfig* mc_s = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig")); 
        RooAbsPdf &pdf = *mc_s->GetPdf();
        RooAbsData *dobs = w->data("data_obs");
        const RooCmdArg &constrainCmdArg = RooFit::Constrain(*mc_s->GetNuisanceParameters());
        std::auto_ptr<RooAbsReal> nll;
        nll.reset(); // first delete the old one, to avoid using more memory, even if temporarily
        nll.reset(pdf.createNLL(*dobs, constrainCmdArg, RooFit::Extended(pdf.canBeExtended()), RooFit::Offset(true))); // make a new nll

	if(dump_) {
		dump_("profdump_nll.out.gz");
	}

	return 0;
}
