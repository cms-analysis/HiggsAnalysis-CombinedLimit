#include "RooWorkspace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TEnv.h"
#include "RooAbsPdf.h"
#include "RooMsgService.h"
#include <dlfcn.h>
#include <RooStats/ModelConfig.h>
#include "../interface/CachingNLL.h"
#include "../interface/ProfilingTools.h"

void (*dump_)(const char *);


int main(int argc, char *argv[]) {
	bool doCombineRuntimes = true;
	if (doCombineRuntimes) {
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
	}
	runtimedef::set("ADDNLL_VERBOSE_CACHING", 0);
	runtimedef::set("fullCloneFunc_VERBOSE", 1);
	runtimedef::set("CACHINGPDF_NOCLONE", 0); // ! THIS ONE IS NOT SAFE YET

	// Uncomment below for more info
	// RooMsgService::instance().addStream(RooFit::INFO);

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
	auto allfuncs = w->allFunctions();
	for (auto *a: allfuncs){
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
        std::unique_ptr<RooAbsReal> nll;
        nll.reset(pdf.createNLL(*dobs, constrainCmdArg, RooFit::Extended(pdf.canBeExtended()), RooFit::Offset(true))); // make a new nll
        nll->getVal();
	if(dump_) {
		dump_("profdump_nll.out.gz");
	}

	return 0;
}
