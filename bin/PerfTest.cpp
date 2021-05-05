#include "RooWorkspace.h"
#include "TFile.h"
#include "TROOT.h"
#include "TEnv.h"
#include "RooAbsPdf.h"
#include <dlfcn.h>

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
			std::cout << rrv->GetName() << " " << rrv->getVal() << "\n";
		}
	}
	if(dump_) {
		dump_("profdump_wsp.out.gz");
	}

	// Load the workspace a second time, there seems to be some overhead
	// in the first load that won't appear in the second...
	// RooWorkspace *w2 = (RooWorkspace*)gDirectory->Get(argv[2]);


	// if(dump_) {
	// 	dump_("profdump_wsp2.out.gz");
	// }

	return 0;
}