#include "../interface/AsimovUtils.h"

#include <memory>
#include <stdexcept>
#include <TIterator.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooProdPdf.h>
#include <RooUniform.h>
#include "../interface/utils.h"
#include "../interface/ToyMCSamplerOpt.h"
#include "../interface/CloseCoutSentry.h"
#include "../interface/CascadeMinimizer.h"
#include "../interface/Logger.h"
#include "../interface/ProfilingTools.h"

RooAbsData *asimovutils::asimovDatasetNominal(RooStats::ModelConfig *mc, double poiValue, int verbose) {
        RooArgSet  poi(*mc->GetParametersOfInterest());
        RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
        r->setConstant(true); r->setVal(poiValue);
        toymcoptutils::SimPdfGenInfo newToyMC(*mc->GetPdf(), *mc->GetObservables(), false); 

	if (verbose>2) {
	    Logger::instance().log(std::string(Form("AsimovUtils.cc: %d -- Parameters after fit for asimov dataset",__LINE__)),Logger::kLogLevelInfo,__func__);
    	    std::unique_ptr<TIterator> iter(mc->GetPdf()->getParameters((const RooArgSet*) 0)->createIterator());
    	    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
	  	TString varstring = utils::printRooArgAsString(a);
	  	Logger::instance().log(std::string(Form("AsimovUtils.cc: %d -- %s",__LINE__,varstring.Data())),Logger::kLogLevelInfo,__func__);
	    }
	}

        RooRealVar *weightVar = 0;
        RooAbsData *asimov = newToyMC.generateAsimov(weightVar,verbose); 
        delete weightVar;
        return asimov;
}

RooAbsData *asimovutils::asimovDatasetWithFit(RooStats::ModelConfig *mc, RooAbsData &realdata, RooAbsCollection &snapshot, bool needsFit, double poiValue, int verbose) {
        RooArgSet  poi(*mc->GetParametersOfInterest());
        RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
        r->setConstant(true); r->setVal(poiValue);
        bool mvg_gen_obs = runtimedef::get("ASIMOV_MULTIVAR_GAUSS_GEN");
        {
            CloseCoutSentry sentry(verbose < 3);
            if (mc->GetNuisanceParameters()) {
                needsFit &= true;
            } else {
                // Do we have free parameters anyway that need fitting?
                std::unique_ptr<RooArgSet> params(mc->GetPdf()->getParameters(realdata));
                std::unique_ptr<TIterator> iter(params->createIterator());
                for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
                    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
                    if ( rrv != 0 && rrv->isConstant() == false ) { needsFit &= true; break; }
                } 
            }
            if (needsFit) {
                //mc->GetPdf()->fitTo(realdata, RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(1), RooFit::Constrain(*mc->GetNuisanceParameters()));
                const RooCmdArg &constrain = (mc->GetNuisanceParameters() ? RooFit::Constrain(*mc->GetNuisanceParameters()) : RooCmdArg());
                std::unique_ptr<RooAbsReal> nll(mc->GetPdf()->createNLL(realdata, constrain, RooFit::Extended(mc->GetPdf()->canBeExtended())));
                CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
                minim.setStrategy(1);
                minim.minimize(verbose-1);
            }
        }
        if (mc->GetNuisanceParameters() && verbose > 1) {
            std::cout << "Nuisance parameters after fit for asimov dataset: " << std::endl;
            mc->GetNuisanceParameters()->Print("V");
        }

	if (verbose>2) { 
	    Logger::instance().log(std::string(Form("AsimovUtils.cc: %d -- Parameters after fit for asimov dataset",__LINE__)),Logger::kLogLevelInfo,__func__);
    	    std::unique_ptr<TIterator> iter(mc->GetPdf()->getParameters(realdata)->createIterator());
    	    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
	  	TString varstring = utils::printRooArgAsString(a);
	  	Logger::instance().log(std::string(Form("AsimovUtils.cc: %d -- %s",__LINE__,varstring.Data())),Logger::kLogLevelInfo,__func__);
	    }
	}

        toymcoptutils::SimPdfGenInfo newToyMC(*mc->GetPdf(), *mc->GetObservables(), false); 
        RooRealVar *weightVar = 0;
        RooAbsData *asimov = newToyMC.generateAsimov(weightVar,verbose); 
        delete weightVar;

        // NOW SNAPSHOT THE GLOBAL OBSERVABLES
        if (mc->GetGlobalObservables() && mc->GetGlobalObservables()->getSize() > 0) {
            RooArgSet gobs(*mc->GetGlobalObservables());

            // snapshot data global observables
            RooArgSet snapGlobalObsData;
            utils::setAllConstant(gobs, true);
            gobs.snapshot(snapGlobalObsData);

            RooArgSet nuis(*mc->GetNuisanceParameters());
            std::unique_ptr<RooAbsPdf> nuispdf(utils::makeNuisancePdf(*mc));
            RooProdPdf *prod = dynamic_cast<RooProdPdf *>(nuispdf.get());
            if (prod == 0) throw std::runtime_error("AsimovUtils: the nuisance pdf is not a RooProdPdf!");
            std::unique_ptr<TIterator> iter(prod->pdfList().createIterator());
            for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
                RooAbsPdf *cterm = dynamic_cast<RooAbsPdf *>(a); 
                if (!cterm) throw std::logic_error("AsimovUtils: a factor of the nuisance pdf is not a Pdf!");
                if (mvg_gen_obs && cterm->InheritsFrom("RooMultiVarGaussian")) {
                    std::unique_ptr<RooArgSet> comps(cterm->getComponents());
                    std::unique_ptr<RooArgSet> vars(cterm->getParameters(poi));
                    auto iter = vars->fwdIterator();
                    for (RooAbsArg *a = iter.next(); a != nullptr; a = iter.next()) {
                        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
                        if (rrv == nullptr) continue;
                        if (verbose > 1)
                            std::cout << "var: " << a->GetName() <<
                                     ", constant? " << rrv->isConstant() <<
                                     ", comp? " << (comps->find(*rrv) != nullptr) << 
                                     ", gobs? " << (gobs.find(*rrv) != nullptr) <<
                                     ", poi? " << (poi.find(*rrv) != nullptr) << std::endl;
                        if (poi.find(*rrv)) continue;
                        std::string vname = rrv->GetName();
                        if (vname.size() > 3 && vname.substr(vname.size() - 3) == "_In") {
                            auto match = comps->find(vname.substr(0, vname.size() - 3).c_str());
                            if (match != nullptr) {
                                if (verbose > 1)
                                    std::cout << "    --> matched to " << match->ClassName() << " " << match->GetName() << std::endl;
                                RooAbsReal *rfunc = dynamic_cast<RooAbsReal *>(match);
                                if (rfunc != nullptr) {
                                    if (rfunc->getVal() == rrv->getVal()) {
                                        if (verbose > 1)
                                            std::cout << "        --> the two already have the same " "value. nothing to do." << std::endl;
                                    } else {
                                        if (verbose > 1)
                                            std::cout << "        --> set " << a->GetName() << ", currently at " << rrv->getVal() << ", to " << rfunc->getVal() << std::endl;
                                        rrv->setVal(rfunc->getVal());
                                    }
                                }
                                if (!gobs.find(*rrv)) {
                                    if (verbose > 1)
                                        std::cout << "    --> added to the global observables: " << std::endl;
                                    gobs.add(*rrv);
                                }
                            }
                        }
                    }
                    continue;
                }
                if (!cterm->dependsOn(nuis)) continue; // dummy constraints
                if (typeid(*cterm) == typeid(RooUniform)) continue;
                std::unique_ptr<RooArgSet> cpars(cterm->getParameters(&gobs));
                std::unique_ptr<RooArgSet> cgobs(cterm->getObservables(&gobs));
                if (cgobs->getSize() != 1) {
                    throw std::runtime_error(Form("AsimovUtils: constraint term %s has multiple global observables", cterm->GetName()));
                }
                RooRealVar &rrv = dynamic_cast<RooRealVar &>(*cgobs->first());

                RooAbsReal *match = 0;
                if (cpars->getSize() == 1) {
                    match = dynamic_cast<RooAbsReal *>(cpars->first());
                } else {
                    std::unique_ptr<TIterator> iter2(cpars->createIterator());
                    for (RooAbsArg *a2 = (RooAbsArg *) iter2->Next(); a2 != 0; a2 = (RooAbsArg *) iter2->Next()) {
                        RooRealVar *rrv2 = dynamic_cast<RooRealVar *>(a2); 
                        if (rrv2 != 0 && !rrv2->isConstant()) {
                            if (match != 0) throw std::runtime_error(Form("AsimovUtils: constraint term %s has multiple floating params", cterm->GetName()));
                            match = rrv2;
                        }
                    }
                }
                if (match == 0) {   
                    std::cerr << "ERROR: AsimovUtils: can't find nuisance for constraint term " << cterm->GetName() << std::endl;
                    std::cerr << "Parameters: " << std::endl;
                    cpars->Print("V");
                    std::cerr << "Observables: " << std::endl;
                    cgobs->Print("V");
                    throw std::runtime_error(Form("AsimovUtils: can't find nuisance for constraint term %s", cterm->GetName()));
                }
                std::string pdfName(cterm->ClassName());
                if (pdfName == "RooGaussian" || pdfName == "SimpleGaussianConstraint"  || pdfName == "RooBifurGauss" || pdfName == "RooPoisson"  || pdfName == "SimplePoissonConstraint" || pdfName == "RooGenericPdf") {
                    // this is easy
                    rrv.setVal(match->getVal());
                } else if (pdfName == "RooGamma") {
                    // notation as in http://en.wikipedia.org/wiki/Gamma_distribution
                    //     nuisance = x
                    //     global obs = kappa ( = observed sideband events + 1)
                    //     scaling    = theta ( = extrapolation from sideband to signal)
                    // we want to set the global obs to a value for which the current value 
                    // of the nuisance is the best fit one.
                    // best fit x = (k-1)*theta    ---->  k = x/theta + 1
                    RooArgList leaves; cterm->leafNodeServerList(&leaves);
                    std::unique_ptr<TIterator> iter2(leaves.createIterator());
                    RooAbsReal *match2 = 0;
                    for (RooAbsArg *a2 = (RooAbsArg *) iter2->Next(); a2 != 0; a2 = (RooAbsArg *) iter2->Next()) {
                        RooAbsReal *rar = dynamic_cast<RooAbsReal *>(a2);
                        if (rar == 0 || rar == match || rar == &rrv) continue;
                        if (!rar->isConstant()) throw std::runtime_error(Form("AsimovUtils: extra floating parameter %s of RooGamma %s.", rar->GetName(), cterm->GetName()));
                        if (rar->getVal() == 0) continue; // this could be mu
                        if (match2 != 0) throw std::runtime_error(Form("AsimovUtils: extra constant non-zero parameter %s of RooGamma %s.", rar->GetName(), cterm->GetName()));
                        match2 = rar;
                    } 
                    if (match2 == 0) throw std::runtime_error(Form("AsimovUtils: could not find the scaling term for  RooGamma %s.", cterm->GetName()));
                    //std::cout << " nuisance "   << match->GetName() << " = x = " << match->getVal() << std::endl;
                    //std::cout << " scaling param "   << match2->GetName() << " = theta = " << match2->getVal() << std::endl;
                    //std::cout << " global obs " << rrv.GetName() << " = kappa = " << rrv.getVal() << std::endl;
                    //std::cout << " new value for global obs " << rrv.GetName() << " = x/theta+1 = " << (match->getVal()/match2->getVal() + 1) << std::endl;
                    rrv.setVal(match->getVal()/match2->getVal() + 1.);
                } else {
                    throw std::runtime_error(Form("AsimovUtils: can't handle constraint term %s of type %s", cterm->GetName(), cterm->ClassName()));
                }
            }

            // snapshot
            snapshot.removeAll();
            utils::setAllConstant(gobs, true);
            gobs.snapshot(snapshot);

            // revert things to normal
            gobs = snapGlobalObsData;
    
            if (verbose > 1) {
                std::cout << "Global observables for data: " << std::endl;
                snapGlobalObsData.Print("V");
                std::cout << "Global observables for asimov: " << std::endl;
                snapshot.Print("V");
            }
        }

        return asimov;
}
