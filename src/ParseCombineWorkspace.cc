#include "CombineHarvester/CombineTools/interface/ParseCombineWorkspace.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistErrorPropagator.h"
#include <iostream>
#include <string>
#include <vector>
#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooStats/ModelConfig.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace ch {
void ParseCombineWorkspacePy(CombineHarvester& cb, RooWorkspace const& ws,
                           std::string const& modelcfg,
                           std::string const& data, bool verbose) {
  ParseCombineWorkspace(cb, const_cast<RooWorkspace&>(ws), modelcfg, data, verbose);
}

void ParseCombineWorkspace(CombineHarvester& cb, RooWorkspace& ws,
                           std::string const& modelcfg_name,
                           std::string const& data_name, bool verbose) {
  bool v = verbose;

  RooStats::ModelConfig* cfg =
      dynamic_cast<RooStats::ModelConfig*>(ws.genobj(modelcfg_name.c_str()));
  if (!cfg) {
    throw std::runtime_error(
        FNERROR("Could not get ModelConfig from workspace"));
  }
  RooSimultaneous* pdf = dynamic_cast<RooSimultaneous*>(cfg->GetPdf());
  if (!pdf) {
    throw std::runtime_error(FNERROR("PDF is not a RooSimultaneous"));
  }

  std::vector<std::string> cats;

  struct ProcInfo {
    std::string bin;
    std::string proc;
    std::string pdf;
    std::string norm;
  };
  std::map<std::string, std::vector<ProcInfo>> proc_infos;

  RooAbsData *data = ws.data(data_name.c_str());
  // Hard-coded the category as "CMS_channel". Could deduce instead...
  std::unique_ptr<TList> split(data->split(RooCategory("CMS_channel", "")));
  for (int i = 0; i < split->GetSize(); ++i) {
    RooAbsData *idat = dynamic_cast<RooAbsData*>(split->At(i));
    cats.push_back(idat->GetName());
    FNLOGC(std::cout, v) << "Found data for category: " << cats.back() << "\n";

    ch::Observation obs;
    obs.set_bin(cats.back());
    cb.InsertObservation(obs);

    bool delete_pdfs = false;

    RooAbsReal *ipdf = FindAddPdf(pdf->getPdf(cats.back().c_str()));
    if (ipdf) {
      std::unique_ptr<RooArgSet> pdf_obs(ipdf->getObservables(data->get()));
      idat = idat->reduce(*pdf_obs);
      ws.import(*idat);
      RooAddPdf *ipdf_add = dynamic_cast<RooAddPdf*>(ipdf);
      RooRealSumPdf *ipdf_sum = dynamic_cast<RooRealSumPdf*>(ipdf);
      RooArgList const* coeffs = nullptr;
      RooArgList const* pdfs = nullptr;
      if (ipdf_add) {
        FNLOGC(std::cout, v) << "Found RooAddPdf: " << ipdf_add->GetName() << "\n";
        coeffs = &(ipdf_add->coefList());
        pdfs = &(ipdf_add->pdfList());
        if (coeffs->getSize() != pdfs->getSize()) {
          throw std::runtime_error(FNERROR("This RooAddPdf is not extended!"));
        }
      }
      if (ipdf_sum) {
        FNLOGC(std::cout, v) << "Found RooRealSumPdf: " << ipdf_sum->GetName() << "\n";
        coeffs = &(ipdf_sum->coefList());
        pdfs = &(ipdf_sum->funcList());
        if (coeffs->getSize() != pdfs->getSize()) {
          throw std::runtime_error(FNERROR("This RooRealSumPdf is not extended!"));
        }
        if (pdfs->getSize() == 1) {
           CMSHistErrorPropagator *err = dynamic_cast<CMSHistErrorPropagator*>(pdfs->at(0));
           if (err) {
             coeffs = &(err->coefList());
             pdfs = new RooArgList(err->wrapperList());
             delete_pdfs = true;
           }
        }
      }
      for (int j = 0; j < coeffs->getSize(); ++j) {
        RooAbsReal *jcoeff = dynamic_cast<RooAbsReal*>(coeffs->at(j));
        RooAbsReal *jpdf = dynamic_cast<RooAbsReal*>(pdfs->at(j));
        FNLOGC(std::cout, v) << "Component " << j << "\t" << jcoeff->GetName()
                             << "\t" << jpdf->GetName() << "\n";
        ch::Process proc;
        proc.set_bin(cats.back());
        // Get the process name & signal flag from the pdf attributes that are
        // set by text2workspace. Should really check they exist first...
        proc.set_process(jpdf->getStringAttribute("combine.process"));
        proc.set_rate(1.);
        proc.set_signal(jcoeff->getAttribute("combine.signal"));
        proc_infos[proc.bin()].push_back(
            {proc.bin(), proc.process(), jpdf->GetName(), jcoeff->GetName()});
        cb.InsertProcess(proc);
      }
      if (delete_pdfs) delete pdfs;
    }
  }
  cb.AddWorkspace(ws);
  cb.ExtractData(ws.GetName(), "$BIN");
  // This loop is slow if we have have a lot of procs
  // Better way would be to filter to each bin, then do all the procs for that
  // bin
  for (auto b : cb.bin_set()) {
    auto cb_bin = cb.cp().bin({b});
      for (auto const& pinfo: proc_infos[b]) {
        cb_bin.cp().process({pinfo.proc}).ExtractPdfs(
            cb, ws.GetName(), pinfo.pdf, pinfo.norm);
      }
  }
}

RooAbsReal* FindAddPdf(RooAbsReal* input) {
  RooAddPdf *as_add = dynamic_cast<RooAddPdf*>(input);
  RooRealSumPdf *as_sum = dynamic_cast<RooRealSumPdf*>(input);
  if (as_add) return as_add;
  if (as_sum) return as_sum;
  RooProdPdf *as_prod = dynamic_cast<RooProdPdf*>(input);
  if (as_prod) {
    RooArgList const& comps = as_prod->pdfList();
    for (int i = 0; i < comps.getSize(); ++i) {
      RooAbsReal* try_add = FindAddPdf(dynamic_cast<RooAbsReal*>(comps.at(i)));
      if (try_add) return try_add;
    }
  }
  return nullptr;
}
}
