#include "CombineTools/interface/ParseCombineWorkspace.h"
#include <iostream>
#include <string>
#include <vector>
#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooStats/ModelConfig.h"
#include "CombineTools/interface/Logging.h"

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
    ws.import(*idat);

    ch::Observation obs;
    obs.set_bin(cats.back());
    cb.InsertObservation(obs);

    RooAddPdf *ipdf = FindAddPdf(pdf->getPdf(cats.back().c_str()));
    if (ipdf) {
      FNLOGC(std::cout, v) << "Found RooAddPdf: " << ipdf->GetName() << "\n";
      RooArgList const& coeffs = ipdf->coefList();
      RooArgList const& pdfs = ipdf->pdfList();
      if (coeffs.getSize() != pdfs.getSize()) {
        throw std::runtime_error(FNERROR("This RooAddPdf is not extended!"));
      }
      for (int j = 0; j < coeffs.getSize(); ++j) {
        RooAbsReal *jcoeff = dynamic_cast<RooAbsReal*>(coeffs.at(j));
        RooAbsPdf *jpdf = dynamic_cast<RooAbsPdf*>(pdfs.at(j));
        FNLOGC(std::cout, v) << "Component " << j << "\t" << jcoeff->GetName()
                             << "\t" << jpdf->GetName() << "\n";
        ch::Process proc;
        proc.set_bin(cats.back());
        // Get the process name & signal flag from the pdf attributes that are
        // set by text2workspace. Should really check they exist first...
        proc.set_process(jpdf->getStringAttribute("combine.process"));
        proc.set_rate(1.);
        proc.set_signal(jpdf->getAttribute("combine.signal"));
        proc_infos[proc.bin()].push_back(
            {proc.bin(), proc.process(), jpdf->GetName(), jcoeff->GetName()});
        cb.InsertProcess(proc);
      }
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

RooAddPdf* FindAddPdf(RooAbsPdf* input) {
  RooAddPdf *as_add = dynamic_cast<RooAddPdf*>(input);
  if (as_add) return as_add;
  RooProdPdf *as_prod = dynamic_cast<RooProdPdf*>(input);
  if (as_prod) {
    RooArgList const& comps = as_prod->pdfList();
    for (int i = 0; i < comps.getSize(); ++i) {
      RooAddPdf* try_add = FindAddPdf(dynamic_cast<RooAbsPdf*>(comps.at(i)));
      if (try_add) return try_add;
    }    
  }
  return nullptr;
}
}
