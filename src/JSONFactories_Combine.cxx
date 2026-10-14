#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFitHS3/JSONIO.h>
#include "../interface/static_execute.h"

#include "RooFit.h"

#include "../interface/ProcessNormalization.h"

#include "RooAddition.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooNumber.h"
#include "RooPoisson.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <charconv>
#include <map>
#include <memory>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using RooFit::Detail::JSONNode;

namespace {

bool parse_double(const std::string &input, double &out)
{
   const char *begin = input.data();
   const char *end = input.data() + input.size();

   double tmp;
   auto [ptr, ec] = std::from_chars(begin, end, tmp);

   if ((ec == std::errc() && ptr == end)) {
      out = tmp;
      return true;
   }
   return false;
}

inline bool approxSymmetricKappas(double kLo, double kHi, double absLogTol = 1e-10, double relLogTol = 1e-10)
{
   // kappas must be positive and finite
   if (!(kLo > 0.0) || !(kHi > 0.0))
      return false;
   if (!std::isfinite(kLo) || !std::isfinite(kHi))
      return false;

   const double a = std::log(kLo);
   const double b = std::log(kHi);

   // symmetric if a ~ -b  <=>  (a + b) ~ 0
   const double diff = std::abs(a + b);
   const double scale = std::max(std::abs(a), std::abs(b)); // 0 if both ~1
   return diff <= std::max(absLogTol, relLogTol * scale);
}

/// Helper: ensure/lookup a parameter (RooRealVar) used by normsys/normfactor modifiers
RooRealVar &getOrMakeParameter(RooJSONFactoryWSTool *tool, std::string const &parName, double init = 1.0,
                               double lo = 0.0, double hi = 10.0)
{
   if (auto *arg = tool->workspace()->arg(parName.c_str())) {
      auto *rrv = dynamic_cast<RooRealVar *>(arg);
      if (!rrv)
         throw std::runtime_error("histfactory_dist importer: parameter '" + parName + "' exists but is not RooRealVar");
      return *rrv;
   }
   auto *rrv = new RooRealVar(parName.c_str(), parName.c_str(), init, lo, hi);
   tool->workspace()->import(*rrv, RooFit::RecycleConflictNodes());
   return *dynamic_cast<RooRealVar *>(tool->workspace()->arg(parName.c_str()));
}


void exportOtherFactorModifiers(RooJSONFactoryWSTool *tool, JSONNode &mods, RooAbsReal *factor)
{
   auto &m = mods.append_child().set_map();
   m["name"] << factor->GetName();
   m["type"] << "normfactor";
   if (auto *rcv = dynamic_cast<RooConstVar *>(factor)) {
      m["parameter"] << rcv->getVal();
   } else {
      m["parameter"] << factor->GetName();
      tool->queueExport(*factor);
   }
}


bool isSingleBinCountingDist(JSONNode const &elem)
{
   auto const *axes = elem.find("axes");
   if (!axes || axes->num_children() != 1)
      return false;
   for (auto const &ax : axes->children()) {
      if (ax.has_child("edges")) {
         if (ax["edges"].num_children() != 2) // two edges == one bin
            return false;
      } else {
         auto const *nb = ax.find("nbins");
         if (!nb || nb->val_int() != 1)
            return false;
         // the count axis n_obs_bin<b> is unbounded above (writeAxisBound: null); a 1-bin shape
         // channel (CMS_TH1obs_<channel>, [0, 1]) is not
         auto const *mx = ax.find("max");
         if (mx && !mx->is_null())
            return false;
      }
   }

   auto const *samples = elem.find("samples");
   if (!samples || samples->num_children() == 0)
      return false;
   for (auto const &s : samples->children()) {
      auto const *data = s.find("data");
      if (!data || !data->has_child("contents") || (*data)["contents"].num_children() != 1)
         return false;
      if (auto const *mods = s.find("modifiers")) {
         for (auto const &m : mods->children()) {
            std::string const mtype = m["type"].val();
            if (mtype != "normsys" && mtype != "normfactor" && mtype != "custom")
               return false;
         }
      }
   }
   return true;
}


void writeAxisBound(JSONNode &node, double value)
{
   if (RooNumber::isInfinite(value)) {
      node.set_null();
   } else {
      node << value;
   }
}

double readAxisBound(JSONNode const &node, const char *key, double defaultValue)
{
   auto const *bound = node.find(key);
   if (!bound || bound->is_null())
      return defaultValue;
   return bound->val_double();
}

/// Write the single-element data.contents array of a sample.
void writeSingleBinContents(JSONNode &sample, double value)
{
   auto &dataNode = sample["data"].set_map();
   RooJSONFactoryWSTool::exportArray(1, &value, dataNode["contents"]);
}

std::string resolveConstraintPdfName(RooJSONFactoryWSTool *tool, RooAbsArg const &theta,
                                     RooAbsCollection const &observables,
                                     std::map<std::string, std::string> &cache)
{
   const std::string thetaName(theta.GetName());
   auto cached = cache.find(thetaName);
   if (cached != cache.end())
      return cached->second;

   RooWorkspace *ws = tool->workspace();
   RooAbsPdf *constraint = nullptr;

   // Fast path: ShapeTools.py names every constraint pdf "<theta>_Pdf".
   if (auto *byName = ws->pdf((thetaName + "_Pdf").c_str())) {
      if (byName->dependsOn(theta))
         constraint = byName;
   }

   if (!constraint) {
      for (auto *arg : ws->allPdfs()) {
         auto *pdf = dynamic_cast<RooAbsPdf *>(arg);
         if (!pdf || !pdf->dependsOn(theta))
            continue;
         if (pdf->dependsOn(observables))
            continue;
         constraint = pdf;
         break;
      }
   }

   std::string name;
   if (constraint) {
      name = constraint->GetName();
      tool->queueExport(*constraint);
   }
   cache[thetaName] = name;
   return name;
}

void exportNormSys(RooJSONFactoryWSTool *tool, JSONNode &mods, RooAbsArg const &theta, double kappaLo, double kappaHi,
                   RooAbsCollection const &observables, std::map<std::string, std::string> &constraintCache)
{
   auto &m = mods.append_child().set_map();
   m["name"] << theta.GetName();
   m["parameter"] << theta.GetName();
   m["type"] << "normsys";
   const std::string constraint = resolveConstraintPdfName(tool, theta, observables, constraintCache);
   if (!constraint.empty())
      m["constraint"] << constraint;
   auto &data = m["data"].set_map();
   data["lo"] << kappaLo;
   data["hi"] << kappaHi;
   tool->queueExport(theta);
}

bool tryExportCountingHistFactory(RooJSONFactoryWSTool *tool, const RooPoisson *pois, JSONNode &elem)
{
   auto const *x = dynamic_cast<RooRealVar const *>(&pois->getX());
   if (!x) {
      return false;
   }
   auto const *add = dynamic_cast<RooAddition const *>(&pois->getMean());
   if (!add) {
      return false;
   }

   RooArgList const &summands = add->list();
   if (summands.getSize() == 0) {
      RooJSONFactoryWSTool::warning("SumOfRooPoisson '" + std::string(pois->GetName()) + "': no processes found!");
      return false;
   }
   for (auto *summand : summands) {
      if (!dynamic_cast<ProcessNormalization const *>(summand) && !dynamic_cast<RooRealVar const *>(summand) &&
          !dynamic_cast<RooConstVar const *>(summand)) {
         RooJSONFactoryWSTool::warning("SumOfRooPoisson '" + std::string(pois->GetName()) + "': summand '" +
                                       summand->GetName() + "' is a " + summand->ClassName() +
                                       ", not a ProcessNormalization or a constant!");
         return false;
      }
   }

   auto &axes = elem["axes"].set_seq();
   auto &ax = axes.append_child().set_map();
   std::string axname(x->GetName());
   RooJSONFactoryWSTool::testValidName(axname, false);
   ax["name"] << axname;
   ax["nbins"] << 1;

   
   writeAxisBound(ax["min"], x->getMin());
   writeAxisBound(ax["max"], x->getMax());

   auto &samples = elem["samples"].set_seq();

   RooArgSet observables;
   observables.add(*x);
   std::map<std::string, std::string> constraintCache;

   for (auto *summand : summands) {
      std::string sampleName(summand->GetName());
      RooJSONFactoryWSTool::testValidName(sampleName, false);

      auto &s = samples.append_child().set_map();
      s["name"] << sampleName;

      if (auto const *pn = dynamic_cast<ProcessNormalization const *>(summand)) {
         writeSingleBinContents(s, pn->nominalValue());
         auto &mods = s["modifiers"].set_seq();

         // symmetric log-normals: kappaLo == 1/kappaHi
         auto const &logKappa = pn->logKappa();
         for (int i = 0; i < static_cast<int>(logKappa.size()); ++i) {
            exportNormSys(tool, mods, pn->thetaList()[i], std::exp(-logKappa[i]), std::exp(logKappa[i]), observables,
                          constraintCache);
         }
         // asymmetric log-normals
         auto const &logAsymmKappa = pn->logAsymmKappa();
         for (int i = 0; i < static_cast<int>(logAsymmKappa.size()); ++i) {
            exportNormSys(tool, mods, pn->asymmThetaList()[i], std::exp(logAsymmKappa[i].first),
                          std::exp(logAsymmKappa[i].second), observables, constraintCache);
         }
         for (auto *of : pn->otherFactorList()) {
            exportOtherFactorModifiers(tool, mods, static_cast<RooAbsReal *>(of));
         }
      } else if (auto const *rrv = dynamic_cast<RooRealVar const *>(summand); rrv && !rrv->isConstant()) {
         // A floating yield: keep it floating instead of freezing its current value.
         writeSingleBinContents(s, 1.0);
         auto &mods = s["modifiers"].set_seq();
         auto &m = mods.append_child().set_map();
         m["name"] << rrv->GetName();
         m["parameter"] << rrv->GetName();
         m["type"] << "normfactor";
         tool->queueExport(*rrv);
      } else {
         writeSingleBinContents(s, static_cast<RooAbsReal const *>(summand)->getVal());
         s["modifiers"].set_seq();
      }
   }

   return true;
}

class SumOfRooPoissonStreamer_RooAddition : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring("histfactory_dist");
      return keystring;
   }

   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {

      auto const *pois = dynamic_cast<const RooPoisson *>(p);
      if (!pois)
         return false;
      if (!tryExportCountingHistFactory(tool, pois, elem))
         return false;
      elem["type"] << key();
           std::cout << "exporting RooPoisson of RooAddition to histfactory" << std::endl;

      return true;
   }
};

class SumOfRooPoissonFactory_RooAddition : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, JSONNode const &elem) const override
   {
      if (!elem.has_child("name"))
         return false;
      // Multi-bin (and any shape-carrying) node belongs to the hist-stack importers.
      if (!isSingleBinCountingDist(elem))
         return false;

      const std::string distName = elem["name"].val();

      std::string axName;
      double axMin = -RooNumber::infinity();
      double axMax = RooNumber::infinity();
      for (auto const &ax : elem["axes"].children()) {
         axName = ax["name"].val();
         axMin = readAxisBound(ax, "min", -RooNumber::infinity());
         axMax = readAxisBound(ax, "max", RooNumber::infinity());
      }
      if (axName.empty())
         throw std::runtime_error("SumOfRooPoisson importer: '" + distName + "' has no usable axis");

      auto *x = tool->workspace()->var(axName);
      std::unique_ptr<RooRealVar> ownedX;
      if (!x) {
         // The value-only constructor leaves the range at (-inf, +inf), matching how combine
         // creates n_obs_bin<b>; only finite bounds have to be applied on top.
         ownedX = std::make_unique<RooRealVar>(axName.c_str(), axName.c_str(), 0.0);
         if (!RooNumber::isInfinite(axMin))
            ownedX->setMin(axMin);
         if (!RooNumber::isInfinite(axMax))
            ownedX->setMax(axMax);
         x = ownedX.get();
      }

      RooArgList summands;

      for (const auto &s : elem["samples"].children()) {
         const std::string sampleName = s["name"].val();

         // Samples are named after the ProcessNormalization they came from, so a process shared
         // between the s and the b-only distribution is reconstructed exactly once.
         if (auto *existing = tool->workspace()->function(sampleName)) {
            summands.add(*existing);
            continue;
         }

         double nominal = 0.;
         for (auto const &c : s["data"]["contents"].children()) {
            nominal = c.val_double();
            break;
         }

         RooArgList normsysPars;
         std::vector<std::pair<double, double>> normsysVals;
         RooArgList normfactorPars;

         if (s.has_child("modifiers")) {
            for (const auto &m : s["modifiers"].children()) {
               const std::string mtype = m["type"].val();
               const std::string mname = m["name"].val();

               if (mtype == "normsys") {
                  auto const *par = m.find("parameter");
                  std::string parname(par ? par->val() : "alpha_" + mname);
                  RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10.0, 10.0);
                  auto const &data = m["data"];
                  normsysPars.add(p);
                  normsysVals.emplace_back(data["lo"].val_double(), data["hi"].val_double());
               } else if (mtype == "normfactor") {
                  auto const *par = m.find("parameter");
                  const std::string parval(par ? par->val() : mname);
                  double literal = 0.;
                  if (parse_double(parval, literal)) {
                     if (auto *existingCv = tool->workspace()->function(mname)) {
                        normfactorPars.add(*existingCv);
                     } else {
                        RooConstVar cv(mname.c_str(), mname.c_str(), literal);
                        tool->workspace()->import(cv, RooFit::RecycleConflictNodes());
                        normfactorPars.add(*tool->workspace()->function(mname));
                     }
                  } else {
                     normfactorPars.add(*tool->request<RooAbsReal>(parval, sampleName));
                  }
               } else if (mtype == "custom") {
                  normfactorPars.add(*tool->request<RooAbsReal>(mname, sampleName));
               } else {
                  throw std::runtime_error("SumOfRooPoisson importer: modifier type '" + mtype + "' on sample '" +
                                           sampleName + "' is not supported by this importer");
               }
            }
         }

         ProcessNormalization pn(sampleName.c_str(), sampleName.c_str(), nominal);
         for (int i = 0; i < normsysPars.getSize(); ++i) {
            auto &par = static_cast<RooAbsReal &>(normsysPars[i]);
            if (approxSymmetricKappas(normsysVals[i].first, normsysVals[i].second)) {
               pn.addLogNormal(normsysVals[i].second, par);
            } else {
               pn.addAsymmLogNormal(normsysVals[i].first, normsysVals[i].second, par);
            }
         }
         for (int i = 0; i < normfactorPars.getSize(); ++i) {
            pn.addOtherFactor(static_cast<RooAbsReal &>(normfactorPars[i]));
         }
         tool->workspace()->import(pn, RooFit::RecycleConflictNodes());

         auto *imported = tool->workspace()->function(sampleName);
         if (!imported)
            throw std::runtime_error("SumOfRooPoisson importer: failed to import ProcessNormalization '" + sampleName +
                                     "'");
         summands.add(*imported);
      }

      // The name of the inner sum is derived from the distribution name. 
      // Sample names and the distribution name itself round-trip exactly.
      const std::string sumName = distName + "_n_exp";
      RooAddition addition(sumName.c_str(), sumName.c_str(), summands);
      tool->workspace()->import(addition, RooFit::RecycleConflictNodes());
      auto *n_exp = tool->workspace()->function(sumName);
      if (!n_exp)
         throw std::runtime_error("SumOfRooPoisson importer: failed to import RooAddition '" + sumName + "'");

      // combine always builds these counting Poissons with noRounding=1.
      RooPoisson pdf(distName.c_str(), distName.c_str(), *x, *n_exp, /*noRounding=*/true);
      tool->workspace()->import(pdf, RooFit::RecycleConflictNodes());
      if (!tool->workspace()->pdf(distName))
         throw std::runtime_error("SumOfRooPoisson importer: failed to import RooPoisson '" + distName + "'");

      return true;
   }
};

}

STATIC_EXECUTE([]() {
   using namespace RooFit::JSONIO;

   registerExporter<SumOfRooPoissonStreamer_RooAddition>("RooPoisson", true);
   registerImporter<SumOfRooPoissonFactory_RooAddition>("histfactory_dist", true);
});
