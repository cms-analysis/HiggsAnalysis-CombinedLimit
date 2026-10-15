#include "../interface/ChannelSumNLL.h"

// ROOT_VERSION: ROOT/RVersion.hxx exists since ROOT 6.30, older versions only have RVersion.h.
#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#include "../interface/CMSHistBarlowBeestonGradientTerm.h"
#include "../interface/CMSHistHelpers.h"
#include "../interface/CombineLogger.h"
#include "../interface/RooMultiPdfCombine.h"
#include "../interface/RooSimultaneousOpt.h"
#include "../interface/SimNLLBase.h"
#include "../interface/SimpleGaussianConstraint.h"
#include "../interface/SimplePoissonConstraint.h"
#include "../interface/utils.h"

#include <RooAbsCategoryLValue.h>
#include <RooCmdArg.h>
#include <RooConstVar.h>
#include <RooGlobalFunc.h>
#include <RooPoisson.h>
#include <RooRealSumPdf.h>
#include <RooUniform.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
#include <RooEvaluatorWrapper.h>
#endif

#include <Math/Util.h>
#include <TMath.h>
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 37, 0)
#include <TList.h>
#include <TObject.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

static std::string subName(RooAbsArg const &owner, const char *suffix) {
  return std::string(owner.GetName()) + "_" + suffix;
}

// The "codegen" backend of RooFit evaluates the likelihood value with the
// Evaluator like the "cpu" backend and uses the generated code only for the
// gradient. The generated function is much cheaper to call, in particular in
// the phases where the minimizer evaluates the function many times with a
// single parameter changed (numerical seed, line searches), so the value is
// taken from the generated code as well where RooFit allows it.
static void useGeneratedFunctionCode(RooAbsReal &nll) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
  if (auto *wrapper = dynamic_cast<RooFit::Experimental::RooEvaluatorWrapper *>(&nll)) {
    wrapper->setUseGeneratedFunctionCode(true);
  }
#else
  (void)nll;
#endif
}

/// The dataset named like the channel, or nullptr if there is none.
static RooAbsData *findDataset(std::vector<std::unique_ptr<RooAbsData>> const &datasets, std::string const &label) {
  for (auto const &dataset : datasets) {
    if (label == dataset->GetName())
      return dataset.get();
  }
  return nullptr;
}

cacheutils::ChannelSumNLL::ChannelSumNLL(const char *name,
                                         const char *title,
                                         RooSimultaneous &simPdf,
                                         RooAbsData &data,
                                         const RooArgSet *constrain,
                                         std::string evalBackend)
    : SimNLLBase(name, title),
      m_simPdf(&simPdf),
      m_data(&data),
      m_hasConstrain(constrain != nullptr),
      m_evalBackend(std::move(evalBackend)),
      m_originalCategories(
          /*inName=*/"originalCategories", "Categories replaced by private copies in the channel pdfs", this),
      m_components(/*inName=*/"components", "Likelihood terms", this) {
  if (constrain)
    m_constrain.add(*constrain);
  setup(data);
}

cacheutils::ChannelSumNLL::ChannelSumNLL(const ChannelSumNLL &other, const char *name)
    : SimNLLBase(other, name),
      m_simPdf(other.m_simPdf),
      m_data(other.m_data),
      m_hasConstrain(other.m_hasConstrain),
      m_constrain(other.m_constrain),
      m_evalBackend(other.m_evalBackend),
      m_originalCategories(
          /*inName=*/"originalCategories", "Categories replaced by private copies in the channel pdfs", this),
      m_components(/*inName=*/"components", "Likelihood terms", this) {
  // The RooFit likelihoods can't be shared, so the copy rebuilds them from the
  // same pdf and data; it starts unmasked and without zero points.
  setup(*m_data);
}

void cacheutils::ChannelSumNLL::setup(RooAbsData &data) {
  // Split off the constraint terms, deduplicated across channels and including
  // the extra constraints of a RooSimultaneousOpt.
  RooAbsPdf *factorized = utils::factorizePdf(*data.get(), *m_simPdf, m_constraints);
  m_factorizedPdf = dynamic_cast<RooSimultaneous *>(factorized);
  if (!m_factorizedPdf) {
    throw std::runtime_error(std::string("ChannelSumNLL: the pdf ") + m_simPdf->GetName() +
                             " did not factorize into a RooSimultaneous");
  }
  if (factorized != m_simPdf)
    m_factorizedPdfOwner.reset(factorized);

  RooAbsCategoryLValue const &indexCat = m_factorizedPdf->indexCat();
  // A RooSimultaneousOpt without channels (the FitDiagnostics prefit fit) has
  // nothing to split: its category has no states and is not in the data.
  if (data.get()->find(indexCat.GetName()))
    m_datasets = splitData(data);

  RooFit::EvalBackend backendArg(m_evalBackend);
  const bool codegen = backendArg.value() == RooFit::EvalBackend::Value::Codegen;

  // The channel masks of a RooSimultaneousOpt are aligned with the bins of
  // the index category, i.e. with the states in insertion order like
  // RooAbsCategoryLValue::setBin(), which is not the iteration order.
  RooArgList const *maskVars = nullptr;
  if (auto const *simOpt = dynamic_cast<RooSimultaneousOpt const *>(m_simPdf)) {
    if (!simOpt->channelMasks().empty())
      maskVars = &simOpt->channelMasks();
  }
  auto const nBins = static_cast<std::size_t>(indexCat.numBins(/*rangeName=*/nullptr));
  if (maskVars && maskVars->size() != nBins) {
    throw std::logic_error(std::string("ChannelSumNLL: the pdf ") + m_simPdf->GetName() + " has " +
                           std::to_string(maskVars->size()) + " channel masks for " + std::to_string(nBins) +
                           " channels");
  }

  // One RooFit likelihood per channel; the channel pdfs hold no constraints anymore.
  for (std::size_t iBin = 0; iBin < nBins; ++iBin) {
    std::string const &label = indexCat.getOrdinal(iBin).first;
    RooAbsPdf *pdf = m_factorizedPdf->getPdf(label.c_str());
    if (!pdf)
      continue;
    RooAbsData *channelData = findDataset(m_datasets, label);
    if (!channelData) {
      throw std::logic_error("ChannelSumNLL: no data for channel " + label);
    }
    // The RooMultiPdfs of the original pdf, for the penalty of the original
    // index, and its histogram nodes: they share their parameters with the
    // copies inside the RooFit likelihood, so the Barlow-Beeston step can run
    // on them.
    std::vector<RooMultiPdf const *> multiPdfs;
    std::vector<RooAbsReal const *> hists;
    std::unique_ptr<RooArgSet> components{pdf->getComponents()};
    for (RooAbsArg const *component : *components) {
      if (auto const *multiPdf = dynamic_cast<RooMultiPdf const *>(component))
        multiPdfs.push_back(multiPdf);
      if (isCMSHist(*component))
        hists.push_back(dynamic_cast<RooAbsReal const *>(component));
    }
    std::vector<CategorySync> syncs;
    RooAbsPdf *nllPdf = &channelPdfForNLL(*pdf, syncs);
    m_channelGradientTerms.emplace_back();
    if (codegen) {
      nllPdf = &addBarlowBeestonGradientTerms(*nllPdf, hists);
    }
    std::unique_ptr<RooAbsReal> nll{nllPdf->createNLL(*channelData, backendArg)};
    if (!nll) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of channel " + label);
    }
    nll->SetName(subName(*this, label.c_str()).c_str());
    if (codegen) {
      useGeneratedFunctionCode(*nll);
    }
    m_channelNames.push_back(label);
    m_channelMaskVars.push_back(maskVars ? dynamic_cast<RooAbsReal const *>(maskVars->at(static_cast<int>(iBin)))
                                         : nullptr);
    m_categorySyncs.push_back(syncs);
    m_channelMultiPdfs.push_back(multiPdfs);
    m_channelHists.push_back(hists);
    setHistogramData(m_channelHists.size() - 1, *channelData);
    RooArgSet params;
    nll->getParameters(/*observables=*/nullptr, params, /*stripDisconnected=*/false);
    replaceClonesByOriginals(params);
    params.remove(m_gradientTermVars, /*silent=*/true);
    m_channelParams.push_back(params);
    m_components.add(*nll);
    m_terms.push_back(std::move(nll));
  }
  m_nChannels = m_terms.size();

  // The simple constraint terms are summed in closed form, the others get a
  // likelihood of their own, see the class comment. Passing the constraint
  // pdfs explicitly also covers constraints on parameters outside of the
  // "constrain" set (e.g. the "constr" directive).
  RooArgSet rooFitConstraints;
  for (RooAbsArg const *arg : m_constraints) {
    auto const &pdf = dynamic_cast<RooAbsPdf const &>(*arg);
    if (!addFastConstraint(pdf))
      rooFitConstraints.add(pdf);
  }
  if (!rooFitConstraints.empty()) {
    m_constraintDummyObs =
        std::make_unique<RooRealVar>(subName(*this, /*suffix=*/"constraintDummyObs").c_str(), "", 0.5, 0., 1.);
    m_constraintPdf =
        std::make_unique<RooUniform>(subName(*this, /*suffix=*/"constraintPdf").c_str(), "", *m_constraintDummyObs);
    m_constraintData = std::make_unique<RooDataSet>(
        subName(*this, /*suffix=*/"constraintData").c_str(), "", RooArgSet(*m_constraintDummyObs));
    m_constraintData->add(RooArgSet(*m_constraintDummyObs));
    RooCmdArg constrainArg = m_hasConstrain ? RooFit::Constrain(m_constrain) : RooCmdArg::none();
    std::unique_ptr<RooAbsReal> nll{m_constraintPdf->createNLL(
        *m_constraintData, backendArg, RooFit::ExternalConstraints(rooFitConstraints), constrainArg)};
    if (!nll) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of the constraint terms");
    }
    nll->SetName(subName(*this, /*suffix=*/"constraints").c_str());
    if (codegen) {
      useGeneratedFunctionCode(*nll);
    }
    m_components.add(*nll);
    m_terms.push_back(std::move(nll));
  }
  m_zeroPoints.assign(m_terms.size(), 0.);

  m_hasGradient = m_nChannels > 0 &&
                  std::all_of(m_terms.begin(), m_terms.end(), [](auto const &term) { return term->hasGradient(); });

  CombineLogger::instance().log("ChannelSumNLL.cc",
                                __LINE__,
                                "Created " + std::string(GetName()) + " with " + std::to_string(m_nChannels) +
                                    " channels and " + std::to_string(m_constraints.size()) + " constraint terms (" +
                                    std::to_string(m_fastConstraintPdfs.size()) +
                                    " in closed form), RooFit backend \"" + m_evalBackend +
                                    "\", analytic gradient: " + (m_hasGradient ? "yes" : "no"),
                                __func__);

  m_gradientMapsDirty = true;
  setValueDirty();
}

RooAbsPdf &cacheutils::ChannelSumNLL::channelPdfForNLL(RooAbsPdf &pdf, std::vector<CategorySync> &syncs) {
  // All parameters, also those of the inactive pdfs of a RooMultiPdf.
  std::unique_ptr<RooArgSet> allParams{utils::getAllParameters(pdf)};
  RooArgSet privateCategories;
  for (RooAbsArg *arg : *allParams) {
    auto *original = dynamic_cast<RooCategory *>(arg);
    if (!original)
      continue;
    // Same name as the original, so that servers and parameters match by name.
    auto clone = std::make_unique<RooCategory>(*original);
    syncs.push_back({original, clone.get()});
    privateCategories.add(*clone);
    m_categoryClones.push_back(std::move(clone));
    if (!m_originalCategories.containsInstance(*original))
      m_originalCategories.add(*original);
  }
  if (privateCategories.empty())
    return pdf;

  // cloneTree() copies everything: reattach the parameters to the shared
  // originals, then replace the categories by the private copies.
  auto *clonedPdf = dynamic_cast<RooAbsPdf *>(pdf.cloneTree());
  if (!clonedPdf) {
    throw std::runtime_error(std::string("ChannelSumNLL: failed to clone the pdf ") + pdf.GetName());
  }
  m_privateChannelPdfs.emplace_back(clonedPdf);
  clonedPdf->recursiveRedirectServers(*allParams);
  clonedPdf->recursiveRedirectServers(privateCategories);
  return *clonedPdf;
}

RooAbsPdf &cacheutils::ChannelSumNLL::addBarlowBeestonGradientTerms(RooAbsPdf &pdf,
                                                                    std::vector<RooAbsReal const *> const &hists) {
  if (hists.empty())
    return pdf;
  // Combine puts the histogram node of a channel directly under an extended
  // RooRealSumPdf with a unit coefficient. The sum gets a gradient term per
  // histogram node, with a unit coefficient as well: it adds nothing to the
  // value and the correction to the normalization integral. The term is fed
  // from the original node, which has the same name and parameters as the
  // copy inside the pdf given to RooFit (see channelPdfForNLL()).
  auto *sumPdf = dynamic_cast<RooRealSumPdf *>(&pdf);
  if (!sumPdf || sumPdf->funcList().size() != sumPdf->coefList().size()) {
    throw std::runtime_error(std::string("ChannelSumNLL: the channel pdf ") + pdf.GetName() +
                             " has histogram nodes with bin-by-bin parameters, but it is not an extended "
                             "RooRealSumPdf, which the analytic gradient of the Barlow-Beeston minimization needs");
  }
  RooArgList funcs(sumPdf->funcList());
  RooArgList coefs(sumPdf->coefList());
  std::vector<std::unique_ptr<CMSHistBarlowBeestonGradientTerm>> terms;
  for (RooAbsReal const *hist : hists) {
    RooAbsReal const &x = cmsHistVisit(*hist, [](auto const &h) -> RooAbsReal const & { return h.getXVar(); });
    terms.push_back(std::make_unique<CMSHistBarlowBeestonGradientTerm>(
        subName(*this, (std::string(hist->GetName()) + "_bbGradient").c_str()).c_str(),
        "",
        const_cast<RooAbsReal &>(x),  // NOLINT(cppcoreguidelines-pro-type-const-cast)
        *hist));
    m_gradientTermVars.add(terms.back()->slopes());
    m_gradientTermVars.add(terms.back()->ratios());
    funcs.add(*terms.back());
    coefs.add(RooFit::RooConst(1.));
  }
  auto newPdf = std::make_unique<RooRealSumPdf>(
      subName(*this, (std::string(pdf.GetName()) + "_bbGradient").c_str()).c_str(), "", funcs, coefs, /*extended=*/true);
  for (auto &term : terms) {
    m_channelGradientTerms.back().push_back(term.get());
    newPdf->addOwnedComponents(std::move(term));
  }
  m_privateChannelPdfs.push_back(std::move(newPdf));
  return *m_privateChannelPdfs.back();
}

void cacheutils::ChannelSumNLL::syncCategories() const {
  for (auto const &syncs : m_categorySyncs) {
    for (CategorySync const &sync : syncs) {
      if (sync.clone->getCurrentIndex() != sync.original->getCurrentIndex()) {
        sync.clone->setIndex(sync.original->getCurrentIndex());
      }
    }
  }
}

void cacheutils::ChannelSumNLL::replaceClonesByOriginals(RooArgSet &params) const {
  for (auto const &syncs : m_categorySyncs) {
    for (CategorySync const &sync : syncs) {
      if (params.containsInstance(*sync.clone)) {
        params.remove(*sync.clone);
        params.add(*sync.original, /*silent=*/true);
      }
    }
  }
}

std::vector<std::unique_ptr<RooAbsData>> cacheutils::ChannelSumNLL::splitData(RooAbsData &data) const {
  // Split by the index category, keeping weights and the zero-weight entries
  // that binned likelihoods need, with pieces also for empty channels.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 37, 0)
  return data.split(*m_factorizedPdf, /*createEmptyDataSets=*/true);
#else
  // Before 6.37, split() returns a non-owning TList of the new datasets.
  std::unique_ptr<TList> list{data.split(*m_factorizedPdf, /*createEmptyDataSets=*/true)};
  std::vector<std::unique_ptr<RooAbsData>> datasets;
  datasets.reserve(list->GetSize());
  for (TObject *obj : *list) {
    datasets.emplace_back(dynamic_cast<RooAbsData *>(obj));
  }
  return datasets;
#endif
}

bool cacheutils::ChannelSumNLL::setData(RooAbsData &data, bool /*cloneData*/) {
  if (m_nChannels == 0) {
    m_data = &data;
    return true;
  }
  auto newDatasets = splitData(data);
  for (std::size_t i = 0; i < m_nChannels; ++i) {
    std::string const &label = m_channelNames[i];
    RooAbsData *channelData = findDataset(newDatasets, label);
    if (!channelData) {
      throw std::logic_error("ChannelSumNLL::setData: no data for channel " + label);
    }
    // The channel likelihoods keep a pointer to the data, we own it.
    m_terms[i]->setData(*channelData, /*cloneData=*/false);
    setHistogramData(i, *channelData);
  }
  // Only release the old datasets after all terms point to the new ones.
  m_datasets = std::move(newDatasets);
  m_data = &data;
  setValueDirty();
  return true;
}

void cacheutils::ChannelSumNLL::setHistogramData(std::size_t iChannel, RooAbsData const &channelData) const {
  for (RooAbsReal const *hist : m_channelHists[iChannel]) {
    cmsHistVisit(*hist, [&](auto const &h) { h.setData(channelData); });
  }
}

void cacheutils::ChannelSumNLL::setAnalyticBarlowBeeston(bool flag) {
  for (std::size_t i = 0; i < m_nChannels; ++i) {
    // Like in CachingSimNLL, a masked channel never has it switched on.
    bool channelFlag = flag && !isMasked(i);
    for (RooAbsReal const *hist : m_channelHists[i]) {
      cmsHistVisit(*hist, [&](auto const &h) { h.setAnalyticBarlowBeeston(channelFlag); });
    }
  }
  m_analyticBarlowBeeston = flag;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::runAnalyticBarlowBeeston() const {
  if (!m_analyticBarlowBeeston)
    return;
  // The best fit values are pushed into the shared bin-by-bin parameters, and
  // the likelihood terms see the change like any other parameter update.
  for (std::size_t i = 0; i < m_nChannels; ++i) {
    if (isMasked(i))
      continue;
    for (RooAbsReal const *hist : m_channelHists[i]) {
      cmsHistVisit(*hist, [](auto const &h) { h.runBarlowBeeston(); });
    }
  }
}

double cacheutils::ChannelSumNLL::evaluate() const {
  syncCategories();
  runAnalyticBarlowBeeston();
  ROOT::Math::KahanSum<double> sum;
  for (std::size_t i = 0; i < m_terms.size(); ++i) {
    if (isMasked(i))
      continue;
    sum += m_terms[i]->getVal() - m_zeroPoints[i];
    if (i < m_nChannels) {
      for (RooMultiPdf const *multiPdf : m_channelMultiPdfs[i]) {
        sum += multiPdf->getCorrection();
      }
    }
  }
  sum += fastConstraintsValue() - m_fastConstraintZeroPoint;
  sum += m_maskingOffset - m_maskingOffsetZero;
  return sum.Sum();
}

bool cacheutils::ChannelSumNLL::getParameters(const RooArgSet *observables,
                                              RooArgSet &outputSet,
                                              bool stripDisconnected) const {
  if (m_internalMasks.empty()) {
    // Collect per term and map the category copies to the originals before
    // merging: the copies of different channels share their names.
    for (auto const &term : m_terms) {
      RooArgSet termParams;
      term->getParameters(observables, termParams, stripDisconnected);
      replaceClonesByOriginals(termParams);
      outputSet.add(termParams, /*silent=*/true);
    }
    for (RooAbsPdf const *pdf : m_fastConstraintPdfs) {
      RooArgSet pdfParams;
      pdf->getParameters(observables, pdfParams, stripDisconnected);
      outputSet.add(pdfParams, /*silent=*/true);
    }
  } else {
    // Only the parameters of active channels: those of masked channels and of the
    // constraint terms stay out, which shrinks the discrete profiling fits.
    outputSet.add(m_activeParameters, /*silent=*/true);
    if (observables) {
      outputSet.remove(*observables, /*silent=*/true, /*matchByNameOnly=*/true);
    }
  }
  outputSet.remove(m_gradientTermVars, /*silent=*/true);
  // Same ordering convention as RooAbsArg::getParameters().
  outputSet.sort();
  return false;
}

void cacheutils::ChannelSumNLL::applyWeightSquared(bool flag) {
  for (std::size_t i = 0; i < m_nChannels; ++i) {
    m_terms[i]->applyWeightSquared(flag);
  }
}

void cacheutils::ChannelSumNLL::updateGradientMaps() const {
  // Each term's gradient follows its own unstripped parameter list (see
  // RooMinimizerFcn); map those positions into our list.
  RooArgSet params;
  getParameters(/*observables=*/nullptr, params, /*stripDisconnected=*/false);
  m_nGradientParams = params.size();
  std::unordered_map<std::string, int> indexByName;
  int idx = 0;
  for (RooAbsArg const *arg : params) {
    indexByName[arg->GetName()] = idx++;
  }

  m_gradientMaps.clear();
  for (auto const &term : m_terms) {
    RooArgSet termParams;
    term->getParameters(/*observables=*/nullptr, termParams, /*stripDisconnected=*/false);
    std::vector<int> map;
    map.reserve(termParams.size());
    for (RooAbsArg const *arg : termParams) {
      auto found = indexByName.find(arg->GetName());
      map.push_back(found == indexByName.end() ? -1 : found->second);
    }
    m_gradientMaps.push_back(std::move(map));
  }
  auto indexOf = [&indexByName](RooAbsArg const *arg) -> int {
    if (!arg)
      return -1;
    auto found = indexByName.find(arg->GetName());
    return found == indexByName.end() ? -1 : found->second;
  };
  m_gaussianGradientIndices.clear();
  for (GaussianConstraint const &c : m_gaussianConstraints) {
    m_gaussianGradientIndices.push_back(
        {indexOf(c.x), indexOf(c.mean), indexOf(dynamic_cast<RooRealVar const *>(c.sigma))});
  }
  m_poissonGradientIndices.clear();
  for (PoissonConstraint const &c : m_poissonConstraints) {
    m_poissonGradientIndices.push_back(indexOf(c.mean));
  }
  m_gradientMapsDirty = false;
}

bool cacheutils::ChannelSumNLL::addFastConstraint(RooAbsPdf const &pdf) {
  if (typeid(pdf) == typeid(SimpleGaussianConstraint)) {
    // Self-normalized, so RooFit doesn't normalize it either.
    auto const &gauss = dynamic_cast<SimpleGaussianConstraint const &>(pdf);
    auto const *x = dynamic_cast<RooRealVar const *>(&gauss.getX());
    auto const *mean = dynamic_cast<RooRealVar const *>(&gauss.getMean());
    if (!x || !mean)
      return false;
    m_gaussianConstraints.push_back({x, mean, &gauss.getSigma()});
  } else if (typeid(pdf) == typeid(RooPoisson) || typeid(pdf) == typeid(SimplePoissonConstraint)) {
    auto const &pois = dynamic_cast<RooPoisson const &>(pdf);
    auto const *mean = dynamic_cast<RooRealVar const *>(&pois.getMean());
    RooAbsReal const &count = pois.getX();
    bool const countIsVar = dynamic_cast<RooRealVar const *>(&count) || dynamic_cast<RooConstVar const *>(&count);
    if (!mean || !countIsVar || !count.isConstant())
      return false;
    // RooFit normalizes a constraint over the constrained parameters among its
    // variables. A count in that set would mean a numeric normalization: leave
    // that to RooFit.
    if (m_hasConstrain && m_constrain.find(count))
      return false;
    bool const normalized = m_hasConstrain && m_constrain.find(*mean);
    m_poissonConstraints.push_back({&pois, &count, mean, normalized});
  } else {
    return false;
  }
  m_fastConstraintPdfs.push_back(&pdf);
  m_components.add(pdf);
  return true;
}

double cacheutils::ChannelSumNLL::fastConstraintsValue() const {
  ROOT::Math::KahanSum<double> sum;
  for (GaussianConstraint const &c : m_gaussianConstraints) {
    double const pull = (c.x->getVal() - c.mean->getVal()) / c.sigma->getVal();
    sum += 0.5 * pull * pull;
  }
  for (PoissonConstraint const &c : m_poissonConstraints) {
    sum += poissonConstraintValue(c);
  }
  return sum.Sum();
}

/// Minus the logarithm of the RooPoisson value, computed like RooFit does but
/// without the exp/log round trip.
double cacheutils::ChannelSumNLL::poissonConstraintValue(PoissonConstraint const &c) {
  RooPoisson const &pois = *c.pdf;
  double const k = pois.getNoRounding() ? c.count->getVal() : std::floor(c.count->getVal());
  double const mu = c.mean->getVal();
  double value = 0.;
  if (mu < 0.) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)
    bool const protectNegative = pois.getProtectNegativeMean();
#else
    bool const protectNegative = true;  // the RooPoisson default, the getter only exists since ROOT 6.36
#endif
    value = protectNegative ? -std::log(1e-3) : std::numeric_limits<double>::quiet_NaN();
  } else if (k < 0.) {
    value = std::numeric_limits<double>::infinity();
  } else if (k == 0.) {
    value = mu;
  } else {
    value = mu + TMath::LnGamma(k + 1.) - (k * std::log(mu));
  }
  if (c.normalized) {
    if (mu < 0.) {
      // Depends on the mean in this case (RooFit's protection against negative means).
      return value + std::log(pois.analyticalIntegral(/*code=*/2));
    }
    if (!c.normValid || c.normCount != k || c.normMin != c.mean->getMin() || c.normMax != c.mean->getMax()) {
      c.normCount = k;
      c.normMin = c.mean->getMin();
      c.normMax = c.mean->getMax();
      c.logNorm = std::log(pois.analyticalIntegral(/*code=*/2));
      c.normValid = true;
    }
    value += c.logNorm;
  }
  return value;
}

void cacheutils::ChannelSumNLL::fastConstraintsGradient(double *out) const {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (std::size_t i = 0; i < m_gaussianConstraints.size(); ++i) {
    GaussianConstraint const &c = m_gaussianConstraints[i];
    std::array<int, 3> const &idx = m_gaussianGradientIndices[i];
    double const sigma = c.sigma->getVal();
    double const pull = (c.x->getVal() - c.mean->getVal()) / sigma;
    if (idx[0] >= 0)
      out[idx[0]] += pull / sigma;
    if (idx[1] >= 0)
      out[idx[1]] -= pull / sigma;
    if (idx[2] >= 0)
      out[idx[2]] -= pull * pull / sigma;
  }
  for (std::size_t i = 0; i < m_poissonConstraints.size(); ++i) {
    int const idx = m_poissonGradientIndices[i];
    PoissonConstraint const &c = m_poissonConstraints[i];
    double const mu = c.mean->getVal();
    if (idx < 0 || mu <= 0.)
      continue;
    double const k = c.pdf->getNoRounding() ? c.count->getVal() : std::floor(c.count->getVal());
    out[idx] += 1. - (k / mu);
  }
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

void cacheutils::ChannelSumNLL::gradient(double *out) const {
  if (!m_hasGradient) {
    throw std::logic_error("ChannelSumNLL::gradient() called, but not all likelihood terms provide a gradient");
  }
  if (m_gradientMapsDirty)
    updateGradientMaps();
  syncCategories();
  // The Barlow-Beeston parameters are constant at their conditional minimum,
  // where the derivative of the likelihood with respect to them vanishes, so
  // the gradient with respect to the other parameters is the full one. Where
  // the minimum is not a stationary point but a kink, the gradient terms of
  // the channel pdfs add the missing piece (CMSHistBarlowBeestonGradientTerm).
  runAnalyticBarlowBeeston();
  if (m_analyticBarlowBeeston) {
    for (std::size_t i = 0; i < m_nChannels; ++i) {
      if (isMasked(i))
        continue;
      for (CMSHistBarlowBeestonGradientTerm *term : m_channelGradientTerms[i]) {
        term->setFromHistogram();
      }
    }
  }
  std::fill_n(out, m_nGradientParams, 0.0);

  for (std::size_t i = 0; i < m_terms.size(); ++i) {
    if (isMasked(i))
      continue;
    std::vector<int> const &map = m_gradientMaps[i];
    m_gradientBuffer.assign(map.size(), 0.0);
    m_terms[i]->gradient(m_gradientBuffer.data());
    for (std::size_t j = 0; j < map.size(); ++j) {
      if (map[j] >= 0)
        out[map[j]] += m_gradientBuffer[j];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
  }
  fastConstraintsGradient(out);
  // Outside of gradient calls the terms must not contribute.
  for (auto const &channelTerms : m_channelGradientTerms) {
    for (CMSHistBarlowBeestonGradientTerm *term : channelTerms) {
      term->clearSlopes();
    }
  }
}

void cacheutils::ChannelSumNLL::setZeroPoint() {
  syncCategories();
  runAnalyticBarlowBeeston();
  for (std::size_t i = 0; i < m_terms.size(); ++i) {
    m_zeroPoints[i] = m_terms[i]->getVal();
  }
  m_fastConstraintZeroPoint = fastConstraintsValue();
  m_maskingOffsetZero = m_maskingOffset;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::clearZeroPoint() {
  m_zeroPoints.assign(m_zeroPoints.size(), 0.0);
  m_fastConstraintZeroPoint = 0;
  m_maskingOffsetZero = 0;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::setMaskNonDiscreteChannels(bool mask) {
  double nllBefore = evaluate();
  m_internalMasks.clear();
  m_activeParameters.removeAll();
  if (mask) {
    m_internalMasks.resize(m_nChannels, false);
    for (std::size_t i = 0; i < m_nChannels; ++i) {
      for (RooAbsArg *param : m_channelParams[i]) {
        auto *cat = dynamic_cast<RooCategory *>(param);
        if (cat && !cat->isConstant()) {
          m_internalMasks[i] = true;
          m_activeParameters.add(m_channelParams[i], /*silent=*/true);
          CombineLogger::instance().log(
              "ChannelSumNLL.cc",
              __LINE__,
              "Enabling channel " + m_channelNames[i] + " that depends on non-constant category " + cat->GetName(),
              __func__);
          break;
        }
      }
    }
  }
  m_gradientMapsDirty = true;
  setValueDirty();
  double nllAfter = evaluate();
  m_maskingOffset += (nllBefore - nllAfter);
  setValueDirty();
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
