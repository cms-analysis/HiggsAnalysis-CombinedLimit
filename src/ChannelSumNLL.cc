#include "../interface/ChannelSumNLL.h"

// ROOT_VERSION: ROOT/RVersion.hxx exists since ROOT 6.30, older versions only have RVersion.h.
#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)

#include "../interface/CombineLogger.h"
#include "../interface/RooMultiPdfCombine.h"
#include "../interface/SimNLLBase.h"
#include "../interface/utils.h"

#include <RooAbsCategoryLValue.h>
#include <RooCmdArg.h>
#include <RooGlobalFunc.h>
#include <RooUniform.h>

#include <Math/Util.h>
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 37, 0)
#include <TList.h>
#include <TObject.h>
#endif

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static std::string subName(RooAbsArg const &owner, const char *suffix) {
  return std::string(owner.GetName()) + "_" + suffix;
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

  // One RooFit likelihood per channel; the channel pdfs hold no constraints anymore.
  for (auto const &state : indexCat) {
    std::string const &label = state.first;
    RooAbsPdf *pdf = m_factorizedPdf->getPdf(label.c_str());
    if (!pdf)
      continue;
    RooAbsData *channelData = findDataset(m_datasets, label);
    if (!channelData) {
      throw std::logic_error("ChannelSumNLL: no data for channel " + label);
    }
    std::vector<CategorySync> syncs;
    RooAbsPdf &nllPdf = channelPdfForNLL(*pdf, syncs);
    std::unique_ptr<RooAbsReal> nll{nllPdf.createNLL(*channelData, backendArg)};
    if (!nll) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of channel " + label);
    }
    nll->SetName(subName(*this, label.c_str()).c_str());
    m_channelNames.push_back(label);
    m_categorySyncs.push_back(syncs);
    // The RooMultiPdfs of the original pdf, for the penalty of the original index.
    std::vector<RooMultiPdf const *> multiPdfs;
    std::unique_ptr<RooArgSet> components{pdf->getComponents()};
    for (RooAbsArg const *component : *components) {
      if (auto const *multiPdf = dynamic_cast<RooMultiPdf const *>(component))
        multiPdfs.push_back(multiPdf);
    }
    m_channelMultiPdfs.push_back(multiPdfs);
    RooArgSet params;
    nll->getParameters(/*observables=*/nullptr, params, /*stripDisconnected=*/false);
    replaceClonesByOriginals(params);
    m_channelParams.push_back(params);
    m_components.add(*nll);
    m_terms.push_back(std::move(nll));
  }
  m_nChannels = m_terms.size();

  // The constraint terms get a likelihood of their own, see the class comment.
  // Passing the constraint pdfs explicitly also covers constraints on parameters
  // outside of the "constrain" set (e.g. the "constr" directive).
  if (!m_constraints.empty()) {
    m_constraintDummyObs =
        std::make_unique<RooRealVar>(subName(*this, /*suffix=*/"constraintDummyObs").c_str(), "", 0.5, 0., 1.);
    m_constraintPdf =
        std::make_unique<RooUniform>(subName(*this, /*suffix=*/"constraintPdf").c_str(), "", *m_constraintDummyObs);
    m_constraintData = std::make_unique<RooDataSet>(
        subName(*this, /*suffix=*/"constraintData").c_str(), "", RooArgSet(*m_constraintDummyObs));
    m_constraintData->add(RooArgSet(*m_constraintDummyObs));
    RooCmdArg constrainArg = m_hasConstrain ? RooFit::Constrain(m_constrain) : RooCmdArg::none();
    std::unique_ptr<RooAbsReal> nll{m_constraintPdf->createNLL(
        *m_constraintData, backendArg, RooFit::ExternalConstraints(RooArgSet(m_constraints)), constrainArg)};
    if (!nll) {
      throw std::runtime_error("ChannelSumNLL: failed to create the likelihood of the constraint terms");
    }
    nll->SetName(subName(*this, /*suffix=*/"constraints").c_str());
    m_components.add(*nll);
    m_terms.push_back(std::move(nll));
  }
  m_zeroPoints.assign(m_terms.size(), 0.);

  m_hasGradient = m_nChannels > 0 &&
                  std::all_of(m_terms.begin(), m_terms.end(), [](auto const &term) { return term->hasGradient(); });

  CombineLogger::instance().log("ChannelSumNLL.cc",
                                __LINE__,
                                "Created " + std::string(GetName()) + " with " + std::to_string(m_nChannels) +
                                    " channels and " + std::to_string(m_constraints.size()) +
                                    " constraint terms, RooFit backend \"" + m_evalBackend +
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
  }
  // Only release the old datasets after all terms point to the new ones.
  m_datasets = std::move(newDatasets);
  m_data = &data;
  setValueDirty();
  return true;
}

double cacheutils::ChannelSumNLL::evaluate() const {
  syncCategories();
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
  } else {
    // Only the parameters of active channels: those of masked channels and of the
    // constraint terms stay out, which shrinks the discrete profiling fits.
    outputSet.add(m_activeParameters, /*silent=*/true);
    if (observables) {
      outputSet.remove(*observables, /*silent=*/true, /*matchByNameOnly=*/true);
    }
  }
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
  m_gradientMapsDirty = false;
}

void cacheutils::ChannelSumNLL::gradient(double *out) const {
  if (!m_hasGradient) {
    throw std::logic_error("ChannelSumNLL::gradient() called, but not all likelihood terms provide a gradient");
  }
  if (m_gradientMapsDirty)
    updateGradientMaps();
  syncCategories();
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
}

void cacheutils::ChannelSumNLL::setZeroPoint() {
  for (std::size_t i = 0; i < m_terms.size(); ++i) {
    m_zeroPoints[i] = m_terms[i]->getVal();
  }
  m_maskingOffsetZero = m_maskingOffset;
  setValueDirty();
}

void cacheutils::ChannelSumNLL::clearZeroPoint() {
  m_zeroPoints.assign(m_zeroPoints.size(), 0.0);
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
