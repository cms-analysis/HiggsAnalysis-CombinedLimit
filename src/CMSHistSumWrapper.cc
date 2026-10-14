#include "../interface/CMSHistSumWrapper.h"

#include <cassert>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>

#include "../interface/CMSHistSum.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RtypesCore.h"
#include "TString.h"

CMSHistSumWrapper::CMSHistSumWrapper(const char* name, const char* title, RooAbsReal& x, CMSHistSum& sum, unsigned idx)
    : RooAbsReal(name, title),
      m_x(/*theName=*/"x", /*desc=*/"", this, x),
      m_sum(/*theName=*/"sum", /*desc=*/"", this, sum),
      m_cache(sum.cache()),
      m_idx(idx),
      m_sentry(TString(name) + "_sentry", /*title=*/"") {
  m_cache.Clear();
}

CMSHistSumWrapper::CMSHistSumWrapper(CMSHistSumWrapper const& other, const char* name)
    : RooAbsReal(other, name),
      m_x(/*theName=*/"x", this, other.m_x),
      m_sum(/*theName=*/"sum", this, other.m_sum),
      m_cache(other.m_cache),
      m_idx(other.m_idx),
      m_sentry(name ? TString(name) + "_sentry" : TString(other.GetName()) + "_sentry", /*title=*/"") {}

void CMSHistSumWrapper::initialize() const {
  if (m_initialized)
    return;
  m_sentry.SetName(TString(this->GetName()) + "_sentry");
  m_psum = dynamic_cast<CMSHistSum const*>(&(m_sum.arg()));
  if (!m_psum) {
    throw std::runtime_error(std::string("CMSHistSumWrapper ") + GetName() +
                             ": the wrapped function is not a CMSHistSum");
  }
  auto sentry_args = m_psum->getSentryArgs();
  for (RooAbsArg* arg : *sentry_args) {
    m_sentry.addArg(*arg);
  }
  m_sentry.setValueDirty();
  m_initialized = true;
}

void CMSHistSumWrapper::updateCache() const {
  initialize();

  // The CMSHistSum sentries will send a dirty flag whenever we need to update the cache
  if (!m_sentry.good()) {
    m_psum->fillProcessShape(m_idx, m_cache);
  }
  m_sentry.reset();
}

Double_t CMSHistSumWrapper::evaluate() const {
  updateCache();
  return m_cache.GetAt(m_x);
}

void CMSHistSumWrapper::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  m_cache.Dump();
  std::cout << ">> Sentry: " << m_sentry.good() << "\n";
  m_sentry.Print(/*options=*/"v");
}

Int_t CMSHistSumWrapper::getAnalyticalIntegral(RooArgSet& allVars,
                                               RooArgSet& analVars,
                                               const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, m_x))
    return 1;
  return 0;
}

Double_t CMSHistSumWrapper::analyticalIntegral(Int_t code, const char* /*rangeName*/) const {
  assert(code == 1);
  updateCache();
  return m_cache.IntegralWidth();
}
