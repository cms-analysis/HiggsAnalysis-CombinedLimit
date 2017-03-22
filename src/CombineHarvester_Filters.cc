#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include <vector>
#include <string>
#include <utility>
#include <set>
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"

namespace ch {

CombineHarvester& CombineHarvester::bin(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::bin), cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::bin), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::bin), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::bin), cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::bin), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::bin), cond);  
  }
  return *this;
}

CombineHarvester& CombineHarvester::bin_id(
    std::vector<int> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::bin_id), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::bin_id), cond);
  FilterContaining(systs_, vec, std::mem_fn(&Systematic::bin_id), cond);
  return *this;
}

CombineHarvester& CombineHarvester::process(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::process), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::process), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::process), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::process), cond);
  }

  return *this;
}

CombineHarvester& CombineHarvester::process_rgx(
    std::vector<std::string> const& vec, bool cond) {
  FilterContainingRgx(procs_, vec, std::mem_fn(&Process::process), cond);
  FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::process), cond);
  return *this;
}

CombineHarvester& CombineHarvester::analysis(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::analysis), cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::analysis), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::analysis), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::analysis), cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::analysis), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::analysis), cond);
  }
  return *this;
}

CombineHarvester& CombineHarvester::era(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::era), cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::era), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::era), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::era), cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::era), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::era), cond);
  }
  return *this;
}

CombineHarvester& CombineHarvester::channel(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::channel), cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::channel), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::channel), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::channel), cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::channel), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::channel), cond);
  }
  return *this;
}

CombineHarvester& CombineHarvester::mass(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::mass), cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::mass), cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::mass), cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::mass), cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::mass), cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::mass), cond);
  }
  return *this;
}

CombineHarvester& CombineHarvester::attr(
    std::vector<std::string> const& vec, std::string attr_label, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(procs_, vec, std::mem_fn(&Process::attribute), attr_label, cond);
    FilterContainingRgx(obs_, vec, std::mem_fn(&Observation::attribute), attr_label, cond);
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::attribute), attr_label, cond);
  } else {
    FilterContaining(procs_, vec, std::mem_fn(&Process::attribute), attr_label, cond);
    FilterContaining(obs_, vec, std::mem_fn(&Observation::attribute), attr_label, cond);
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::attribute), attr_label, cond);  
  }
  return *this;
}


CombineHarvester& CombineHarvester::syst_name(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::name), cond);
  } else {
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::name), cond);
  }
  return *this;
}

CombineHarvester& CombineHarvester::syst_type(
    std::vector<std::string> const& vec, bool cond) {
  if (GetFlag("filters-use-regex")) {
    FilterContainingRgx(systs_, vec, std::mem_fn(&Systematic::type), cond);
  } else {
    FilterContaining(systs_, vec, std::mem_fn(&Systematic::type), cond);
  }
  return *this;
}

CombineHarvester & CombineHarvester::signals() {
  ch::erase_if(systs_, [&] (std::shared_ptr<Systematic> val) {
    return !val->signal();
  });
  ch::erase_if(procs_, [&] (std::shared_ptr<Process> val) {
    return !val->signal();
  });
  return *this;
}

CombineHarvester & CombineHarvester::backgrounds() {
  ch::erase_if(systs_, [&] (std::shared_ptr<Systematic> val) {
    return val->signal();
  });
  ch::erase_if(procs_, [&] (std::shared_ptr<Process> val) {
    return val->signal();
  });
  return *this;
}

CombineHarvester & CombineHarvester::histograms() {
  ch::erase_if(obs_, [&] (std::shared_ptr<Observation> val) {
    return val->shape() == nullptr;
  });
  ch::erase_if(procs_, [&] (std::shared_ptr<Process> val) {
    return val->shape() == nullptr;
  });
  return *this;
}

CombineHarvester & CombineHarvester::pdfs() {
  ch::erase_if(procs_, [&] (std::shared_ptr<Process> val) {
    return val->pdf() == nullptr;
  });
  return *this;
}

CombineHarvester & CombineHarvester::data() {
  ch::erase_if(obs_, [&] (std::shared_ptr<Observation> val) {
    return val->data() == nullptr;
  });
  return *this;
}

std::set<std::string> CombineHarvester::bin_set() {
  std::set<std::string> result =
      this->SetFromObs(std::mem_fn(&ch::Observation::bin));
  std::set<std::string> result2 =
      this->SetFromProcs(std::mem_fn(&ch::Process::bin));
  std::set<std::string> result3 =
      this->SetFromSysts(std::mem_fn(&ch::Systematic::bin));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<int> CombineHarvester::bin_id_set() {
  std::set<int> result =
      this->SetFromObs(std::mem_fn(&ch::Observation::bin_id));
  std::set<int> result2 =
      this->SetFromProcs(std::mem_fn(&ch::Process::bin_id));
  std::set<int> result3 =
      this->SetFromSysts(std::mem_fn(&ch::Systematic::bin_id));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::process_set() {
  std::set<std::string> result = this->SetFromProcs(
      std::mem_fn(&ch::Process::process));
  std::set<std::string> result2 = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::process));
  result.insert(result2.begin(), result2.end());
  return result;
}

std::set<std::string> CombineHarvester::analysis_set() {
  std::set<std::string> result = this->SetFromObs(
      std::mem_fn(&ch::Observation::analysis));
  std::set<std::string> result2 = this->SetFromProcs(
      std::mem_fn(&ch::Process::analysis));
  std::set<std::string> result3 = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::analysis));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::era_set() {
  std::set<std::string> result =
      this->SetFromObs(std::mem_fn(&ch::Observation::era));
  std::set<std::string> result2 =
      this->SetFromProcs(std::mem_fn(&ch::Process::era));
  std::set<std::string> result3 =
      this->SetFromSysts(std::mem_fn(&ch::Systematic::era));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::channel_set() {
  std::set<std::string> result = this->SetFromObs(
      std::mem_fn(&ch::Observation::channel));
  std::set<std::string> result2 = this->SetFromProcs(
      std::mem_fn(&ch::Process::channel));
  std::set<std::string> result3 = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::channel));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::mass_set() {
  std::set<std::string> result = this->SetFromObs(
      std::mem_fn(&ch::Observation::mass));
  std::set<std::string> result2 = this->SetFromProcs(
      std::mem_fn(&ch::Process::mass));
  std::set<std::string> result3 = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::mass));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::syst_name_set() {
  std::set<std::string> result = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::name));
  return result;
}

std::set<std::string> CombineHarvester::syst_type_set() {
  std::set<std::string> result = this->SetFromSysts(
      std::mem_fn(&ch::Systematic::type));
  return result;
}
}
