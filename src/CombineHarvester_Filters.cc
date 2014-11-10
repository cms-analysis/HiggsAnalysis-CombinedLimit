#include "CombineTools/interface/CombineHarvester.h"
#include <vector>
#include <string>
#include <utility>
#include <set>
#include "CombineTools/interface/Observation.h"
#include "CombineTools/interface/Process.h"
#include "CombineTools/interface/Nuisance.h"
#include "CombineTools/interface/Algorithm.h"

namespace ch {

CombineHarvester& CombineHarvester::bin(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::bin), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::bin), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::bin), cond);
  return *this;
}

CombineHarvester& CombineHarvester::bin_id(
    std::vector<int> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::bin_id), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::bin_id), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::bin_id), cond);
  return *this;
}

CombineHarvester& CombineHarvester::process(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::process), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::process), cond);
  return *this;
}

CombineHarvester& CombineHarvester::process_rgx(
    std::vector<std::string> const& vec, bool cond) {
  FilterContainingRgx(procs_, vec, std::mem_fn(&Process::process), cond);
  FilterContainingRgx(nus_, vec, std::mem_fn(&Nuisance::process), cond);
  return *this;
}

CombineHarvester& CombineHarvester::analysis(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::analysis), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::analysis), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::analysis), cond);
  return *this;
}

CombineHarvester& CombineHarvester::era(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::era), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::era), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::era), cond);
  return *this;
}

CombineHarvester& CombineHarvester::channel(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::channel), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::channel), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::channel), cond);
  return *this;
}

CombineHarvester& CombineHarvester::mass(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(procs_, vec, std::mem_fn(&Process::mass), cond);
  FilterContaining(obs_, vec, std::mem_fn(&Observation::mass), cond);
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::mass), cond);
  return *this;
}

CombineHarvester& CombineHarvester::nus_name(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::name), cond);
  return *this;
}

CombineHarvester& CombineHarvester::nus_type(
    std::vector<std::string> const& vec, bool cond) {
  FilterContaining(nus_, vec, std::mem_fn(&Nuisance::type), cond);
  return *this;
}

CombineHarvester & CombineHarvester::signals() {
  ch::erase_if(nus_, [&] (std::shared_ptr<Nuisance> val) {
    return !val->signal();
  });
  ch::erase_if(procs_, [&] (std::shared_ptr<Process> val) {
    return !val->signal();
  });
  return *this;
}

CombineHarvester & CombineHarvester::backgrounds() {
  ch::erase_if(nus_, [&] (std::shared_ptr<Nuisance> val) {
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
      this->GenerateSetFromObs<std::string>(std::mem_fn(&ch::Observation::bin));
  std::set<std::string> result2 =
      this->GenerateSetFromProcs<std::string>(std::mem_fn(&ch::Process::bin));
  std::set<std::string> result3 =
      this->GenerateSetFromNus<std::string>(std::mem_fn(&ch::Nuisance::bin));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<int> CombineHarvester::bin_id_set() {
  std::set<int> result =
      this->GenerateSetFromObs<int>(std::mem_fn(&ch::Observation::bin_id));
  std::set<int> result2 =
      this->GenerateSetFromProcs<int>(std::mem_fn(&ch::Process::bin_id));
  std::set<int> result3 =
      this->GenerateSetFromNus<int>(std::mem_fn(&ch::Nuisance::bin_id));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::process_set() {
  std::set<std::string> result = this->GenerateSetFromProcs<std::string>(
      std::mem_fn(&ch::Process::process));
  std::set<std::string> result2 = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::process));
  result.insert(result2.begin(), result2.end());
  return result;
}

std::set<std::string> CombineHarvester::analysis_set() {
  std::set<std::string> result = this->GenerateSetFromObs<std::string>(
      std::mem_fn(&ch::Observation::analysis));
  std::set<std::string> result2 = this->GenerateSetFromProcs<std::string>(
      std::mem_fn(&ch::Process::analysis));
  std::set<std::string> result3 = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::analysis));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::era_set() {
  std::set<std::string> result =
      this->GenerateSetFromObs<std::string>(std::mem_fn(&ch::Observation::era));
  std::set<std::string> result2 =
      this->GenerateSetFromProcs<std::string>(std::mem_fn(&ch::Process::era));
  std::set<std::string> result3 =
      this->GenerateSetFromNus<std::string>(std::mem_fn(&ch::Nuisance::era));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::channel_set() {
  std::set<std::string> result = this->GenerateSetFromObs<std::string>(
      std::mem_fn(&ch::Observation::channel));
  std::set<std::string> result2 = this->GenerateSetFromProcs<std::string>(
      std::mem_fn(&ch::Process::channel));
  std::set<std::string> result3 = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::channel));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::mass_set() {
  std::set<std::string> result = this->GenerateSetFromObs<std::string>(
      std::mem_fn(&ch::Observation::mass));
  std::set<std::string> result2 = this->GenerateSetFromProcs<std::string>(
      std::mem_fn(&ch::Process::mass));
  std::set<std::string> result3 = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::mass));
  result.insert(result2.begin(), result2.end());
  result.insert(result3.begin(), result3.end());
  return result;
}

std::set<std::string> CombineHarvester::nus_name_set() {
  std::set<std::string> result = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::name));
  return result;
}

std::set<std::string> CombineHarvester::nus_type_set() {
  std::set<std::string> result = this->GenerateSetFromNus<std::string>(
      std::mem_fn(&ch::Nuisance::type));
  return result;
}
}
