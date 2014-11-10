#include "CombineTools/interface/Nuisance.h"
#include <iostream>
#include "boost/format.hpp"

namespace ch {

Nuisance::Nuisance()
    : bin_(""),
      process_(""),
      signal_(false),
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_(""),
      shape_u_(),
      shape_d_() {
  }

Nuisance::~Nuisance() { }

void swap(Nuisance& first, Nuisance& second) {
  using std::swap;
  swap(first.bin_, second.bin_);
  swap(first.process_, second.process_);
  swap(first.signal_, second.signal_);
  swap(first.name_, second.name_);
  swap(first.type_, second.type_);
  swap(first.value_u_, second.value_u_);
  swap(first.value_d_, second.value_d_);
  swap(first.scale_, second.scale_);
  swap(first.asymm_, second.asymm_);
  swap(first.analysis_, second.analysis_);
  swap(first.era_, second.era_);
  swap(first.channel_, second.channel_);
  swap(first.bin_id_, second.bin_id_);
  swap(first.mass_, second.mass_);
  swap(first.shape_u_, second.shape_u_);
  swap(first.shape_d_, second.shape_d_);
}

Nuisance::Nuisance(Nuisance const& other)
    : bin_(other.bin_),
      process_(other.process_),
      signal_(other.signal_),
      name_(other.name_),
      type_(other.type_),
      value_u_(other.value_u_),
      value_d_(other.value_d_),
      scale_(other.scale_),
      asymm_(other.asymm_),
      analysis_(other.analysis_),
      era_(other.era_),
      channel_(other.channel_),
      bin_id_(other.bin_id_),
      mass_(other.mass_) {
  TH1 *h_u = nullptr;
  if (other.shape_u_) {
    h_u = dynamic_cast<TH1*>(other.shape_u_->Clone());
    h_u->SetDirectory(0);
  }
  shape_u_ = std::unique_ptr<TH1>(h_u);
  TH1 *h_d = nullptr;
  if (other.shape_d_) {
    h_d = dynamic_cast<TH1*>(other.shape_d_->Clone());
    h_d->SetDirectory(0);
  }
  shape_d_ = std::unique_ptr<TH1>(h_d);
}

Nuisance::Nuisance(Nuisance&& other)
    : bin_(""),
      process_(""),
      signal_(false),
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_(""),
      shape_u_(),
      shape_d_() {
  swap(*this, other);
}

Nuisance& Nuisance::operator=(Nuisance other) {
  swap(*this, other);
  return (*this);
}

std::ostream& Nuisance::PrintHeader(std::ostream &out) {
  std::string line =
   (boost::format("%-6s %-9s %-6s %-8s %-28s %-3i"
    " %-22s %-4i %-40s %-8s %-9s %-10i %-10i")
    % "mass" % "analysis" % "era" % "channel" % "bin" % "id" % "process" % "sig"
    % "nuisance" % "type" % "value" % "shape_d" % "shape_u").str();
  std::string div(line.length(), '-');
  out << div  << std::endl;
  out << line << std::endl;
  out << div  << std::endl;
  return out;
}

std::ostream& operator<< (std::ostream &out, Nuisance &val) {
  std::string value_fmt;
  if (val.asymm()) {
    value_fmt = (boost::format("%-4.4g/%-4.4g")
      % val.value_d() % val.value_u()).str();
  } else {
    value_fmt = (boost::format("%-9.4g") % val.value_u()).str();
  }
  out << boost::format("%-6s %-9s %-6s %-8s %-28s %-3i"
    " %-22s %-4i %-40s %-8s %-9s %-10i %-10i")
  % val.mass()
  % val.analysis()
  % val.era()
  % val.channel()
  % val.bin()
  % val.bin_id()
  % val.process()
  % val.signal()
  % val.name()
  % val.type()
  % value_fmt
  % bool(val.shape_d())
  % bool(val.shape_u());
  return out;
}
}
