#include "CombineTools/interface/Systematic.h"
#include <iostream>
#include "boost/format.hpp"
#include "CombineTools/interface/Logging.h"

namespace ch {

Systematic::Systematic()
    : Object(),
      signal_(false),
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      shape_u_(),
      shape_d_(),
      data_u_(nullptr),
      data_d_(nullptr) {
  }

Systematic::~Systematic() { }

void swap(Systematic& first, Systematic& second) {
  using std::swap;
  swap(static_cast<Object&>(first), static_cast<Object&>(second));
  swap(first.signal_, second.signal_);
  swap(first.name_, second.name_);
  swap(first.type_, second.type_);
  swap(first.value_u_, second.value_u_);
  swap(first.value_d_, second.value_d_);
  swap(first.scale_, second.scale_);
  swap(first.asymm_, second.asymm_);
  swap(first.shape_u_, second.shape_u_);
  swap(first.shape_d_, second.shape_d_);
  swap(first.data_u_, second.data_u_);
  swap(first.data_d_, second.data_d_);
}

Systematic::Systematic(Systematic const& other)
    : Object(other),
      signal_(other.signal_),
      name_(other.name_),
      type_(other.type_),
      value_u_(other.value_u_),
      value_d_(other.value_d_),
      scale_(other.scale_),
      asymm_(other.asymm_),
      data_u_(other.data_u_),
      data_d_(other.data_d_) {
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

Systematic::Systematic(Systematic&& other)
    : Object(),
      signal_(false),
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      shape_u_(),
      shape_d_(),
      data_u_(nullptr),
      data_d_(nullptr) {
  swap(*this, other);
}

Systematic& Systematic::operator=(Systematic other) {
  swap(*this, other);
  return (*this);
}

void Systematic::set_shapes(std::unique_ptr<TH1> shape_u,
                            std::unique_ptr<TH1> shape_d, TH1 const* nominal) {
  // Check that the inputs make sense
  if (bool(shape_u) != bool(shape_d)) {
    throw std::runtime_error(
        "shape_u and shape_d must be either both valid or both null");
  }
  if (!shape_u && !shape_d) {
    shape_u_ = nullptr;
    shape_d_ = nullptr;
    return;
  }

  // for (int i = 1; i <= shape_u->GetNbinsX(); ++i) {
  //   if (shape_u->GetBinContent(i) < 0.) {
  //     throw std::runtime_error(
  //         FNERROR("shape_u TH1 has a bin with content < 0"));
  //   }
  // }
  // for (int i = 1; i <= shape_d->GetNbinsX(); ++i) {
  //   if (shape_d->GetBinContent(i) < 0.) {
  //     throw std::runtime_error(
  //         FNERROR("shape_d TH1 has a bin with content < 0"));
  //   }
  // }

  // if (nominal) {
  //   for (int i = 1; i <= nominal->GetNbinsX(); ++i) {
  //     if (nominal->GetBinContent(i) < 0.) {
  //       throw std::runtime_error(
  //           FNERROR("nominal TH1 has a bin with content < 0"));
  //     }
  //   }
  // }

  shape_u_ = std::move(shape_u);
  shape_d_ = std::move(shape_d);
  shape_u_->SetDirectory(0);
  shape_d_->SetDirectory(0);

  if (nominal && nominal->Integral() > 0.) {
    this->set_value_u(shape_u_->Integral() / nominal->Integral());
    this->set_value_d(shape_d_->Integral() / nominal->Integral());
  }

  if (shape_u_->Integral() > 0.) shape_u_->Scale(1. / shape_u_->Integral());
  if (shape_d_->Integral() > 0.) shape_d_->Scale(1. / shape_d_->Integral());
}

void Systematic::set_data(RooDataHist* data_u, RooDataHist* data_d,
                          RooDataHist const* nominal) {
  if (nominal && nominal->sumEntries() > 0.) {
    this->set_value_u(data_u->sumEntries() / nominal->sumEntries());
    this->set_value_d(data_d->sumEntries() / nominal->sumEntries());
  }
  data_u_ = data_u;
  data_d_ = data_d;
}


std::unique_ptr<TH1> Systematic::ClonedShapeU() const {
  if (!shape_u_) return std::unique_ptr<TH1>();
  std::unique_ptr<TH1> res(static_cast<TH1 *>(shape_u_->Clone()));
  res->SetDirectory(0);
  return res;
}

std::unique_ptr<TH1> Systematic::ClonedShapeD() const {
  if (!shape_d_) return std::unique_ptr<TH1>();
  std::unique_ptr<TH1> res(static_cast<TH1 *>(shape_d_->Clone()));
  res->SetDirectory(0);
  return res;
}

std::ostream& Systematic::PrintHeader(std::ostream &out) {
  std::string line =
   (boost::format("%-6s %-9s %-6s %-8s %-28s %-3i"
    " %-16s %-4i %-45s %-8s %-13s %-4i %-4i")
    % "mass" % "analysis" % "era" % "channel" % "bin" % "id" % "process" % "sig"
    % "nuisance" % "type" % "value" % "sh_d" % "sh_u").str();
  std::string div(line.length(), '-');
  out << div  << "\n";
  out << line << "\n";
  out << div  << "\n";
  return out;
}

std::ostream& operator<< (std::ostream &out, Systematic const& val) {
  std::string value_fmt;
  if (val.asymm()) {
    value_fmt = (boost::format("%.4g/%.4g")
      % val.value_d() % val.value_u()).str();
  } else {
    value_fmt = (boost::format("%.4g") % val.value_u()).str();
  }
  out << boost::format("%-6s %-9s %-6s %-8s %-28s %-3i"
    " %-16s %-4i %-45s %-8s %-13s %-4i %-4i")
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
  % (bool(val.shape_d()) || bool(val.data_d()))
  % (bool(val.shape_u()) || bool(val.data_u()));
  return out;
}
}
