#include "CombineHarvester/CombineTools/interface/Systematic.h"
#include <iostream>
#include "boost/format.hpp"
#include "CombineHarvester/CombineTools/interface/Logging.h"
#include <regex>

namespace {
auto format_syst(const ch::Systematic& val) {
  std::string value_fmt;
  if (val.asymm()) {
    value_fmt = (boost::format("%.4g/%.4g")
      % val.value_d() % val.value_u()).str();
  } else {
    value_fmt = (boost::format("%.4g") % val.value_u()).str();
  }
  return boost::format("%-6s %-9s %-6s %-8s %-28s %-3i"
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
  % (bool(val.shape_d()) || bool(val.data_d()) || bool(val.pdf_d()))
  % (bool(val.shape_u()) || bool(val.data_u()) || bool(val.pdf_u()));
}
}

namespace ch {

Systematic::Systematic()
    : Object(),
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      shape_u_(),
      shape_d_(),
      pdf_u_(nullptr),
      pdf_d_(nullptr),
      data_u_(nullptr),
      data_d_(nullptr) {
  }

Systematic::~Systematic() { }

void Systematic::set_name(std::string const& name) { 
//test = std::regex_replace(test, std::regex("def"), "klm");
    if (data_u_) data_u_->SetName(std::regex_replace(data_u_->GetName(),std::regex(name_),name).c_str());
    if (data_d_) data_d_->SetName(std::regex_replace(data_d_->GetName(),std::regex(name_),name).c_str());
    if (pdf_u_) pdf_u_->SetName(std::regex_replace(pdf_u_->GetName(),std::regex(name_),name).c_str());
    if (pdf_d_) pdf_d_->SetName(std::regex_replace(pdf_d_->GetName(),std::regex(name_),name).c_str());
    name_ = name; 
}

void swap(Systematic& first, Systematic& second) {
  using std::swap;
  swap(static_cast<Object&>(first), static_cast<Object&>(second));
  swap(first.name_, second.name_);
  swap(first.type_, second.type_);
  swap(first.value_u_, second.value_u_);
  swap(first.value_d_, second.value_d_);
  swap(first.scale_, second.scale_);
  swap(first.asymm_, second.asymm_);
  swap(first.shape_u_, second.shape_u_);
  swap(first.shape_d_, second.shape_d_);
  swap(first.pdf_u_, second.pdf_u_);
  swap(first.pdf_d_, second.pdf_d_);
  swap(first.data_u_, second.data_u_);
  swap(first.data_d_, second.data_d_);
}

Systematic::Systematic(Systematic const& other)
    : Object(other),
      name_(other.name_),
      type_(other.type_),
      value_u_(other.value_u_),
      value_d_(other.value_d_),
      scale_(other.scale_),
      asymm_(other.asymm_),
      pdf_u_(other.pdf_u_),
      pdf_d_(other.pdf_d_),
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
      name_(""),
      type_(""),
      value_u_(0.0),
      value_d_(0.0),
      scale_(1.0),
      asymm_(false),
      shape_u_(),
      shape_d_(),
      pdf_u_(nullptr),
      pdf_d_(nullptr),
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

void Systematic::set_shapes(TH1 const& shape_u, TH1 const& shape_d,
                            TH1 const& nominal) {
  set_shapes(std::unique_ptr<TH1>(static_cast<TH1*>(shape_u.Clone())),
             std::unique_ptr<TH1>(static_cast<TH1*>(shape_d.Clone())),
             &nominal);
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

void Systematic::set_pdf(RooAbsReal* pdf_u, RooAbsReal* pdf_d,
                          RooAbsReal const* nominal) {
  pdf_u_ = pdf_u;
  pdf_d_ = pdf_d;
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

TH1F Systematic::ShapeUAsTH1F() const {
  TH1F res;
  if (this->shape_u()) {
    // Need to get the shape as a concrete type (TH1F or TH1D)
    // A nice way to do this is just to use TH1D::Copy into a fresh TH1F
    TH1F const* test_f = dynamic_cast<TH1F const*>(this->shape_u());
    TH1D const* test_d = dynamic_cast<TH1D const*>(this->shape_u());
    if (test_f) {
      test_f->Copy(res);
    } else if (test_d) {
      test_d->Copy(res);
    } else {
      throw std::runtime_error(FNERROR("TH1 shape is not a TH1F or a TH1D"));
    }
  }  
 return res;
}

TH1F Systematic::ShapeDAsTH1F() const {
  TH1F res;
  if (this->shape_d()) {
    // Need to get the shape as a concrete type (TH1F or TH1D)
    // A nice way to do this is just to use TH1D::Copy into a fresh TH1F
    TH1F const* test_f = dynamic_cast<TH1F const*>(this->shape_d());
    TH1D const* test_d = dynamic_cast<TH1D const*>(this->shape_d());
    if (test_f) {
      test_f->Copy(res);
    } else if (test_d) {
      test_d->Copy(res);
    } else {
      throw std::runtime_error(FNERROR("TH1 shape is not a TH1F or a TH1D"));
    }
  }  
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

std::string Systematic::to_string() const {
  return ::format_syst(*this).str();
}

std::ostream& operator<< (std::ostream &out, Systematic const& val) {
  out << ::format_syst(val);
  return out;
}

void Systematic::SwapUpAndDown() {
  if (not asymm()) 
      value_u_=1./value_u_;
  else{
      double tmp = value_u_;
      value_u_ = value_d_;
      value_d_ = tmp;
  }
  shape_u_.swap(shape_d_);
}
}
