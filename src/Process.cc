#include "CombineTools/interface/Process.h"
#include <iostream>
#include <string>
#include "boost/format.hpp"
#include "CombineTools/interface/Logging.h"

namespace ch {

Process::Process()
    : bin_(""),
      rate_(0.0),
      process_(""),
      signal_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_(""),
      shape_(),
      pdf_(nullptr),
      norm_(nullptr) {
  }

Process::~Process() { }

void swap(Process& first, Process& second) {
  using std::swap;
  swap(first.bin_, second.bin_);
  swap(first.rate_, second.rate_);
  swap(first.process_, second.process_);
  swap(first.signal_, second.signal_);
  swap(first.analysis_, second.analysis_);
  swap(first.era_, second.era_);
  swap(first.channel_, second.channel_);
  swap(first.bin_id_, second.bin_id_);
  swap(first.mass_, second.mass_);
  swap(first.shape_, second.shape_);
  swap(first.pdf_, second.pdf_);
  swap(first.norm_, second.norm_);
}

Process::Process(Process const& other)
    : bin_(other.bin_),
      rate_(other.rate_),
      process_(other.process_),
      signal_(other.signal_),
      analysis_(other.analysis_),
      era_(other.era_),
      channel_(other.channel_),
      bin_id_(other.bin_id_),
      mass_(other.mass_),
      pdf_(other.pdf_),
      norm_(other.norm_) {
  TH1 *h = nullptr;
  if (other.shape_) {
    h = static_cast<TH1*>(other.shape_->Clone());
    h->SetDirectory(0);
  }
  shape_ = std::unique_ptr<TH1>(h);
}

Process::Process(Process&& other)
    : bin_(""),
      rate_(0.0),
      process_(""),
      signal_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_(""),
      shape_(),
      pdf_(nullptr),
      norm_(nullptr) {
  swap(*this, other);
}

Process& Process::operator=(Process other) {
  swap(*this, other);
  return (*this);
}

void Process::set_shape(std::unique_ptr<TH1> shape, bool set_rate) {
  // We were given a nullptr - this is fine, and so we're done
  if (!shape) {
    // This will safely release any existing TH1 held by shape_
    shape_ = nullptr;
    return;
  }
  // Otherwise we validate this new hist. First is to check that all bins have
  // +ve values, otherwise the interpretation as a pdf is tricky
  for (int i = 1; i <= shape->GetNbinsX(); ++i) {
    if (shape->GetBinContent(i) < 0.) {
      throw std::runtime_error(FNERROR("TH1 has a bin with content < 0"));
    }
  }
  // At this point we can safely move the shape in and take ownership
  shape_ = std::move(shape);
  // Ensure that root will not try and clean this up
  shape_->SetDirectory(0);
  if (set_rate) {
    this->set_rate(shape_->Integral());
  }
  if (shape_->Integral() > 0.) shape_->Scale(1. / shape_->Integral());
}

std::unique_ptr<TH1> Process::ClonedShape() const {
  if (!shape_) return std::unique_ptr<TH1>();
  std::unique_ptr<TH1> res(static_cast<TH1 *>(shape_->Clone()));
  res->SetDirectory(0);
  return res;
}

std::unique_ptr<TH1> Process::ClonedScaledShape() const {
  if (!shape_) return std::unique_ptr<TH1>();
  std::unique_ptr<TH1> res(ClonedShape());
  res->Scale(this->no_norm_rate());
  return res;
}

TH1F Process::ShapeAsTH1F() const {
  if (!shape_) {
    throw std::runtime_error(
        FNERROR("Process object does not contain a shape"));
  }
  TH1F res;
  // Need to get the shape as a concrete type (TH1F or TH1D)
  // A nice way to do this is just to use TH1D::Copy into a fresh TH1F
  TH1F const* test_f = dynamic_cast<TH1F const*>(this->shape());
  TH1D const* test_d = dynamic_cast<TH1D const*>(this->shape());
  if (test_f) {
    test_f->Copy(res);
  } else if (test_d) {
    test_d->Copy(res);
  } else {
    throw std::runtime_error(FNERROR("TH1 shape is not a TH1F or a TH1D"));
  }
  return res;
}


std::ostream& Process::PrintHeader(std::ostream& out) {
  std::string line =
      (boost::format(
           "%-6s %-9s %-6s %-8s %-28s %-3i %-22s %-4i %-10.5g %-10i") %
       "mass" % "analysis" % "era" % "channel" % "bin" % "id" % "process" %
       "sig" % "rate" % "shape/pdf").str();
  std::string div(line.length(), '-');
  out << div  << std::endl;
  out << line << std::endl;
  out << div  << std::endl;
  return out;
}

std::ostream& operator<< (std::ostream &out, Process &val) {
  out << boost::format(
             "%-6s %-9s %-6s %-8s %-28s %-3i %-22s %-4i %-10.5g %-10i") %
             val.mass() % val.analysis() % val.era() % val.channel() %
             val.bin() % val.bin_id() % val.process() % val.signal() %
             val.rate() % (val.shape() != nullptr || val.pdf() != nullptr);
  return out;
}
}
