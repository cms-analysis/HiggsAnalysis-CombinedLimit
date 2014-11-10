#include "CombineTools/interface/Process.h"
#include <iostream>
#include <string>
#include "boost/format.hpp"

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

void Process::SetNormShape(std::unique_ptr<TH1> shape) {
  if (!shape) return;
  if (shape->Integral() > 0.0) shape->Scale(1.0 / shape->Integral());
  this->set_shape(std::move(shape));
}

void Process::SetNormShapeAndRate(std::unique_ptr<TH1> shape) {
  if (!shape) return;
  this->set_rate(shape->Integral());
  SetNormShape(std::move(shape));
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
