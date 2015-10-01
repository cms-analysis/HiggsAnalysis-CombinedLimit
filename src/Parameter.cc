#include "CombineHarvester/CombineTools/interface/Parameter.h"
#include <iostream>
#include <string>
#include <limits>
#include "boost/format.hpp"

namespace ch {

Parameter::Parameter()
    : name_(""),
      val_(0.0),
      err_u_(1.0),
      err_d_(-1.0),
      range_u_(std::numeric_limits<double>::max()),
      range_d_(std::numeric_limits<double>::lowest()) {
  }

Parameter::~Parameter() { }

void swap(Parameter& first, Parameter& second) {
  using std::swap;
  swap(first.name_, second.name_);
  swap(first.val_, second.val_);
  swap(first.err_u_, second.err_u_);
  swap(first.err_d_, second.err_d_);
  swap(first.range_u_, second.range_u_);
  swap(first.range_d_, second.range_d_);
  swap(first.vars_, second.vars_);
}

Parameter::Parameter(Parameter const& other)
    : name_(other.name_),
      val_(other.val_),
      err_u_(other.err_u_),
      err_d_(other.err_d_),
      range_u_(other.range_u_),
      range_d_(other.range_d_),
      vars_(other.vars_) {
}

Parameter::Parameter(Parameter&& other)
    : name_(""),
      val_(0.0),
      err_u_(1.0),
      err_d_(-1.0),
      range_u_(std::numeric_limits<double>::max()),
      range_d_(std::numeric_limits<double>::lowest()) {
  swap(*this, other);
}

Parameter& Parameter::operator=(Parameter other) {
  swap(*this, other);
  return (*this);
}

std::ostream& Parameter::PrintHeader(std::ostream &out) {
  std::string line =
    (boost::format("%-70s %-10.4f %-9.4f %-7.4f")
    % "name" % "value" % "error_d" % "error_u").str();
  std::string div(line.length(), '-');
  out << div << std::endl;
  out << line << std::endl;
  out << div << std::endl;
  return out;
}

std::ostream& operator<< (std::ostream &out, Parameter &val) {
  out << boost::format("%-70s %-10.4f %-9.4f %-7.4f")
  % val.name()
  % val.val()
  % val.err_d()
  % val.err_u();
  return out;
}
}
