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
      range_d_(std::numeric_limits<double>::lowest()),
      frozen_(false) {
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
  swap(first.frozen_, second.frozen_);
  swap(first.vars_, second.vars_);
  swap(first.groups_, second.groups_);
}

Parameter::Parameter(Parameter const& other)
    : name_(other.name_),
      val_(other.val_),
      err_u_(other.err_u_),
      err_d_(other.err_d_),
      range_u_(other.range_u_),
      range_d_(other.range_d_),
      frozen_(other.frozen_),
      vars_(other.vars_),
      groups_(other.groups_) {
}

Parameter::Parameter(Parameter&& other)
    : name_(""),
      val_(0.0),
      err_u_(1.0),
      err_d_(-1.0),
      range_u_(std::numeric_limits<double>::max()),
      range_d_(std::numeric_limits<double>::lowest()),
      frozen_(false) {
  swap(*this, other);
}

Parameter& Parameter::operator=(Parameter other) {
  swap(*this, other);
  return (*this);
}

std::ostream& Parameter::PrintHeader(std::ostream &out) {
  std::string line =
    (boost::format("%-70s %-10.4f %-9.4f %-7.4f %-6s %s")
    % "name" % "value" % "error_d" % "error_u" % "frozen" % "groups").str();
  std::string div(line.length(), '-');
  out << div << std::endl;
  out << line << std::endl;
  out << div << std::endl;
  return out;
}

std::ostream& operator<< (std::ostream &out, Parameter &val) {
  std::string grps = "";
  unsigned i = 0;
  unsigned n = val.groups().size();
  for (auto it = val.groups().begin(); i < n; ++it, ++i) {
    if (i < (n-1)) {
      grps += (*it + ",");
    } else {
      grps += (*it);
    }
  }
  out << boost::format("%-70s %-10.4f %-9.4f %-7.4f %-6i %s")
  % val.name()
  % val.val()
  % val.err_d()
  % val.err_u()
  % val.frozen()
  % grps;
  return out;
}
}
