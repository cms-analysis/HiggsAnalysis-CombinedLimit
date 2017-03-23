#include "CombineHarvester/CombineTools/interface/Object.h"
#include <iostream>
namespace ch {

Object::Object()
    : bin_(""),
      process_(""),
      signal_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_("") {
  }

Object::~Object() { }

void swap(Object& first, Object& second) {
  using std::swap;
  swap(first.bin_, second.bin_);
  swap(first.process_, second.process_);
  swap(first.signal_, second.signal_);
  swap(first.analysis_, second.analysis_);
  swap(first.era_, second.era_);
  swap(first.channel_, second.channel_);
  swap(first.bin_id_, second.bin_id_);
  swap(first.mass_, second.mass_);
  swap(first.attributes_, second.attributes_);
}

Object::Object(Object const& other)
    : bin_(other.bin_),
      process_(other.process_),
      signal_(other.signal_),
      analysis_(other.analysis_),
      era_(other.era_),
      channel_(other.channel_),
      bin_id_(other.bin_id_),
      mass_(other.mass_),
      attributes_(other.attributes_) {
}

Object::Object(Object&& other)
    : bin_(""),
      process_(""),
      signal_(false),
      analysis_(""),
      era_(""),
      channel_(""),
      bin_id_(0),
      mass_("") {
  swap(*this, other);
}

void Object::set_attribute(std::string const& attr_label, std::string const& attr_value){
    if(attributes_.count(attr_label)>0){
        attributes_[attr_label]=attr_value;
    } else {
        attributes_.insert(std::pair<std::string,std::string>(attr_label,attr_value));
    }
}

Object& Object::operator=(Object other) {
  swap(*this, other);
  return (*this);
}
}
