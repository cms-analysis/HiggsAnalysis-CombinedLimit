#ifndef CombineTools_Object_h
#define CombineTools_Object_h
#include <string>
#include <map>

namespace ch {

class Object {
 public:
  Object();
  virtual ~Object();
  Object(Object const& other);
  Object(Object&& other);
  Object& operator=(Object other);

  virtual void set_bin(std::string const& bin) { bin_ = bin; }
  virtual std::string const& bin() const { return bin_; }

  virtual void set_process(std::string const& process) { process_ = process; }
  virtual std::string const& process() const { return process_; }

  void set_signal(bool const& signal) { signal_ = signal; }
  bool signal() const { return signal_; }

  virtual void set_analysis(std::string const& analysis) { analysis_ = analysis; }
  virtual std::string const& analysis() const { return analysis_; }

  virtual void set_era(std::string const& era) { era_ = era; }
  virtual std::string const& era() const { return era_; }

  virtual void set_channel(std::string const& channel) { channel_ = channel; }
  virtual std::string const& channel() const { return channel_; }

  virtual void set_bin_id(int const& bin_id) { bin_id_ = bin_id; }
  virtual int bin_id() const { return bin_id_; }

  virtual void set_mass(std::string const& mass) { mass_ = mass; }
  virtual std::string const& mass() const { return mass_; }

  virtual void set_attribute(std::string const& attr_label, std::string const& attr_value);
  virtual void delete_attribute(std::string const& attr_label) { attributes_.erase(attr_label); }
  virtual void set_all_attributes(std::map<std::string,std::string> const& attrs_) { attributes_ = attrs_;} 
  virtual std::map<std::string,std::string> const& all_attributes() const { return attributes_;}
  virtual std::string const attribute(std::string const& attr_label) const { return attributes_.count(attr_label) >0 ? attributes_.at(attr_label) : "" ; }

 private:
  std::string bin_;
  std::string process_;
  bool signal_;
  std::string analysis_;
  std::string era_;
  std::string channel_;
  int bin_id_;
  std::string mass_;
  std::map<std::string,std::string> attributes_;
  friend void swap(Object& first, Object& second);
};
}

#endif
