#ifndef CombineTools_Object_h
#define CombineTools_Object_h
#include <string>

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

 private:
  std::string bin_;
  std::string process_;
  std::string analysis_;
  std::string era_;
  std::string channel_;
  int bin_id_;
  std::string mass_;
  friend void swap(Object& first, Object& second);
};
}

#endif
