#ifndef CombineTools_Observation_h
#define CombineTools_Observation_h
#include <memory>
#include <string>
#include "TH1.h"
#include "RooAbsData.h"
#include "CombineTools/interface/MakeUnique.h"

namespace ch {

class Observation {
 public:
  Observation();
  ~Observation();
  Observation(Observation const& other);
  Observation(Observation&& other);
  Observation& operator=(Observation other);

  void set_bin(std::string const& bin) { bin_ = bin; }
  std::string const& bin() const { return bin_; }

  void set_rate(double const& rate) { rate_ = rate; }
  double rate() const { return rate_; }

  void set_analysis(std::string const& analysis) { analysis_ = analysis; }
  std::string const& analysis() const { return analysis_; }

  void set_era(std::string const& era) { era_ = era; }
  std::string const& era() const { return era_; }

  void set_channel(std::string const& channel) { channel_ = channel; }
  std::string const& channel() const { return channel_; }

  void set_bin_id(int const& bin_id) { bin_id_ = bin_id; }
  int bin_id() const { return bin_id_; }

  void set_mass(std::string const& mass) { mass_ = mass; }
  std::string const& mass() const { return mass_; }

  void set_shape(std::unique_ptr<TH1> shape, bool set_rate);
  TH1 const* shape() const { return shape_.get(); }

  std::unique_ptr<TH1> ClonedShape() const;
  std::unique_ptr<TH1> ClonedScaledShape() const;

  TH1F ShapeAsTH1F() const;

  void set_data(RooAbsData* data) { data_ = data; }
  RooAbsData const* data() const { return data_; }

  std::string process() const { return "data_obs"; }

  friend std::ostream& operator<< (std::ostream &out, Observation &val);
  static std::ostream& PrintHeader(std::ostream &out);

 private:
  std::string bin_;
  double rate_;
  std::string analysis_;
  std::string era_;
  std::string channel_;
  int bin_id_;
  std::string mass_;
  std::unique_ptr<TH1> shape_;
  RooAbsData* data_;

  friend void swap(Observation& first, Observation& second);
};
}

#endif
