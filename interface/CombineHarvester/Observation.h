#ifndef CombineTools_Observation_h
#define CombineTools_Observation_h
#include <memory>
#include <string>
#include "TH1.h"
#include "RooAbsData.h"
#include "CombineHarvester/CombineTools/interface/MakeUnique.h"
#include "CombineHarvester/CombineTools/interface/Object.h"

namespace ch {

class Observation : public Object {
 public:
  Observation();
  ~Observation();
  Observation(Observation const& other);
  Observation(Observation&& other);
  Observation& operator=(Observation other);

  void set_rate(double const& rate) { rate_ = rate; }
  double rate() const { return rate_; }

  void set_shape(std::unique_ptr<TH1> shape, bool set_rate);

  void set_shape(TH1 const& shape, bool set_rate);
  
  TH1 const* shape() const { return shape_.get(); }

  std::unique_ptr<TH1> ClonedShape() const;
  std::unique_ptr<TH1> ClonedScaledShape() const;

  TH1F ShapeAsTH1F() const;

  void set_data(RooAbsData* data) { data_ = data; }
  RooAbsData const* data() const { return data_; }

  std::string to_string() const;
  friend std::ostream& operator<< (std::ostream &out, Observation const& val);
  static std::ostream& PrintHeader(std::ostream &out);

 private:
  double rate_;
  std::unique_ptr<TH1> shape_;
  RooAbsData* data_;

  friend void swap(Observation& first, Observation& second);
};
}

#endif
