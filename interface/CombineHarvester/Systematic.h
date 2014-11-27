#ifndef CombineTools_Systematic_h
#define CombineTools_Systematic_h
#include <memory>
#include <string>
#include "TH1.h"
#include "CombineTools/interface/MakeUnique.h"

namespace ch {

class Systematic {
 public:
  Systematic();
  ~Systematic();
  Systematic(Systematic const& other);
  Systematic(Systematic&& other);
  Systematic& operator=(Systematic other);

  void set_bin(std::string const& bin) { bin_ = bin; }
  std::string const& bin() const { return bin_; }

  void set_process(std::string const& process) { process_ = process; }
  std::string const& process() const { return process_; }

  void set_signal(bool const& signal) { signal_ = signal; }
  bool signal() const { return signal_; }

  void set_name(std::string const& name) { name_ = name; }
  std::string const& name() const { return name_; }

  void set_type(std::string const& type) { type_ = type; }
  std::string const& type() const { return type_; }

  void set_value_u(double const& value_u) { value_u_ = value_u; }
  double value_u() const { return value_u_; }

  void set_value_d(double const& value_d) { value_d_ = value_d; }
  double value_d() const { return value_d_; }

  void set_scale(double const& scale) { scale_ = scale; }
  double scale() const { return scale_; }

  void set_asymm(bool const& asymm) { asymm_ = asymm; }
  bool asymm() const { return asymm_; }

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

  void set_shape_u(std::unique_ptr<TH1> shape_u) { shape_u_ = std::move(shape_u); }
  TH1 const* shape_u() const { return shape_u_.get(); }

  void set_shape_d(std::unique_ptr<TH1> shape_d) { shape_d_ = std::move(shape_d); }
  TH1 const* shape_d() const { return shape_d_.get(); }

  void SetShapesAndVals(std::unique_ptr<TH1> shape_u,
                        std::unique_ptr<TH1> shape_d, TH1 const* nominal);

  friend std::ostream& operator<< (std::ostream &out, Systematic const& val);
  static std::ostream& PrintHeader(std::ostream &out);

 private:
  std::string bin_;
  std::string process_;
  bool signal_;
  std::string name_;
  std::string type_;
  double value_u_;
  double value_d_;
  double scale_;
  bool asymm_;
  std::string analysis_;
  std::string era_;
  std::string channel_;
  int bin_id_;
  std::string mass_;
  std::unique_ptr<TH1> shape_u_;
  std::unique_ptr<TH1> shape_d_;

  friend void swap(Systematic& first, Systematic& second);
};
}

#endif
