#ifndef CombineTools_Systematic_h
#define CombineTools_Systematic_h
#include <memory>
#include <string>
#include "TH1.h"
#include "RooDataHist.h"
#include "CombineHarvester/CombineTools/interface/MakeUnique.h"
#include "CombineHarvester/CombineTools/interface/Object.h"

namespace ch {

class Systematic : public Object {
 public:
  Systematic();
  ~Systematic();
  Systematic(Systematic const& other);
  Systematic(Systematic&& other);
  Systematic& operator=(Systematic other);

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

  TH1 const* shape_u() const { return shape_u_.get(); }

  std::unique_ptr<TH1> ClonedShapeU() const;
  std::unique_ptr<TH1> ClonedShapeD() const;

  TH1F ShapeUAsTH1F() const;
  TH1F ShapeDAsTH1F() const;

  TH1 const* shape_d() const { return shape_d_.get(); }

  RooDataHist const* data_u() const { return data_u_; }

  RooDataHist const* data_d() const { return data_d_; }

  void set_data(RooDataHist* data_u, RooDataHist* data_d,
                RooDataHist const* nominal);

  void set_shapes(std::unique_ptr<TH1> shape_u, std::unique_ptr<TH1> shape_d,
                  TH1 const* nominal);

  void set_shapes(TH1 const& shape_u, TH1 const& shape_d,
                  TH1 const& nominal);

  friend std::ostream& operator<< (std::ostream &out, Systematic const& val);
  static std::ostream& PrintHeader(std::ostream &out);

  void SwapUpAndDown();

 private:
  std::string name_;
  std::string type_;
  double value_u_;
  double value_d_;
  double scale_;
  bool asymm_;
  std::unique_ptr<TH1> shape_u_;
  std::unique_ptr<TH1> shape_d_;
  RooDataHist * data_u_;
  RooDataHist * data_d_;

  friend void swap(Systematic& first, Systematic& second);
};
}

#endif
