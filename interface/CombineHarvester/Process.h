#ifndef CombineTools_Process_h
#define CombineTools_Process_h
#include <memory>
#include <string>
#include "TH1.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "CombineTools/interface/MakeUnique.h"

namespace ch {

class Process {
 public:
  Process();
  ~Process();
  Process(Process const& other);
  Process(Process&& other);
  Process& operator=(Process other);

  void set_bin(std::string const& bin) { bin_ = bin; }
  std::string const& bin() const { return bin_; }

  void set_rate(double const& rate) { rate_ = rate; }
  double rate() const { return norm_ ? norm_->getVal() * rate_ : rate_; }
  double no_norm_rate() const { return rate_; }

  void set_process(std::string const& process) { process_ = process; }
  std::string const& process() const { return process_; }

  void set_signal(bool const& signal) { signal_ = signal; }
  bool signal() const { return signal_; }

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

  void set_shape(std::unique_ptr<TH1> shape) { shape_ = std::move(shape); }
  TH1 const* shape() const { return shape_.get(); }

  void set_pdf(RooAbsPdf* pdf) { pdf_ = pdf; }
  RooAbsPdf const* pdf() const { return pdf_; }

  void set_norm(RooAbsReal* norm) { norm_ = norm; }
  RooAbsReal const* norm() const { return norm_; }

  void SetNormShape(std::unique_ptr<TH1> shape);

  void SetNormShapeAndRate(std::unique_ptr<TH1> shape);

  friend std::ostream& operator<< (std::ostream &out, Process &val);
  static std::ostream& PrintHeader(std::ostream &out);


 private:
  std::string bin_;
  double rate_;
  std::string process_;
  bool signal_;
  std::string analysis_;
  std::string era_;
  std::string channel_;
  int bin_id_;
  std::string mass_;
  std::unique_ptr<TH1> shape_;
  RooAbsPdf* pdf_;
  RooAbsReal* norm_;



  friend void swap(Process& first, Process& second);
};
}

#endif
