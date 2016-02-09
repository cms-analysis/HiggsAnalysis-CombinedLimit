#ifndef CombineTools_Systematics_h
#define CombineTools_Systematics_h
#include <vector>
#include <string>
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Logging.h"

namespace ch {
namespace syst {

  struct bin {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->bin(); }
  };

  struct analysis {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->analysis(); }
  };

  struct era {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->era(); }
  };

  struct channel {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->channel(); }
  };

  struct mass {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->mass(); }
  };

  struct process {
    typedef std::string type;
    inline static type get(ch::Process *p) { return p->process(); }
  };

  class bin_id {
   public:
    typedef int type;
    inline static type get(ch::Process *p) { return p->bin_id(); }
  };

  namespace detail {
    template <typename F>
    inline void cross_imp(F f) {
      f();
    }
    template <typename F, typename H, typename... Ts>
    inline void cross_imp(F f, std::vector<H> const &h,
                          std::vector<Ts> const &... t) {
      for (H const &he : h)
        cross_imp([&](Ts const &... ts) { f(he, ts...); }, t...);
    }

    template <typename... Ts>
    std::vector<std::tuple<Ts...>> cross(std::vector<Ts> const &... in) {
      std::vector<std::tuple<Ts...>> res;
      cross_imp([&](Ts const &... ts) {
        res.emplace_back(ts...);
      }, in...);
      return res;
    }
  }

  template<class... T>
  class SystMap {
   private:
    std::map<std::tuple<typename T::type...>, double> tmap_;

   public:
    SystMap& operator()(std::vector<typename T::type>... input, double val) {
      auto res = ch::syst::detail::cross(input...);
      for (auto const& a : res) {
        tmap_.insert(std::make_pair(a, val));
      }
      return *this;
    }

    bool Contains(ch::Process *p) const {
      if (p) {
        return tmap_.count(std::make_tuple(T::get(p)...));
      } else {
        return false;
      }
    }

    double ValU(ch::Process *p) const {
      if (p) {
        return tmap_.find(std::make_tuple(T::get(p)...))->second;
      } else {
        return 0.0;
      }
    }

    double ValD(ch::Process * /*p*/) const {
      return 0.0;
    }

    std::string Formula(ch::Process * /*p*/) const {
      return std::string("");
    }

    std::string Args(ch::Process * /*p*/) const {
      return std::string("");
    }

    static SystMap<T...> init(std::vector<typename T::type>... input,
                              double val) {
      SystMap<T...> x;
      return x(input..., val);
    }

    bool IsAsymm() const { return false; }

    std::set<std::tuple<typename T::type...>> GetTupleSet() const {
      std::set<std::tuple<typename T::type...>> res;
      for (auto const& x : tmap_) res.insert(x.first);
      return res;
    }

    auto GetTuple(ch::Process *p) const -> decltype( std::make_tuple(T::get(p)...) ) {
      if (p) {
        return std::make_tuple(T::get(p)...);
      } else {
        throw std::runtime_error(FNERROR("Supplied pointer is null"));
      }
    }
  };

  template<class... T>
  class SystMapAsymm {
   private:
    std::map<std::tuple<typename T::type...>, std::pair<double, double>> tmap_;

   public:
    SystMapAsymm &operator()(std::vector<typename T::type>... input,
                             double val_d, double val_u) {
      auto res = ch::syst::detail::cross(input...);
      for (auto const& a : res)
        tmap_.insert(std::make_pair(a, std::make_pair(val_d, val_u)));
      return *this;
    }

    bool Contains(ch::Process *p) const {
      if (p) {
        return tmap_.count(std::make_tuple(T::get(p)...));
      } else {
        return false;
      }
    }

    double ValD(ch::Process *p) const {
      if (p) {
        return tmap_.find(std::make_tuple(T::get(p)...))->second.first;
      } else {
        return 0.0;
      }
    }
    double ValU(ch::Process *p) const {
      if (p) {
        return tmap_.find(std::make_tuple(T::get(p)...))->second.second;
      } else {
        return 0.0;
      }
    }

    std::string Formula(ch::Process * /*p*/) const {
      return std::string("");
    }

    std::string Args(ch::Process * /*p*/) const {
      return std::string("");
    }

    static SystMapAsymm<T...> init(std::vector<typename T::type>... input,
                                   double val_d, double val_u) {
      SystMapAsymm<T...> x;
      return x(input..., val_d, val_u);
    }

    bool IsAsymm() const { return true; }

    std::set<std::tuple<typename T::type...>> GetTupleSet() const {
      std::set<std::tuple<typename T::type...>> res;
      for (auto const& x : tmap_) res.insert(x.first);
      return res;
    }

    auto GetTuple(ch::Process *p) const -> decltype( std::make_tuple(T::get(p)...) ) {
      if (p) {
        return std::make_tuple(T::get(p)...);
      } else {
        throw std::runtime_error(FNERROR("Supplied pointer is null"));
      }
    }
  };

  template<class... T>
  class SystMapFunc {
   private:
    std::map<std::tuple<typename T::type...>, std::pair<std::string, std::string>> tmap_;

   public:
    SystMapFunc &operator()(std::vector<typename T::type>... input,
                             std::string formula, std::string args) {
      auto res = ch::syst::detail::cross(input...);
      for (auto const& a : res)
        tmap_.insert(std::make_pair(a, std::make_pair(formula, args)));
      return *this;
    }

    bool Contains(ch::Process *p) const {
      if (p) {
        return tmap_.count(std::make_tuple(T::get(p)...));
      } else {
        return false;
      }
    }

    double ValD(ch::Process * /*p*/) const {
      return 0.0;
    }
    double ValU(ch::Process * /*p*/) const {
      return 0.0;
    }

    std::string Formula(ch::Process *p) const {
      if (p) {
        return tmap_.find(std::make_tuple(T::get(p)...))->second.first;
      } else {
        return std::string("");
      }
    }
    std::string Args(ch::Process *p) const {
      if (p) {
        return tmap_.find(std::make_tuple(T::get(p)...))->second.second;
      } else {
        return std::string("");
      }
    }

    static SystMapFunc<T...> init(std::vector<typename T::type>... input,
                                   std::string formula, std::string args) {
      SystMapFunc<T...> x;
      return x(input..., formula, args);
    }

    bool IsAsymm() const { return false; }

    std::set<std::tuple<typename T::type...>> GetTupleSet() const {
      std::set<std::tuple<typename T::type...>> res;
      for (auto const& x : tmap_) res.insert(x.first);
      return res;
    }

    auto GetTuple(ch::Process *p) const -> decltype( std::make_tuple(T::get(p)...) ) {
      if (p) {
        return std::make_tuple(T::get(p)...);
      } else {
        throw std::runtime_error(FNERROR("Supplied pointer is null"));
      }
    }
  };
}
}

#endif
