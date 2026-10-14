#ifndef CombineTools_Algorithm_h
#define CombineTools_Algorithm_h
#include <algorithm>
#include <vector>
#include "boost/range/begin.hpp"
#include "boost/range/end.hpp"
#include "boost/regex.hpp"
#include "boost/range/algorithm_ext/erase.hpp"

namespace ch {
template <typename Range, typename Predicate>
void erase_if(Range & r, Predicate p) {
  r.erase(std::remove_if(r.begin(), r.end(), p), r.end());
}

template<typename Range, typename Predicate>
bool any_of(const Range &r, Predicate p) {
  return std::any_of(boost::begin(r), boost::end(r), p);
}

template<typename Range, typename T>
bool contains(const Range &r, T p) {
  return std::find(boost::begin(r), boost::end(r), p) != boost::end(r);
}

// This version needed for cases like: contains({"A", "B"}, "B"),
// which are not implicitly convertible in the above generic version
template<typename R, typename T>
bool contains(const std::initializer_list<R> &r, T p) {
  return std::find(boost::begin(r), boost::end(r), p) != boost::end(r);
}

template <typename T>
bool contains_rgx(const std::vector<boost::regex>& r, T p) {
  for (auto const& rgx : r)
    if (regex_match(p, rgx)) return true;
  return false;
}

template <typename Input, typename Filter, typename Converter>
void FilterContaining(Input& in, Filter const& filter, Converter fn,
                      bool cond) {
  boost::remove_erase_if(in, [&](typename Input::value_type const& p) {
    return cond != ch::contains(filter, fn(p));
  });
}

template <typename Input, typename Filter, typename Converter>
void FilterContainingRgx(Input& in, Filter const& filter, Converter fn,
                         bool cond) {
  std::vector<boost::regex> rgx;
  for (auto const& ele : filter) rgx.emplace_back(ele);
  boost::remove_erase_if(in, [&](typename Input::value_type const& p) {
    return cond != ch::contains_rgx(rgx, fn(p));
  });
}

template <typename Input, typename Filter, typename Converter, typename Funcarg>
void FilterContaining(Input& in, Filter const& filter, Converter fn, Funcarg arg,
                      bool cond) {
  boost::remove_erase_if(in, [&](typename Input::value_type const& p) {
    return cond != ch::contains(filter, fn(p,arg));
  });
}

template <typename Input, typename Filter, typename Converter, typename Funcarg>
void FilterContainingRgx(Input& in, Filter const& filter, Converter fn, Funcarg arg,
                         bool cond) {
  std::vector<boost::regex> rgx;
  for (auto const& ele : filter) rgx.emplace_back(ele);
  boost::remove_erase_if(in, [&](typename Input::value_type const& p) {
    return cond != ch::contains_rgx(rgx, fn(p,arg));
  });
}
}

#endif
