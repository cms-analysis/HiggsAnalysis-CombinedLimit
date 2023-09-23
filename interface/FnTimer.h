#ifndef FnTimer_h
#define FnTimer_h
#include <string>
#include <iostream>
#include <chrono>


/**
 * \brief Extracts the fully-qualified function name from a complete function
 *signature
 * \details The input str will typically be the compiler variable
 * `__PRETTY_FUNCTION__`. The return type and arguments are then removed, e.g.:
 *
 *      std::string GetQualififedName(std::string const& str)
 *
 * would be converted to
 *
 *      GetQualififedName
 */
std::string GetQualififedName(std::string const& str);

/**
 * Conveniently initialise a FnTimer instance
 *
 * This macro should be placed at the start of a function, e.g.:
 *
 *     void MyFunction() {
 *      LAUNCH_FUNCTION_TIMER(__timer__, __token__)
 *     }
 *
 * The arguments are the names of two objects (a FnTimer and a
 * FnTimer::Token) that will be created by this macro. Note that the
 * FnTimer will be assigned the current function name automatically.
 */
#define LAUNCH_FUNCTION_TIMER(x, y)                             \
  static FnTimer x(GetQualififedName(__PRETTY_FUNCTION__)); \
  auto y = x.Inc();

/**
 * Determine the total amount of time spent in a function
 *
 * An FnTimer instance should typically be declared as a static variable at
 * the beginning of a function, follwed by a call to the Inc() method, which
 * will increment the counter. The Inc() method also returns an FnTimer::Token
 * object that records the time at which it is constructed and then destroyed,
 * the latter occurring automatically at the end of the function. At the end
 * of the program the FnTimer destructor will write a message to the screen
 * summarising the number of calls and the time information.
 *
 *  \note A simple way of using this class is via the LAUNCH_FUNCTION_TIMER(x,y)
 *  macro
 */
class FnTimer {
 public:
  class Token {
    public:
      explicit Token(FnTimer *src);
      ~Token();
    private:
      FnTimer *src_;
  };

  explicit FnTimer(std::string name);
  ~FnTimer();
  Token Inc();
  inline void StartTimer() {
    start_ = std::chrono::high_resolution_clock::now();
  }
  inline void StopTimer() {
    end_ = std::chrono::high_resolution_clock::now();
    elapsed_ += std::chrono::duration<double>(end_ - start_).count();
    start_ = std::chrono::high_resolution_clock::now();
    end_ = std::chrono::high_resolution_clock::now();
    elapsed_overhead_ += std::chrono::duration<double>(end_ - start_).count();

  }

 private:
  std::string name_;
  unsigned calls_;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
  std::chrono::time_point<std::chrono::high_resolution_clock> end_;
  double elapsed_;
  double elapsed_overhead_;
};


#endif
