#ifndef CombineTools_Logging_h
#define CombineTools_Logging_h
#include <string>
#include <iostream>
#include <chrono>

namespace ch {

#define FNERROR(x) ch::FnError(x, __FILE__, __LINE__, __PRETTY_FUNCTION__)

#define LOGLINE(x, y) LogLine(x, __func__, y)

#define FNLOG(x) x << "[" << __func__ << "] "
#define FNLOGC(x, y) if (y) x << "[" << __func__ << "] "
/**
 * \brief Writes a logging message to a given ostream
 * \details Should be used with the macro LOGLINE which will call this function
 * and automatically insert the name of the calling function, e.g.:
 *
 *      void MyFunc() {
 *        LOGLINE(std::cout, "Some message");
 *      }
 *
 * produces
 *
 *      [MyFunc] Some message
 */
void LogLine(std::ostream& stream, std::string const& func,
             std::string const& message);

/**
 * \brief Generates an error message which includes the name of the calling
 * function and the filename and line number where the error occured
 * \details Should be used via the macro FNERROR which calls this function and
 * inserts the file, line and fn arguments automatically.
 */
std::string FnError(std::string const& message, std::string const& file,
                    unsigned line, std::string const& fn);

/**
 * \brief Extracts the fully-qualified function name from a complete function
 *signature
 * \details The input str will typically be the compiler variable
 * `__PRETTY_FUNCTION__`. The return type and arguments are then removed, e.g.:
 *
 *      std::string ch::GetQualififedName(std::string const& str)
 *
 * would be converted to
 *
 *      ch::GetQualififedName
 */
std::string GetQualififedName(std::string const& str);

/**
 * Conveniently initialise a ch::FnTimer instance
 *
 * This macro should be placed at the start of a function, e.g.:
 *
 *     void MyFunction() {
 *      LAUNCH_FUNCTION_TIMER(__timer__, __token__)
 *     }
 *
 * The arguments are the names of two objects (a ch::FnTimer and a
 * ch::FnTimer::Token) that will be created by this macro. Note that the
 * ch::FnTimer will be assigned the current function name automatically.
 */
#define LAUNCH_FUNCTION_TIMER(x, y)                             \
  static FnTimer x(ch::GetQualififedName(__PRETTY_FUNCTION__)); \
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
  void StartTimer();
  void StopTimer();

 private:
  std::string name_;
  unsigned calls_;
  std::chrono::time_point<std::chrono::system_clock> start_;
  std::chrono::time_point<std::chrono::system_clock> end_;
  double elapsed_;
};
}

#endif
