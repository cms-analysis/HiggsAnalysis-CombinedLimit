#ifndef CombineTools_Logging_h
#define CombineTools_Logging_h
#include <string>
#include <iostream>

namespace ch {

#define FNERROR(x) FnError(x, __FILE__, __LINE__, __PRETTY_FUNCTION__)

#define LOGLINE(x, y) LogLine(x, __func__, y)

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
}

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

#endif
