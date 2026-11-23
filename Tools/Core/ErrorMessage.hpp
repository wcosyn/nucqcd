/*! \file ErrorMessage.hpp
 * \brief Builds exception messages with stream syntax
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * \details The libraries report errors by throwing the standard exceptions (until 2026-11 they
 * printed to std::cerr and called exit(1)):
 * - std::invalid_argument: an invalid option, name, format string or quantum number;
 * - std::out_of_range: a value outside the range of a grid or of a parametrization;
 * - std::logic_error: a call in the wrong order, e.g. reading a grid before it is filled;
 * - std::runtime_error: input that cannot be read, or a calculation that fails.
 *
 * An uncaught exception ends the program with its message, like the old exit(1) did, but a caller
 * (or a test) can now catch it. ErrorMessage builds the message from the same stream expressions
 * that went to std::cerr:
 * \code
 * throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range: " << r);
 * \endcode
 * Exceptions cannot pass through the C and Fortran integration routines (Cuba, cubature, GSL);
 * thrown inside an integrand called from there, they still end the program.
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef ERRORMESSAGE_HPP
#define ERRORMESSAGE_HPP

#include <sstream>
#include <stdexcept>
#include <string>

namespace nucqcd {

/*! Collects a message with operator<<; converts to std::string, so it can be passed to the
 * constructor of a standard exception. */
class ErrorMessage {
 public:
  /// appends a value to the message
  template <class T>
  ErrorMessage &operator<<(const T &value) {
    stream_ << value;
    return *this;
  }
  operator std::string() const { return stream_.str(); } ///< the message

 private:
  std::ostringstream stream_; ///< the message being built
};

}  // namespace nucqcd

/** @} */
#endif
