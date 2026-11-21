/*! \file FortranString.hpp
 * \brief Bounded copies of C++ strings into fixed-size character buffers shared with Fortran
 * \author Wim Cosyn
 * \date 21/11/2026
 *
 * \details Fortran code reads file names from fixed-size buffers: CHARACTER*N variables in COMMON
 * blocks (Alekhin, GRSV) or a bind(C) C_CHAR array (He3 wave function). Copying into them with
 * strcpy/strcat overflowed silently for long paths. These functions check the length and throw
 * std::length_error instead.
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef FORTRANSTRING_HPP
#define FORTRANSTRING_HPP

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>

namespace nucqcd {

/*! Copies a string into a Fortran CHARACTER*N variable (for example in a COMMON block).
 * Fortran strings are padded with blanks and have no terminating NUL.
 * \param dest the CHARACTER*N buffer, declared as char[N] on the C++ side
 * \param s string to copy, at most N characters
 * \throw std::length_error if s is longer than N characters
 */
template <std::size_t N>
void to_fortran_string(char (&dest)[N], const std::string &s) {
  if (s.size() > N)
    throw std::length_error("to_fortran_string: '" + s + "' is longer than the " + std::to_string(N) +
                            " characters of the Fortran variable");
  std::memcpy(dest, s.data(), s.size());
  std::memset(dest + s.size(), ' ', N - s.size());
}

/*! Copies a string into a NUL-terminated C buffer (for example a bind(C) C_CHAR array read by Fortran).
 * \param dest the buffer, char[N]
 * \param s string to copy, at most N-1 characters
 * \throw std::length_error if s does not fit with its terminating NUL
 */
template <std::size_t N>
void to_c_string(char (&dest)[N], const std::string &s) {
  if (s.size() >= N)
    throw std::length_error("to_c_string: '" + s + "' does not fit in a buffer of " + std::to_string(N) +
                            " characters");
  std::memcpy(dest, s.data(), s.size());
  dest[s.size()] = '\0';
}

}  // namespace nucqcd

/** @} */
#endif
