// Tools/FortranString.hpp: bounded copies into the fixed-size character buffers shared with Fortran
// (COMMON blocks of the Alekhin and GRSV code, the He3 wave-function file name), and the length
// checks in the classes that use them. Before 2026-11 these were strcpy/strcat copies that overflowed
// silently for long paths.
#include <test_support.hpp>

#include <FortranString.hpp>
#include <He3wf.hpp>
#include <TKinematics2to2.h>
#include <constants.hpp>

#include <stdexcept>
#include <string>

TEST_CASE("to_fortran_string pads with blanks and has no NUL", "[fortran]") {
  char buf[8];
  nucqcd::to_fortran_string(buf, "abc");
  CHECK(std::string(buf, 8) == "abc     ");
  nucqcd::to_fortran_string(buf, "abcdefgh");  // exactly N characters fit
  CHECK(std::string(buf, 8) == "abcdefgh");
  nucqcd::to_fortran_string(buf, "");
  CHECK(std::string(buf, 8) == "        ");
  CHECK_THROWS_AS(nucqcd::to_fortran_string(buf, "abcdefghi"), std::length_error);
}

TEST_CASE("to_c_string NUL-terminates and leaves room for the NUL", "[fortran]") {
  char buf[8];
  nucqcd::to_c_string(buf, "abcdefg");  // N-1 characters fit
  CHECK(std::string(buf) == "abcdefg");
  nucqcd::to_c_string(buf, "ab");
  CHECK(std::string(buf) == "ab");
  CHECK_THROWS_AS(nucqcd::to_c_string(buf, "abcdefgh"), std::length_error);
}

TEST_CASE("He3wf rejects a data directory too long for the Fortran file name buffer", "[fortran]") {
  const std::string longdir = "/" + std::string(300, 'x');
  CHECK_THROWS_AS(He3wf("AV18", longdir), std::length_error);
}

TEST_CASE("TKinematics rejects a format string longer than its 100-character buffer", "[fortran]") {
  std::string format = "qsquared:wlab:pklab";
  while (format.size() < 100) format += ":pklab";
  CHECK_THROWS_AS(TKinematics2to2("", "", nucqcd::mass_d, nucqcd::mass_p, nucqcd::mass_n, format.c_str(), 1.e06, 600., 150.),
                  std::length_error);
}
