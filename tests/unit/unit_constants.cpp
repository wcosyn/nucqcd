// Tools/constants.hpp(.in): physical constants. Known bugs are documented as [!mayfail] tests
// (they are reported but do not fail the suite); remove the tag when the bug is fixed.
#include <test_support.hpp>

#include <constants.hpp>

#include <cmath>
#include <complex>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Unit conversions and basic constants", "[constants]") {
  CHECK_THAT(PI, WithinRel(std::acos(-1.), 1e-15));
  CHECK_THAT(HBARC * INVHBARC, WithinRel(1., 1e-5));
  CHECK_THAT(DEGRTORAD * RADTODEGR, WithinRel(1., 1e-8));
  CHECK_THAT(DEGRTORAD, WithinRel(PI / 180., 1e-8));
  CHECK_THAT(1. / ALPHA, WithinRel(137.036, 1e-5));
  CHECK(I_UNIT * I_UNIT == std::complex<double>(-1., 0.));
}

TEST_CASE("Hadron masses: MeV and GeV versions agree, proton lighter than neutron", "[constants]") {
  CHECK(MASSP < MASSN);
  CHECK_THAT(MASSP_G * 1e3, WithinRel(MASSP, 1e-12));
  CHECK_THAT(MASSN_G * 1e3, WithinRel(MASSN, 1e-12));
  CHECK_THAT(MASSD_G * 1e3, WithinRel(MASSD, 1e-12));
  CHECK_THAT(MASSn, WithinRel(0.5 * (MASSP + MASSN), 1e-15));
  // the second mass set agrees with the first to the quoted digits
  CHECK_THAT(M_P, WithinRel(MASSP, 1e-5));
  CHECK_THAT(M_N, WithinRel(MASSN, 1e-5));
  CHECK_THAT(M_NUCL, WithinRel(MASSn, 1e-5));
  // deuteron binding energy ~2.2 MeV
  CHECK_THAT(MASSP + MASSN - MASSD, WithinAbs(2.22, 0.05));
}

TEST_CASE("Known bug: MASSn_G is in MeV instead of GeV", "[constants][known-bug][!mayfail]") {
  CHECK_THAT(MASSn_G, WithinRel(0.5 * (MASSP_G + MASSN_G), 1e-12));
}

TEST_CASE("Known bug: M_Z is 10x too large (used in Knockout/WeakQECross.cpp)", "[constants][known-bug][!mayfail]") {
  CHECK_THAT(M_Z, WithinRel(91187.6, 1e-4));
}

TEST_CASE("Known bug: INVHBARC = 1/197.328 is inconsistent with HBARC = 197.327 (5e-6 relative)",
          "[constants][known-bug][!mayfail]") {
  CHECK_THAT(HBARC * INVHBARC, WithinRel(1., 1e-9));
}

TEST_CASE("Electroweak constants", "[constants]") {
  CHECK_THAT(M_W, WithinRel(80385., 1e-4));
  CHECK_THAT(G_FERMI, WithinRel(1.16637e-11, 1e-5));  // MeV^-2
  CHECK(SIN2W > 0.22);
  CHECK(SIN2W < 0.24);
}
