// Tools/constants.hpp(.in): physical constants.
#include <test_support.hpp>

#include <constants.hpp>

#include <cmath>
#include <complex>
#include <string>

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
  // the second set of mass names are aliases of the first
  CHECK(M_P == MASSP);
  CHECK(M_N == MASSN);
  CHECK(M_NUCL == MASSn);
  // deuteron binding energy ~2.2 MeV
  CHECK_THAT(MASSP + MASSN - MASSD, WithinAbs(2.22, 0.05));
}

TEST_CASE("MASSn_G is the average nucleon mass in GeV", "[constants]") {
  CHECK_THAT(MASSn_G, WithinRel(0.5 * (MASSP_G + MASSN_G), 1e-12));
}

TEST_CASE("Z boson mass", "[constants]") {
  CHECK_THAT(M_Z, WithinRel(91187.6, 1e-4));
}

TEST_CASE("INVHBARC is the inverse of HBARC", "[constants]") {
  CHECK_THAT(HBARC * INVHBARC, WithinRel(1., 1e-9));
}

TEST_CASE("Electroweak constants", "[constants]") {
  CHECK_THAT(M_W, WithinRel(80385., 1e-4));
  CHECK_THAT(G_FERMI, WithinRel(1.16637e-11, 1e-5));  // MeV^-2
  CHECK(SIN2W > 0.22);
  CHECK(SIN2W < 0.24);
}

// The constants are inline constexpr variables in namespace nucqcd; the macros are aliases (step 5).
static_assert(nucqcd::mass_nucleon == 0.5 * (nucqcd::mass_n + nucqcd::mass_p), "constexpr average nucleon mass");
static_assert(nucqcd::constants::hbarc == nucqcd::hbarc, "inline namespace nucqcd::constants");

TEST_CASE("Macros are aliases of the nucqcd constants", "[constants]") {
  CHECK(PI == nucqcd::pi);  // PI is a literal (shared with Integrator/numint/macros.hpp)
  CHECK(MASSP == nucqcd::mass_p);
  CHECK(M_P == nucqcd::mass_p);
  CHECK(MASSn == nucqcd::mass_nucleon);
  CHECK(M_PIP == MASSPI);
  CHECK(M_RHO == MASSRHO);  // one rho mass since 2026-11 (was 775.49 and 775.5)
  CHECK(DEGRTORAD * 180. == PI);
  CHECK(I_UNIT == std::complex<double>(0., 1.));
  CHECK(std::string(HOMEDIR) == nucqcd::share_dir);
}
