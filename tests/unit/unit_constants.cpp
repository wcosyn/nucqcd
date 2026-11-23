// Tools/Core/constants.hpp(.in): physical constants.
#include <test_support.hpp>

#include <Tools/Core/constants.hpp>

#include <cmath>
#include <complex>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Unit conversions and basic constants", "[constants]") {
  CHECK_THAT(PI, WithinRel(std::acos(-1.), 1e-15));
  CHECK_THAT(nucqcd::hbarc * nucqcd::inv_hbarc, WithinRel(1., 1e-5));
  CHECK_THAT(nucqcd::deg_to_rad * nucqcd::rad_to_deg, WithinRel(1., 1e-8));
  CHECK_THAT(nucqcd::deg_to_rad, WithinRel(PI / 180., 1e-8));
  CHECK_THAT(1. / nucqcd::alpha_em, WithinRel(137.036, 1e-5));
  CHECK(nucqcd::i_unit * nucqcd::i_unit == std::complex<double>(-1., 0.));
}

TEST_CASE("Hadron masses: MeV and GeV versions agree, proton lighter than neutron", "[constants]") {
  CHECK(nucqcd::mass_p < nucqcd::mass_n);
  CHECK_THAT(nucqcd::mass_p_GeV * 1e3, WithinRel(nucqcd::mass_p, 1e-12));
  CHECK_THAT(nucqcd::mass_n_GeV * 1e3, WithinRel(nucqcd::mass_n, 1e-12));
  CHECK_THAT(nucqcd::mass_d_GeV * 1e3, WithinRel(nucqcd::mass_d, 1e-12));
  CHECK_THAT(nucqcd::mass_nucleon, WithinRel(0.5 * (nucqcd::mass_p + nucqcd::mass_n), 1e-15));
  // deuteron binding energy ~2.2 MeV
  CHECK_THAT(nucqcd::mass_p + nucqcd::mass_n - nucqcd::mass_d, WithinAbs(2.22, 0.05));
}

TEST_CASE("mass_nucleon_GeV is the average nucleon mass in GeV", "[constants]") {
  CHECK_THAT(nucqcd::mass_nucleon_GeV, WithinRel(0.5 * (nucqcd::mass_p_GeV + nucqcd::mass_n_GeV), 1e-12));
}

TEST_CASE("Z boson mass", "[constants]") {
  CHECK_THAT(nucqcd::mass_z, WithinRel(91187.6, 1e-4));
}

TEST_CASE("inv_hbarc is the inverse of hbarc", "[constants]") {
  CHECK_THAT(nucqcd::hbarc * nucqcd::inv_hbarc, WithinRel(1., 1e-9));
}

TEST_CASE("Electroweak constants", "[constants]") {
  CHECK_THAT(nucqcd::mass_w, WithinRel(80385., 1e-4));
  CHECK_THAT(nucqcd::g_fermi, WithinRel(1.16637e-11, 1e-5));  // MeV^-2
  CHECK(nucqcd::sin2_theta_w > 0.22);
  CHECK(nucqcd::sin2_theta_w < 0.24);
}

// The constants are inline constexpr variables in namespace nucqcd. The macro aliases (MASSP, HBARC,
// M_P, HOMEDIR, ...) were removed in 2026-11 (roadmap step 6); only PI is left, a literal shared with
// Integrator/numint/macros.hpp.
static_assert(nucqcd::mass_nucleon == 0.5 * (nucqcd::mass_n + nucqcd::mass_p), "constexpr average nucleon mass");
static_assert(nucqcd::constants::hbarc == nucqcd::hbarc, "inline namespace nucqcd::constants");
#if defined(MASSP) || defined(MASSN) || defined(HBARC) || defined(INVHBARC) || defined(MASSD) || defined(I_UNIT) || \
    defined(DEGRTORAD) || defined(RADTODEGR) || defined(STRANGEUFLOW) || defined(M_P) || defined(M_W) || defined(HOMEDIR)
#error "constants.hpp should no longer define the old macro aliases"
#endif

TEST_CASE("PI is the literal pi, consistent with nucqcd::pi", "[constants]") {
  CHECK(PI == nucqcd::pi);
  CHECK(nucqcd::deg_to_rad * 180. == PI);
  CHECK(nucqcd::i_unit == std::complex<double>(0., 1.));
  CHECK(std::string(nucqcd::share_dir).size() > 0);
}
