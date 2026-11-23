// MePhys/DIS: He3 wave function (Tools/WaveFunctions/He3wf, a wrapper around Fortran code by Misak Sargsian, data from the
// Bochum group): normalization, rotation invariance and antisymmetry.
// Until 2026-11 He3wf::getWF fixed ms3 = mA - ms1 - ms2, which left out the D- and P-wave components
// whose nucleon spins do not add up to mA: the norm was 0.971 for AV18 and the spin-summed density was
// not rotation invariant (ROADMAP step 4, "He3 normalization").
#include <test_support.hpp>

#include <DIS/He3Cross.hpp>
#include <Tools/WaveFunctions/He3wf.hpp>
#include <Tools/Core/constants.hpp>

#include <Tools/Core/Vector3.hpp>
#include <gsl/gsl_integration.h>

#include <cmath>
#include <complex>
#include <random>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
// sum over all spins (mA = +1) and isospins of |Psi|^2 at (p1, p2) [MeV]; mt3 = +1 (proton) and -1 (neutron) separately
void density(He3wf &wf, const nucqcd::Vector3 &p1, const nucqcd::Vector3 &p2, double &dens_p3, double &dens_n3) {
  dens_p3 = dens_n3 = 0.;
  for (int ms1 : {-1, 1})
    for (int ms2 : {-1, 1})
      for (int ms3 : {-1, 1})
        for (int mt1 : {-1, 1})
          for (int mt2 : {-1, 1}) {
            const int mt3 = 1 - mt1 - mt2;
            if (std::abs(mt3) != 1) continue;
            (mt3 == 1 ? dens_p3 : dens_n3) += std::norm(wf.getWF(p1, p2, ms1, ms2, ms3, mt1, mt2, 1));
          }
}
}  // namespace

TEST_CASE("He3wf AV18: normalized to 1 in Jacobi coordinates, 2/3 and 1/3 for the third nucleon a proton or neutron",
          "[he3]") {
  // Jacobi momenta p12 = (p1-p2)/2, p3 (Jacobian 1); the spin-summed density is rotation invariant, so
  // p3 along z (4 pi) and p12 in the xz plane (2 pi). Gauss-Legendre with p = a t/(1-t).
  He3wf wf("AV18", nucqcd_test::share_dir() + "/he3");
  const int np = 24, nc = 12;
  const double a = 300.;  // MeV
  gsl_integration_glfixed_table *tp = gsl_integration_glfixed_table_alloc(np);
  gsl_integration_glfixed_table *tc = gsl_integration_glfixed_table_alloc(nc);
  double norm_p3 = 0., norm_n3 = 0.;
  for (int i = 0; i < np; i++) {
    double t, wt;
    gsl_integration_glfixed_point(0., 1., i, &t, &wt, tp);
    const double p12 = a * t / (1. - t), w12 = wt * a / ((1. - t) * (1. - t)) * p12 * p12;
    for (int j = 0; j < np; j++) {
      double u, wu;
      gsl_integration_glfixed_point(0., 1., j, &u, &wu, tp);
      const double p3 = a * u / (1. - u), w3 = wu * a / ((1. - u) * (1. - u)) * p3 * p3;
      for (int k = 0; k < nc; k++) {
        double c, wc;
        gsl_integration_glfixed_point(-1., 1., k, &c, &wc, tc);
        const nucqcd::Vector3 q12(p12 * std::sqrt(1. - c * c), 0., p12 * c), q3(0., 0., p3);
        double dp, dn;
        density(wf, q12 - 0.5 * q3, -q12 - 0.5 * q3, dp, dn);
        norm_p3 += 8. * PI * PI * w12 * w3 * wc * dp;
        norm_n3 += 8. * PI * PI * w12 * w3 * wc * dn;
      }
    }
  }
  gsl_integration_glfixed_table_free(tp);
  gsl_integration_glfixed_table_free(tc);
  // converged value 0.99985 (np = 48); 0.971 until 2026-11
  CHECK_THAT(norm_p3 + norm_n3, WithinAbs(1., 3e-3));
  CHECK_THAT(norm_p3, WithinRel(2. / 3., 3e-3));
  CHECK_THAT(norm_n3, WithinRel(1. / 3., 3e-3));
}

TEST_CASE("He3wf AV18: spin-summed density rotation invariant, antisymmetric under exchange", "[he3]") {
  He3wf wf("AV18", nucqcd_test::share_dir() + "/he3");
  std::mt19937 gen(7);
  std::normal_distribution<double> mom(0., 150.);  // MeV
  std::uniform_real_distribution<double> angle(0., 2. * PI);
  for (int it = 0; it < 10; it++) {
    const nucqcd::Vector3 p1(mom(gen), mom(gen), mom(gen)), p2(mom(gen), mom(gen), mom(gen));
    nucqcd::Vector3 axis(mom(gen), mom(gen), mom(gen)), r1 = p1, r2 = p2;
    const double phi = angle(gen);
    r1.Rotate(phi, axis);
    r2.Rotate(phi, axis);
    double dp, dn, rp, rn;
    density(wf, p1, p2, dp, dn);
    density(wf, r1, r2, rp, rn);
    CHECK_THAT(rp + rn, WithinRel(dp + dn, 1e-5));  // until 2026-11 up to 30% off
    // exchanging nucleons 1 and 2 (momentum, spin, isospin) changes the sign
    for (int ms1 : {-1, 1})
      for (int ms3 : {-1, 1}) {
        const std::complex<double> psi = wf.getWF(p1, p2, ms1, -1, ms3, 1, -1, 1);
        const std::complex<double> swapped = wf.getWF(p2, p1, -1, ms1, ms3, -1, 1, 1);
        CHECK(std::abs(psi + swapped) <= 1e-6 * std::sqrt(dp + dn));
      }
  }
  CHECK_THROWS_AS(wf.getWF(nucqcd::Vector3(), nucqcd::Vector3(), 1, 1, 3, 1, 1, 1), std::invalid_argument);
  CHECK_THROWS_AS(He3wf("AV19", nucqcd_test::share_dir() + "/he3"), std::invalid_argument);
}

TEST_CASE("He3Cross: spectator momentum distribution integrates to Z = 2 (struck proton) and N = 1 (struck neutron)",
          "[he3]") {
  // the factor A = 3 between the He3wf normalization (over all isospin arrangements) and the note's
  // (per struck nucleon) was missing until 2026-11, giving 2/3 and 1/3
  for (bool proton : {true, false}) {
    He3Cross he3("AV18", nucqcd_test::share_dir(), "SLAC", proton);
    const int np = 24, nc = 12;
    const double a = 300.;
    gsl_integration_glfixed_table *tp = gsl_integration_glfixed_table_alloc(np);
    gsl_integration_glfixed_table *tc = gsl_integration_glfixed_table_alloc(nc);
    double integral = 0.;
    for (int i = 0; i < np; i++) {
      double t, wt;
      gsl_integration_glfixed_point(0., 1., i, &t, &wt, tp);
      const double p12 = a * t / (1. - t), w12 = wt * a / ((1. - t) * (1. - t)) * p12 * p12;
      for (int j = 0; j < np; j++) {
        double u, wu;
        gsl_integration_glfixed_point(0., 1., j, &u, &wu, tp);
        const double p3 = a * u / (1. - u), w3 = wu * a / ((1. - u) * (1. - u)) * p3 * p3;
        for (int k = 0; k < nc; k++) {
          double c, wc;
          gsl_integration_glfixed_point(-1., 1., k, &c, &wc, tc);
          const nucqcd::Vector3 q12(p12 * std::sqrt(1. - c * c), 0., p12 * c), q3(0., 0., p3);
          integral += 8. * PI * PI * w12 * w3 * wc * he3.getMomentumDistribution(q12 - 0.5 * q3, -q12 - 0.5 * q3);
        }
      }
    }
    gsl_integration_glfixed_table_free(tp);
    gsl_integration_glfixed_table_free(tc);
    CAPTURE(proton);
    CHECK_THAT(integral, WithinRel(proton ? 2. : 1., 3e-3));
  }
}
