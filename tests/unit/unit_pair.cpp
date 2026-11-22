// MePhys: Pair (Tools/pair.h), the wave function of a nucleon pair in the harmonic-oscillator shell
// model, optionally with central and tensor correlations in the relative wave function. Unfinished
// code (used only by DoubleNModel); these tests cover the bugs fixed in 2026-11.
#include <test_support.hpp>

#include <pair.h>

#include <cmath>
#include <complex>

using Catch::Matchers::WithinAbs;

namespace {
struct Spherical {
  double r, costh, phi;
};
Spherical spherical(double x, double y, double z) {
  const double r = std::sqrt(x * x + y * y + z * z);
  return {r, r == 0. ? 1. : z / r, std::atan2(y, x)};
}
// two 0p3/2 nucleons of carbon: a proton (m_j = 1/2, m_s = 1/2) and a neutron (m_j = -1/2, m_s = -1/2)
Pair carbon_pn() { return Pair("paris", nucqcd_test::share_dir(), 12, 0, 1, 1, 3, 1, 1, 0, 1, -1, 3, -1, -1); }
}  // namespace

TEST_CASE("Pair: the relative wave function depends only on r1 - r2", "[pair]") {
  // until 2026-11 the azimuth of the relative coordinate was computed from r1 + r2
  const Pair pair = carbon_pn();
  // same r1 - r2 = (1.2, -0.4, 0.3), different c.m.
  const Spherical a1 = spherical(0.9, 0.5, 0.6), a2 = spherical(-0.3, 0.9, 0.3);
  const Spherical b1 = spherical(1.4, -1.1, -0.2), b2 = spherical(0.2, -0.7, -0.5);
  // relative state n=0, l=1, la=1, S=0, j=1, m_j=1, T=1: depends on the azimuth through m_l = 1
  const std::complex<double> a = pair.getrelwf(a1.r, a1.costh, a1.phi, a2.r, a2.costh, a2.phi, 0, 1, 1, 0, 1, 1, 1);
  const std::complex<double> b = pair.getrelwf(b1.r, b1.costh, b1.phi, b2.r, b2.costh, b2.phi, 0, 1, 1, 0, 1, 1, 1);
  REQUIRE(std::abs(a) > 1e-3);
  CHECK_THAT(std::abs(a - b), WithinAbs(0., 1e-12));
}

TEST_CASE("Pair: the correlated wave function is finite at r1 = r2 and r1 = -r2", "[pair]") {
  // the coordinate transformations divided 0 by 0 there, and u(r)/r of the deuteron at r = 0
  const Pair pair = carbon_pn();
  const double pts[][6] = {{1., 0.3, 1., 1., 0.3, 1.},   // r1 = r2
                           {1., 1., 0., 1., 1., 0.},     // r1 = r2 on the z axis
                           {1., -1., 0., 1., 1., 0.},    // r1 = -r2
                           {0., 1., 0., 0., 1., 0.}};    // both at the origin
  for (const auto &q : pts) {
    for (bool corr : {false, true}) {
      const std::complex<double> wf = pair.getwf(corr, q[0], q[1], q[2], q[3], q[4], q[5]);
      CHECK(std::isfinite(wf.real()));
      CHECK(std::isfinite(wf.imag()));
    }
  }
  // (no continuity check through r1 = r2: the Paris u(r)/r and w(r)/r are numerically noisy below
  // about 0.01 fm, see docs/ROADMAP.md)
}
