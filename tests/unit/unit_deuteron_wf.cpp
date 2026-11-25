// MePhys: deuteron wave functions (TDeuteron, OldDeuteron). Normalization, D-state probability,
// quadrupole moment, spin structure; logic from progs/poldeut/wftest.cpp.
#include <test_support.hpp>

#include <Tools/WaveFunctions/OldDeuteron.hpp>
#include <Tools/WaveFunctions/TDeuteron.hpp>
#include <Tools/WaveFunctions/TInterpolatingWavefunction.hpp>
#include <Tools/Core/constants.hpp>

#include <Tools/Core/Vector3.hpp>

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
using namespace nucqcd;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
// composite Simpson rule on [a,b] with n (even) intervals
double simpson(const std::function<double(double)> &f, double a, double b, int n) {
  const double h = (b - a) / n;
  double s = f(a) + f(b);
  for (int i = 1; i < n; i++) s += f(a + i * h) * (i % 2 ? 4. : 2.);
  return s * h / 3.;
}

std::unique_ptr<TDeuteron::Wavefunction> make_wf(const std::string &name) {
  return std::unique_ptr<TDeuteron::Wavefunction>(TDeuteron::Wavefunction::CreateWavefunction(name));
}

struct WfData {
  const char *name;
  double PD;  // D-state probability [%]
  double Qd;  // quadrupole moment [fm^2]
};
// published values: AV18 (Wiringa et al. PRC51, 38), CD-Bonn (Machleidt PRC63, 024001), Paris (Lacombe et al.)
const WfData wfs[] = {{"AV18", 5.76, 0.270}, {"CDBonn", 4.85, 0.270}, {"Paris", 5.77, 0.279}};
}  // namespace

TEST_CASE("TDeuteron wave functions are normalized in r and p space", "[deuteron]") {
  for (const auto &d : wfs) {
    INFO("wave function " << d.name);
    auto wf = make_wf(d.name);
    REQUIRE(wf);
    const double pnorm = simpson([&](double k) { return k * k * (std::pow(wf->GetUp(k), 2) + std::pow(wf->GetWp(k), 2)); },
                                 0., 6000., 6000);
    const double rnorm = simpson([&](double r) { return std::pow(wf->GetUr(r), 2) + std::pow(wf->GetWr(r), 2); },
                                 0., 40., 8000);
    CHECK_THAT(pnorm, WithinAbs(1., 2e-3));
    CHECK_THAT(rnorm, WithinAbs(1., 2e-3));
  }
}

TEST_CASE("D-state probability and quadrupole moment match the published values", "[deuteron]") {
  for (const auto &d : wfs) {
    INFO("wave function " << d.name);
    auto wf = make_wf(d.name);
    const double PD = simpson([&](double k) { return k * k * std::pow(wf->GetWp(k), 2); }, 0., 6000., 6000);
    const double PDr = simpson([&](double r) { return std::pow(wf->GetWr(r), 2); }, 0., 40., 8000);
    // Q_d = 1/20 int r^2 w (sqrt(8) u - w) dr  (Machleidt nucl-th/0006014 Eq. D16)
    const double Qd = simpson([&](double r) {
      return r * r * (std::sqrt(8.) * wf->GetUr(r) * wf->GetWr(r) - std::pow(wf->GetWr(r), 2));
    }, 0., 40., 8000) / 20.;
    CHECK_THAT(100. * PD, WithinAbs(d.PD, 0.02));
    CHECK_THAT(100. * PDr, WithinAbs(d.PD, 0.02));
    CHECK_THAT(Qd, WithinAbs(d.Qd, 0.003));
  }
}

TEST_CASE("OldDeuteron wave functions are normalized", "[deuteron]") {
  for (const char *name : {"AV18", "AV18b", "CDBonn", "Paris"}) {
    INFO("wave function " << name);
    OldDeuteron wf(name);
    const double norm = 4. * PI * simpson([&](double k) { return k * k * (std::pow(wf.U(k), 2) + std::pow(wf.W(k), 2)); },
                                          0., 1000., 5000);
    CHECK_THAT(norm, WithinAbs(1., 5e-3));
  }
}

TEST_CASE("DeuteronPState: spin structure against closed forms", "[deuteron]") {
  auto wf = make_wf("AV18");
  const double k = 100., theta = 1., phi = 2.;
  const nucqcd::Vector3 p(k * std::sin(theta) * std::cos(phi), k * std::sin(theta) * std::sin(phi), k * std::cos(theta));
  // M=0, both nucleon spins up: pure D-wave (m_l=-1), -W(k) 3/2 cos sin e^{-i phi}/sqrt(4 pi)   [wftest.cpp]
  const std::complex<double> expected =
      -wf->GetWp(k) * 1.5 * std::cos(theta) * std::sin(theta) * std::exp(std::complex<double>(0., -phi)) /
      std::sqrt(4. * PI);
  const std::complex<double> got = wf->DeuteronPState(0, 1, 1, p);
  CHECK_THAT(got.real(), WithinAbs(expected.real(), 1e-12 * std::abs(expected) + 1e-300));
  CHECK_THAT(got.imag(), WithinAbs(expected.imag(), 1e-12 * std::abs(expected) + 1e-300));

  // momentum density summed over nucleon spins, averaged over M: (U^2+W^2)/(4 pi)
  double unpol = 0., tensor = 0.;
  for (int M = -2; M <= 2; M += 2)
    for (int s1 = -1; s1 <= 1; s1 += 2)
      for (int s2 = -1; s2 <= 1; s2 += 2) {
        const double n = std::norm(wf->DeuteronPState(M, s2, s1, p));
        unpol += n / 3.;
        tensor += (M == 0 ? -2. : 1.) * n / 3.;
      }
  const double U = wf->GetUp(k), W = wf->GetWp(k);
  CHECK_THAT(unpol, WithinRel((U * U + W * W) / (4. * PI), 1e-10));
  // tensor-polarized density [wftest.cpp]
  CHECK_THAT(tensor, WithinRel(-(U * W / std::sqrt(2.) + W * W / 4.) * (3. * std::cos(2. * theta) + 1.) / (4. * PI), 1e-10));
}

TEST_CASE("Clebsch-Gordan coefficients and spherical harmonics", "[deuteron]") {
  using WF = TDeuteron::Wavefunction;
  CHECK_THAT(WF::ClebschGordan(4, 4, 1, 1, 5, 5), WithinAbs(1., 1e-14));
  CHECK_THAT(WF::ClebschGordan(1, 1, 1, -1, 0, 0), WithinAbs(1. / std::sqrt(2.), 1e-14));
  CHECK_THAT(WF::ClebschGordan(1, -1, 1, 1, 0, 0), WithinAbs(-1. / std::sqrt(2.), 1e-14));
  // orthonormality: sum_{m1,m2} <j1 m1 j2 m2|J M><j1 m1 j2 m2|J' M> = delta_{J J'}, j1=2, j2=1/2
  for (int twoM = -3; twoM <= 3; twoM += 2)
    for (int twoJ : {3, 5})
      for (int twoJp : {3, 5}) {
        if (std::abs(twoM) > std::min(twoJ, twoJp)) continue;
        double sum = 0.;
        for (int m2 = -1; m2 <= 1; m2 += 2) {
          const int m1 = twoM - m2;
          if (std::abs(m1) > 4) continue;
          sum += WF::ClebschGordan(4, m1, 1, m2, twoJ, twoM) * WF::ClebschGordan(4, m1, 1, m2, twoJp, twoM);
        }
        CHECK_THAT(sum, WithinAbs(twoJ == twoJp ? 1. : 0., 1e-13));
      }
  const double th = 0.8, ph = 1.3;
  CHECK_THAT(std::abs(WF::SphericalHarmonic(0, 0, th, ph)), WithinRel(1. / std::sqrt(4. * PI), 1e-13));
  CHECK_THAT(WF::SphericalHarmonic(1, 0, th, ph).real(), WithinRel(std::sqrt(3. / (4. * PI)) * std::cos(th), 1e-13));
  CHECK_THAT(std::abs(WF::SphericalHarmonic(2, 2, th, ph)),
             WithinRel(0.25 * std::sqrt(15. / (2. * PI)) * std::pow(std::sin(th), 2), 1e-13));
}

TEST_CASE("D-wave conventions: sign of w(r) and W(p), and the Bessel transform between them", "[deuteron]") {
  // documented in TDeuteron::Wavefunction: u and w have the same sign at large r, W(p) > 0 at low p,
  // and U, W are the plain Bessel transforms (no factor (-i)^L, which DeuteronPState adds itself)
  for (const std::string name : {"Paris", "AV18", "CDBonn"}) {
    auto wf = make_wf(name);
    CHECK(wf->GetUr(3.) > 0.);
    CHECK(wf->GetWr(3.) > 0.);
    CHECK(wf->GetWp(100.) > 0.);
    for (double p : {50., 150., 300.}) {
      const double k = p / nucqcd::hbarc;
      const double f = std::sqrt(2. / PI) / std::pow(nucqcd::hbarc, 1.5);
      const double U = f * simpson([&](double r) { return wf->GetUr(r) * std::sph_bessel(0, k * r) * r; }, 0., 40., 40000);
      const double W = f * simpson([&](double r) { return wf->GetWr(r) * std::sph_bessel(2, k * r) * r; }, 0., 40., 40000);
      CHECK_THAT(U, WithinRel(wf->GetUp(p), 2e-3));
      CHECK_THAT(W, WithinRel(wf->GetWp(p), 2e-3));
    }
    // DeuteronPState: the D wave enters with a minus sign, DeuteronRState with a plus sign
    // (M=0, both spins up is pure D wave with m_L=-1)
    const nucqcd::Vector3 dir(std::sin(1.) * std::cos(2.), std::sin(1.) * std::sin(2.), std::cos(1.));
    const std::complex<double> Y2m1 = TDeuteron::Wavefunction::SphericalHarmonicCos(2, -1, std::cos(1.), 2.);
    const double cg = TDeuteron::Wavefunction::ClebschGordan(4, -2, 2, 2, 2, 0);
    CHECK(std::abs(wf->DeuteronPState(0, 1, 1, 150. * dir) - (-wf->GetWp(150.) * cg * Y2m1)) < 1e-12 * wf->GetWp(150.));
    CHECK(std::abs(wf->DeuteronRState(0, 1, 1, 2. * dir) - (wf->GetWr(2.) * cg * Y2m1)) < 1e-12 * wf->GetWr(2.));
  }
}

TEST_CASE("TDeuteron::Wavefunction: errors are exceptions", "[deuteron][errors]") {
  // until 2026-11 an unknown name returned a null pointer (with a warning), and getResidu() of a
  // non-Yukawa wave function had no return value when assertions were off
  CHECK_THROWS_AS(TDeuteron::Wavefunction::CreateWavefunction("Nijmegen"), std::invalid_argument);
  CHECK_NOTHROW(make_wf("nijm3"));
  CHECK_THROWS_AS(make_wf("nijm3")->getResidu(), std::logic_error);
  CHECK(make_wf("AV18")->getResidu() > 0.);
}

TEST_CASE("TInterpolatingWavefunction states have the TDeuteron conventions", "[deuteron]") {
  // Deut_Conv_EMTLF and Poldeut use an interpolated copy of a TDeuteron wave function, and
  // TInterpolatingWavefunction has its own DeuteronPState/DeuteronRState: same signs and couplings
  auto wf = make_wf("AV18");
  TInterpolatingWavefunction interp;
  for (int i = 0; i <= 1000; i++) {
    interp.AddUp(i, wf->GetUp(i));
    interp.AddWp(i, wf->GetWp(i));
    interp.AddUr(i * 0.02, wf->GetUr(i * 0.02));
    interp.AddWr(i * 0.02, wf->GetWr(i * 0.02));
  }
  const nucqcd::Vector3 dir(std::sin(1.) * std::cos(2.), std::sin(1.) * std::sin(2.), std::cos(1.));
  const nucqcd::Vector3 p = 150. * dir, r = 2. * dir;  // on grid points of the tables, so no interpolation error
  for (int M = -2; M <= 2; M += 2)
    for (int s1 = -1; s1 <= 1; s1 += 2)
      for (int s2 = -1; s2 <= 1; s2 += 2) {
        CHECK(std::abs(interp.DeuteronPState(M, s2, s1, p) - wf->DeuteronPState(M, s2, s1, p)) < 1e-12 * wf->GetUp(150.));
        CHECK(std::abs(interp.DeuteronRState(M, s2, s1, r) - wf->DeuteronRState(M, s2, s1, r)) < 1e-12 * wf->GetUr(2.));
      }
}
