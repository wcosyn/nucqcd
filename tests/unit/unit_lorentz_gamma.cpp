// MePhys: FourVector algebra and Dirac (gamma-matrix) algebra; logic from progs/attic/gammatest/gamma.cpp.
#include <test_support.hpp>

#include <Tools/Core/FourVector.hpp>
#include <Tools/Core/GammaStructure.hpp>
#include <Tools/Core/Matrix.hpp>

#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/LorentzVector.hpp>

#include <cmath>
#include <complex>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using cd = std::complex<double>;

namespace {
const FourVector<GammaStructure> gamma_mu(GammaStructure(0., 0., 1.), GammaStructure(0., 0., 0., 1.),
                                          GammaStructure(0., 0., 0., 0., 1.), GammaStructure(0., 0., 0., 0., 0., 1.));
const GammaStructure gamma_5(0., 1.);
const GammaStructure unity(1.);

double metric(int mu, int nu) { return mu != nu ? 0. : (mu == 0 ? 1. : -1.); }

// max |A_ij - B_ij|
double max_diff(const Matrix<4, 4> &a, const Matrix<4, 4> &b) {
  double d = 0.;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) d = std::max(d, std::abs(a(i, j) - b(i, j)));
  return d;
}
}  // namespace

TEST_CASE("FourVector Minkowski product uses the (+,-,-,-) metric", "[fourvector]") {
  const double m = 938.272, px = 120., py = -45., pz = 800.;
  const double E = std::sqrt(m * m + px * px + py * py + pz * pz);
  FourVector<double> p(E, px, py, pz);
  CHECK_THAT(p * p, WithinRel(m * m, 1e-12));
  FourVector<double> q(3., 1., 2., 0.5);
  CHECK_THAT(p * q, WithinRel(3. * E - px - 2. * py - 0.5 * pz, 1e-14));
  CHECK_THAT((p + q) * (p + q), WithinRel(p * p + 2. * (p * q) + q * q, 1e-12));
}

TEST_CASE("FourVector products are Lorentz invariant", "[fourvector]") {
  FourVector<double> p(1500., 300., -200., 900.), q(700., -100., 50., 400.);
  nucqcd::LorentzRotation L;
  L.Boost(0.3, -0.2, 0.5);
  L.RotateY(0.7);
  const FourVector<double> pp = L * p, qq = L * q;
  CHECK_THAT(pp * qq, WithinRel(p * q, 1e-12));
  CHECK_THAT(pp * pp, WithinRel(p * p, 1e-12));
  // round trip through the ROOT conversion helpers
  const FourVector<double> back = ToFourVector(ToLorentzVector(p));
  CHECK(back[0] == p[0]);
  CHECK(back[3] == p[3]);
}

TEST_CASE("Clifford algebra {gamma^mu, gamma^nu} = 2 g^{mu nu}", "[gamma]") {
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++) {
      const GammaStructure anti = gamma_mu[mu] * gamma_mu[nu] + gamma_mu[nu] * gamma_mu[mu];
      CHECK(max_diff(anti.value(), (unity * (2. * metric(mu, nu))).value()) < 1e-14);
    }
}

TEST_CASE("gamma5 identities and traces", "[gamma]") {
  // gamma5 = i gamma0 gamma1 gamma2 gamma3, gamma5^2 = 1
  const GammaStructure prod = gamma_mu[0] * gamma_mu[1] * gamma_mu[2] * gamma_mu[3];
  CHECK(max_diff((cd(0., 1.) * prod).value(), gamma_5.value()) < 1e-14);
  CHECK(max_diff((gamma_5 * gamma_5).value(), unity.value()) < 1e-14);
  // Tr(gamma5 gamma0 gamma1 gamma2 gamma3) = -4i
  const cd tr5 = Trace((gamma_5 * prod).value());
  CHECK_THAT(tr5.real(), WithinAbs(0., 1e-14));
  CHECK_THAT(tr5.imag(), WithinAbs(-4., 1e-14));
  // traceless gamma^mu, gamma5 anticommutes with gamma^mu
  for (int mu = 0; mu < 4; mu++) {
    CHECK(std::abs(Trace(gamma_mu[mu].value())) < 1e-14);
    CHECK(max_diff((gamma_5 * gamma_mu[mu] + gamma_mu[mu] * gamma_5).value(), GammaStructure(0.).value()) < 1e-14);
  }
  // Tr(gamma^mu gamma^nu) = 4 g^{mu nu}
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      CHECK_THAT(Trace((gamma_mu[mu] * gamma_mu[nu]).value()).real(), WithinAbs(4. * metric(mu, nu), 1e-14));
}

TEST_CASE("pslash pslash = p^2 and Tr(pslash qslash) = 4 p.q", "[gamma]") {
  FourVector<double> p(1500., 300., -200., 900.), q(700., -100., 50., 400.);
  const GammaStructure pslash = p * gamma_mu, qslash = q * gamma_mu;
  CHECK(max_diff((pslash * pslash).value(), (unity * (p * p)).value()) < 1e-8 * (p * p));
  CHECK_THAT(Trace((pslash * qslash).value()).real(), WithinRel(4. * (p * q), 1e-12));
  // Tr(a b c d) = 4[(ab)(cd) - (ac)(bd) + (ad)(bc)]
  FourVector<double> c(10., 1., 2., 3.), d(20., -3., 0.5, 4.);
  const GammaStructure cslash = c * gamma_mu, dslash = d * gamma_mu;
  const double expected = 4. * ((p * q) * (c * d) - (p * c) * (q * d) + (p * d) * (q * c));
  CHECK_THAT(Trace((pslash * qslash * cslash * dslash).value()).real(), WithinRel(expected, 1e-10));
}

TEST_CASE("Photon polarization vectors are complete: sum = -g + q q / Q^2", "[gamma]") {
  // virtual photon along z, from progs/attic/gammatest/gamma.cpp
  const double Q2 = 1.e06, nu = 800.;
  const double qvec = std::sqrt(Q2 + nu * nu);
  FourVector<cd> q(nu, 0., 0., qvec);
  FourVector<cd> eps_plus(0., -1. / std::sqrt(2.), cd(0., -1. / std::sqrt(2.)), 0.);
  FourVector<cd> eps_min(0., 1. / std::sqrt(2.), cd(0., -1. / std::sqrt(2.)), 0.);
  FourVector<cd> eps_0(qvec / std::sqrt(Q2), 0., 0., nu / std::sqrt(Q2));
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      const cd lhs = -eps_min[j] * std::conj(eps_min[i]) + std::conj(eps_0[i]) * eps_0[j] -
                     eps_plus[j] * std::conj(eps_plus[i]);
      const cd rhs = metric(i, j) + q[i] * q[j] / Q2;
      CHECK(std::abs(lhs - rhs) < 1e-12);
    }
}
