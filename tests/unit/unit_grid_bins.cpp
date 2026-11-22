// FastParticle::gridBinCentre: the FSI grids are computed for the particle moved to the centre of
// its grid-cache bin (AbstractFsiGrid::addParticle), so a cached grid does not depend on which
// kinematics in its bin came first. Bins: 10 MeV in p, 0.1 in cos(theta) and phi, 0.1 GeV^2 in the
// hard scale, 0.01 in lc_mod and nkt_mod.
#include <test_support.hpp>

#include <FastParticle.hpp>
#include <constants.hpp>

#include <cmath>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
// FastParticle(type, inc, p, theta, phi, hard_scale, Gamma, lc_mod, nkt_mod, dir)
FastParticle proton(double p, double theta, double phi, double hs = 8., double lc_mod = 1., double nkt_mod = 1.) {
  return FastParticle(0, 0, p, theta, phi, hs, 0., lc_mod, nkt_mod, nucqcd_test::share_dir());
}
}  // namespace

TEST_CASE("gridBinCentre leaves a particle at a bin centre unchanged", "[gridbins]") {
  const FastParticle a = proton(1500., std::acos(0.9), 0.);
  const FastParticle c = a.gridBinCentre();
  CHECK(c.getP() == 1500.);
  CHECK_THAT(c.getCosTheta(), WithinAbs(0.9, 1e-15));
  CHECK(c.getPhi() == 0.);
  CHECK(c.getHardScale() == 8.);
  CHECK(c.get_lc_mod() == 1.);
  CHECK(c.get_nkt_mod() == 1.);
  CHECK(c.getSigmap() == a.getSigmap());
  CHECK(c.getBeta2p() == a.getBeta2p());
  CHECK(c.getEpsilonn() == a.getEpsilonn());
  CHECK(c.getParticletype() == 0);
  CHECK_FALSE(c.getIncoming());
  CHECK_FALSE(c.getUserset());
}

TEST_CASE("gridBinCentre moves the kinematics to the bin centre", "[gridbins]") {
  CHECK(proton(1503., 0., 0.).gridBinCentre().getP() == 1500.);
  CHECK(proton(1506., 0., 0.).gridBinCentre().getP() == 1510.);
  CHECK(proton(3., 0., 0.).gridBinCentre().getP() == 10.);  // at least 10 MeV
  // the scattering parameters follow from the binned momentum
  CHECK(proton(1506., 0., 0.).gridBinCentre().getSigmap() == proton(1510., 0., 0.).getSigmap());
  // cos(theta) to 0.1: theta = 0.3 rad (cos 0.955) is computed along the z axis
  // exactly 1, so the faster single-particle code for cos(theta) == 1 is used
  CHECK(proton(1500., 0.3, 0.).gridBinCentre().getCosTheta() == 1.);
  CHECK(proton(1500., 0.3, 0.).gridBinCentre().getTheta() == 0.);
  CHECK(proton(1500., 0., 0.).gridBinCentre().getCosTheta() == 1.);
  CHECK(proton(1500., PI - 0.2, 0.).gridBinCentre().getCosTheta() == -1.);
  // along the z axis phi is set to 0 (in-plane symmetry, one cache file)
  CHECK(proton(1500., 0.3, 0.73).gridBinCentre().getPhi() == 0.);
  CHECK(proton(1500., PI - 0.2, -2.).gridBinCentre().getPhi() == 0.);
  CHECK_THAT(proton(1500., 0.47, 0.).gridBinCentre().getCosTheta(), WithinAbs(0.9, 1e-15));
  // phi to 0.1, taken in [-pi,pi]; the edge bin |phi| > 3.05 goes to pi (the phi=pi plane)
  CHECK_THAT(proton(1500., 1., 0.73).gridBinCentre().getPhi(), WithinAbs(0.7, 1e-15));
  CHECK_THAT(proton(1500., 1., -0.73).gridBinCentre().getPhi(), WithinAbs(-0.7, 1e-15));
  CHECK_THAT(proton(1500., 1., 2. * PI + 0.21).gridBinCentre().getPhi(), WithinAbs(0.2, 1e-12));
  CHECK(proton(1500., 1., 3.12).gridBinCentre().getPhi() == PI);
  CHECK(proton(1500., 1., -3.12).gridBinCentre().getPhi() == PI);
  CHECK(proton(1500., 1., 0.02).gridBinCentre().getPhi() == 0.);
  // CT parameters
  CHECK_THAT(proton(1500., 0., 0., 7.96).gridBinCentre().getHardScale(), WithinAbs(8., 1e-15));
  CHECK_THAT(proton(1500., 0., 0., 8., 0.996, 1.234).gridBinCentre().get_lc_mod(), WithinAbs(1., 1e-15));
  CHECK_THAT(proton(1500., 0., 0., 8., 0.996, 1.234).gridBinCentre().get_nkt_mod(), WithinAbs(1.23, 1e-15));
}

TEST_CASE("gridBinCentre: two particles in the same bin give the same particle", "[gridbins]") {
  const FastParticle a = proton(1497., 0.44, 0.68, 7.98).gridBinCentre();
  const FastParticle b = proton(1504., 0.47, 0.72, 8.03).gridBinCentre();
  CHECK(a.getP() == b.getP());
  CHECK(a.getCosTheta() == b.getCosTheta());
  CHECK(a.getPhi() == b.getPhi());
  CHECK(a.getHardScale() == b.getHardScale());
  CHECK(a.getSigmap() == b.getSigmap());
  CHECK(a.getLc() == b.getLc());
}

TEST_CASE("gridBinCentre bins user-set scattering parameters", "[gridbins]") {
  FastParticle a = proton(1500., 0., 0.);
  a.setScatter(40.3, 0.5, -0.27, 1500.);  // sigma [mb], beta [GeV^-2], eps
  const FastParticle c = a.gridBinCentre();
  CHECK(c.getUserset());
  CHECK_THAT(c.getSigmap(), WithinAbs(4.0, 1e-15));  // fm^2, to 0.1 fm^2 = 1 mb
  CHECK_THAT(c.getSigman(), WithinAbs(4.0, 1e-15));
  CHECK_THAT(c.getBeta2p(), WithinAbs(std::round(a.getBeta2p() * 100.) / 100., 1e-15));
  CHECK_THAT(c.getEpsilonp(), WithinAbs(-0.27, 1e-15));
}

TEST_CASE("FastParticle copies keep the decay width", "[gridbins]") {
  // the copy constructor and assignment did not copy Gamma before 2026-11
  const FastParticle rho(4, 0, 3000., 0., 0., 2., 149., 1., 1., nucqcd_test::share_dir());
  REQUIRE(rho.getDecay_dil() > 0.);
  FastParticle copy(rho);
  FastParticle assigned = proton(1000., 0., 0.);
  assigned = rho;
  // decay_dil is recomputed from Gamma at the binned momentum (3000 MeV is a bin centre), to 1 MeV
  CHECK(copy.gridBinCentre().getDecay_dil() == std::round(rho.getDecay_dil()));
  CHECK(assigned.gridBinCentre().getDecay_dil() == std::round(rho.getDecay_dil()));
  CHECK_THAT(rho.gridBinCentre().getDecay_dil(), WithinRel(rho.getDecay_dil(), 0.5 / rho.getDecay_dil()));
}
