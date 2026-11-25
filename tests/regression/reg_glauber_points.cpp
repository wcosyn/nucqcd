// Regression (fast): Glauber FSI factors at a few points. The Glauber phase is only computed on
// grid nodes (GlauberGridThick::calcGlauberphasesBoth is private), so a tiny grid is used: its
// handful of nodes are the "points". Proton along z (in-plane symmetry), carbon.
// The full grid machinery (fill, interpolation, cache file round trip) is in reg_glauber_grid [slow].
#include <test_support.hpp>

#include <Glauber/FastParticle.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>

#include <cmath>
#include <complex>
#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_glauber_points");
}

TEST_CASE("Glauber and CT FSI factors for a 1.5 GeV proton in carbon (tiny grid)", "[glauber]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  // FastParticle(type, inc, p, theta, phi, hard_scale, Gamma, lc_mod, nkt_mod, dir); type 0 = proton
  FastParticle proton(0, 0, 1500., 0., 0., 8., 0., 1., 1., dir);
  CHECK_REF(refs, "proton_sigmap", proton.getSigmap(), 1e-10);
  CHECK_REF(refs, "proton_sigman", proton.getSigman(), 1e-10);
  GlauberGridThick grid(4, 4, 1, &carbon, 1.E-03, 2, dir);
  grid.addParticle(proton);
  grid.fillGrids();
  const double R = carbon.getRange();
  for (double r : {0., 0.25 * R, 0.5 * R})
    for (double cth : {-1., 0., 1.}) {
      const std::complex<double> fsi = grid.getFsiGridFull_interp3(r, cth, 0.);
      const std::complex<double> fsict = grid.getFsiCtGridFull_interp3(r, cth, 0.);
      CHECK(std::abs(fsi) <= 1. + 1e-6);         // absorption only
      CHECK(std::abs(fsict) >= std::abs(fsi) - 1e-6);  // color transparency reduces the absorption
      std::ostringstream key;
      key << "C_r" << r / R << "R_cth" << cth;
      CHECK_REF(refs, key.str() + "_fsi_re", fsi.real(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsi_im", fsi.imag(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsict_re", fsict.real(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsict_im", fsict.imag(), 1e-6, 1e-9);
    }
  // a proton created at the edge of the nucleus and leaving it forward hardly interacts
  CHECK(std::abs(grid.getFsiGridFull_interp3(R, 1., 0.)) > 0.99);
}
