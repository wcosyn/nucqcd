// Regression: the plane-wave momentum distribution grid of DistMomDistrGrid (a 1D radial integral for each
// momentum). Until 2026-11 it was computed with the nested Romberg routine rombergerN; the values were recorded
// before replacing it with a 1D numint::cube_adaptive call (ROADMAP step 6), which moved them by less than the
// tolerance below.
// The tolerance (1e-6 relative) is the precision of the radial integral.

#include <test_support.hpp>

#include <Glauber/DistMomDistrGrid.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Tools/Core/constants.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_distmomdistr_pw");
}

TEST_CASE("DistMomDistrGrid: plane-wave momentum distribution grid of carbon s and p shells", "[distmomdistr]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  GlauberGridThick fsi(4, 4, 1, &carbon, 1.E-03, IntegrationMethod::Adaptive, dir);  // only its number of grids is used by the constructor
  for (int shell : {0, 1}) {
    DistMomDistrGrid dist(shell, 400., 8, 1, 1, &fsi, 1.E-03, IntegrationMethod::Adaptive, 20000, 0., dir);
    const std::vector<double> &pw = dist.getRhopwgrid();
    REQUIRE(pw.size() == 9);
    for (int i = 0; i <= 8; i++) {
      std::ostringstream key;
      key << "C_shell" << shell << "_p" << i * 50;
      CHECK(pw[i] >= 0.);
      CHECK_REF(refs, key.str(), pw[i], 1e-6, 1e-30);
    }
    if (shell == 0) CHECK(pw[0] > pw[8]);  // the s wave falls with p (the p wave is zero at p = 0)
  }
}
