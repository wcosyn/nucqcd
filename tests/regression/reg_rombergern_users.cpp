// Regression: the classes that still integrate with the nested 1D Romberg routine rombergerN
// (ROADMAP step 6, "Remove the integrators rombergerN and cube_romb?"), with the values it gives, so that
// replacing it can be measured:
// - the plane-wave momentum distribution grid of DistMomDistrGrid (1D radial integral),
// - RhoDeuteron: the coherent rho cross section (nested complex 3D integral over the spectator momentum).
// (InclusiveCrossRes, the third user, is not built: it does not compile, see the ROADMAP.)
// The nested routine has no accuracy control that can be stated (see docs/integrator_benchmark.md), so
// the tolerances are loose where the integrals are expected to move: they must stay below the
// precision the classes ask for.
#include <test_support.hpp>

#include <Glauber/DistMomDistrGrid.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Knockout/RhoDeuteron.hpp>
#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/Core/constants.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_rombergern_users");
}

TEST_CASE("DistMomDistrGrid: plane-wave momentum distribution grid of carbon s and p shells", "[rombergern][distmomdistr]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  GlauberGridThick fsi(4, 4, 1, &carbon, 1.E-03, 2, dir);  // only its number of grids is used by the constructor
  for (int shell : {0, 1}) {
    DistMomDistrGrid dist(shell, 400., 8, 1, 1, &fsi, 1.E-03, 2, 20000, 0., dir);
    const double *pw = dist.getRhopwgrid();
    REQUIRE(pw != nullptr);
    for (int i = 0; i <= 8; i++) {
      std::ostringstream key;
      key << "C_shell" << shell << "_p" << i * 50;
      CHECK(pw[i] >= 0.);
      CHECK_REF(refs, key.str(), pw[i], 1e-6, 1e-30);
    }
    if (shell == 0) CHECK(pw[0] > pw[8]);  // the s wave falls with p (the p wave is zero at p = 0)
  }
}

TEST_CASE("RhoDeuteron: coherent rho cross section at fixed z and at fixed t", "[rombergern][rhodeuteron]") {
  // Paris wave function, no experimental cuts, pm up to 400 MeV, no FSI (the coherent part does not use it)
  RhoDeuteron rho("paris", 400., true, 2, 1.E-03, 20000, false);
  const double Ebeam = 5., Q2 = 2.;
  int nonzero = 0;
  for (double nu : {3.5, 4.0}) {
    for (double z : {0.97, 0.98, 0.99}) {
      double result = 0.;
      rho.getCrossz_coh(&result, Ebeam, Q2, nu, z);
      std::ostringstream key;
      key << "z_nu" << nu << "_z" << z;
      if (result != 0.) nonzero++;
      CHECK_REF(refs, key.str(), result, 1e-4, 1e-30);
    }
    for (double t : {-0.15, -0.25, -0.35}) {
      double result = 0.;
      rho.getCrosst_coh(&result, Ebeam, Q2, nu, t);
      std::ostringstream key;
      key << "t_nu" << nu << "_t" << t;
      if (result != 0.) nonzero++;
      CHECK_REF(refs, key.str(), result, 1e-4, 1e-30);
    }
  }
  CHECK(nonzero >= 4);  // the kinematic points must give something to compare
}
