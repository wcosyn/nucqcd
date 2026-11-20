// Regression [slow]: the full Glauber grid machinery (Glauber/AbstractFsiGrid*, GlauberGridThick):
// filling a grid for an off-axis proton, interpolation, writing the grid cache file and reading it
// back, and RMSGA A(e,e'p) through Knockout/Cross with its default 120x36x5 grids. Grids are
// computed from scratch in the sandbox share dir (its grids/ is emptied before every ctest run).
#include <test_support.hpp>

#include <Cross.hpp>
#include <FastParticle.hpp>
#include <GlauberGridThick.hpp>
#include <MeanFieldNucleusThick.hpp>
#include <TElectronKinematics.h>
#include <TKinematics2to2.h>
#include <constants.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

namespace {
nucqcd_test::RefStore refs("reg_glauber_grid");
}

TEST_CASE("Glauber grid: fill, interpolate, cache round trip (off-axis proton in carbon)", "[glauber]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  // off-axis momentum, in the phi=0 plane
  FastParticle proton(0, 0, 1500., 0.3, 0., 8., 0., 1., 1., dir);
  const double pts[][3] = {{0.5, -0.8, 0.4}, {2., 0., 3.}, {4., 0.9, 1.}, {1., 0.5, 5.5}};

  std::complex<double> computed[4][2];
  {
    GlauberGridThick grid(12, 6, 3, &carbon, 1.E-03, 2, dir);
    grid.addParticle(proton);
    grid.fillGrids();  // computes and writes the cache file
    for (int i = 0; i < 4; i++) {
      computed[i][0] = grid.getFsiGridFull_interp3(pts[i][0], pts[i][1], pts[i][2]);
      computed[i][1] = grid.getFsiCtGridFull_interp3(pts[i][0], pts[i][1], pts[i][2]);
      CHECK(std::abs(computed[i][0]) <= 1. + 1e-6);
      std::ostringstream key;
      key << "C_r" << pts[i][0] << "_cth" << pts[i][1] << "_phi" << pts[i][2];
      CHECK_REF(refs, key.str() + "_fsi_re", computed[i][0].real(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsi_im", computed[i][0].imag(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsict_re", computed[i][1].real(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_fsict_im", computed[i][1].imag(), 1e-6, 1e-9);
    }
  }
  // the grid files were written to the sandbox
  int files = 0;
  for (const auto &entry : std::filesystem::directory_iterator(std::filesystem::path(dir) / "grids"))
    if (entry.path().filename().string().rfind("Gl.Thick.C.", 0) == 0) files++;
  CHECK(files >= 1);
  // a new grid object with the same parameters reads the cache and reproduces the values exactly
  GlauberGridThick reread(12, 6, 3, &carbon, 1.E-03, 2, dir);
  reread.addParticle(proton);
  reread.fillGrids();
  for (int i = 0; i < 4; i++) {
    CHECK(reread.getFsiGridFull_interp3(pts[i][0], pts[i][1], pts[i][2]) == computed[i][0]);
    CHECK(reread.getFsiCtGridFull_interp3(pts[i][0], pts[i][1], pts[i][2]) == computed[i][1]);
  }
}

TEST_CASE("RMSGA A(e,e'p) cross section for carbon (default 120x36x5 grids)", "[glauber][eep]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  std::unique_ptr<TElectronKinematics> elec(TElectronKinematics::CreateWithBeamEnergy(3245.));
  Cross obs(*elec, &carbon, 1.E-03, 2, dir, false, 1.);
  TKinematics2to2 kin("", "", carbon.getMassA(), carbon.getMassA_min_proton(), MASSP, "qsquared:wlab:pklab", 1.8E06,
                      1000., 150.);
  for (int shell = 0; shell < carbon.getPLevels(); shell++) {
    std::ostringstream key;
    key << "C_RMSGA_shell" << shell;
    const double pw = obs.getDiffCross(kin, 2, 0, 0, 0, 1, shell, 0., 20000, true, false, 1., 1.);
    const double rmsga = obs.getDiffCross(kin, 2, 1, 0, 0, 0, shell, 0., 20000, true, false, 1., 1.);
    CHECK(rmsga < pw);  // FSI attenuates
    CHECK_REF(refs, key.str(), rmsga, 1e-5, 1e-30);
    CHECK_REF(refs, key.str() + "_transparency", rmsga / pw, 1e-5);
  }
}
