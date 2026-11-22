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
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

namespace {
nucqcd_test::RefStore refs("reg_glauber_grid");
}

TEST_CASE("Glauber grid: fill, interpolate, cache round trip (off-axis proton in carbon)", "[glauber]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  // off-axis momentum, in the phi=0 plane, at a bin centre of the grid cache key (cos theta = 0.9;
  // until 2026-11 this was theta = 0.3, which is now computed at theta = 0, see FastParticle::gridBinCentre)
  FastParticle proton(0, 0, 1500., std::acos(0.9), 0., 8., 0., 1., 1., dir);
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

TEST_CASE("Glauber grid cache: NUCQCD_GRID_DIR and recovery from a truncated file", "[glauber][gridcache]") {
  namespace fs = std::filesystem;
  const std::string dir = nucqcd_test::share_dir();
  const fs::path cache = fs::path(dir) / "grid_cache_test" / "glauber_env";
  fs::remove_all(cache);
  setenv("NUCQCD_GRID_DIR", cache.c_str(), 1);
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  FastParticle proton(0, 0, 1200., 0.5, 0., 8., 0., 1., 1., dir);  // not used by the other tests
  std::complex<double> computed;
  std::string fsi_file, fsict_file;
  {
    GlauberGridThick grid(8, 4, 2, &carbon, 1.E-03, 2, dir);
    grid.addParticle(proton);
    grid.fillGrids();
    computed = grid.getFsiCtGridFull_interp3(1.5, 0.3, 1.);
    fsi_file = grid.getFsi_Filename();
    fsict_file = grid.getFsi_Ct_Filename();
  }
  // both files are in the directory from NUCQCD_GRID_DIR, without temporary files, and not in share/grids
  CHECK(fs::path(fsi_file).parent_path() == cache);
  CHECK(fs::exists(fsi_file));
  CHECK(fs::exists(fsict_file));
  int files = 0;
  for ([[maybe_unused]] const auto &e : fs::directory_iterator(cache)) files++;
  CHECK(files == 2);
  CHECK_FALSE(fs::exists(fs::path(dir) / "grids" / fs::path(fsi_file).filename()));

  // a truncated cache file (e.g. from a crash under the old non-atomic writes) is recomputed
  fs::resize_file(fsict_file, fs::file_size(fsict_file) / 2);
  {
    GlauberGridThick grid(8, 4, 2, &carbon, 1.E-03, 2, dir);
    grid.addParticle(proton);
    grid.fillGrids();
    // constructCtGrid() integrates the CT grid on its own instead of together with the FSI grid
    // (constructAllGrids), so the adaptive integrator converges differently: equal within the
    // integration precision 1e-3 (the difference is about 6e-4)
    CHECK(std::abs(grid.getFsiCtGridFull_interp3(1.5, 0.3, 1.) - computed) <= 2e-3 * std::abs(computed));
  }
  CHECK(fs::file_size(fsict_file) == fs::file_size(fsi_file));
  unsetenv("NUCQCD_GRID_DIR");
  fs::remove_all(cache);
}

TEST_CASE("Glauber grid cache: grids are computed at the bin centre, independent of history", "[glauber][gridcache]") {
  namespace fs = std::filesystem;
  const std::string dir = nucqcd_test::share_dir();
  const fs::path cache1 = fs::path(dir) / "grid_cache_test" / "bins1";
  const fs::path cache2 = fs::path(dir) / "grid_cache_test" / "bins2";
  fs::remove_all(cache1);
  fs::remove_all(cache2);
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  // two protons in the same bin: p 1500 MeV, cos theta 0.9, phi 0.7, hard scale 8 GeV^2
  FastParticle first(0, 0, 1497., 0.44, 0.68, 7.98, 0., 1., 1., dir);
  FastParticle second(0, 0, 1504., 0.47, 0.72, 8.03, 0., 1., 1., dir);
  auto compute = [&](FastParticle &particle, const fs::path &cache, std::string &file) {
    setenv("NUCQCD_GRID_DIR", cache.c_str(), 1);
    GlauberGridThick grid(8, 4, 2, &carbon, 1.E-03, 2, dir);
    grid.addParticle(particle);
    grid.fillGrids();
    file = fs::path(grid.getFsi_Ct_Filename()).filename().string();
    return std::make_pair(grid.getFsiGridFull_interp3(1.5, 0.3, 1.), grid.getFsiCtGridFull_interp3(1.5, 0.3, 1.));
  };
  std::string file1, file2, file3;
  const auto v1 = compute(first, cache1, file1);   // computed
  const auto v2 = compute(second, cache2, file2);  // computed in an empty cache
  const auto v3 = compute(second, cache1, file3);  // read from the grid of the first particle
  unsetenv("NUCQCD_GRID_DIR");
  CHECK(file1 == file2);
  CHECK(file1 == file3);
  CHECK(file1.find(".key2.") != std::string::npos);
  CHECK(v1 == v2);  // the same grid whichever particle of the bin computed it
  CHECK(v1 == v3);
  fs::remove_all(cache1);
  fs::remove_all(cache2);
}
