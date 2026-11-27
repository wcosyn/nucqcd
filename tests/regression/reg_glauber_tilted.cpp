// Regression [slow]: the integrals of the thick Glauber grids for a single particle off the z axis are done in the
// coordinates around its momentum (GlauberThickGridCommon::integrateTilted, the default) instead of in the spherical coordinates
// of the lab. Both are the same integral: the grids must agree within the integration precision (the lab version at the
// precision asked for here has errors of a few 1e-4 on these tiny grids and hits the cap on the number of evaluations on
// larger ones). The tilted one is also much faster: it needs a second or less here.
#include <test_support.hpp>

#include <Glauber/FastParticle.hpp>
#include <Glauber/GlauberDecayGridThick.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>

#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <complex>
#include <filesystem>
#include <string>
#include <vector>
using namespace nucqcd;

namespace {
template <class Grid>
std::vector<std::complex<double> > values(Grid &grid, int ngrids, double R) {
  const double pts[][3] = {{0.1, -0.9, 1.}, {0.3, -0.5, 0.4}, {0.5, 0.2, 0.}, {0.8, 0.9, 2.5}};
  std::vector<std::complex<double> > v;
  for (const auto &pt : pts)
    for (int g = 0; g < ngrids; g++) v.push_back(grid.getFsiGridN_interp3(g, pt[0] * R, pt[1], pt[2]));
  return v;
}
void compare(const std::vector<std::complex<double> > &a, const std::vector<std::complex<double> > &b, double tol) {
  REQUIRE(a.size() == b.size());
  for (size_t i = 0; i < a.size(); i++)
    if (std::abs(b[i]) > 1e-2) CHECK(std::abs(a[i] - b[i]) <= tol * std::abs(b[i]));
}
// the file name of the grid cache does not say which coordinates were used: every variant gets its own empty cache directory
// (otherwise the second one would read the grid of the first and the comparison would be empty)
std::string fresh_cache(const std::string &dir, const std::string &name) {
  const std::filesystem::path cache = std::filesystem::path(dir) / "grid_cache_test" / name;
  std::filesystem::remove_all(cache);
  setenv("NUCQCD_GRID_DIR", cache.c_str(), 1);
  return cache.string();
}
void remove_cache(const std::string &cache) {
  unsetenv("NUCQCD_GRID_DIR");
  std::filesystem::remove_all(cache);
}
}  // namespace

TEST_CASE("Tilted frame and lab frame integrals agree: GlauberGridThick, proton", "[glauber][tilted]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  std::vector<std::complex<double> > v[2];
  for (int tilted = 0; tilted < 2; tilted++) {
    const std::string cache = fresh_cache(dir, "tilted_proton" + std::to_string(tilted));
    FastParticle proton(0, 0, 1500., 0.5, 0., 8., 0., 1., 1., dir);
    GlauberGridThick grid(4, 4, 3, &carbon, 1.E-03, IntegrationMethod::Adaptive, dir);
    grid.setUseTiltedFrame(tilted);
    grid.addParticle(proton);
    grid.fillGrids();
    v[tilted] = values(grid, 4, carbon.getRange());
    remove_cache(cache);
  }
  compare(v[1], v[0], 2e-3);
  // computed separately, not identical (this would show that one read the other's cache)
  bool differ = false;
  for (size_t i = 0; i < v[0].size(); i++) differ = differ || v[0][i] != v[1][i];
  CHECK(differ);
}

TEST_CASE("Tilted frame and lab frame integrals agree: GlauberDecayGridThick, rho and its pion pair", "[glauber][tilted]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  for (int type : {4, 7}) {
    std::vector<std::complex<double> > v[2];
    for (int tilted = 0; tilted < 2; tilted++) {
      const std::string cache = fresh_cache(dir, "tilted_rho" + std::to_string(type) + std::to_string(tilted));
      FastParticle rho(type, 0, 2000., 0.5, 0., 2., 145., 1., 1., dir);
      GlauberDecayGridThick grid(4, 4, 3, &carbon, 1.E-03, IntegrationMethod::Adaptive, dir);
      grid.setUseTiltedFrame(tilted);
      grid.addParticle(rho);
      grid.fillGrids();
      v[tilted] = values(grid, 8, carbon.getRange());
      remove_cache(cache);
    }
    compare(v[1], v[0], 2e-3);
  }
}

TEST_CASE("Tilted frame: an incoming particle", "[glauber][tilted]") {
  // an incoming particle scatters on the nucleons it has passed (z < hit z): the integral over the other half of the z range
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  std::vector<std::complex<double> > v[2];
  for (int tilted = 0; tilted < 2; tilted++) {
    const std::string cache = fresh_cache(dir, "tilted_incoming" + std::to_string(tilted));
    FastParticle proton(0, 1, 1500., 0.5, 0., 8., 0., 1., 1., dir);
    GlauberGridThick grid(4, 4, 3, &carbon, 1.E-03, IntegrationMethod::Adaptive, dir);
    grid.setUseTiltedFrame(tilted);
    grid.addParticle(proton);
    grid.fillGrids();
    v[tilted] = values(grid, 4, carbon.getRange());
    remove_cache(cache);
  }
  compare(v[1], v[0], 2e-3);
}
