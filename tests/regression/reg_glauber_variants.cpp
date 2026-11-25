// Regression [slow]: the Glauber grid classes that reg_glauber_grid and reg_glauber_points do not
// cover: GlauberGrid (plain thin, off axis only), GlauberGridThick with the SRC grids, GlauberDecayGridThick (all
// eight grids; the rho and the rho-decay pion pair), OneGlauberGrid and ROMEAGrid. Tiny grids (a
// handful of nodes), on the z axis and off axis. Safety net for the deduplication of the grid
// hierarchy (ROADMAP step 6): the values must not change.
#include <test_support.hpp>

#include <Glauber/FastParticle.hpp>
#include <Glauber/GlauberDecayGridThick.hpp>
#include <Glauber/GlauberGrid.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Glauber/OneGlauberGrid.hpp>
#include <Glauber/ROMEAGrid.hpp>

#include <complex>
#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_glauber_variants");

// records grids 0..ngrids-1 of a filled grid at a few points
template <class Grid>
void record(Grid &grid, const std::string &name, int ngrids, double R) {
  const double pts[][3] = {{0.1, -0.9, 1.}, {0.3, -0.5, 0.4}, {0.5, 0.2, 0.}, {0.8, 0.9, 2.5}};
  for (const auto &pt : pts)
    for (int g = 0; g < ngrids; g++) {
      const std::complex<double> v = grid.getFsiGridN_interp3(g, pt[0] * R, pt[1], pt[2]);
      std::ostringstream key;
      key << name << "_g" << g << "_r" << pt[0] << "_cth" << pt[1] << "_phi" << pt[2];
      CHECK_REF(refs, key.str() + "_re", v.real(), 1e-6, 1e-9);
      CHECK_REF(refs, key.str() + "_im", v.imag(), 1e-6, 1e-9);
    }
}
}  // namespace

TEST_CASE("GlauberGrid (thin nucleus), carbon, off axis", "[glauber][variants]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleus carbon(MeanFieldNucleus::C, dir);
  // GlauberGrid only handles particles off the z axis (it asks for OneGlauberGrid otherwise)
  FastParticle proton(0, 0, 1500., 0.5, 0., 8., 0., 1., 1., dir);
  GlauberGrid grid(4, 4, 3, &carbon, 1.E-03, 2, dir);
  grid.addParticle(proton);
  grid.fillGrids();
  record(grid, "GG_off", 2, carbon.getRange());
}

TEST_CASE("GlauberGridThick with SRC grids, carbon, off axis", "[glauber][variants]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  FastParticle proton(0, 0, 1500., 0.5, 0., 8., 0., 1., 1., dir);
  GlauberGridThick grid(4, 4, 3, &carbon, 1.E-03, 2, dir);
  grid.addParticle(proton);
  grid.fillGrids();
  record(grid, "GGT_off", 4, carbon.getRange());
}

TEST_CASE("GlauberDecayGridThick: rho and rho->pi pi along z and off axis", "[glauber][variants]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  // type 4: rho0, type 7: the pion pair of a rho decay; Gamma [MeV] as in RhoTCross
  for (int type : {4, 7})
    for (double theta : {0., 0.5}) {
      FastParticle rho(type, 0, 2000., theta, 0., 2., 145., 1., 1., dir);
      GlauberDecayGridThick grid(4, 4, 3, &carbon, 1.E-03, 2, dir);
      grid.addParticle(rho);
      grid.fillGrids();
      std::ostringstream name;
      name << "GDT_type" << type << (theta == 0. ? "_z" : "_off");
      record(grid, name.str(), 8, carbon.getRange());
    }
}

TEST_CASE("GlauberDecayGridThick: proton (no decay) reproduces the plain thick grid", "[glauber][variants]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  FastParticle proton(0, 0, 1500., 0., 0., 8., 0., 1., 1., dir);
  GlauberDecayGridThick dgrid(4, 4, 3, &carbon, 1.E-03, 2, dir);
  dgrid.addParticle(proton);
  dgrid.fillGrids();
  record(dgrid, "GDT_proton", 8, carbon.getRange());
  GlauberGridThick grid(4, 4, 3, &carbon, 1.E-03, 2, dir);
  grid.addParticle(proton);
  grid.fillGrids();
  const double R = carbon.getRange();
  for (int g = 0; g < 4; g++) {
    // without decay width the decay grids equal the others, and grid g of the thick grid is
    // grid g of the decay grid (indices 0-3: no decay)
    const std::complex<double> a = dgrid.getFsiGridN_interp3(g, 0.5 * R, 0.2, 0.);
    const std::complex<double> b = grid.getFsiGridN_interp3(g, 0.5 * R, 0.2, 0.);
    CHECK(std::abs(a - b) <= 2e-3 * std::abs(b));  // different class, same physics (integration precision)
  }
}

TEST_CASE("OneGlauberGrid and ROMEAGrid, carbon, along z", "[glauber][variants]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleus carbon(MeanFieldNucleus::C, dir);
  FastParticle proton(0, 0, 1500., 0., 0., 8., 0., 1., 1., dir);
  {
    OneGlauberGrid grid(4, 4, &carbon, 1.E-03, 2, dir);
    grid.addParticle(proton);
    grid.fillGrids();
    record(grid, "One_z", 2, carbon.getRange());
  }
  {
    ROMEAGrid grid(4, 4, 3, &carbon, 1.E-03, 2, ROMEAGrid::EDAIC12, dir);
    grid.addParticle(proton);
    grid.fillGrids();
    record(grid, "ROMEA_z", 1, carbon.getRange());
  }
}
