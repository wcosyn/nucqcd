// Regression: nucleon structure functions (Tools/StructureFunctions/NucleonStructure, mostly Fortran parametrizations).
// Each parametrization is its own TEST_CASE, i.e. its own process (Fortran COMMON/SAVE state).
// Tools/StructureFunctions/NuclStructure is a second wrapper around the same parametrizations; it is the one used by
// DeuteronStructure (DeuteronCross, InclusiveCross, He3Cross), so the two are checked against each other.
#include <test_support.hpp>

#include <Tools/StructureFunctions/LHAPDFStructure.hpp>
#include <Tools/StructureFunctions/NuclStructure.hpp>
#include <Tools/StructureFunctions/NucleonStructure.hpp>
#include <Tools/Core/constants.hpp>

#include <sstream>
#include <string>
#include <vector>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_nucleon_structure");

// "MSTW" is the LHAPDF set MSTW2008lo68cl since 2026-11 (re-baselined then: LHAPDF interpolates the
// MSTW grid differently from the former bundled MSTW code, up to ~1e-4 relative at these points)
void require_mstw() {
  if (!nucqcd::LHAPDFStructure::available() || !nucqcd::LHAPDFStructure::installedSet("MSTW2008lo68cl"))
    SKIP("MSTW needs nucqcd built with LHAPDF and the set MSTW2008lo68cl");
}

// default grid of (x, Q^2[MeV^2]) points; Christy-Bosted only covers W < ~3 GeV (and throws std::out_of_range
// outside its range), so it gets its own grid
const std::vector<double> xs_default = {0.1, 0.3, 0.6}, Q2s_default = {2.e06, 5.e06};
const std::vector<double> xs_CB = {0.3, 0.5, 0.7}, Q2s_CB = {1.e06, 2.e06};

void check_param(const std::string &name, const std::vector<double> &xs = xs_default,
                 const std::vector<double> &Q2s = Q2s_default) {
  NucleonStructure sf(name);
  for (int proton = 0; proton <= 1; proton++)
    for (double x : xs)
      for (double Q2 : Q2s) {
        double F1 = 0., F2 = 0.;
        sf.getF_xQ(F1, F2, proton, x, Q2);
        std::ostringstream key;
        key << name << (proton ? "_p" : "_n") << "_x" << x << "_Q2_" << Q2 * 1e-6;
        CHECK_REF(refs, key.str() + "_F1", F1);
        CHECK_REF(refs, key.str() + "_F2", F2);
      }
}

// NuclStructure (used by the deuteron code) and NucleonStructure must give the same F1 and F2. Since
// 2026-11 NuclStructure is a wrapper around NucleonStructure; before, it had its own copy of the
// parametrizations (with F1 = inf for SLAC and CTEQ, and until 2026 Q^2 instead of Q for MSTW).
void check_same(const std::string &name, const std::vector<double> &xs = xs_default,
                const std::vector<double> &Q2s = Q2s_default) {
  NucleonStructure sf(name);
  for (int proton = 0; proton <= 1; proton++)
    for (double x : xs)
      for (double Q2 : Q2s) {
        CAPTURE(name, proton, x, Q2);
        double F1 = 0., F2 = 0.;
        sf.getF_xQ(F1, F2, proton, x, Q2);
        NuclStructure nsf(proton, Q2, x, 0, name);
        CHECK_THAT(nsf.getF2(), Catch::Matchers::WithinRel(F2, 1e-12));
        CHECK_THAT(nsf.getF1(), Catch::Matchers::WithinRel(F1, 1e-12));
      }
}
}  // namespace

TEST_CASE("NuclStructure and NucleonStructure agree: SLAC", "[structure]") { check_same("SLAC"); }
TEST_CASE("NuclStructure and NucleonStructure agree: Christy-Bosted", "[structure]") { check_same("CB", xs_CB, Q2s_CB); }
TEST_CASE("NuclStructure and NucleonStructure agree: Alekhin", "[structure]") { check_same("Alekhin"); }
TEST_CASE("NuclStructure and NucleonStructure agree: CTEQ", "[structure]") { check_same("CTEQ"); }
TEST_CASE("NuclStructure and NucleonStructure agree: MSTW", "[structure]") {
  require_mstw();
  check_same("MSTW");
}

TEST_CASE("NuclStructure off-shell kinematics: W^2 independent of x", "[structure]") {
  // DeuteronStructure passes an off-shell W^2; Christy-Bosted and SLAC (Bodek) use it, the others x.
  // Checked against the former stand-alone implementation of NuclStructure on 216 points in 2026-11.
  const double Q2 = 2.e06, x = 0.5, mass = nucqcd::mass_p;
  const double Wsq_on = Q2 * (1. / x - 1.) + mass * mass, Wsq_off = Wsq_on - 0.3e06;
  for (const std::string name : {"CB", "SLAC"}) {
    CAPTURE(name);
    NucleonStructure sf(name);
    double F1 = 0., F2 = 0., F1on = 0., F2on = 0.;
    sf.getF_xQW(F1, F2, true, x, Q2, Wsq_off);
    sf.getF_xQ(F1on, F2on, true, x, Q2);
    NuclStructure off(true, Q2, x, Wsq_off, name);
    CHECK(off.getF2() == F2);
    CHECK(off.getF1() == F1);
    CHECK(F2 != F2on);  // W^2 matters for these two
  }
  // above W^2 = 9 GeV^2 Christy-Bosted switches to SLAC (DeuteronStructure asks isF2Fit; until 2026-11
  // getName() returned "SLAC" there and DeuteronStructure compared names)
  CHECK(NuclStructure(true, 4.e06, 0.2, 0, "CB").isF2Fit());
  CHECK_FALSE(NuclStructure(true, 2.e06, 0.5, 0, "CB").isF2Fit());
  CHECK(NuclStructure(true, 4.e06, 0.2, 0, "CB").getName() == "CB");
}

TEST_CASE("NuclStructure MSTW copy and assignment", "[structure]") {
  require_mstw();
  NuclStructure a(true, 5.e06, 0.3, 0, "MSTW");
  const double F2 = a.getF2();
  NuclStructure b(a);
  CHECK(b.getF2() == F2);
  NuclStructure c(false, 2.e06, 0.1, 0, "SLAC");
  c = a;
  CHECK(c.getF2() == F2);
}

TEST_CASE("NucleonStructure SLAC", "[structure]") { check_param("SLAC"); }
TEST_CASE("NucleonStructure Christy-Bosted", "[structure]") { check_param("CB", xs_CB, Q2s_CB); }
TEST_CASE("NucleonStructure Alekhin", "[structure]") { check_param("Alekhin"); }
TEST_CASE("NucleonStructure CTEQ", "[structure]") { check_param("CTEQ"); }
TEST_CASE("NucleonStructure MSTW", "[structure]") {
  require_mstw();
  check_param("MSTW");
}

TEST_CASE("Polarized structure functions g1, g1+g2 (GRSV2000)", "[structure]") {
  for (int proton = 0; proton <= 1; proton++)
    for (double x : {0.1, 0.3, 0.6}) {
      std::ostringstream key;
      key << "GRSV" << (proton ? "_p" : "_n") << "_x" << x << "_Q2_4";
      CHECK_REF(refs, key.str() + "_g1", NucleonStructure::getG1_grsv2000(proton, x, 4.e06));
      CHECK_REF(refs, key.str() + "_g1plusg2", NucleonStructure::getG1plusg2_grsv2000(proton, x, 4.e06), 1e-6);
    }
}
