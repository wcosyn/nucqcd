// Regression: nucleon structure functions (Tools/NucleonStructure, mostly Fortran parametrizations).
// Each parametrization is its own TEST_CASE, i.e. its own process (Fortran COMMON/SAVE state).
// Tools/NuclStructure is a second wrapper around the same parametrizations; it is the one used by
// DeuteronStructure (DeuteronCross, InclusiveCross, He3Cross), so the two are checked against each other.
#include <test_support.hpp>

#include <NuclStructure.hpp>
#include <NucleonStructure.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace {
nucqcd_test::RefStore refs("reg_nucleon_structure");

// default grid of (x, Q^2[MeV^2]) points; Christy-Bosted only covers W < ~3 GeV (and calls exit()
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

// NuclStructure (used by the deuteron code) and NucleonStructure must give the same F2 (and F1,
// where NuclStructure defines it: for SLAC and CTEQ it returns inf, and DeuteronStructure builds
// F1 from F2 itself for SLAC). Until 2026 NuclStructure passed Q^2 instead of Q to MSTW.
void check_same(const std::string &name, bool with_F1, const std::vector<double> &xs = xs_default,
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
        if (with_F1) CHECK_THAT(nsf.getF1(), Catch::Matchers::WithinRel(F1, 1e-12));
      }
}
}  // namespace

TEST_CASE("NuclStructure and NucleonStructure agree: SLAC", "[structure]") { check_same("SLAC", false); }
TEST_CASE("NuclStructure and NucleonStructure agree: Christy-Bosted", "[structure]") {
  check_same("CB", true, xs_CB, Q2s_CB);
}
TEST_CASE("NuclStructure and NucleonStructure agree: Alekhin", "[structure]") { check_same("Alekhin", true); }
TEST_CASE("NuclStructure and NucleonStructure agree: CTEQ", "[structure]") { check_same("CTEQ", false); }
TEST_CASE("NuclStructure and NucleonStructure agree: MSTW", "[structure]") { check_same("MSTW", true); }

TEST_CASE("NuclStructure MSTW copy and assignment", "[structure]") {
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
TEST_CASE("NucleonStructure MSTW", "[structure]") { check_param("MSTW"); }

TEST_CASE("Polarized structure functions g1, g1+g2 (GRSV2000)", "[structure]") {
  for (int proton = 0; proton <= 1; proton++)
    for (double x : {0.1, 0.3, 0.6}) {
      std::ostringstream key;
      key << "GRSV" << (proton ? "_p" : "_n") << "_x" << x << "_Q2_4";
      CHECK_REF(refs, key.str() + "_g1", NucleonStructure::getG1_grsv2000(proton, x, 4.e06));
      CHECK_REF(refs, key.str() + "_g1plusg2", NucleonStructure::getG1plusg2_grsv2000(proton, x, 4.e06), 1e-6);
    }
}
