// Regression: nucleon structure functions (Tools/NucleonStructure, mostly Fortran parametrizations).
// Each parametrization is its own TEST_CASE, i.e. its own process (Fortran COMMON/SAVE state).
#include <test_support.hpp>

#include <NucleonStructure.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace {
physcode_test::RefStore refs("reg_nucleon_structure");

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
}  // namespace

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
