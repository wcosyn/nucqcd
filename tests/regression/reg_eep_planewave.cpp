// Regression: A(e,e'p) in the plane-wave approximation (Knockout/Cross); free electron-proton
// cross section and bound-nucleon cross sections per shell. No Glauber grids are computed.
#include <test_support.hpp>

#include <Cross.hpp>
#include <MeanFieldNucleusThick.hpp>
#include <TElectronKinematics.h>
#include <TKinematics2to2.h>
#include <constants.hpp>

#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

namespace {
physcode_test::RefStore refs("reg_eep_planewave");

// kinematics of progs/onenucl/A_QE_cross.cpp (Dutta data, Ebeam 3245 MeV)
void check_nucleus(const char *name) {
  const std::string dir = physcode_test::share_dir();
  MeanFieldNucleusThick nucl(MeanFieldNucleus::TypeNames.at(name), dir);
  std::unique_ptr<TElectronKinematics> elec(TElectronKinematics::CreateWithBeamEnergy(3245.));
  Cross obs(*elec, &nucl, 1.E-05, 2, dir, false, 1.);
  for (double pm : {50., 150., 250.}) {
    TKinematics2to2 kin("", "", nucl.getMassA(), nucl.getMassA_min_proton(), MASSP, "qsquared:wlab:pklab", 1.8E06,
                        1000., pm);
    std::ostringstream base;
    base << name << "_pm" << pm;
    for (int current : {1, 2}) {
      std::ostringstream key;
      key << base.str() << "_cc" << current;
      CHECK_REF(refs, key.str() + "_free", obs.getElCross(kin, current, 0., 2000) * HBARC * HBARC, 1e-8, 1e-30);
    }
    for (int shell = 0; shell < nucl.getPLevels(); shell++)
      for (double phi : {0., PI}) {
        std::ostringstream key;
        key << base.str() << "_shell" << shell << "_phi" << phi;
        // getDiffCross(kin, current, thick, SRC, CT, pw, shell, phi, maxEval, lab, phiavg, lc_mod, nkt_mod)
        CHECK_REF(refs, key.str() + "_pw", obs.getDiffCross(kin, 2, 0, 0, 0, 1, shell, phi, 20000, true, false, 1., 1.),
                  1e-7, 1e-30);
      }
  }
}
}  // namespace

TEST_CASE("A(e,e'p) plane wave: carbon", "[eep]") { check_nucleus("C"); }
TEST_CASE("A(e,e'p) plane wave: iron", "[eep]") { check_nucleus("Fe"); }

TEST_CASE("Plane-wave calculations do not write Glauber grids", "[eep]") {
  check_nucleus("O");
  const std::filesystem::path grids = std::filesystem::path(physcode_test::share_dir()) / "grids";
  for (const auto &entry : std::filesystem::directory_iterator(grids))
    CHECK(entry.path().filename().string().rfind("O.", 0) != 0);  // no grid files for oxygen
}
