// Regression: (anti)neutrino quasi-elastic scattering (Knockout/WeakQECross): elementary nucleon
// cross sections (charged current only; the elementary routine returns 0 for neutral current) and
// plane-wave nuclear cross sections for charged and neutral current.
// NOTE: the neutral-current nuclear prefactor (Knockout/WeakQECross.cpp:637-638) is suspect:
//   pow(E_beam - omega*G_FERMI/(Q2/M_Z/M_Z+1.)/PI, 2.)  -- probably meant (E_beam-omega)*G_FERMI/...
// (compare the charged-current branch, line 157), which makes the NC values ~1e23 too large, and
// M_Z in constants.hpp is 10x too large (see unit_constants). Fixing either moves the "NC" values.
#include <test_support.hpp>

#include <MeanFieldNucleusThick.hpp>
#include <TElectronKinematics.h>
#include <TKinematics2to2.h>
#include <TLeptonKinematics.h>
#include <WeakQECross.hpp>
#include <constants.hpp>

#include <memory>
#include <sstream>
#include <string>

namespace {
nucqcd_test::RefStore refs("reg_weakqe");
const double M_A = 1.03E03, mmu = 105.6583715;
}  // namespace

TEST_CASE("Elementary charged-current QE cross sections", "[weakqe]") {
  for (int proton = 0; proton <= 1; proton++)
    for (double Q2 : {0.2e06, 1.0e06})
      for (double E : {1000., 3000.}) {
        std::ostringstream key;
        key << "CC" << (proton ? "_p" : "_n") << "_Q2_" << Q2 * 1e-6 << "_E" << E;
        for (auto type : {WeakQECross::realQ2, WeakQECross::costheta_lab}) {
          const double xs = WeakQECross::getElWeakQECross(Q2, E, proton, true, M_A, mmu * mmu, type);
          CHECK(xs > 0.);
          CHECK_REF(refs, key.str() + (type == WeakQECross::realQ2 ? "_dQ2" : "_dOmega"), xs, 1e-8, 1e-40);
        }
      }
}

TEST_CASE("Nuclear plane-wave charged-current QE on carbon", "[weakqe]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  std::unique_ptr<TLeptonKinematics> muon(TLeptonKinematics::CreateWithBeamEnergy(TLeptonKinematics::muon, 1500.));
  WeakQECross cc(muon.get(), &carbon, 1.E-05, 2, dir, true, M_A, false, false);
  for (double pm : {50., 200.}) {
    // proton shells: antineutrino CC knocks out a neutron
    TKinematics2to2 kin("", "", carbon.getMassA(), carbon.getMassA_min_proton(), MASSN, "qsquared:wlab:pklab", 0.5E06,
                        350., pm);
    for (int shell = 0; shell < carbon.getPLevels(); shell++) {
      std::ostringstream key;
      key << "C_pm" << pm << "_shell" << shell << "_CC_pw";
      // getDiffWeakQECross(kin, current, thick, SRC, CT, pw, shell, phi, maxEval, d_type, phi_int)
      const double xs = cc.getDiffWeakQECross(kin, 2, 0, 0, 0, 1, shell, 0., 20000, WeakQECross::ctheta_lab, false);
      CHECK(xs > 0.);
      CHECK_REF(refs, key.str(), xs, 1e-7, 1e-40);
    }
  }
}

// Hidden known bug (run by hand: reg_weakqe "[known-bug]"): the neutral-current constructor
// WeakQECross(TElectronKinematics*, ...) never initializes the member `lepton`, but
// getDiffWeakQECross dereferences lepton->GetBeamEnergy(kin) for every process, so NC nuclear
// cross sections segfault (or return garbage). See also the NC prefactor issue at the top of this
// file. Once fixed: remove [.] and [known-bug], and capture reference values.
TEST_CASE("Nuclear plane-wave neutral-current QE on carbon", "[.][weakqe][known-bug]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  std::unique_ptr<TElectronKinematics> elec(TElectronKinematics::CreateWithBeamEnergy(1500.));
  WeakQECross nc(elec.get(), &carbon, 1.E-05, 2, dir, false, M_A, false, false);
  TKinematics2to2 kin("", "", carbon.getMassA(), carbon.getMassA_min_proton(), MASSP, "qsquared:wlab:pklab", 0.5E06,
                      350., 50.);
  for (int shell = 0; shell < carbon.getPLevels(); shell++) {
    std::ostringstream key;
    key << "C_pm50_shell" << shell << "_NC_pw";
    CHECK_REF(refs, key.str(), nc.getDiffWeakQECross(kin, 2, 0, 0, 0, 1, shell, 0., 20000, WeakQECross::ctheta_lab, false),
              1e-7, 1e-40);
  }
}
