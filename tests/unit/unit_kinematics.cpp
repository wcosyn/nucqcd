// MePhys: kinematics classes (TKinematics2to2, TElectronKinematics) obey the exact relations of
// energy-momentum conservation.
#include <test_support.hpp>

#include <TElectronKinematics.h>
#include <TKinematics2to2.h>
#include <constants.hpp>

#include <cmath>
#include <memory>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("TKinematics2to2 for A(e,e'p): conservation laws", "[kinematics]") {
  // 12C(e,e'p) 11B: target A, residual A-1 ("meson"), ejected proton ("hyperon")
  const double mA = 11177.93, mA1 = 10255.10, Q2 = 1.0e06, omega = 600., pm = 150.;
  TKinematics2to2 kin("", "", mA, mA1, MASSP, "qsquared:wlab:pklab", Q2, omega, pm);
  CHECK_THAT(kin.GetQsquared(), WithinRel(Q2, 1e-12));
  CHECK_THAT(kin.GetWlab(), WithinRel(omega, 1e-12));
  CHECK_THAT(kin.GetPklab(), WithinRel(pm, 1e-10));
  const double k = kin.GetKlab();
  CHECK_THAT(k * k - omega * omega, WithinRel(Q2, 1e-12));
  // invariant mass
  CHECK_THAT(kin.GetW() * kin.GetW(), WithinRel((mA + omega) * (mA + omega) - k * k, 1e-12));
  // energy conservation in the lab
  const double EY = std::sqrt(MASSP * MASSP + std::pow(kin.GetPYlab(), 2));
  const double Ek = std::sqrt(mA1 * mA1 + pm * pm);
  CHECK_THAT(EY + Ek, WithinRel(mA + omega, 1e-12));
  // momentum conservation along and perpendicular to q
  const double cY = kin.GetCosthYlab(), ck = kin.GetCosthklab();
  CHECK_THAT(kin.GetPYlab() * cY + pm * ck, WithinRel(k, 1e-9));
  CHECK_THAT(kin.GetPYlab() * std::sqrt(1. - cY * cY), WithinRel(pm * std::sqrt(1. - ck * ck), 1e-8));
  // s + t + u = sum of the masses squared (virtual photon: -Q^2)
  CHECK_THAT(kin.GetS() + kin.GetT() + kin.GetU(), WithinRel(mA * mA + mA1 * mA1 + MASSP * MASSP - Q2, 1e-10));
}

TEST_CASE("TElectronKinematics: scattering angle reproduces Q^2", "[kinematics]") {
  const double Ebeam = 3245.;
  std::unique_ptr<TElectronKinematics> elec(TElectronKinematics::CreateWithBeamEnergy(Ebeam));
  for (double Q2 : {0.5e06, 1.0e06, 2.0e06}) {
    const double nu = 700.;
    const double cth = elec->GetCosScatterAngle(Q2, nu);
    CHECK_THAT(elec->GetBeamEnergy(Q2, nu), WithinRel(Ebeam, 1e-12));
    CHECK_THAT(2. * Ebeam * (Ebeam - nu) * (1. - cth), WithinRel(Q2, 1e-6));  // massless electron
  }
}
