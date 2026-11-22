// MePhys: kinematics classes (TKinematics2to2, TKinematics2to3WithLabAngles, TElectronKinematics) obey the exact relations of
// energy-momentum conservation.
#include <test_support.hpp>

#include <TElectronKinematics.h>
#include <TKinematics2to2.h>
#include <TKinematics2to3WithLabAngles.h>
#include <constants.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>

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

// Until 2026-11 the check on the number of variables read (n<2) && (n>4), which is never true.
TEST_CASE("TKinematics rejects format strings without 2, 3 or 4 variables", "[kinematics]") {
  const double mA = 11177.93, mA1 = 10255.10;
  CHECK_THROWS_AS(TKinematics2to2("", "", mA, mA1, MASSP, "qsquared", 1.e06, 600., 150.), std::invalid_argument);
  CHECK_THROWS_AS(TKinematics2to2("", "", mA, mA1, MASSP, "qsquared:wlab:pklab:costhkcm:phi", 1.e06, 600., 150.),
                  std::invalid_argument);
  CHECK_NOTHROW(TKinematics2to2("", "", mA, mA1, MASSP, "qsquared:wlab:pklab", 1.e06, 600., 150.));
}

TEST_CASE("TKinematics2to3WithLabAngles for A(e,e'pn): solutions and conservation laws", "[kinematics]") {
  // 12C(e,e'pn) as in progs/attic/twonucl: the proton ("N") at 1000 MeV and theta_N, the neutron
  // ("K") at theta_K, phi_K = pi; the residual A-2 nucleus is "Y" (isospin channel 7: user masses).
  // Until 2026-11 the format string of each solution was built with a sprintf whose destination was
  // also a source (undefined behaviour); with gcc -O3 constructing the solutions failed with "Format
  // string should provide 6 independent variables".
  const double mA = 11177.93, mA2 = mA - MASSP - MASSN + 30.;
  const double Q2 = 1631389.67942, omega = 580., pn = 1000., thn = 10. * DEGRTORAD, thk = 40. * DEGRTORAD;
  TKinematics2to3WithLabAngles kin("", "", 7, TKinematics2to3::kN, "qsquared:wlab:pn:thn:thk:phik", Q2, omega, pn,
                                   thn, thk, PI, mA, MASSP, MASSN, mA2);
  REQUIRE(kin.GetNumberOfSolutions() >= 1);
  CHECK(kin.IsPhysical());
  CHECK_THAT(kin.GetPn(), WithinRel(pn, 1e-9));
  CHECK_THAT(kin.GetN3().Theta(), WithinAbs(thn, 1e-9));
  CHECK_THAT(kin.GetK3().Theta(), WithinAbs(thk, 1e-9));
  const TLorentzVector in = kin.GetG4() + kin.GetD4();
  const TLorentzVector out = kin.GetN4() + kin.GetK4() + kin.GetY4();
  CHECK_THAT(out.E(), WithinRel(in.E(), 1e-10));
  CHECK_THAT((out - in).Vect().Mag(), WithinAbs(0., 1e-6));
  for (int i = 0; i < kin.GetNumberOfSolutions(); i++) {
    // the solutions were built from a valid format string, and with the masses (until 2026-11 they
    // got none in isospin channel 7, so their vectors were NaN)
    const TKinematics2to3 &sol = kin.GetSolution(i);
    CHECK_THAT(sol.GetQsquared(), WithinRel(Q2, 1e-9));
    CHECK_THAT(sol.GetWlab(), WithinRel(omega, 1e-9));
    CHECK(sol.IsPhysical());
    CHECK_THAT(sol.GetD4().M(), WithinRel(mA, 1e-9));
    CHECK_THAT(sol.GetY4().M(), WithinRel(mA2, 1e-9));
    CHECK_THAT(sol.GetPn(), WithinRel(pn, 1e-9));
    CHECK_THAT(sol.GetN3().Theta(), WithinAbs(thn, 1e-9));
    CHECK_THAT(sol.GetK3().Theta(), WithinAbs(thk, 1e-9));
    const TLorentzVector sol_in = sol.GetG4() + sol.GetD4();
    const TLorentzVector sol_out = sol.GetN4() + sol.GetK4() + sol.GetY4();
    CHECK_THAT(sol_out.E(), WithinRel(sol_in.E(), 1e-10));
    CHECK_THAT((sol_out - sol_in).Vect().Mag(), WithinAbs(0., 1e-6));
  }
}
