// MePhys: kinematics classes (Kinematics2to2, Kinematics2to3WithLabAngles, ElectronKinematics) obey the exact relations of
// energy-momentum conservation.
#include <test_support.hpp>

#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Tools/Kinematics/Kinematics2to3WithLabAngles.hpp>
#include <Tools/Core/constants.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
using namespace nucqcd;

#include <catch2/matchers/catch_matchers_string.hpp>
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Kinematics2to2 for A(e,e'p): conservation laws", "[kinematics]") {
  // 12C(e,e'p) 11B: target A, residual A-1 ("meson"), ejected proton ("hyperon")
  const double mA = 11177.93, mA1 = 10255.10, Q2 = 1.0e06, omega = 600., pm = 150.;
  Kinematics2to2 kin(mA, mA1, nucqcd::mass_p, "qsquared:wlab:pklab", Q2, omega, pm);
  CHECK_THAT(kin.GetQsquared(), WithinRel(Q2, 1e-12));
  CHECK_THAT(kin.GetWlab(), WithinRel(omega, 1e-12));
  CHECK_THAT(kin.GetPklab(), WithinRel(pm, 1e-10));
  const double k = kin.GetKlab();
  CHECK_THAT(k * k - omega * omega, WithinRel(Q2, 1e-12));
  // invariant mass
  CHECK_THAT(kin.GetW() * kin.GetW(), WithinRel((mA + omega) * (mA + omega) - k * k, 1e-12));
  // energy conservation in the lab
  const double EY = std::sqrt(nucqcd::mass_p * nucqcd::mass_p + std::pow(kin.GetPYlab(), 2));
  const double Ek = std::sqrt(mA1 * mA1 + pm * pm);
  CHECK_THAT(EY + Ek, WithinRel(mA + omega, 1e-12));
  // momentum conservation along and perpendicular to q
  const double cY = kin.GetCosthYlab(), ck = kin.GetCosthklab();
  CHECK_THAT(kin.GetPYlab() * cY + pm * ck, WithinRel(k, 1e-9));
  CHECK_THAT(kin.GetPYlab() * std::sqrt(1. - cY * cY), WithinRel(pm * std::sqrt(1. - ck * ck), 1e-8));
  // s + t + u = sum of the masses squared (virtual photon: -Q^2)
  CHECK_THAT(kin.GetS() + kin.GetT() + kin.GetU(), WithinRel(mA * mA + mA1 * mA1 + nucqcd::mass_p * nucqcd::mass_p - Q2, 1e-10));
}

TEST_CASE("ElectronKinematics: scattering angle reproduces Q^2", "[kinematics]") {
  const double Ebeam = 3245.;
  std::unique_ptr<ElectronKinematics> elec(ElectronKinematics::CreateWithBeamEnergy(Ebeam));
  for (double Q2 : {0.5e06, 1.0e06, 2.0e06}) {
    const double nu = 700.;
    const double cth = elec->GetCosScatterAngle(Q2, nu);
    CHECK_THAT(elec->GetBeamEnergy(Q2, nu), WithinRel(Ebeam, 1e-12));
    CHECK_THAT(2. * Ebeam * (Ebeam - nu) * (1. - cth), WithinRel(Q2, 1e-6));  // massless electron
  }
}

// Until 2026-11 the check on the number of variables read (n<2) && (n>4), which is never true.
TEST_CASE("Kinematics rejects format strings without 2, 3 or 4 variables", "[kinematics]") {
  const double mA = 11177.93, mA1 = 10255.10;
  CHECK_THROWS_AS(Kinematics2to2(mA, mA1, nucqcd::mass_p, "qsquared", 1.e06, 600., 150.), std::invalid_argument);
  CHECK_THROWS_AS(Kinematics2to2(mA, mA1, nucqcd::mass_p, "qsquared:wlab:pklab:costhkcm:phi", 1.e06, 600., 150.),
                  std::invalid_argument);
  CHECK_NOTHROW(Kinematics2to2(mA, mA1, nucqcd::mass_p, "qsquared:wlab:pklab", 1.e06, 600., 150.));
}

TEST_CASE("Kinematics2to3WithLabAngles for A(e,e'pn): solutions and conservation laws", "[kinematics]") {
  // 12C(e,e'pn) as in progs/attic/twonucl: the proton ("N") at 1000 MeV and theta_N, the neutron
  // ("K") at theta_K, phi_K = pi; the residual A-2 nucleus is "Y" (isospin channel 7: user masses).
  // Until 2026-11 the format string of each solution was built with a sprintf whose destination was
  // also a source (undefined behaviour); with gcc -O3 constructing the solutions failed with "Format
  // string should provide 6 independent variables".
  const double mA = 11177.93, mA2 = mA - nucqcd::mass_p - nucqcd::mass_n + 30.;
  const double Q2 = 1631389.67942, omega = 580., pn = 1000., thn = 10. * nucqcd::deg_to_rad, thk = 40. * nucqcd::deg_to_rad;
  Kinematics2to3WithLabAngles kin(7, Kinematics2to3::kN, "qsquared:wlab:pn:thn:thk:phik", Q2, omega, pn,
                                   thn, thk, PI, mA, nucqcd::mass_p, nucqcd::mass_n, mA2);
  REQUIRE(kin.GetNumberOfSolutions() >= 1);
  CHECK(kin.IsPhysical());
  CHECK_THAT(kin.GetPn(), WithinRel(pn, 1e-9));
  CHECK_THAT(kin.GetN3().Theta(), WithinAbs(thn, 1e-9));
  CHECK_THAT(kin.GetK3().Theta(), WithinAbs(thk, 1e-9));
  const nucqcd::LorentzVector in = kin.GetG4() + kin.GetD4();
  const nucqcd::LorentzVector out = kin.GetN4() + kin.GetK4() + kin.GetY4();
  CHECK_THAT(out.E(), WithinRel(in.E(), 1e-10));
  CHECK_THAT((out - in).Vect().Mag(), WithinAbs(0., 1e-6));
  for (int i = 0; i < kin.GetNumberOfSolutions(); i++) {
    // the solutions were built from a valid format string, and with the masses (until 2026-11 they
    // got none in isospin channel 7, so their vectors were NaN)
    const Kinematics2to3 &sol = kin.GetSolution(i);
    CHECK_THAT(sol.GetQsquared(), WithinRel(Q2, 1e-9));
    CHECK_THAT(sol.GetWlab(), WithinRel(omega, 1e-9));
    CHECK(sol.IsPhysical());
    CHECK_THAT(sol.GetD4().M(), WithinRel(mA, 1e-9));
    CHECK_THAT(sol.GetY4().M(), WithinRel(mA2, 1e-9));
    CHECK_THAT(sol.GetPn(), WithinRel(pn, 1e-9));
    CHECK_THAT(sol.GetN3().Theta(), WithinAbs(thn, 1e-9));
    CHECK_THAT(sol.GetK3().Theta(), WithinAbs(thk, 1e-9));
    const nucqcd::LorentzVector sol_in = sol.GetG4() + sol.GetD4();
    const nucqcd::LorentzVector sol_out = sol.GetN4() + sol.GetK4() + sol.GetY4();
    CHECK_THAT(sol_out.E(), WithinRel(sol_in.E(), 1e-10));
    CHECK_THAT((sol_out - sol_in).Vect().Mag(), WithinAbs(0., 1e-6));
  }
}

TEST_CASE("Kinematics2to3: azimuthal angles with the standard axis convention", "[kinematics]") {
  // GetPhik(kS) lacked a return until 2026-11 and gave the kY result whatever the convention (kY
  // itself needs a hyperon-first format, which these A(e,e'pN) kinematics do not use)
  const double mA = 11177.93, mA2 = mA - nucqcd::mass_p - nucqcd::mass_n + 30.;
  // the first particle of the format fixes the axes: with kK that is the kaon (phik is then zero by
  // definition and phin is the input)
  for (auto ac : {Kinematics2to3::kN, Kinematics2to3::kK}) {
    Kinematics2to3WithLabAngles kin(7, ac,
                                     ac == Kinematics2to3::kK ? "qsquared:wlab:pk:thk:thn:phin" : "qsquared:wlab:pn:thn:thk:phik",
                                     1631389.67942, 580., 1000., 10. * nucqcd::deg_to_rad, 40. * nucqcd::deg_to_rad, 2.5, mA, nucqcd::mass_p, nucqcd::mass_n, mA2);
    REQUIRE(kin.GetNumberOfSolutions() >= 1);
    CHECK(kin.GetPhik() == kin.GetPhik(ac));
    CHECK(kin.GetPhiy() == kin.GetPhiy(ac));
    CHECK(kin.GetPhin() == kin.GetPhin(ac));
    CHECK(kin.GetVar("phik") == kin.GetPhik(ac));
  }
}

TEST_CASE("Format strings and variable names are case insensitive, empty tokens are left out", "[kinematics]") {
  // Since 2026-11 the parsing uses std::string instead of ROOT's TString (roadmap step 9, phase 3);
  // the behaviour of TString::Tokenize (no empty tokens) and of the case-insensitive comparisons is kept.
  const double mA = 11177.93, mA1 = 10255.10, Q2 = 1.0e06, omega = 600., pm = 150.;
  Kinematics2to2 ref(mA, mA1, nucqcd::mass_p, "qsquared:wlab:pklab", Q2, omega, pm);
  Kinematics2to2 upper(mA, mA1, nucqcd::mass_p, "QSquared::WLAB:pkLab:", Q2, omega, pm);
  CHECK(upper.GetPklab() == ref.GetPklab());
  CHECK(upper.GetCosthklab() == ref.GetCosthklab());
  CHECK(upper.GetVarName(2) == "wlab");
  CHECK_THROWS_AS(Kinematics2to2(mA, mA1, nucqcd::mass_p, "qsquared:wlab:nonsense", Q2, omega, pm),
                  std::invalid_argument);

  const double mA2 = mA - nucqcd::mass_p - nucqcd::mass_n + 30.;
  Kinematics2to3WithLabAngles k3(7, Kinematics2to3::kN, "QSQUARED:wlab:Pn:thn:THK:phik", 1631389.67942, 580.,
                                  1000., 10. * nucqcd::deg_to_rad, 40. * nucqcd::deg_to_rad, PI, mA, nucqcd::mass_p,
                                  nucqcd::mass_n, mA2);
  REQUIRE(k3.GetNumberOfSolutions() >= 1);
  const Kinematics2to3 &sol = k3.GetSolution(0);
  CHECK(sol.GetVar("PN") == sol.GetPn());
  CHECK(sol.GetVar("thk") == sol.GetThk());
  CHECK(sol.GetVarName(3) == "pn");
  CHECK_THROWS_AS(sol.GetVar("nonsense"), std::invalid_argument);
}

TEST_CASE("Kinematics2to3: variables by name (GetVar, GetVarArray, GetPhysicalVarArray)", "[kinematics]") {
  // the name -> getter table (Kinematics2to3_VarNames.cpp) replaced four copies of an if-chain in 2026-11;
  // it was checked against the old code on all 152 names, including GetVarArray and GetPhysicalVarArray
  const double mA = 11177.93, mA2 = mA - nucqcd::mass_p - nucqcd::mass_n + 30.;
  const double Q2 = 1631389.67942, omega = 580., pn = 1000., thn = 10. * nucqcd::deg_to_rad, thk = 40. * nucqcd::deg_to_rad;
  Kinematics2to3WithLabAngles lab(7, Kinematics2to3::kN, "qsquared:wlab:pn:thn:thk:phik", Q2, omega, pn, thn,
                                   thk, PI, mA, nucqcd::mass_p, nucqcd::mass_n, mA2);
  REQUIRE(lab.GetNumberOfSolutions() >= 1);
  Kinematics2to3 k(lab.GetSolution(0));
  CHECK(k.GetVar("pk") == k.GetPk());
  CHECK(k.GetVar("PK") == k.GetPk());                 // case insensitive
  CHECK(k.GetVar("thk_deg") == k.GetThk_deg());
  CHECK(k.GetVar("ynphincm") == k.GetVar("nyphincm"));  // alias
  // a variable that is not in the format string: kyphikcm_deg, the example of the documentation, was
  // misspelled "kyhpikcm_deg" in the table (so it was unknown) until 2026-11; the old name is an alias
  CHECK(k.GetVar("kyphikcm_deg") == k.GetKYsystem_phikcm_deg());
  CHECK(k.GetVar("kyhpikcm_deg") == k.GetVar("kyphikcm_deg"));
  std::unique_ptr<double[]> all(k.GetVarArray("pn"));
  CHECK(all[0] == k.GetPn());
  std::unique_ptr<double[]> physical(k.GetPhysicalVarArray("pn"));
  CHECK(physical[0] == k.GetPn());
  CHECK_THROWS_WITH(k.GetVar("nonsense"), ContainsSubstring("GetVar(const std::string&)") && ContainsSubstring("\"nonsense\" unknown"));
  CHECK_THROWS_WITH(k.GetVarArray("nonsense"), ContainsSubstring("GetVarArray") && ContainsSubstring("unknown"));
  CHECK_THROWS_WITH(k.GetPhysicalVarArray("nonsense"), ContainsSubstring("GetPhysicalVarArray") && ContainsSubstring("unknown"));
}

TEST_CASE("Kinematics2to3: format strings of the three scenarios reproduce the same event", "[kinematics]") {
  // ReadFormatString and UpdateKinematics were split into helpers and steps in 2026-11 (checked against the
  // old code on 900k random format strings, bitwise but for 14 angles within 1e-8 of 0 or pi): an event is
  // built from lab angles, then rebuilt from other sets of variables: an energy and angles in the CM of the
  // (2+3) subsystem (scenario 1), and mandelstam variables (scenario 3)
  Kinematics2to3WithLabAngles lab(1, Kinematics2to3::kN, "qsquared:wlab:pn:thn:thk:phik", 2.e5, 1800., 300., 1.0,
                                   0.7, 1.0);
  REQUIRE(lab.GetNumberOfSolutions() >= 1);
  const Kinematics2to3 &ev = lab.GetSolution(0);
  REQUIRE(ev.IsPhysical());
  auto same_event = [&](const Kinematics2to3 &k) {
    CHECK(k.IsPhysical());
    for (const char *v : {"pk", "py", "pn", "thk", "thy", "thn", "wlab", "stot", "sky", "skn", "syn"})
      CHECK_THAT(k.GetVar(v), WithinRel(ev.GetVar(v), 1e-6));
  };
  // scenario 1: wlab, pn and thn fix the first particle, the nucleon (the axis convention kN), the kaon's
  // angles are in the KY CM
  same_event(Kinematics2to3(1, Kinematics2to3::kN, "qsquared:wlab:pn:thn:kythkcm:kyphikcm", ev.GetVar("qsquared"),
                             ev.GetVar("wlab"), ev.GetVar("pn"), ev.GetVar("thn"), ev.GetVar("kythkcm"), ev.GetVar("kyphikcm")));
  // scenario 3: s_tot, t_gN, t_gK, t_dY, s_KY
  same_event(Kinematics2to3(1, Kinematics2to3::kN, "qsquared:stot:tgn:tgk:tdy:sky", ev.GetVar("qsquared"),
                             ev.GetVar("stot"), ev.GetVar("tgn"), ev.GetVar("tgk"), ev.GetVar("tdy"), ev.GetVar("sky")));
}

TEST_CASE("Kinematics2to3: format strings that make no sense are refused", "[kinematics]") {
  auto make = [](Kinematics2to3::EAxisConvention ac, const char *format) {
    return Kinematics2to3(1, ac, format, 2.e5, 1800., 300., 0.5, 0.5, 0.5);
  };
  // the same energy twice, a particle's energy and its momentum, an unknown name, no Q^2
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "qsquared:wlab:klab:pn:thn:phik"), ContainsSubstring("Format string is not sane!"));
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "qsquared:wlab:pn:en:thn:thk"), ContainsSubstring("Format string is not sane!"));
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "qsquared:wlab:pn:thn:thk:nonsense"), ContainsSubstring("\"nonsense\" unknown"));
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "wlab:pn:thn:thk:phik:ek"), ContainsSubstring("Q^2 needs to be one of 6"));
  // the azimuth of the particle on the x axis is zero by definition (named in the message)
  CHECK_THROWS_WITH(make(Kinematics2to3::kK, "phik:qsquared:wlab:thn:thk:ek"), ContainsSubstring("phik is zero by definition"));
  CHECK_THROWS_WITH(make(Kinematics2to3::kY, "phiy:qsquared:wlab:thn:thk:ey"), ContainsSubstring("phiy is zero by definition"));
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "phin:qsquared:wlab:thy:thk:en"), ContainsSubstring("phin is zero by definition"));
  // the variables of the second particle given in the lab frame together with a CM angle
  CHECK_THROWS_WITH(make(Kinematics2to3::kN, "qsquared:wlab:pn:thn:thk:kyphikcm"), ContainsSubstring("Format string is not sane!"));
}

TEST_CASE("Kinematics2to3: a format string without all variables of its scenario is refused (was a segmentation fault)", "[kinematics]") {
  // until 2026-11 these were accepted and UpdateKinematics dereferenced a null pointer: Q^2 twice, two
  // energies and no angle, ... (found by trying 900k format strings; each crashed)
  auto make = [](Kinematics2to3::EAxisConvention ac, const char *format) {
    return Kinematics2to3(1, ac, format, 2.e5, 1800., 300., 0.5, 0.5, 0.5);
  };
  for (const auto &c : {std::make_pair(Kinematics2to3::kN, "qsquared:ey:ykcosthycm:wlab:ykphiycm:en"),
                        std::make_pair(Kinematics2to3::kY, "qsquared:qsquared:ey:wkn:nkphikcm:tgk"),
                        std::make_pair(Kinematics2to3::kN, "costhn:kyphiycm:qsquared:qsquared:tgy:wtot")})
    CHECK_THROWS_WITH(make(c.first, c.second), ContainsSubstring("Format string is not sane!") && ContainsSubstring("Variables are missing"));
}
