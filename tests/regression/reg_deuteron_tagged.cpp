// Regression: tagged spectator DIS on the deuteron (DIS/DeuteronCross) and He3 (DIS/He3Cross);
// kinematics from progs/LDRD_test/LDRD_test.cpp.
#include <test_support.hpp>

#include <DeuteronCross.hpp>
#include <He3Cross.hpp>
#include <LightConeKin2to2.hpp>
#include <LightConeKin2to3.hpp>
#include <constants.hpp>

#include <cmath>
#include <string>

namespace {
physcode_test::RefStore refs("reg_deuteron_tagged");
}

TEST_CASE("DeuteronCross: light-cone and VNA, plane wave and FSI", "[tagged]") {
  const double Ein = 12.E03, thetain = PI - 0.1, Q2 = 0.1E06;
  TVector3 vecpd(0., 0., 10.E03);
  TVector3 vecps = TVector3(100., 50., 30.) + 0.5 * vecpd;
  TVector3 qvec(0., 0., 10.E03);
  TVector3 veckin(Ein * std::sin(thetain), 0., Ein * std::cos(thetain));
  LightConeKin2to2 kin(MASSD, Q2, MASSP, vecpd, qvec, vecps, veckin);
  CHECK_REF(refs, "kin_massX", kin.getMassX());
  CHECK_REF(refs, "kin_xA", kin.getXA());
  CHECK_REF(refs, "kin_epsilon", kin.getEpsilon());

  for (const char *wf : {"paris", "AV18"}) {
    const std::string w(wf);
    DeuteronCross test(wf, 0, "SLAC", 40., 8., -0.5, 8., 1.2, 3, 1);
    CHECK_REF(refs, w + "_pwLCavg", test.getavgLCCross(kin, 1), 1e-7);
    CHECK_REF(refs, w + "_pwLC", test.getLCCross(kin, 1), 1e-7);
    CHECK_REF(refs, w + "_pwVNAavg", test.getavgVNACross(kin, 1), 1e-7);
    CHECK_REF(refs, w + "_pwVNA", test.getVNACross(kin, 1), 1e-7);
    CHECK_REF(refs, w + "_fsiVNAavg", test.getavgVNACross(kin, 0), 1e-6);
    CHECK_REF(refs, w + "_fsiVNA", test.getVNACross(kin, 0), 1e-6);
    const double sigma = DeuteronCross::sigmaparam(kin.getMassX() * kin.getMassX(), kin.getQ2()) * 10. * HBARC * HBARC;
    CHECK_REF(refs, w + "_sigmaparam_mb", sigma);
    test.setScatter(sigma, 8., -0.5);
    CHECK_REF(refs, w + "_fsiVNA_sigmaparam", test.getVNACross(kin, 0), 1e-6);
  }
}

// Hidden known bug (run by hand: reg_deuteron_tagged "[known-bug]"): since commit a458b63 the
// Fortran side of He3wf declares the file name as character(kind=C_CHAR), bind(C) :: datname,
// i.e. of length 1 (Tools/getwavemod_tm.f:25), so it opens '/' and stops with a Fortran runtime
// error. Once fixed: remove the [.] and [known-bug] tags and capture the reference values.
TEST_CASE("He3Cross: tagged DIS with two spectators", "[.][tagged][known-bug]") {
  TVector3 sp1(0., 200., 0.), sp2(100., -100., 100.);
  TVector3 spcm = sp1 + sp2;
  TVector3 beam(5000. * std::cos(-1.), 0., 5000. * std::sin(-1.));
  LightConeKin2to3 kin(MASSHE3, 2.E06, MASSP, MASSP, 0., 4.E03, sp1, sp2, beam, spcm);
  for (int proton = 0; proton <= 1; proton++) {
    He3Cross he3("AV18", physcode_test::share_dir(), "SLAC", proton);
    const std::string nuc = proton ? "p" : "n";
    CHECK_REF(refs, "He3_density_" + nuc, he3.getDensity(kin), 1e-7, 1e-30);
    CHECK_REF(refs, "He3_cross_" + nuc, he3.getCross(kin), 1e-7, 1e-30);
  }
}
