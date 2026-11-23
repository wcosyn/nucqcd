// MePhys: LHAPDF structure functions (Tools/LHAPDFStructure, NucleonStructure "LHAPDF:<set>") and the
// combinations of parametrizations in NucleonStructure ("<low>+<high>@<W2>"), roadmap step 10.
// The LHAPDF tests follow tests/test_lhapdf.py of the poldeut repository; they need nucqcd built with
// LHAPDF and the sets CT18LO and NNPDFpol11_100 (lhapdf install CT18LO NNPDFpol11_100), and are skipped
// otherwise.
#include <test_support.hpp>

#include <LHAPDFStructure.hpp>
#include <NuclStructure.hpp>
#include <NucleonStructure.hpp>
#include <constants.hpp>

#if __has_include(<LHAPDF/LHAPDF.h>)
#include <LHAPDF/LHAPDF.h>
#define NUCQCD_TEST_LHAPDF_HEADER 1
#endif

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const char *const UNPOL = "CT18LO";
const char *const POL = "NNPDFpol11_100";
const double Q2_GeV2 = 10.;
const double Q2 = Q2_GeV2 * 1.E06;  // MeV^2
const std::vector<double> X = {0.01, 0.05, 0.1, 0.3, 0.5, 0.7};

void require_sets() {
  if (!nucqcd::LHAPDFStructure::available()) SKIP("nucqcd built without LHAPDF");
  if (!nucqcd::LHAPDFStructure::installedSet(UNPOL) || !nucqcd::LHAPDFStructure::installedSet(POL))
    SKIP("LHAPDF sets CT18LO and NNPDFpol11_100 not installed");
}

#ifdef NUCQCD_TEST_LHAPDF_HEADER
/* sum e_q^2 x (q + qbar) straight from LHAPDF, as in poldeut's test */
double manual(const LHAPDF::PDF &pdf, double x, double q2, bool proton, int nf) {
  const double e2[6] = {0., 1. / 9., 4. / 9., 1. / 9., 4. / 9., 1. / 9.};
  double out = 0.;
  for (int f = 1; f <= nf; f++) {
    const int fl = (proton || f > 2) ? f : 3 - f;
    out += e2[f] * (pdf.xfxQ2(fl, x, q2) + pdf.xfxQ2(-fl, x, q2));
  }
  return out;
}
#endif
}  // namespace

TEST_CASE("LHAPDF: F2 and g1 are the LO parton-model sums, the neutron by isospin", "[lhapdf]") {
  require_sets();
#ifndef NUCQCD_TEST_LHAPDF_HEADER
  SKIP("LHAPDF headers not available to the tests");
#else
  LHAPDF::setVerbosity(0);
  std::unique_ptr<LHAPDF::PDF> u(LHAPDF::mkPDF(UNPOL, 0)), p(LHAPDF::mkPDF(POL, 0));
  NucleonStructure lh(std::string("LHAPDF:") + UNPOL + ":" + POL);
  for (bool proton : {true, false}) {
    for (double x : X) {
      double F1, F2;
      lh.getF_xQ(F1, F2, proton, x, Q2);
      CHECK_THAT(F2, WithinRel(manual(*u, x, Q2_GeV2, proton, 5), 1e-12));
      CHECK_THAT(lh.getG1(proton, x, Q2), WithinRel(0.5 * manual(*p, x, Q2_GeV2, proton, 3) / x, 1e-12));
      // F1 from F2 with R1998 and the target-mass factor, as for MSTW
      const double mass = proton ? nucqcd::mass_p : nucqcd::mass_n;
      const double R = NucleonStructure::getr1998(x, Q2);
      CHECK_THAT(F1, WithinRel(F2 / (2. * x * (1. + R)) * (1. + 4. * mass * mass * x * x / Q2), 1e-12));
    }
  }
  // the proton and neutron differ (u <-> d)
  CHECK(std::abs(lh.getF2_xQ(true, 0.3, Q2) - lh.getF2_xQ(false, 0.3, Q2)) > 0.01);
#endif
}

TEST_CASE("LHAPDF: options for F1 (constant R, no target-mass factor) and the number of flavours", "[lhapdf]") {
  require_sets();
  nucqcd::LHAPDFOptions options;
  options.unpolarized = UNPOL;
  options.R = 0.18;
  options.target_mass = false;
  NucleonStructure lh(options);
  CHECK(lh.getName() == std::string("LHAPDF:") + UNPOL);
  for (double x : X) {
    double F1, F2;
    lh.getF_xQ(F1, F2, true, x, Q2);
    CHECK_THAT(F1, WithinRel(F2 / (2. * x * 1.18), 1e-12));  // as poldeut's default
  }
  options.nf_unpol = 3;
  NucleonStructure lh3(options);
  CHECK(lh3.getF2_xQ(true, 0.01, Q2) < lh.getF2_xQ(true, 0.01, Q2));  // no charm and bottom
}

TEST_CASE("LHAPDF: zero outside 0 < x < 1, one set only, malformed names", "[lhapdf]") {
  require_sets();
  NucleonStructure lh(std::string("LHAPDF:") + UNPOL + ":" + POL);
  for (double x : {0., 1., 1.2}) {
    double F1, F2;
    lh.getF_xQ(F1, F2, true, x, Q2);
    CHECK(F1 == 0.);
    CHECK(F2 == 0.);
    CHECK(lh.getG1(true, x, Q2) == 0.);
    CHECK(lh.getG1plusg2(true, x, Q2) == 0.);
  }
  NucleonStructure unpol(std::string("LHAPDF:") + UNPOL);
  CHECK(std::isfinite(unpol.getF2_xQ(true, 0.1, Q2)));
  CHECK_THROWS_AS(unpol.getG1(true, 0.1, Q2), std::logic_error);
  NucleonStructure polonly(std::string("LHAPDF::") + POL);
  CHECK(std::isfinite(polonly.getG1(true, 0.1, Q2)));
  CHECK_THROWS_AS(polonly.getF2_xQ(true, 0.1, Q2), std::logic_error);
  CHECK_THROWS_AS(NucleonStructure("LHAPDF:"), std::invalid_argument);
  CHECK_THROWS_AS(NucleonStructure("LHAPDF:a:b:c"), std::invalid_argument);
  CHECK_THROWS_AS(NucleonStructure("LHAPDF:NoSuchSet123"), std::invalid_argument);
  // copies share the loaded sets and give the same values
  NucleonStructure copy(lh);
  CHECK(copy.getF2_xQ(false, 0.2, Q2) == lh.getF2_xQ(false, 0.2, Q2));
  CHECK(copy.getG1plusg2(false, 0.2, Q2) == lh.getG1plusg2(false, 0.2, Q2));
}

TEST_CASE("LHAPDF: the Wandzura-Wilczek table against a direct integration", "[lhapdf]") {
  require_sets();
  nucqcd::LHAPDFOptions options;
  options.polarized = POL;
  const auto lh = nucqcd::LHAPDFStructure::get(options);
  for (bool proton : {true, false})
    for (double x : {0.02, 0.1, 0.3, 0.6}) {
      const double ref = lh->getG1plusg2_direct(proton, x, Q2);
      // poldeut's tolerance: rel = 2e-3, abs = 2e-5
      CHECK_THAT(lh->getG1plusg2(proton, x, Q2), WithinAbs(ref, 2e-5 + 2e-3 * std::abs(ref)));
    }
}

TEST_CASE("LHAPDF: close to SLAC F2 and GRSV2000 g1 of the proton", "[lhapdf]") {
  // rough check, as in poldeut: different fits, orders and higher-twist content (up to ~14% at Q^2 = 10)
  require_sets();
  NucleonStructure lh(std::string("LHAPDF:") + UNPOL + ":" + POL), slac("SLAC");
  for (double x : {0.1, 0.2, 0.3, 0.5}) {
    CHECK_THAT(lh.getF2_xQ(true, x, Q2), WithinRel(slac.getF2_xQ(true, x, Q2), 0.2));
    CHECK_THAT(lh.getG1(true, x, Q2), WithinRel(NucleonStructure::getG1_grsv2000(true, x, Q2), 0.2));
  }
}

TEST_CASE("LHAPDF: \"MSTW\" is the set MSTW2008lo68cl", "[lhapdf]") {
  if (!nucqcd::LHAPDFStructure::available() || !nucqcd::LHAPDFStructure::installedSet("MSTW2008lo68cl"))
    SKIP("needs LHAPDF and the set MSTW2008lo68cl");
  NucleonStructure mstw("MSTW"), lh("LHAPDF:MSTW2008lo68cl");
  CHECK(mstw.getName() == "MSTW");
  for (bool proton : {true, false})
    for (double x : X) {
      double F1, F2, F1l, F2l;
      mstw.getF_xQ(F1, F2, proton, x, Q2);
      lh.getF_xQ(F1l, F2l, proton, x, Q2);
      CHECK(F1 == F1l);
      CHECK(F2 == F2l);
    }
  // values of the former bundled MSTW code (share/mstw_grids/mstw2008lo.00.dat) at a grid node
  // (Q^2 = 100 GeV^2), where both interpolations give the grid value: p, x = 0.3 and 0.6
  CHECK_THAT(mstw.getF2_xQ(true, 0.3, 100.E06), WithinRel(0.2326534444444445, 1e-12));
  CHECK_THAT(mstw.getF2_xQ(true, 0.6, 100.E06), WithinRel(0.044726789422222225, 1e-12));
}

TEST_CASE("Combinations: CB+SLAC is CB up to the W^2 boundary and SLAC above", "[structure]") {
  NucleonStructure cb("CB"), slac("SLAC"), comb("CB+SLAC"), comb4("CB+SLAC@4"), explicit4("CB", "SLAC", 4.E06);
  CHECK(comb.getName() == "CB+SLAC");
  CHECK(explicit4.getName() == "CB+SLAC@4");
  const double x = 0.5;
  for (double Wsq : {2.E06, 3.9E06, 4.1E06, 8.9E06, 9.1E06, 16.E06}) {
    double F1, F2, F1c, F2c;
    (Wsq > 9.E06 ? slac : cb).getF_xQW(F1, F2, true, x, Q2, Wsq);
    comb.getF_xQW(F1c, F2c, true, x, Q2, Wsq);
    CHECK(F1c == F1);
    CHECK(F2c == F2);
    CHECK(comb.isF2Fit(Wsq) == (Wsq > 9.E06));
    CHECK(comb.getName(Wsq) == (Wsq > 9.E06 ? "SLAC" : "CB"));
    (Wsq > 4.E06 ? slac : cb).getF_xQW(F1, F2, false, x, Q2, Wsq);
    comb4.getF_xQW(F1c, F2c, false, x, Q2, Wsq);
    CHECK(F2c == F2);
    explicit4.getF_xQW(F1c, F2c, false, x, Q2, Wsq);
    CHECK(F2c == F2);
  }
  // copies and assignment keep the combination
  NucleonStructure copy(comb4), assigned;
  assigned = comb4;
  CHECK(copy.getName() == "CB+SLAC@4");
  CHECK(assigned.getF2_WQ(true, 5.E06, Q2) == comb4.getF2_WQ(true, 5.E06, Q2));
  CHECK_FALSE(NucleonStructure("CB").isF2Fit(16.E06));
  CHECK(NucleonStructure("SLAC").isF2Fit(1.E06));
  CHECK_THROWS_AS(NucleonStructure("+SLAC"), std::invalid_argument);
  CHECK_THROWS_AS(NucleonStructure("CB+"), std::invalid_argument);
  CHECK_THROWS_AS(NucleonStructure("CB+SLAC@"), std::invalid_argument);
  CHECK_THROWS_AS(NucleonStructure("CB+SLAC@4GeV"), std::invalid_argument);
}

TEST_CASE("NuclStructure \"CB\" is CB up to W^2 = 9 GeV^2 and SLAC above, and keeps its name", "[structure]") {
  NucleonStructure cb("CB"), slac("SLAC");
  for (double Wsq : {4.E06, 8.9E06, 9.1E06, 16.E06}) {
    NuclStructure nucl(true, Q2, 0.4, Wsq, "CB");
    double F1, F2, F1n, F2n;
    (Wsq > 9.E06 ? slac : cb).getF_xQW(F1, F2, true, 0.4, Q2, Wsq);
    nucl.getF(F1n, F2n);
    CHECK(F1n == F1);
    CHECK(F2n == F2);
    CHECK(nucl.getName() == "CB");  // until 2026-11 "SLAC" above 9 GeV^2
    CHECK(nucl.isF2Fit() == (Wsq > 9.E06));
  }
}
