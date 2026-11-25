// Regression: deuteron wave function values (Deuteron, OldDeuteron) in r and p space, all
// parametrizations; also pins the sign convention of the D wave.
#include <test_support.hpp>

#include <Tools/WaveFunctions/OldDeuteron.hpp>
#include <Tools/WaveFunctions/Deuteron.hpp>

#include <memory>
#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_deuteron_wf");
}

TEST_CASE("Deuteron wave functions at selected r and p", "[deuteron]") {
  for (const char *name : {"AV18", "AV18b", "CDBonn", "Paris", "CDBonn87", "nijm3", "grossIIb"}) {
    INFO("wave function " << name);
    std::unique_ptr<Deuteron::Wavefunction> wf(Deuteron::Wavefunction::CreateWavefunction(name));
    REQUIRE(wf);
    for (double k : {20., 100., 300., 700.}) {
      std::ostringstream key;
      key << name << "_p" << k;
      CHECK_REF(refs, key.str() + "_U", wf->GetUp(k));
      CHECK_REF(refs, key.str() + "_W", wf->GetWp(k));
    }
    if (wf->IsRelativistic()) continue;
    for (double r : {0.5, 1., 2., 5.}) {
      std::ostringstream key;
      key << name << "_r" << r;
      CHECK_REF(refs, key.str() + "_u", wf->GetUr(r));
      CHECK_REF(refs, key.str() + "_w", wf->GetWr(r));
    }
  }
}

TEST_CASE("Sign convention of the D wave (U*W at fixed p and r)", "[deuteron]") {
  // pins the relative sign of S and D waves; the LF densities in Poldeut assume a convention
  for (const char *name : {"AV18", "CDBonn", "Paris"}) {
    std::unique_ptr<Deuteron::Wavefunction> wf(Deuteron::Wavefunction::CreateWavefunction(name));
    CHECK_REF(refs, std::string(name) + "_sign_UW_p100", wf->GetUp(100.) * wf->GetWp(100.) > 0 ? 1. : -1.);
    CHECK_REF(refs, std::string(name) + "_sign_uw_r2", wf->GetUr(2.) * wf->GetWr(2.) > 0 ? 1. : -1.);
  }
}

TEST_CASE("OldDeuteron wave functions", "[deuteron]") {
  for (const char *name : {"AV18", "AV18b", "CDBonn", "Paris"}) {
    OldDeuteron wf(name);
    for (double k : {20., 100., 300.}) {
      std::ostringstream key;
      key << "Old_" << name << "_p" << k;
      CHECK_REF(refs, key.str() + "_U", wf.U(k));
      CHECK_REF(refs, key.str() + "_W", wf.W(k));
    }
    const auto psi = wf.deuteronwf(0, 1, 1, -1, 150., 0.7, 0.4);
    CHECK_REF(refs, std::string("Old_") + name + "_psi_re", psi.real(), 1e-8, 1e-14);
    CHECK_REF(refs, std::string("Old_") + name + "_psi_im", psi.imag(), 1e-8, 1e-14);
  }
}
