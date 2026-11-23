// Glauber: relativistic mean-field nuclei (MeanFieldNucleus[Thick]) read from the share/ input.
// Shell wave functions and densities must be normalized, shell occupations must match Z and N.
#include <test_support.hpp>

#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Tools/Core/constants.hpp>

#include <gsl/gsl_sf_bessel.h>

#include <cmath>
#include <functional>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
double simpson(const std::function<double(double)> &f, double a, double b, int n) {
  const double h = (b - a) / n;
  double s = f(a) + f(b);
  for (int i = 1; i < n; i++) s += f(a + i * h) * (i % 2 ? 4. : 2.);
  return s * h / 3.;
}
}  // namespace

TEST_CASE("Mean-field nuclei: quantum numbers, normalization of shells and densities", "[meanfield]") {
  const std::string dir = nucqcd_test::share_dir();
  for (const char *name : {"He", "C", "O", "Ca40", "Fe", "Pb"}) {
    INFO("nucleus " << name);
    MeanFieldNucleusThick nucl(MeanFieldNucleus::TypeNames.at(name), dir);
    CHECK(nucl.getNucleusName() == name);
    CHECK(nucl.getA() == nucl.getZ() + nucl.getN());
    CHECK(nucl.getTotalLevels() == nucl.getPLevels() + nucl.getNLevels());
    // occupation: (2j+1) per shell (J_array holds 2j); the last shell may be partially filled
    int protons = 0, neutrons = 0;
    for (int i = 0; i < nucl.getPLevels(); i++) protons += nucl.getJ_array()[i] + 1;
    for (int i = nucl.getPLevels(); i < nucl.getTotalLevels(); i++) neutrons += nucl.getJ_array()[i] + 1;
    CHECK(protons >= nucl.getZ());
    CHECK(protons - (nucl.getJ_array()[nucl.getPLevels() - 1] + 1) < nucl.getZ());
    CHECK(neutrons >= nucl.getN());
    CHECK(neutrons - (nucl.getJ_array()[nucl.getTotalLevels() - 1] + 1) < nucl.getN());
    // radial Dirac wave functions: int (G^2+F^2) dr = 1
    const double R = nucl.getRange();
    for (int shell = 0; shell < nucl.getTotalLevels(); shell++) {
      INFO("shell " << shell);
      const double norm = simpson([&](double r) {
        return std::pow(nucl.getWave_G(shell, r), 2) + std::pow(nucl.getWave_F(shell, r), 2);
      }, 0., R, 2000);
      CHECK_THAT(norm, WithinAbs(1., 5e-3));
    }
    // densities include r^2 and are normalized as 4 pi int dr r^2 rho = 1
    CHECK_THAT(4. * PI * simpson([&](double r) { return nucl.getProtonDensity(r); }, 0., R, 2000), WithinAbs(1., 5e-3));
    CHECK_THAT(4. * PI * simpson([&](double r) { return nucl.getNeutronDensity(r); }, 0., R, 2000), WithinAbs(1., 5e-3));
    // binding: separation energies positive
    CHECK(nucl.getMassA() < nucl.getMassA_min_proton() + nucqcd::mass_p);
    CHECK(nucl.getMassA() < nucl.getMassA_min_neutron() + nucqcd::mass_n);
  }
}

TEST_CASE("Mean-field momentum distributions are normalized (Parseval)", "[meanfield]") {
  // (2/pi) int dk k^2 [ (int dr r G j_l(kr))^2 + (int dr r F j_lbar(kr))^2 ] = int (G^2+F^2) dr = 1
  // (note: progs/momdistr/momdistr.cpp uses l instead of lbar for the lower component F)
  MeanFieldNucleusThick nucl(MeanFieldNucleus::TypeNames.at("C"), nucqcd_test::share_dir());
  const double R = nucl.getRange();
  for (int shell = 0; shell < nucl.getTotalLevels(); shell++) {
    INFO("shell " << shell);
    const int l = nucl.getL_array()[shell], lbar = nucl.getLbar_array()[shell];
    auto density = [&](double k) {  // k in fm^-1
      const double up = simpson([&](double r) { return r * nucl.getWave_G(shell, r) * gsl_sf_bessel_jl(l, k * r); }, 0., R, 400);
      const double low = simpson([&](double r) { return r * nucl.getWave_F(shell, r) * gsl_sf_bessel_jl(lbar, k * r); }, 0., R, 400);
      return k * k * (up * up + low * low);
    };
    CHECK_THAT(2. / PI * simpson(density, 0., 8., 400), WithinAbs(1., 5e-3));
  }
}
