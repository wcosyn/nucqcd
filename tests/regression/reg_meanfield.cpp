// Regression: relativistic mean-field nuclei, wave functions and densities (Glauber/MeanFieldNucleus*).
#include <test_support.hpp>

#include <Glauber/MeanFieldNucleusThick.hpp>

#include <sstream>
#include <string>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_meanfield");
}

TEST_CASE("Mean-field wave functions, densities and masses", "[meanfield]") {
  for (const char *name : {"C", "O", "Fe", "Pb"}) {
    MeanFieldNucleusThick nucl(MeanFieldNucleus::TypeNames.at(name), nucqcd_test::share_dir());
    const std::string n(name);
    CHECK_REF(refs, n + "_massA", nucl.getMassA());
    CHECK_REF(refs, n + "_massA_min_p", nucl.getMassA_min_proton());
    CHECK_REF(refs, n + "_massA_min_n", nucl.getMassA_min_neutron());
    CHECK_REF(refs, n + "_range", nucl.getRange());
    for (int shell = 0; shell < nucl.getTotalLevels(); shell++)
      for (double r : {0.5, 2., 4.}) {
        std::ostringstream key;
        key << n << "_shell" << shell << "_r" << r;
        CHECK_REF(refs, key.str() + "_G", nucl.getWave_G(shell, r), 1e-8, 1e-14);
        CHECK_REF(refs, key.str() + "_F", nucl.getWave_F(shell, r), 1e-8, 1e-14);
      }
    for (double r : {0.5, 2., 4., 6.}) {
      std::ostringstream key;
      key << n << "_r" << r;
      CHECK_REF(refs, key.str() + "_rho_p", nucl.getProtonDensity(r), 1e-8, 1e-16);
      CHECK_REF(refs, key.str() + "_rho_n", nucl.getNeutronDensity(r), 1e-8, 1e-16);
    }
  }
}
