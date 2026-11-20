// Regression: polarized deuteron light-front spectator tagging (DIS/Poldeut).
#include <test_support.hpp>

#include <Poldeut.hpp>

#include <sstream>
#include <string>

namespace {
physcode_test::RefStore refs("reg_poldeut");
}

TEST_CASE("Poldeut light-front distributions", "[poldeut]") {
  for (const char *wf : {"AV18", "CDBonn", "Paris"}) {
    Poldeut pd(wf, "SLAC");
    const std::string w(wf);
    for (double alpha : {0.9, 1.0, 1.15})
      for (double pt : {20., 150.}) {
        double S[6];
        pd.getLFdistributions(alpha, pt, S[0], S[1], S[2], S[3], S[4], S[5]);
        const char *names[6] = {"S_L", "S_T", "deltaS_L", "deltaS_T", "deltaT_S_L", "deltaT_S_T"};
        for (int i = 0; i < 6; i++) {
          std::ostringstream key;
          key << w << "_LFold_a" << alpha << "_pt" << pt << "_" << names[i];
          CHECK_REF(refs, key.str(), S[i], 1e-8, 1e-20);
        }
        double lfratio = 0., nrratio = 0.;
        pd.Tensor_Compare_nonrel(lfratio, nrratio, alpha, pt);
        std::ostringstream key;
        key << w << "_tensor_a" << alpha << "_pt" << pt;
        CHECK_REF(refs, key.str() + "_lf", lfratio, 1e-8, 1e-14);
        CHECK_REF(refs, key.str() + "_nr", nrratio, 1e-8, 1e-14);
      }
    for (double alpha : {0.9, 1.0, 1.15})
      for (double pt : {20., 150.}) {
        double P[9];
        pd.getLFdistributions_new(alpha, pt, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
        const char *names[9] = {"P_U", "P_TLL", "P_TLT", "P_TTT", "P_SLSL", "P_STSL", "P_SLST", "P_STST_par", "P_STST_perp"};
        for (int i = 0; i < 9; i++) {
          std::ostringstream key;
          key << w << "_LF_a" << alpha << "_pt" << pt << "_" << names[i];
          CHECK_REF(refs, key.str(), P[i], 1e-8, 1e-20);
        }
      }
  }
}

TEST_CASE("Poldeut S-only, D-only and no-Melosh variants", "[poldeut]") {
  Poldeut pd("AV18", "SLAC");
  double P[9];
  pd.set_D(false);
  pd.getLFdistributions_new(1.05, 100., P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
  CHECK_REF(refs, "AV18_Sonly_P_U", P[0], 1e-8, 1e-20);
  CHECK_REF(refs, "AV18_Sonly_P_SLSL", P[4], 1e-8, 1e-20);
  pd.set_D(true);
  pd.set_S(false);
  pd.getLFdistributions_new(1.05, 100., P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
  CHECK_REF(refs, "AV18_Donly_P_U", P[0], 1e-8, 1e-20);
  CHECK_REF(refs, "AV18_Donly_P_TLL", P[1], 1e-8, 1e-20);
  pd.set_S(true);
  pd.set_melosh(false);
  pd.getLFdistributions_new(1.05, 100., P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
  CHECK_REF(refs, "AV18_noMelosh_P_SLSL", P[4], 1e-8, 1e-20);
  CHECK_REF(refs, "AV18_noMelosh_P_STST_par", P[7], 1e-8, 1e-20);
}

TEST_CASE("Poldeut observables: double-spin asymmetry and unpolarized tagged cross section", "[poldeut]") {
  Poldeut pd("AV18", "SLAC");
  for (int proton = 0; proton <= 1; proton++) {
    const std::string nuc = proton ? "p" : "n";
    double F_LSL = 0., F_U = 0.;
    pd.calc_Double_Asymm(0.3, 4.e06, 0.5, 1.02, 60., proton, F_LSL, F_U);
    CHECK_REF(refs, "Asymm_" + nuc + "_F_LSL", F_LSL, 1e-8, 1e-20);
    CHECK_REF(refs, "Asymm_" + nuc + "_F_U", F_U, 1e-8, 1e-20);
    // s = 4 E_e E_d for 10 GeV electrons on 100 GeV deuterons (in MeV^2); t' = -0.02 GeV^2
    CHECK_REF(refs, "dsigma_unpol_" + nuc, pd.calc_dsigma_unpol(0.1, 10.e06, 4.e10, 1.02, -2.e04, proton), 1e-8, 1e-20);
  }
}
