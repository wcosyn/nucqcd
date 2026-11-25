// Regression: inclusive deuteron DIS with FSI and tensor asymmetry (DIS/InclusiveCross);
// kinematics as in progs/incldis/disincres.cpp.
#include <test_support.hpp>

#include <DIS/InclusiveCross.hpp>
#include <Tools/Kinematics/ElectronKinematics.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <vector>
using namespace nucqcd;

namespace {
nucqcd_test::RefStore refs("reg_inclusive_dis");
}

TEST_CASE("Inclusive F2D: plane wave and FSI (resonance model)", "[inclusive]") {
  std::unique_ptr<ElectronKinematics> elec(ElectronKinematics::CreateWithBeamEnergy(5000.));
  std::vector<double> res = {1.232E03, 1.5E03, 1.75E03};
  const double Q2 = 2.e06;
  for (int proton = 0; proton <= 1; proton++) {
    InclusiveCross dinc(proton, "SLAC", "AV18", *elec, res, 3);
    for (double x : {0.3, 0.6}) {
      std::ostringstream key;
      key << (proton ? "p" : "n") << "_x" << x;
      CHECK_REF(refs, key.str() + "_F2D_pw", dinc.calc_F2Dinc(Q2, x), 1e-7);
      double fsi1 = 0., fsi2 = 0.;
      dinc.calc_F2DincFSI(fsi1, fsi2, Q2, x);
      CHECK_REF(refs, key.str() + "_F2D_fsi1", fsi1, 1e-6, 1e-12);
      CHECK_REF(refs, key.str() + "_F2D_fsi2", fsi2, 1e-6, 1e-12);
      double off1 = 0., off2 = 0.;
      dinc.calc_F2DincFSI_off(off1, off2, Q2, x);
      CHECK_REF(refs, key.str() + "_F2D_fsioff1", off1, 1e-6, 1e-12);
      CHECK_REF(refs, key.str() + "_F2D_fsioff2", off2, 1e-6, 1e-12);
    }
  }
}

TEST_CASE("Inclusive tensor asymmetry Azz (plane wave)", "[inclusive]") {
  std::unique_ptr<ElectronKinematics> elec(ElectronKinematics::CreateWithBeamEnergy(11000.));
  std::vector<double> res = {1.232E03, 1.5E03, 1.75E03};
  for (int proton = 0; proton <= 1; proton++) {
    InclusiveCross dinc(proton, "CB", "AV18", *elec, res, 3);
    for (double x : {0.4, 0.8}) {
      double Azz = 0., cross = 0.;
      dinc.calc_Azzinc(Azz, cross, 2.e06, x, 0.);
      std::ostringstream key;
      key << (proton ? "p" : "n") << "_x" << x;
      CHECK_REF(refs, key.str() + "_Azz", Azz, 1e-7, 1e-14);
      CHECK_REF(refs, key.str() + "_Azz_cross", cross, 1e-7, 1e-14);
    }
  }
}
