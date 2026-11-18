// Tabulates the neutron structure functions used by the tagged-DIS plot scripts
// (~/Calculations/poldeut/results/neutron_SF*.txt).
// Replaces the commented-out block at the top of poldeut.cpp, which used MASSn_G
// (defined in MeV in constants.hpp despite the _G suffix).
//
// run like  > neutron_sf [nx]     (nx = number of x points per Q2, x = 0.002*j, j=1..nx)
//   nx=100 -> neutron_SF.txt, nx=300 -> neutron_SF_more.txt
//
// Columns: Q2[GeV^2] | x | F1 | F2 | FT | FL | g1 | g1+g2 (WW)
//   F2: SLAC (Bodek) fit; F1 = (1+gamma^2) F2 / (2x(1+R)) with constant R = 0.18,
//   gamma^2 = 4 m^2 x^2/Q^2; FT = 2F1, FL = (1+gamma^2) F2/x - 2F1;
//   g1: GRSV2000 LO; g1+g2 from the Wandzura-Wilczek relation.
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <constants.hpp>
#include <NucleonStructure.hpp>

using namespace std;

int main(int argc, char *argv[])
{
  int nx = (argc > 1) ? atoi(argv[1]) : 100;
  const bool neutron = 0;
  const double R = 0.18;
  const double mN = MASSn * 1.E-03;   // average nucleon mass [GeV]

  NucleonStructure nucl("SLAC");
  cout << "output is for neutron SF: Q2[GeV^2]|x|F1|F2|FT|FL|g1|g1+g2 (WW approx)" << endl;
  for (int i = 0; i < 10; i++) {
    double Q2 = 2. + i;                     // GeV^2
    for (int j = 1; j <= nx; j++) {
      double x = 0.002 * j;
      double F1, F2;
      nucl.getF_xQ(F1, F2, neutron, x, Q2 * 1.E06);
      double gamma2 = pow(2. * x * mN, 2.) / Q2;
      F1 = (1. + gamma2) / (1. + R) / (2. * x) * F2;
      double FT = 2. * F1;
      double FL = (1. + gamma2) * F2 / x - 2. * F1;
      double g1 = NucleonStructure::getG1_grsv2000(neutron, x, Q2 * 1.E06);
      double g1pg2 = NucleonStructure::getG1plusg2_grsv2000(neutron, x, Q2 * 1.E06);
      cout << Q2 << " " << x << " " << F1 << " " << F2 << " " << FT << " " << FL << " "
           << g1 << " " << g1pg2 << endl;
    }
  }
  return 0;
}
