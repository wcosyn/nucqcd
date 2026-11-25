//Simple program to output proton and neutron F2 values

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include <Tools/Core/constants.hpp>
#include <Tools/Core/RunLog.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/StructureFunctions/NuclStructure.hpp>
#include <DIS/InclusiveCross.hpp>
using namespace nucqcd;

// make f2n
// run ../../bin/f2n Q2[GeV^2] nucleonstruct
// fi ../../bin/f2n 5. SLAC

int main(int argc, char *argv[])
{
  nucqcd::RunLog runlog(argc, argv, {"Q2 [GeV^2]", "structure function"});
    double Ein=5000.;  
    double Q2=atof(argv[1])*1.E06;
    string strucname = argv[2];
    for(int i=1;i<50;i++){      
      double x=0.02*i;
      cout << x << " " << sqrt(Q2*(1./x-1.)+nucqcd::mass_n*nucqcd::mass_n)*1.E-03 << " ";
      double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
      for(int proton=0;proton<=1;++proton){
	NuclStructure nucleon(proton, Q2, x, 0, strucname);	
	cout << nucleon.getF2() << " ";
      }
      cout << endl;
    }
}
