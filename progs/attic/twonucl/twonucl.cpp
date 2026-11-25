//testcode for the two nucleon code.  Never used seriously as far as I can remember...
//apparently uses things from maarten's code for the relative wf of the 2 struck nucleons...

#include <iostream>
#include <cstdlib>

using namespace std;

#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Glauber/FsiCorrelator.hpp>
#include <Glauber/FastParticle.hpp>
#include <Tools/Core/constants.hpp>
#include <Glauber/GlauberGrid.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Tools/Kinematics/TKinematics2to3.hpp>
#include <Tools/Kinematics/TElectronKinematics.hpp>
#include "DoubleNCross.hpp"
using namespace nucqcd;


int main(int argc, char *argv[])
{
  
  string homedir=nucqcd::share_dir;
  int shellindex1=0;
  int shellindex2=2;
  bool corr=atoi(argv[1]);
  cout << corr << endl;
  MeanFieldNucleusThick Carbon(1,homedir);
  TKinematics2to3WithLabAngles kin("","",7,TKinematics2to3::kN,"qsquared:wlab:pn:thn:thk:phik",1631389.67942,580,1000,15*nucqcd::deg_to_rad,31.17*nucqcd::deg_to_rad,PI,
		      Carbon.getMassA(), nucqcd::mass_p, nucqcd::mass_n, Carbon.getMassA_min_pn()+Carbon.getExcitation()[shellindex1]+Carbon.getExcitation()[shellindex2]);
  
  cout << kin.GetNumberOfSolutions() << " " << kin.GetMd() << " " << kin.GetMn() << " " << kin.GetMk() << " " << kin.GetMy() << endl;
  cout << kin.GetNumberOfSolutions() << " " << kin.GetPn() << " " << kin.GetPk() << " " << kin.GetPy() << endl;
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(4627.);
  DoubleNCross obs(*elec,&Carbon,1.E-03,2,homedir);
  double cross=obs.getDiffCross(kin,0,0,0,corr,shellindex1,shellindex2,0.);
  delete elec;
  return 0;
}


