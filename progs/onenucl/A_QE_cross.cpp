//program computes A(e,e'p) cross section.  Used for several calculations (Monaghan's paper for instance)

#include <Tools/Core/RunLog.hpp>
#include <iostream>
#include <cstdlib>

using namespace std;

#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Glauber/FsiCorrelator.hpp>
#include <Glauber/FastParticle.hpp>
#include <Tools/Core/constants.hpp>
#include <Glauber/GlauberGrid.hpp>
#include <Glauber/OneGlauberGrid.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Knockout/EMQECross.hpp>
using namespace nucqcd;


//run ./A_QE_cross [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{
  nucqcd::RunLog runlog(argc, argv, {"Q2", "omega", "pm"});
  double Q2=atof(argv[1]);
  double omega=atof(argv[2]);
  double pm=atof(argv[3]);
  bool screening=0;//atoi(argv[4]);
  double scr=1.;//atof(argv[5]);
  int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  
  string homedir=nucqcd::share_dir;

  MeanFieldNucleusThick Carbon(nucleus,homedir);
  Kinematics2to2 kin(Carbon.getMassA(),Carbon.getMassA_min_proton(),nucqcd::mass_p,"qsquared:wlab:pklab",Q2,omega,pm);
//   ElectronKinematics *elec = ElectronKinematics::CreateWithBeamEnergy(4627.);  //Monaghan Data
  ElectronKinematics *elec = ElectronKinematics::CreateWithBeamEnergy(3245.); //Dutta Data
  EMQECross obs(*elec,&Carbon,prec,integrator,homedir,screening,scr);
  double free=obs.getElCross(kin,2,0.,2E03)*nucqcd::hbarc*nucqcd::hbarc;
  vector<double> cross;

// comparison with pascal's code (input_compare_onenucl.dat input file) 
//  ../../bin/onenucl 1696000.0 851.5 190.0 with dutta beam inq
//   double crossp=obs.getDiffCross(kin, 2, 0, 0, 0, 1, 1, 0.,2E05,1);
//   double crosss=obs.getDiffCross(kin, 2, 0, 0, 0, 1,0,0.,2E04,1);
//   cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << 
//   acos(kin.GetCosthYlab())*nucqcd::rad_to_deg << " " << kin.GetPklab() << " " << acos(kin.GetCosthklab())*nucqcd::rad_to_deg << endl;
//   cout << sqrt(Q2+omega*omega) << " " << crossp << " " << crosss << endl;

  MeanFieldNucleusThick Fe(3,homedir);
  Kinematics2to2 kin2(Fe.getMassA(),Fe.getMassA_min_proton(),nucqcd::mass_p,"qsquared:wlab:pklab",Q2,omega,pm);
  EMQECross obs2(*elec,&Fe,prec,integrator,homedir,screening,scr);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 0, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 1, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 2, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 3, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 4, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 5, 0.,2E04,1,0,1.,1.);
  obs2.getDiffCross(kin2, 2, 0, 0, 0, 1, 6, 0.,2E04,1,0,1.,1.);
  
  exit(1);
  

 /* obs.getAllDiffCross(cross,kin,2,1,thick,0.,maxEval,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg << " " << kin.GetPklab() << " " <<  
      cross[0] << " " << cross[1] << " " << cross[4] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl;
  
      
  free=obs.getElCross(kin,2,PI,2E03)*nucqcd::hbarc*nucqcd::hbarc;
  obs.getAllDiffCross(cross,kin,2,1,thick,PI,maxEval,1);
  cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg << " " << -kin.GetPklab() << " " <<  
      cross[0] << " " << cross[1] << " " << cross[4] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl;
 */     
      
  vector<double> observ;
//   obs.getAllObs_tnl(observ,kin,2,1,thick,0,1.,2E04,1);
//   cout << cross[4] << " " << observ[8*4] << endl;
//   for(int i=0;i<8;i++) cout << observ[i] << endl;
//   for(int i=0;i<8;i++) cout << observ[32+i] << endl;
  //cout << observ[0] << " " << observ[8*(thick?4:2)] << " " << observ[3] << " " << observ[8*(thick?4:2)+3] << endl;
//   for(int i=0;i<5;i++) {
//     for(int j=0;j<8;j++) cout << observ[i*8+j] << " ";
//     cout << endl;
//   }
      //   double crossp=obs.getDiffCross(kin, 2, 0, 0, 0, 0, 1, 0.);
// // //   double crosss=obs.getDiffCross(kin, 0, 0, 0, 0, 0.);
//   double crosspsrc=obs.getDiffCross(kin, 2, 0, 1, 0, 0, 1, 0.);
// // //   double crosspct=obs.getDiffCross(kin, 0, 1, 0, 1, 0.);
// // //   double crosspsrcct=obs.getDiffCross(kin, 1, 1, 0, 1, 0.);
// // //   double crossppw=obs.getDiffCross(kin, 0, 0, 1, 1, 0.);
// // //   cout << free/nucqcd::hbarc/nucqcd::hbarc << endl;
//   cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg << " " << kin.GetPklab() << " " <<  
//       crossp << " " << crosspsrc << " " << free << " " << crossp/free << " " << crosspsrc/free << endl;
// 
//       //   cout << crossp << " " << crosss << " " << crosspsrc << " " << crosspct << " " << crosspsrcct << " " << crossppw << endl;
//   //cout << obs.getElCross(kin,0.) << " " << obs.getElCross(kin,PI) <<endl;
// //   EMQEHadronCurrent test(&Carbon,0,homedir);
// //   cout << test.getMatrixEl(kin,-1,0,0,-1,1,1) << endl;
// //   
  delete elec;
  return 0;
}


