//test program for quasi-elastic polarized observables
//used for the paper with eli for Mainz medium modification of FF

#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <Utilfunctions.hpp>
#include <NucleonEMOperator.hpp>

//run ./observables [Q2 [MeV^2]] [omega] [missing momentum]
int main(int argc, char *argv[])
{
//   double Ein=atof(argv[1]);
//   double Q2=atof(argv[2]);
//   double omega=atof(argv[3]);
//   int shell=atoi(argv[4]);
//   double pm=atof(argv[5]);
//   bool screening=0;//atoi(argv[4]);
//   double scr=1.;//atof(argv[5]);
//   int nucleus=1;//atoi(argv[6]);
//   double prec=1.E-05;//atof(argv[7]);
//   int integrator=2;//atoi(argv[8]);
//   int thick=1;//atoi(argv[9]);
//   int maxEval=20000;//atoi(argv[10]);
//   
//   string homedir="/home/wim/Code/share";
// 
//   MeanFieldNucleusThick Carbon(nucleus,homedir);
//   TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton(),nucqcd::mass_p,"qsquared:wlab:pklab",Q2,omega,pm);
//   TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
//   Cross obs(*elec,&Carbon,prec,integrator,homedir,screening,scr);
//   double free=obs.getElCross(kin,2,PI)*nucqcd::hbarc*nucqcd::hbarc;
//   vector<double> cross;
// 
// //   obs.getAllDiffCross(cross,kin,2,1,thick,0.,maxEval,1);
//   cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg << " " << kin.GetPklab() << " "  
//       /*<< cross[0] << " " << cross[1] << " " << cross[thick?4:2] << " " << free << " " << cross[0]/free << " " << cross[1]/free << endl*/;
// 
//   vector<double> observ;
//   obs.getAllObs(observ,kin,2,shell,thick,PI,maxEval,1);
//   //cout << observ[0] << " " << observ[8*(thick?4:2)] << " " << observ[3] << " " << observ[8*(thick?4:2)+3] << endl;
//   for(int i=0;i<40;i++) cout << observ[i] << " ";
//   cout << free*1.E-09 << endl;
//   delete elec;
//   return 0;

// FF_mod printout

// for(int j=0;j<60;j++){
//   for(int i=0;i<15;i++){
//     cout << j*0.05 << " " << i*0.016 << " " << NucleonEMOperator::QMCGE[j][i] << " " <<
// 	NucleonEMOperator::QMCGM[j][i] << " " <<
// 	NucleonEMOperator::CQSMGE[j][i] << " " << 
// 	NucleonEMOperator::CQSMGM[j][i] << endl;
//   }
//   cout << endl << endl << endl;
// }
// exit(1);

  double Ein=600.;
  double Eout=atof(argv[3]);
  int shell=atoi(argv[1]);
  double thetae=atof(argv[2])*nucqcd::deg_to_rad;
  //double thetap=atof(argv[3])*nucqcd::deg_to_rad;
  double omega=Ein-Eout;
  int medium=atoi(argv[4]);
  
  bool screening=0;//atoi(argv[4]);
  double scr=1.;//atof(argv[5]);
  int nucleus=1;//atoi(argv[6]);
  double prec=1.E-05;//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  
  string homedir=nucqcd::share_dir;
  MeanFieldNucleusThick Carbon(nucleus,homedir);
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
  
  double qvec=sqrt(Ein*Ein+Eout*Eout-2.*Ein*Eout*cos(thetae));
  double Q2=qvec*qvec-omega*omega;
  double thetaq=atan2(-sin(thetae)*Eout,Ein-Eout*cos(thetae));
  //thetap-=thetaq;
//   double A=omega+Carbon.getMassA();
//   double B=pow(Carbon.getMassA_min_proton()+Carbon.getExcitation()[shell],2.)+qvec*qvec-A*A-nucqcd::mass_p*nucqcd::mass_p;
//   double D=2.*qvec*cos(thetap);
//   double a=4.*A*A-D*D;
//   double b=2.*B*D;
//   double c=4.*A*A*nucqcd::mass_p*nucqcd::mass_p-B*B;
//   double discr=b*b-4.*a*c;
//   if(discr<0.) cout << "discr <0" << endl;
//   double Pp1=(-b+sqrt(discr))/(2.*a);
//   double Pp2=(-b-sqrt(discr))/(2.*a);
//   if((SIGN(B-D*Pp1)==-1)){
//     double pmx=-Pp1*sin(thetap);
//     double pmz=qvec-Pp1*cos(thetap);
//     double EAmin1=sqrt(pow(Carbon.getMassA_min_proton()+Carbon.getExcitation()[shell],2.)+pmx*pmx+pmz*pmz);
//     double pm=sqrt(pmx*pmx+pmz*pmz);
//     double Ep=sqrt(nucqcd::mass_p*nucqcd::mass_p+Pp1*Pp1);
//     cout << thetae*nucqcd::rad_to_deg << " " << thetaq*nucqcd::rad_to_deg << " " << thetap*nucqcd::rad_to_deg << " " << Pp1 << " " << pm << " " << atan2(pmx,pmz)*nucqcd::rad_to_deg << " " << A << " " << Ep+EAmin1 << endl;
//   }
//   if((SIGN(B-D*Pp2)==-1)){
//     double pmx=-Pp2*sin(thetap);
//     double pmz=qvec-Pp2*cos(thetap);
//     double EAmin1=sqrt(pow(Carbon.getMassA_min_proton()+Carbon.getExcitation()[shell],2.)+pmx*pmx+pmz*pmz);
//     double pm=sqrt(pmx*pmx+pmz*pmz);
//     double Ep=sqrt(nucqcd::mass_p*nucqcd::mass_p+Pp2*Pp2);
//     cout << thetae*nucqcd::rad_to_deg << " " << thetaq*nucqcd::rad_to_deg << " " << thetap*nucqcd::rad_to_deg << " " << Pp2 << " " << pm << " " << atan2(pmx,pmz)*nucqcd::rad_to_deg << " " << A << " " << Ep+EAmin1 << endl;
//   }
  
  
  for(int i=0;i<=20;++i){
    double pm=i*5.;
  
  
    TKinematics2to2 kin("","",Carbon.getMassA(),Carbon.getMassA_min_proton()+Carbon.getExcitation()[shell]
		,nucqcd::mass_p,"qsquared:wlab:pklab",Q2,omega,pm);
    if(kin.IsPhysical()&&kin.GetPYlab()<801.6&&kin.GetPYlab()>534.4){
        Cross obs(*elec,&Carbon,prec,integrator,homedir,screening,scr);
      double free0=obs.getElCross(kin,2,0.,maxEval)*nucqcd::hbarc*nucqcd::hbarc;
//       obs.printDensity_profile(kin,shell,thick,maxEval);
      
      vector<double> observ;
      obs.getAllObs_xyz(observ,kin,2,shell,thick,medium,0.,maxEval,1,1.,1.);
      
     cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg 
      << " " << kin.GetPklab() << " " << acos(kin.GetCosthklab())*nucqcd::rad_to_deg << " ";
      for(int i=0;i<5;i+=4) cout << observ[i*8] << " " << observ[i*8+3] << " " << observ[i*8+5] << " " << observ[i*8+7]<< " ";
      cout << free0*1.E-09 << endl;
      
      cout << kin.GetKlab() << " " << kin.GetWlab() << " " << kin.GetPYlab() << " " << acos(kin.GetCosthYlab())*nucqcd::rad_to_deg 
      << " " << -kin.GetPklab() << " " << acos(kin.GetCosthklab())*nucqcd::rad_to_deg << " ";
      
      double free180=obs.getElCross(kin,2,PI,maxEval)*nucqcd::hbarc*nucqcd::hbarc;
      obs.getAllObs_xyz(observ,kin,2,shell,thick,medium,PI,maxEval,1,1.,1.);
      for(int i=0;i<5;i+=4) cout << observ[i*8] << " " << observ[i*8+3] << " " << observ[i*8+5] << " " << observ[i*8+7]<< " ";
      cout << free180*1.E-09 << endl;
    }

  }  
    
  delete elec;
  return 0;




}


