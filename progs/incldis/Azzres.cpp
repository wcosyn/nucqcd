//program that calculates FSI for resonance region for Azz observable in DIS

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

// run ../bin/inclusive beam[GeV] Q2[GeV] paris Alekhin polangle

int main(int argc, char *argv[])
{
  nucqcd::RunLog runlog(argc, argv, {"Ein [GeV]", "Q2 [GeV^2]", "deuteron wf", "structure function", "polarized beam"});
    double Ein=atof(argv[1])*1.E03;  
    double Q2=atof(argv[2])*1.E06;
    string strucname = argv[4];
    string wf = argv[3];
    bool polbeam=bool(atoi(argv[5]));
    double thetapol=0.;
    
    ElectronKinematics *elec = ElectronKinematics::CreateWithBeamEnergy(Ein);
    std::vector<double> resonances;
     resonances.push_back(1.232E03);
     resonances.push_back(1.5E03);
    resonances.push_back(1.75E03);
//     for(int i=1;i<100;i++){      
//       double x=0.01*i;
//       cout << x << " " << sqrt(Q2*(1./x-1.)+nucqcd::mass_n*nucqcd::mass_n)*1.E-03 << " " ;
//       //Kinematics2to2 kin("","",nucqcd::mass_d,nucqcd::mass_p,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
//       double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=0.;
//       for(int proton=0;proton<=1;++proton){
// 	InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,symm,4);
// 	NuclStructure nucleon(proton, Q2, x, 0, strucname);
// // 	teller+=Dinc.calc_F2Dinc(Q2,x);
// 	noemer+=nucleon.getF2();
// 	cout << nucleon.getF2() << " ";
//       }
//       cout << endl;
//      }
//     exit(1);
    for(int i=(Q2>4.E05?14:8);i<50;i++){      
      double x=0.02*i;
      //Kinematics2to2 kin("","",nucqcd::mass_d,nucqcd::mass_p,Wprime,"qsquared:wlab:pklab",1.8E06,nu,pr);
      double teller=0.,fsi1=0.,fsi2=0.,fsioff1=0.,fsioff2=0.,fsioffsuppr1=0.,fsioffsuppr2=0.,noemer=1.;
      double azz=0.,azzfsi1=0.,azzfsi2=0.,azzfsioff1=0.,azzfsioff2=0.,azzfsioffsuppr1=0.,azzfsioffsuppr2=0.;
      for(int proton=0;proton<=1;++proton){
	if(polbeam){
	  double nu=Q2/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n)*x);
	  double qvec=sqrt(Q2+nu*nu);
	  double Eout=Ein-nu;
	  double thetae=2.*asin(sqrt(Q2/(4.*Ein*Eout)));
	  double qz=Ein-Eout*cos(thetae);
	  double qx=-Eout*sin(thetae);
	  double thetaq=atan2(qx,qz); //angle between beam and qvec
	  thetapol=thetaq;
// 	  cout << thetaq << endl;
	}
	
	InclusiveCross Dinc(proton,strucname,wf,*elec,resonances,4);
	NuclStructure nucleon(proton, Q2, x, 0, strucname);
	Dinc.calc_Azzinc(azz,teller,Q2,x,thetapol);
	double f1,f2,azzf1,azzf2;
  	Dinc.calc_AzzincFSI(azzf1,azzf2,f1,f2,Q2,x,thetapol);
	//cout << "bla " << f1 << " " << f2 << endl;
	azzfsi1+=azzf1;
	azzfsi2+=azzf2;
	fsi1+=f1;
	fsi2+=f2;
	if(!Dinc.getDeutwf()->IsRelativistic()){
// 	  Dinc.calc_AzzincFSI_off(azzf1, azzf2,f1,f2,Q2,x,0.);
  // // 	cout << "bla " << f1 << " " << f2 << endl;
// 	  azzfsioff1+=azzf1;
// 	  azzfsioff2+=azzf2;
// 	  fsioff1+=f1;
// 	  fsioff2+=f2;
	  //suppressed
	  Dinc.setOffshell(3);
	  Dinc.calc_AzzincFSI_off(azzf1,azzf2,f1,f2,Q2,x,thetapol);
	  azzfsioffsuppr1+=azzf1;
	  azzfsioffsuppr2+=azzf2;
	  fsioffsuppr1+=f1;
	  fsioffsuppr2+=f2;
	}
	
      }
      cout << x << " " << sqrt(Q2*(1./x-1.)+nucqcd::mass_n*nucqcd::mass_n)*1.E-03 << " " << thetapol << " " << azz << " " << azzfsi1 << " " << 
      azzfsi2 << " " << azzfsioffsuppr1 << " " << azzfsioffsuppr2 << " " <<  teller << " " << 
      fsi1 << " " << fsi2 << " " << fsioffsuppr1 << " " << fsioffsuppr2 << endl;
      //cout << x << " " << teller << " " << noemer << " " << teller << endl;
    }
    delete elec;
}
