#include <Tools/Core/RunLog.hpp>
#include <iostream>
#include <cstdlib>

using namespace std;

#include <Knockout/RhoTCross.hpp>
#include <Integrator/numint/cube.hpp>
#include <Tools/Core/IntegratorFlag.hpp>
#include <string>
using namespace nucqcd;

void klaas_t(numint::vector_d & res, double nu, double t, RhoTCross & cross, double Q2, double Ebeam){
  res=numint::vector_d(cross.getNrofcross(),0.);
  double results[cross.getNrofcross()];
  cross.getCrosst(results,Ebeam,Q2,nu,t);
  for(int i=0;i<cross.getNrofcross();++i) res[i]=results[i];

//   double qvec=sqrt(Q2+nu*nu);
//   double massA = cross.getNucleusthick().getMassA();
//   massA*=1.E-03;
//   //determine kinematics
//   double A = (t+Q2-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06)*0.5;
//   double B=(2.*massA*nu+nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-Q2)/2.+A*(massA+nu)/nu;
//   double prz=B/((massA/nu)*qvec);
//   double Erho=(prz*qvec-A)/nu;
//   double prx=sqrt(Erho*Erho-prz*prz-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);
//   double EA=sqrt(massA*massA+prx*prx+(qvec-prz)*(qvec-prz));
//   double z=Erho/nu;
//   if((!std::isnan(prx))&&(z>0.9&&z<1.)){
// //     cout << EA+Erho << " " << massA+nu << endl;
//     for(int i=0;i<cross.getNrofcross();++i) res[i]=0.5*Erho/qvec;    
//   }
//   else for(int i=0;i<cross.getNrofcross();++i) res[i]=0.;
  
}
void klaas_z(numint::vector_d & res, double nu, double z, RhoTCross & cross, double Q2, double Ebeam){
  res=numint::vector_d(cross.getNrofcross(),0.);
  double results[cross.getNrofcross()];
  cross.getCrossz(results,Ebeam,Q2,nu,z);
  for(int i=0;i<cross.getNrofcross();++i) res[i]=results[i];

  
  //   double qvec=sqrt(Q2+nu*nu);
//   double massA = cross.getNucleusthick().getMassA();
//   massA*=1.E-03;
//   double Erho=z*nu;
//   double EA=nu+massA-Erho;
//   double prho=sqrt(Erho*Erho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);
//   double pA=sqrt(EA*EA-massA*massA);
//   double prhoz=(prho*prho+qvec*qvec-pA*pA)/(2.*qvec);
//   double pAz=(-prho*prho+qvec*qvec+pA*pA)/(2.*qvec);
//   double prhox=sqrt(prho*prho-prhoz*prhoz);
//   double t = -Q2 + nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2*Erho*nu+2.*qvec*prhoz;
//   if((!std::isnan(prhox))&&(t<-0.1&&t>-0.4)) {
//     for(int i=0;i<cross.getNrofcross();++i) res[i]=nu*massA*Erho/qvec;
//     if(std::isnan(prhox)) cout << "Blaaa" << endl;
// //   cout << pA << " " << sqrt(prhox*prhox+pAz*pAz) << endl;
//   }
//   else for(int i=0;i<cross.getNrofcross();++i) res[i]=0.;
}

int main(int argc, char *argv[])
{
  nucqcd::RunLog runlog(argc, argv, {"torz", "nucleus", "Q2", "nu_min", "nu_max", "p_dil", "nocuts", "usersigma", "integrator", "maxEval", "prec", "ptype"});

  int torz = atoi(argv[1]);
  int nucleus = atoi(argv[2]);
  double Q2 = atof(argv[3]);
  double nu_min = atof(argv[4]);
  double nu_max = atof(argv[5]);
  double p_dil = atof(argv[6]);
  int nocuts = atoi(argv[7]);
  double usersigma = atof(argv[8]);
  IntegrationMethod integrator = integration_method(atoi(argv[9]));  // 1: cube_romb, 2: cube_adaptive (0, the old nested rombergerN, was removed)
  int maxEval = atoi(argv[10]);
  double prec=atof(argv[11]);
  int ptype=atoi(argv[12]);
//   int maxEval2 = atoi(argv[11]);
  string homedir=nucqcd::share_dir;
  double Ebeam = 4.023;
  if(!nocuts) RhoTCross::checkWcut(Q2,nu_min);  // the analysis cuts W > 2 GeV (and z > 0.9, t in -0.4..-0.1, applied in RhoTCross)
  

  RhoTCross test = RhoTCross(nucleus,330,homedir,nocuts,1,usersigma,prec,IntegrationMethod::Adaptive,maxEval,p_dil,ptype);  
  
  {
    /*! struct that is used for integrators (clean ones)*/
    struct Ftor_tz {

      /*! integrandum function */
      static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
	Ftor_tz &p = * (Ftor_tz *) param;
	p.f(ret,x[0],x[1],*p.cross,p.Q2,p.Ebeam);
      }
      RhoTCross *cross;/*!< pointer to the grid where the integration is performed */
      double Q2;/*!< [GeV^2] Qsquared */
      double Ebeam; /*!< [GeV] rho momentum */
      /*! integrandum 
      * \param res results
      * \param pm first integration variable
      * \param costheta second integration variable
      * \param phi third integration variable
      * \param cross the RhoTCross instance
      */
      void (*f)(numint::vector_d & res, double nu, double other, RhoTCross & cross, double Q2, double Ebeam);
    };
    int res=90;
    unsigned count=0;
    numint::vector_d ret(test.getNrofcross(),0.);
    numint::array<double,2> lower = {{nu_min,torz?-0.4:0.9}};
    numint::array<double,2> upper = {{nu_max,torz?-0.1:1.}};
    
    Ftor_tz F;
    F.cross = &test;
    F.Q2=Q2;
    F.Ebeam=Ebeam;
    numint::mdfunction<numint::vector_d,2> mdf;
    mdf.func = &Ftor_tz::exec;
    mdf.param = &F;
    F.f=torz?klaas_t:klaas_z;
    res = numint::cube({integrator,1.E-12,1.E-03,100,20000},mdf,lower,upper,ret,count,0);
    cout << Q2 << " ";
    for (int l=0;l<test.getNrofcross();l++) cout << ret[l] << " ";
    cout << count << " " << res << endl;
  }
  
  
}
