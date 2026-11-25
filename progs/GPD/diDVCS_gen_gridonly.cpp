#include <string>
#include <Tools/Core/RunLog.hpp>
#include "event_random.hpp"
#include <Tools/Core/constants.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include <complex>



using namespace std;

#define PI 3.1415926535897932385 /*!< \def defines constant pi */
#define ALPHA 0.00729735253 /*!< \def defines finestructure constant [] */


#include <Tools/Core/Vector3.hpp>
#include <Tools/Core/LorentzVector.hpp>
#include <iostream>
#include <fstream>

double me=511e-6;
double Mp=0.938272;
double mpi=0.13957;
double MassRho=0.77545;
double GammaRho=0.145;

double s1_min;
double s1_max;
double th_gen;
double mrho;
double psf;


nucqcd::LorentzVector BeamL;
nucqcd::LorentzVector BeamH;
nucqcd::LorentzVector Rho0,eplus,eminus,piplus,piminus;
nucqcd::LorentzVector QPrime;
nucqcd::LorentzVector Proton;
nucqcd::LorentzVector ScattL;
nucqcd::LorentzVector q;
nucqcd::LorentzVector cms;
nucqcd::LorentzVector S2;


double s;


double theta_e,Ep, t_min;
double s2_min;

double xb,q2,s2,trho,tN,qprime2,y,q2p, trho_min;


double Ecm,pcm; //Energy/momentum in center-of-mass of (rho+q')
double Ecm_pr,pcm_pr;  //Energy/momentum in center-of-mass of (p')
double Ecm_rho,pcm_rho; //Energy/momentum in center-of-mass of (rho)

nucqcd::Vector3 dir_pr;
nucqcd::Vector3 dir_qp;
nucqcd::Vector3 rot_dir;

//To generate kinematics
double theta_max=177.9*nucqcd::deg_to_rad;
double W2min=4;//(q+p)^2 in GeV^2
double q2min=1;
double q2max=10;
double q2pmin=2.;
double trange=1;

//interpolation of the CFF grid
void CFF_interpolate(complex<double> &CFF_H, complex<double> &CFF_E, double Q2out_over_Q2in, double xi, double tN, complex<double> ****grid){
    double index_Q=0.,index_xi=0.,index_t=0.;
    double frac_Q=0.,frac_xi=0.,frac_t=0. ;
    frac_xi=modf(xi*100-1,&index_xi);
    frac_t=modf(-tN*20.-1,&index_t);
    frac_Q=modf(4.*log10(Q2out_over_Q2in)+12.,&index_Q);  //grid is logarithmic 
    //cout << index_Q << " " << index_t << " " << index_xi << " " << Q2out_over_Q2in << " " << xi << " " << tN << endl;

    CFF_H=grid[int(index_Q)][int(index_t)][int(index_xi)][0]*(1.-frac_Q)*(1.-frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)][int(index_xi)][0]*(frac_Q)*(1.-frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)][int(index_t)+1][int(index_xi)][0]*(1.-frac_Q)*(frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)+1][int(index_xi)][0]*(frac_Q)*(frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)][int(index_t)][int(index_xi)+1][0]*(1.-frac_Q)*(1.-frac_t)*(frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)][int(index_xi)+1][0]*(frac_Q)*(1.-frac_t)*(frac_xi)
                                        +grid[int(index_Q)][int(index_t)+1][int(index_xi)+1][0]*(1.-frac_Q)*(frac_t)*(frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)+1][int(index_xi)+1][0]*(frac_Q)*(frac_t)*(frac_xi)
                                        ;
    CFF_E=grid[int(index_Q)][int(index_t)][int(index_xi)][1]*(1.-frac_Q)*(1.-frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)][int(index_xi)][1]*(frac_Q)*(1.-frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)][int(index_t)+1][int(index_xi)][1]*(1.-frac_Q)*(frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)+1][int(index_xi)][1]*(frac_Q)*(frac_t)*(1-frac_xi)
                                        +grid[int(index_Q)][int(index_t)][int(index_xi)+1][1]*(1.-frac_Q)*(1.-frac_t)*(frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)][int(index_xi)+1][1]*(frac_Q)*(1.-frac_t)*(frac_xi)
                                        +grid[int(index_Q)][int(index_t)+1][int(index_xi)+1][1]*(1.-frac_Q)*(frac_t)*(frac_xi)
                                        +grid[int(index_Q)+1][int(index_t)+1][int(index_xi)+1][1]*(frac_Q)*(frac_t)*(frac_xi)
                                        ;
    //cout << "interp " << CFF_H << " "<< CFF_E <<endl;

}


//result is electroproduction cross section in nb/GeV^10, interpolated from CFF grid.
// invariants input in GeV^2
double getCross_DiDVCS_interp(const double scale, const double t_rho, const double t_N, const double Q2in, const double Q2out, 
                          const double s2, const double s_eN, double y,complex<double> ****grid){

    double frho0 = 0.216; //rho0 decay constant [GeV]
    double alpha_s = 0.306039; //scale is 2 GeV^2
    int Nc=3; //number of colors
    double CF=(Nc*Nc-1)/2./Nc;

    //kinematics
    double s_gN = s_eN*y;
    double s1=s_gN*Q2out/s2;
    double xi = (Q2in +s_gN*(Q2out/s2))/(2.*s_gN+Q2in -s_gN*(Q2out/s2));
//    cout << "xi =" << xi << endl;

    double K = pow(2,7.)*PI*PI*ALPHA*CF*CF*alpha_s*alpha_s/8.; //[dimensionless]
    
    complex<double> CFF_H, CFF_E;
    CFF_interpolate(CFF_H, CFF_E,Q2out/Q2in,xi,t_N,grid);
    //cout << "CFF " << CFF_H << " " << CFF_E << endl;
    double CFF_asympt = pow(3.*sqrt(2.)*K*frho0/Q2in/sqrt(Q2in*Q2out),2.)*((1-xi*xi)*norm(CFF_H)
                        -2.*xi*xi*((conj(CFF_H)*CFF_E).real())-(xi*xi+t_N/4./Mp/Mp)*norm(CFF_E)); //[GeV^-6]

    double result = CFF_asympt*ALPHA/pow(1+xi,2.)/(48*pow(2.*PI,4.)*Q2out*(s2-t_rho)); //[GeV-10]
    result *=0.389379E06; // photoprod cross section [nb GeV-8]
    double epsilon = (1-y)/(1-y+y*y/2.);
    result *= ALPHA/2./PI*y/Q2in/(1-epsilon); //electroproduction [nb GeV-10]
    return result;
}



//Method to compute the lower and upper bound of the squared momentum transfer to the recoil baryon (Proton, Nucleus, Delta...). It depends on q2, xb, the mass of the recoil baryon mrecoil and the mass of the produced particle mpp 
double tmin(double q, double x, double mrecoil, double mpp){
 double W=std::sqrt(q*(1/x-1)+Mp*Mp);
 double E1cm=(W*W+q+Mp*Mp)/2./W;
 double P1cm=std::sqrt(E1cm*E1cm-Mp*Mp);
 double E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 double P3cm=std::sqrt(E3cm*E3cm-mrecoil*mrecoil);
 double result=(std::pow(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-std::pow(P1cm-P3cm,2));
 return result;
}

double tmax(double q, double x, double mrecoil, double mpp){
 double W=std::sqrt(q*(1/x-1)+Mp*Mp);
 double E1cm=(W*W+q+Mp*Mp)/2./W;
 double P1cm=std::sqrt(E1cm*E1cm-Mp*Mp);
 double E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 double P3cm=std::sqrt(E3cm*E3cm-mrecoil*mrecoil);
 double result=(std::pow(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-std::pow(P1cm+P3cm,2));
 return result;
}

//y the inelasticity which is more convenient to use than xb.
double ymax(double q2){
  return std::min(1-q2/4./BeamL.E()/BeamL.E()/std::sin(theta_max/2.)/std::sin(theta_max/2.),0.95);
}

double ymin(double q2){
  return std::max((W2min-Mp*Mp+q2)/s,0.01);
}


//THE method which computes all LorentzVector. Q2 and xb has already been set and are kept in memory. Then remains t and phi to define the entire final set
//In truth... it must be validated on the recoil baryon
void ComputeLeptScatt(double q2, double y){
  theta_e=2*std::atan(std::sqrt(q2/(1-y)/4./BeamL.E()/BeamL.E()));
  Ep=q2/4./BeamL.E()/std::sin(theta_e/2.)/std::sin(theta_e/2.);
 ScattL.SetPxPyPzE(Ep*std::sin(nucqcd::pi+theta_e),0,Ep*std::cos(nucqcd::pi+theta_e),Ep);
 
  q=BeamL-ScattL;
  cms=q+BeamH;
  //q2=sxy is approximative... need to readjust the value of xb for consistency
  xb=-q.Mag2()/2./(BeamH*q);
   //The collision is easier to understand in the center-of-mass frame
   ScattL.Boost(-cms.BoostVector());
   BeamH.Boost(-cms.BoostVector());
   BeamL.Boost(-cms.BoostVector());
   q.Boost(-cms.BoostVector());
}

double tmin_proton(double s2, double q2prime){
  double Es2_pr=(s2-q2prime+Mp*Mp)/2./std::sqrt(s2);//Energy in cms of s2 for the proton
  BeamH.Boost(-S2.BoostVector()); //Let's check incoming proton momentum in S2 rest frame
  double tN_min=2*Mp*Mp-2*(Es2_pr*BeamH.E()-BeamH.P()*std::sqrt(Es2_pr*Es2_pr-Mp*Mp)); //t=(p-p')^2 minimal when p and p' are collinear 
  BeamH.Boost(S2.BoostVector()); //Let's put Incoming proton back in gamma*N cms.
  //Tmin is obtained when incoming and recoil proton have parallel momentum.
  return tN_min;
}


void ComputePrQprime(double s2, double q2prime, double t_N, double phi_trento){
  double Es2_pr=(s2-q2prime+Mp*Mp)/2./std::sqrt(s2);//Energy in cms of s2 for the proton
  double Es2_qp=(s2+q2prime-Mp*Mp)/2./std::sqrt(s2);//Energy in cms of s2 for the q2'
  double ps2=std::sqrt(((s2-pow(Mp+sqrt(q2prime),2))*(s2-pow(Mp-sqrt(q2prime),2))))/2./std::sqrt(s2);//momentum of both particle in particle s2
  
  //Let's go to S2 cms
  BeamH.Boost(-S2.BoostVector());

  double theta=acos(((t_N-2*Mp*Mp)/2.+BeamH.E()*Es2_pr)/BeamH.Vect().Mag()/ps2);
 
  nucqcd::Vector3 unit=BeamH.Vect().Unit();
  unit.Rotate(theta,unit.Orthogonal());
  unit.Rotate(phi_trento,BeamH.Vect());
  Proton.SetPxPyPzE(ps2*unit.X(),ps2*unit.Y(),ps2*unit.Z(),Es2_pr);
  QPrime.SetPxPyPzE(-Proton.Px(),-Proton.Py(),-Proton.Pz(),Es2_qp);

  //Back to gamma*p cms frame
  BeamH.Boost(S2.BoostVector());
  Proton.Boost(S2.BoostVector());
  QPrime.Boost(S2.BoostVector());
  
  //Back to Lab frame
  BeamH.Boost(cms.BoostVector());
  Proton.Boost(cms.BoostVector());
  QPrime.Boost(cms.BoostVector());
  Rho0.Boost(cms.BoostVector());
  ScattL.Boost(cms.BoostVector());
  BeamL.Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());
 }

void ComputeRho(double mpp, double mrecoil, double t_rho, double phi_trento){
  Ecm=(cms.Mag2()-mrecoil*mrecoil+mpp*mpp)/2./cms.Mag();
  pcm=std::sqrt(Ecm*Ecm-mpp*mpp);
  
  Ecm_pr=(cms.Mag2()+mrecoil*mrecoil-mpp*mpp)/2./cms.Mag();
  pcm_pr=std::sqrt(Ecm_pr*Ecm_pr-mrecoil*mrecoil);


   //t=2q_cm q'_cm(1-cos (theta_gamma*/meson))
  
  double theta_cm=acos(1-(tmin(-q.Mag2(), xb, mrecoil, mpp)-t_rho)/2./pcm/q.Vect().Mag());
   //Phi is a rotation of all four momentum around the axis defined by the virtual photon.
  nucqcd::Vector3 dir_pr=q.Vect().Unit();
  nucqcd::Vector3 dir_qp=q.Vect().Unit();
  nucqcd::Vector3 rot_dir=q.Vect().Orthogonal().Unit();
   dir_pr.Rotate(theta_cm,rot_dir);
   dir_pr.Rotate(phi_trento,q.Vect().Unit());
    Rho0.SetPxPyPzE(pcm*dir_pr.X(),pcm*dir_pr.Y(),pcm*dir_pr.Z(),Ecm);

    S2=q+BeamH-Rho0;
}

void GammaPrimeDecay(){
  double th_gen=event_random::sin2_angle();
  double pel=std::sqrt(QPrime.Mag2()/4.-me*me);
  nucqcd::Vector3 unit_dir=QPrime.Vect().Unit();
  unit_dir.Rotate(th_gen,unit_dir.Orthogonal()); //theta angle with photon direction
  unit_dir.Rotate(event_random::uniform(0,2*nucqcd::pi),QPrime.Vect()); //rotation around the photon direction
  eminus.SetPxPyPzE(pel*unit_dir.X(),pel*unit_dir.Y(),pel*unit_dir.Z(),QPrime.Mag()/2.);
  eplus.SetPxPyPzE(-pel*unit_dir.X(),-pel*unit_dir.Y(),-pel*unit_dir.Z(),QPrime.Mag()/2.);

  //BackToLabFrame
  eminus.Boost(QPrime.BoostVector());
   eplus.Boost(QPrime.BoostVector());
  
}

void RhoDecay(){
  double th_gen_rho=event_random::cos2_angle();
  double ppi=std::sqrt(mrho*mrho/4.-mpi*mpi);
  nucqcd::Vector3 unit_dir=Rho0.Vect().Unit();
  unit_dir.Rotate(th_gen_rho,unit_dir.Orthogonal()); //theta angle with rho direction
  unit_dir.Rotate(event_random::uniform(0,2*nucqcd::pi),Rho0.Vect()); //rotation around the rho direction
  piminus.SetPxPyPzE(ppi*unit_dir.X(),ppi*unit_dir.Y(),ppi*unit_dir.Z(),mrho/2.);
  piplus.SetPxPyPzE(-ppi*unit_dir.X(),-ppi*unit_dir.Y(),-ppi*unit_dir.Z(),mrho/2.);

  //BackToLabFrame
  piminus.Boost(Rho0.BoostVector());
   piplus.Boost(Rho0.BoostVector());
}

//This function is the main function...
//ScattL is the scattered lepton
//Proton is the recoiled proton
//eplus and eminus are the decay product of q'
//piplus piminus are the decay product of rho0
void GetParticles(double EneL, double EneH, double q2, double y, double q2prime, double s2, double t_N, double t_rho){
  //Set Beam particle
//   BeamL.SetPxPyPzE(0,0,-EneL,EneL);
//   BeamH.SetPxPyPzE(0,0,EneH,std::sqrt(Mp*Mp+EneH*EneH));
//   s=(BeamL+BeamH).Mag2();
  ComputeLeptScatt(q2,y);//Compute ScattL and xb... all 4 momenta are exprressed in gamma*p cms.
  
//   while (mrho<0.4||mrho>1.2) mrho=event_random::breit_wigner(MassRho,GammaRho); 
  
  ComputeRho(mrho,std::sqrt(s2),t_rho,0); //Give rho 4-momentum and last argument is the angle around the virtual photon in gamma*p cms.
  
  ComputePrQprime(s2, q2prime, t_N, 3.141); //Give rho 4-momenta of p' and gamma'. last argument is the angle of p' gamma' around the incoming hadron in gamma*p cms. ALL THE 4-MOMENTA ARE NOW EXPRESSED in LAB FRAME.

  GammaPrimeDecay();
  RhoDecay();

  /*cout<<t_N<<" "<<(Proton-BeamH).Mag2()<<endl;
  cout<<t_rho<<" "<<(Rho0-q).Mag2()<<endl;
  cout<<q2prime<<" "<<QPrime.Mag2()<<" "<<(eplus+eminus).Mag2()<<endl;
  cout<<s2<<" "<<(QPrime+Proton).Mag2()<<endl;
  cout<<q2<<" "<<-(BeamL-ScattL).Mag2()<<endl;
  ScattL.Print();
  eplus.Print();
  eminus.Print();
  piplus.Print();
  piminus.Print();
  Proton.Print();*/
 

  return;
}

//num_events is the number of points in q2/xb/t_rho/t_n... ndec is the number of decays for the rho0 and q' photon we perform for each q2/xb/t_rho/t_n. In the end, the total number of particle is ndec*num_events
void GenerateEvents(double EneL, double EneH, int nevents, int ndec, complex<double> ****grid){
  //Set Beam particle
  BeamL.SetPxPyPzE(0,0,-EneL,EneL);
  BeamH.SetPxPyPzE(0,0,EneH,std::sqrt(Mp*Mp+EneH*EneH));
  s=(BeamL+BeamH).Mag2();

  event_random::engine();  // seeds the generator (new seed every run, or NUCQCD_SEED), as SetSeed(0) did
  
  for (int ev=0;ev<nevents;ev++){
    q2=event_random::uniform(1,10); //Generate Q^2 between 1 and 10
    y=event_random::uniform(ymin(q2),ymax(q2)); //Generate y between ymin and ymax
    ComputeLeptScatt(q2,y);//Compute ScattL and xb... all 4 momenta are exprressed in gamma*p cms.
    
    //Now that gamma* is generated, we focus on the rho (mrho,trho)
    while (mrho<0.4||mrho>1.2) mrho=event_random::breit_wigner(MassRho,GammaRho);
    s1_min=10*std::pow(std::sqrt(q2pmin)+Mp,2);//Here it is 10*s2_min with q2p_min=2
    s1_max=std::pow(std::sqrt(q2*(1/xb-1)+Mp*Mp)-Mp,2);
    if (s1_max>s1_min){
      s2_min=std::pow(std::sqrt(q2pmin)+Mp,2);
      s2=event_random::uniform(s2_min,s1_max/10.);
      trho_min=tmin(q2, xb, sqrt(s2), mrho);
      trho=event_random::uniform(trho_min-trange,trho_min);
      ComputeRho(mrho,std::sqrt(s2),trho,0); //Give rho 4-momentum and last argument is the angle around the virtual photon in gamma*p cms.

      //Now that we have s2, we need to generate t_N to have the final state
      q2p=event_random::uniform(q2pmin,min(std::pow(std::sqrt(s2)-Mp,2),20.)/*std::pow(std::sqrt(s2)-Mp,2)*/); //Knowing the s2, q2' and MP, we can derive the momentum of the proton in the cms
      tN=event_random::uniform(tmin_proton(s2,q2p)-trange,tmin_proton(s2,q2p));
        ComputePrQprime(s2, q2p, tN, 3.141); //Last is phi trento
      if(tN>=-0.55){   //limit tNrange...
        for (int dec=0;dec<ndec;dec++){
            RhoDecay();
            GammaPrimeDecay();
            psf=(q2max-q2min)*(ymax(q2)-ymin(q2))*(s1_max/10.-s2_min)*1*1*(std::pow(std::sqrt(s2)-Mp,2)-q2pmin);//psf= deltaQ2 * delta Y * delta s2 * delta t_rho * delta t_N * delta Q'2 () => Be careful about angle around beams (ScattL), virtual photon (Rho0/S2 direction) and phi trento(angle between Pomeron*N and N' Q' plane) 
            // cout<<"//////////////////////"<<endl;
            // ScattL.Print();
            // Proton.Print();
            // QPrime.Print();
            // Rho0.Print();
            // (BeamL+BeamH-eplus-eminus-piplus-piminus-ScattL-Proton).Print();
            double s_gN = y*s;
            double xi = (q2 +s_gN*(q2p/s2))/(2.*s_gN+q2 -s_gN*(q2p/s2));
            double s1 = s_gN*q2p/s2;

            //cross section electroproduction [nb/GeV^-10], *integrated* over all phi angles
            double result_interp = getCross_DiDVCS_interp(2.,trho,tN,q2,q2p,s2,s,y,grid)*exp(-2.*(trho_min-trho));  
            //////////////////OUTPUT to text file to write here
            //===> Need to output all 4 monenta in text file here with a weight=xsec*psf/num_events/ndec and psf= deltaQ2 * delta Y * delta s2 * delta t_rho * delta t_N * delta Q'2 ()
            cout << "ScattL//Proton//Qprime//Rho0" << endl;  
            ScattL.Print();        
            Proton.Print();
            QPrime.Print();
            Rho0.Print();
           cout << "weight " << result_interp*psf/nevents/ndec << endl;
           cout << "=========================" << endl;
        }
      }
    }
  }
}


/*
 * Main function.
 */
int main(int argc, char *argv[]) {
    nucqcd::RunLog runlog(argc, argv, {"electron beam momentum [GeV]", "proton beam momentum [GeV]", "nevents", "ndec"});


    //double xi=atof(argv[1]);
    //cout << argv[1] << endl;
    double El_beam = atof(argv[1]);  //electron momentum : 18 GeV
    double p_beam = atof(argv[2]);    //proton momentum: 275 GeV
    int nevents = atoi(argv[3]);
    int ndec = atoi(argv[4]);

    complex<double> ****CFF_grid = new complex<double> ***[25];
    ifstream input("gridtest_doubletmin");  //CHANGE PATH IF NECESSARY
    double dummy,real,imag;
    if(input.is_open()){
        for(int i=0;i<=24;i++){
            CFF_grid[i] = new complex<double> **[11];
            for(int j=0;j<=10;j++){
                CFF_grid[i][j] = new complex<double> *[99];
                for(int k=0;k<99;k++){
                    CFF_grid[i][j][k] = new complex<double>[2];
                    input >> dummy >> dummy >> dummy >> real >> imag;
                    CFF_grid[i][j][k][0]=complex<double>(real,imag);
                    input >> real >> imag;
                    CFF_grid[i][j][k][1]=complex<double>(real,imag);
                    //cout << i << " " << j << " " << k << " " << CFF_grid[i][j][k][0] << " " << CFF_grid[i][j][k][1] << endl;
                }
            }
        }
    }
    else cerr << "file not found" << endl;

    GenerateEvents(El_beam, p_beam, nevents, ndec, CFF_grid);
    for(int i=0;i<=24;i++){            
        for(int j=0;j<=10;j++){                
            for(int k=0;k<99;k++){
                delete [] CFF_grid[i][j][k];
            }
            delete [] CFF_grid[i][j];
        }
        delete [] CFF_grid[i];
    }
        

        // GetParticles(El_beam,p_beam,Q2in,y,Q2out,s2,t_N,t_rho);

        // cout << results[0] << " " << results[1] << endl;
        // cout << "Rho "; Rho0.Print();
        // cout << "Qprime "; QPrime.Print();
        // cout << "eprime "; ScattL.Print();
        // cout << "proton "; Proton.Print();
        


        

    return 0;
}
