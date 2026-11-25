#include <Tools/WaveFunctions/He3wf.hpp>
#include <Tools/Core/FortranString.hpp>
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/constants.hpp>

namespace nucqcd {


using namespace std;


/// Fortran COMPLEX return value: real part dr, imaginary part di
typedef struct{double dr; /*!< real part */ double di; /*!< imaginary part */} complex_;

//fortran imports
extern"C"{
    //reads inputfile
    //DATNAME: import file
    //mbjset: spin state of he3 nucleus
    //usecdep: include PV parts or not (use FALSE) 
    void readhe3wave(int *mbjset, int *usecdep);
}

extern"C"{
    //resets he3 spin state
    //mbjset: spin state of he3 nucleus
    //usecdep: include PV parts or not (use FALSE) 
    void he3reset(int *mbjset, int *usecdep);
}

extern"C"{
    //returns wf state
    //input momenta are fm^-1 !!!!!!!!!!
     complex_ gethe3wave(double *p1x,double *p1y,double *p1z,int *ms1,int *mt1,
		      double *p2x,double *p2y,double *p2z,int *ms2,int *mt2,int *ms3,int *mt3);
}

extern"C"{
    char he3filename[256];
}




namespace {
/* m_A (twice) for which the Fortran code was last prepared (he3reset): it holds one wave function and
 * one m_A for the whole program, shared by all He3wf objects */
int fortran_mA = 0;

void prepare(int mA) {
  if (mA != fortran_mA) {
    int usecdep = 0;
    he3reset(&mA, &usecdep);
    fortran_mA = mA;
  }
}
}  // namespace

He3wf::He3wf(std::string wf_name,std::string dir){
  string file;
  if(!wf_name.compare("AV18")) file="/he_3.av18.form";
  else if(!wf_name.compare("AV18-UIX")) file="/he_3.av18+urb-ix.jmax=6.esuche.cdep_app.fad.c.form";
  else if(!wf_name.compare("AV18-TML")) file="/he_3.av18+tml5156.jmax=6.esuche.cdep_app.fad.c.form";
  else if(!wf_name.compare("CDBonn")) file="/he_3.cdbonn.jmax=6.esuche.cdep_app.fad.form";
  else if(!wf_name.compare("CDBonn-TML")) file="/he_3.cdbonn+tml4784.jmax=6.esuche.cdep_app.fad.form";
  else throw std::invalid_argument(nucqcd::ErrorMessage() << "He3wf: invalid name for the he3 wave function: " << wf_name);
  nucqcd::to_c_string(he3filename,dir+file);
  int mA=-1, usecdep=0;
  readhe3wave(&mA,&usecdep);
  he3reset(&mA,&usecdep);
  fortran_mA=mA;
}

He3wf::He3wf(){
}

He3wf::He3wf(const He3wf&){
}

He3wf& He3wf::operator=(const He3wf&){
  return *this;
}

complex<double> He3wf::getWF(const nucqcd::Vector3 &p1, const nucqcd::Vector3 &p2, int ms1, int ms2, int ms3, int mt1, int mt2, int mA){
  const double q1[3]={p1.X(),p1.Y(),p1.Z()}, q2[3]={p2.X(),p2.Y(),p2.Z()};
  return getWF(q1,q2,ms1,ms2,ms3,mt1,mt2,mA);
}

complex<double> He3wf::getWF(const double p1[3], const double p2[3], int ms1, int ms2, int ms3, int mt1, int mt2, int mA){
  if(abs(mA)!=1||abs(ms1)!=1||abs(ms2)!=1||abs(ms3)!=1||abs(mt1)!=1||abs(mt2)!=1)
    throw std::invalid_argument(nucqcd::ErrorMessage() << "He3wf::getWF: spin and isospin projections must be +-1 (twice the projection)");
  int mt3=1-mt1-mt2;  // total isospin projection of He3 is +1/2
  if(abs(mt3)!=1) return 0.;
  prepare(mA);
  // the Fortran code takes momenta in fm^-1 and returns fm^3
  double p1x=p1[0]*nucqcd::inv_hbarc;
  double p1y=p1[1]*nucqcd::inv_hbarc;
  double p1z=p1[2]*nucqcd::inv_hbarc;
  double p2x=p2[0]*nucqcd::inv_hbarc;
  double p2y=p2[1]*nucqcd::inv_hbarc;
  double p2z=p2[2]*nucqcd::inv_hbarc;  
  double p3x=-p1x-p2x;
  double p3y=-p1y-p2y;
  double p3z=-p1z-p2z;
  // the three Faddeev components: pair 12 + spectator 3, and the cyclic permutations
  complex_ wfvalue1=gethe3wave(&p1x,&p1y,&p1z,&ms1,&mt1,
			    &p2x,&p2y,&p2z,&ms2,&mt2,&ms3,&mt3);
  complex_ wfvalue2=gethe3wave(&p2x,&p2y,&p2z,&ms2,&mt2,
			    &p3x,&p3y,&p3z,&ms3,&mt3,&ms1,&mt1);
  complex_ wfvalue3=gethe3wave(&p3x,&p3y,&p3z,&ms3,&mt3,
			    &p1x,&p1y,&p1z,&ms1,&mt1,&ms2,&mt2);
  return complex<double>((wfvalue1.dr+wfvalue2.dr+wfvalue3.dr)*pow(nucqcd::inv_hbarc,3.)
			  ,(wfvalue1.di+wfvalue2.di+wfvalue3.di)*pow(nucqcd::inv_hbarc,3.));
}

}  // namespace nucqcd

