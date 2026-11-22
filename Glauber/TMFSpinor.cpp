#include "TMFSpinor.hpp"
#include <constants.hpp>
//#include <numtoa.h>
// #include <TString.h>
#include <FourVector.h>
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
// #include <climits>
#include <iostream>
#include <complex>
#include <TDeuteron.h>

using TMath::Sign;
using namespace std;


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi)
: fComponent(new Matrix<4,1>){

  double Gr=nucleus.getWave_G(shellindex,r);
  double Fr=nucleus.getWave_F(shellindex,r);
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi)
: fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, const double *Spher_harm, double phi)
: fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*Spher_harm[0];
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getL_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*Spher_harm[1];
  (*fComponent)(2,0)=-Fr*spinup*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m-1,1,1, nucleus.getJ_array()[shellindex],m)*Spher_harm[2];
  (*fComponent)(3,0)=-Fr*spindown*TDeuteron::Wavefunction::ClebschGordan(2*nucleus.getLbar_array()[shellindex],m+1,1,-1, nucleus.getJ_array()[shellindex],m)*Spher_harm[3];
  
}

TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double r, double costheta, double phi, Matrix<4,1> CG)
: fComponent(new Matrix<4,1>){

  double Gr=nucleus.getWave_G(shellindex,r);
  double Fr=nucleus.getWave_F(shellindex,r);
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*CG(0,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*CG(1,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*CG(2,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*CG(3,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
  
}


TMFSpinor::TMFSpinor(const MeanFieldNucleus &nucleus, int shellindex, int m, double Gr, double Fr, double costheta, double phi, Matrix<4,1> CG)
: fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*CG(0,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*CG(1,0)*
		      nucleus.Spher_Harm(nucleus.getL_array()[shellindex], (m+1)/2, costheta);
  (*fComponent)(2,0)=-Fr*spinup*CG(2,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m-1)/2, costheta);
  (*fComponent)(3,0)=-Fr*spindown*CG(3,0)*
		      nucleus.Spher_Harm(nucleus.getLbar_array()[shellindex], (m+1)/2, costheta);
}


TMFSpinor::TMFSpinor(double Gr, double Fr, int m, const double *Spher_harm, double phi, Matrix<4,1> CG)
: fComponent(new Matrix<4,1>){
  complex<double> spinup=exp(nucqcd::i_unit*0.5*phi*double(m-1));
  complex<double> spindown=exp(nucqcd::i_unit*0.5*phi*double(m+1));
  (*fComponent)(0,0)=nucqcd::i_unit*Gr*spinup*Spher_harm[0]*CG(0,0);
  (*fComponent)(1,0)=nucqcd::i_unit*Gr*spindown*Spher_harm[1]*CG(1,0);
  (*fComponent)(0,0)=-Fr*spinup*Spher_harm[2]*CG(2,0);
  (*fComponent)(1,0)=-Fr*spindown*Spher_harm[3]*CG(3,0);
}


//_____________________________________________________________________


TMFSpinor::TMFSpinor(const TMFSpinor& rhs)
  : fComponent(new Matrix<4,1>(*rhs.fComponent))
{
  // Copy constructor

  assert(fComponent); // check allocation
}

//_____________________________________________________________________
TMFSpinor::~TMFSpinor()
{
  // Destructor
  delete fComponent;
}

//_____________________________________________________________________
TMFSpinor& TMFSpinor::operator=(const TMFSpinor& rhs)
{
  // Assignment
  if( this!=&rhs) { // avoid self-assignment
    *fComponent=*rhs.fComponent;
  }

  return *this;
}

//_____________________________________________________________________
TMFSpinor* TMFSpinor::Clone(const char*) const
{
  // Virtual copy constructor
  return new TMFSpinor(*this);
}

//______________________________________________________________________________


Matrix<1,4> TMFSpinor::H()
{
  return (((Matrix<4,1>)*this).H());
}

//______________________________________________________________________________
ostream& operator<<(ostream& output, const TMFSpinor& rhs) 
{
  return output << static_cast<const Matrix<4,1>&>(rhs);
}




//______________________________________________________________________________
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TMFSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
complex<double> operator*(const Matrix<1,4>& lhs, const TMFSpinor& rhs)
{ 
  return lhs * (Matrix<4,1>)rhs; 
}

//______________________________________________________________________________
Matrix<4,1> operator*(const GammaStructure& lhs, const TMFSpinor& rhs)
{
  return (Matrix<4,4>)lhs * (Matrix<4,1>)rhs;
}

