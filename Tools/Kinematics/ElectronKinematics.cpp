/*
 * ElectronKinematics.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: March,24 2010
 *
 */

#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Tools/Kinematics/Kinematics2to3.hpp>
// #include <Structures.h>
#include <Tools/Core/constants.hpp>
#include <iostream>
#include <cassert>
#include <cmath>

namespace nucqcd {

using std::cerr;

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// ElectronKinematics                                                   //
//                                                                       //
// In electroproduction reactions we separate the electron scattering
// vertex from the virtual photon - deuteron system.
// This class deals with the kinematics of the former. In conjunction
// with a Kinematics2to3 object the kinematics are completely fixed by
// * the incoming electron's beam energy
// * or the outgoing electron's scattering angle
// * or the virtual photon's transverse linear polarization.
//
// The virtual photon's transverse linear polarization is defined as
//BEGIN_LATEX #varepsilon = #left( 1+2#frac{|p_{#gamma}|^{2}}{Q^{2}}tan^{2}#frac{#theta_{e}}{2} #right)^{-1}  END_LATEX
// and can be obtained or set with GetEpsilon and SetEpsilon.
//                                                                       //
///////////////////////////////////////////////////////////////////////////


//_____________________________________________________________________
ElectronKinematics::ElectronKinematics()
: fInput(kEpsilon), fBeamEnergy(0.), fEpsilon(0.), fScatterAngle(0.)
{
  // Protected Default Constructor
  //
  // Instances of ElectronKinematics should be created with
  // the named constructors.
}

//_____________________________________________________________________

//_____________________________________________________________________
ElectronKinematics::ElectronKinematics(const ElectronKinematics& rhs)
  : fInput(rhs.fInput), fBeamEnergy(rhs.fBeamEnergy), 
    fEpsilon(rhs.fEpsilon), fScatterAngle(rhs.fScatterAngle)
{
  // Copy Constructor
}

//_____________________________________________________________________
ElectronKinematics& ElectronKinematics::operator=(const ElectronKinematics& rhs)
{
  // Assignment
  if( this!=&rhs ) {
    fInput = rhs.fInput;
    fBeamEnergy = rhs.fBeamEnergy;
    fEpsilon = rhs.fEpsilon;
    fScatterAngle = rhs.fScatterAngle;
  }
  return *this;
}

//_____________________________________________________________________
ElectronKinematics* ElectronKinematics::CreateWithBeamEnergy(double beamEnergy)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  //
  // The user is responsible for deleting the returned object.
  ElectronKinematics *kinematics = new ElectronKinematics();
  kinematics->SetBeamEnergy(beamEnergy);
  return kinematics;
}

//_____________________________________________________________________
ElectronKinematics* ElectronKinematics::CreateWithEpsilon(double epsilon)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the 
  // virtual photon's transverse linear polarization.
  //
  // The user is responsible for deleting the returned object.
  ElectronKinematics *kinematics = new ElectronKinematics();
  kinematics->SetEpsilon(epsilon);
  return kinematics;
}

//_____________________________________________________________________
ElectronKinematics* ElectronKinematics::CreateWithCosScatterAngle(double scatterAngle)
{
  // Named constructor
  //
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  //
  // The user is responsible for deleting the returned object.
  ElectronKinematics *kinematics = new ElectronKinematics();
  kinematics->SetCosScatterAngle(scatterAngle);
  return kinematics;
}

//_____________________________________________________________________
void ElectronKinematics::SetBeamEnergy(double beamEnergy)
{
  // Fix the electron kinematics by specifying the energy
  // of the incoming electron beam in [MeV].
  fInput = kBeamEnergy;
  fBeamEnergy = beamEnergy;
}

//_____________________________________________________________________
void ElectronKinematics::SetEpsilon(double epsilon)
{
  // Fix the electron kinematics by specifying the
  // virtual photon's transverse linear polarization.
  fInput = kEpsilon;
  fEpsilon = epsilon;
}

//_____________________________________________________________________
void ElectronKinematics::SetCosScatterAngle(double scatterAngle)
{
  // Fix the electron kinematics by specifying the cosine of the
  // electron's scattering angle in ]-1,1].
  fInput = kScatterAngle;
  fScatterAngle = scatterAngle;
}

//_____________________________________________________________________
double ElectronKinematics::GetBeamEnergy(const Kinematics2to2& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetBeamEnergy(Kinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}

double ElectronKinematics::GetBeamEnergy(const double Q2, const double nu) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(Q2,nu) ) {
    cerr << "ERROR in ElectronKinematics::GetBeamEnergy(const double Q2, const double nu) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}

double ElectronKinematics::GetBeamEnergy(const Kinematics2to3& tk) const
{
  // Returns the incoming electron's energy in [MeV]
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetBeamEnergy(Kinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fBeamEnergy;
}


//_____________________________________________________________________
double ElectronKinematics::GetEpsilon(const Kinematics2to2& tk) const
{
  // Returns the virtual photon's transverse linear polarization
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetEpsilon(Kinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fEpsilon;
}

double ElectronKinematics::GetEpsilon(const double Q2, const double nu) const
{
  // Returns the virtual photon's transverse linear polarization
  // given the kinematics in the g*D system.
  if( !SolveKinematics(Q2,nu) ) {
    cerr << "ERROR in ElectronKinematics::GetEpsilon(const double Q2, const double nu) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fEpsilon;
}

double ElectronKinematics::GetEpsilon(const Kinematics2to3& tk) const
{
  // Returns the virtual photon's transverse linear polarization
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetEpsilon(Kinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fEpsilon;
}

//_____________________________________________________________________
double ElectronKinematics::GetCosScatterAngle(const Kinematics2to2& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(Kinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

//_____________________________________________________________________
double ElectronKinematics::GetCosScatterAngle(const double Q2, const double nu) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(Q2,nu) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(const double Q2, const double nu) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

double ElectronKinematics::GetCosScatterAngle(const Kinematics2to3& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(Kinematics2to3&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return fScatterAngle; 
}

//_____________________________________________________________________
double ElectronKinematics::GetTan2HalfAngle(const Kinematics2to2& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(Kinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return ftan2HalfAngle; 
}


//_____________________________________________________________________
double ElectronKinematics::GetTan2HalfAngle(const double Q2, const double nu) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(Q2,nu) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(const double Q2, const double nu) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return ftan2HalfAngle; 
}

//_____________________________________________________________________
double ElectronKinematics::GetTan2HalfAngle(const Kinematics2to3& tk) const
{
  // Returns the cosine of the electron's scattering angle
  // given the kinematics in the g*D system.
  if( !SolveKinematics(tk) ) {
    cerr << "ERROR in ElectronKinematics::GetCosScatterAngle(Kinematics2to2&) "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  return ftan2HalfAngle; 
}

//_____________________________________________________________________
bool ElectronKinematics::SolveKinematics(const Kinematics2to2& tk) const
{
  return SolveKinematics(tk.GetQsquared(),tk.GetWlab());

}

//_____________________________________________________________________
bool ElectronKinematics::SolveKinematics(const Kinematics2to3& tk) const
{

  return SolveKinematics(tk.GetQsquared(),tk.GetWlab());
}

//_____________________________________________________________________
bool ElectronKinematics::SolveKinematics(const double Q2, const double nu) const
{
  // Solve the electron scattering system kinematics and check whether 
  // the electron-system kinematics are compatible with those of the 
  // virtual photon - deuteron system.
  //
  // Returns 'false' when the kinematics are unphysical.

  double tmp;

  switch( fInput ) {
  case kBeamEnergy:
    // Solve the kinematics starting with the beam energy.
    if( fBeamEnergy <= nu ) return false;
    
    // Find the scattering angle
    fScatterAngle = 1. - Q2/2./fBeamEnergy/(fBeamEnergy-nu);
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < nucqcd::strange_underflow ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Finally calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( Q2<nucqcd::strange_underflow ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -nucqcd::strange_underflow ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > nucqcd::strange_underflow ) return false;
      else fEpsilon = 1.;
    }
    break;

  case kEpsilon:
    // Solve the kinematics starting with epsilon.
    if( fEpsilon<0. || fEpsilon>1. ) return false;

    // Find the scattering angle
    if( Q2<nucqcd::strange_underflow ) fScatterAngle = 1.;
    else {
      tmp = (1./fEpsilon-1.)/2./(1.+nu*nu/Q2);
      fScatterAngle = (1.-tmp)/(1.+tmp);
    }
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < nucqcd::strange_underflow ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Calculate the beam energy
    if( fScatterAngle > 1. ) return false;
    else if( fScatterAngle == 1. ) fBeamEnergy = nu + nucqcd::strange_underflow;
    else
      fBeamEnergy = nu/2. + sqrt(nu*nu/4.
					   +Q2/2./(1.-fScatterAngle));
    if( fBeamEnergy <= nu ) return false;
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;

  case kScatterAngle:
    // Solve the kinematics starting with the scattering angle
    if( fabs(fScatterAngle) > 1. ) {
      if( fabs(fScatterAngle)-1.0 < nucqcd::strange_underflow ) 
     	fScatterAngle = ( fScatterAngle>0.0 ? 1.0 : -1.0 );
      else
	return false;
    }

    // Find the beam energy
    if( fScatterAngle == 1. ) fBeamEnergy = nu + nucqcd::strange_underflow;
    else
      fBeamEnergy = nu/2. + sqrt( nu*nu/4.
					    +Q2/2./(1.-fScatterAngle) );
    if( fBeamEnergy <= nu ) return false;

    // Calculate epsilon
    if( fScatterAngle<= -1. ) fEpsilon = 0.;
    else if( Q2<nucqcd::strange_underflow ) fEpsilon = 1.;
    else {
      ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
      fEpsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *ftan2HalfAngle);
    }
    if( fEpsilon < 0. ) {
      if( fEpsilon < -nucqcd::strange_underflow ) return false;
      else fEpsilon = 0.;
    } else if( fEpsilon > 1. ) {
      if( (fEpsilon-1.) > nucqcd::strange_underflow ) return false;
      else fEpsilon = 1.;
    }
    ftan2HalfAngle = (1.-fScatterAngle)/(1.+fScatterAngle);
    break;
    
  default:
    cerr << "ERROR in ElectronKinematics::SolveKinematics(Kinematics2to2&): "
	 << "Invalid input variable.\n";
    assert(1==0);
  }
  
  return true;
}

void ElectronKinematics::GetLeptonVectors(const Kinematics2to2& kin, FourVector<double> &k_in, 
					 FourVector<double> &k_out){
  
  if( !SolveKinematics(kin) ) {
    cerr << "ERROR in LeptonKinematics::GetLeptonVectors "
	 << "impossible kinematics.\n";
    assert(1==0);
  }
  
  double Eout=fBeamEnergy-kin.GetWlab();
  double pout=Eout;
  double kx=fBeamEnergy*pout*sqrt(1.-fScatterAngle*fScatterAngle)/kin.GetKlab();
  double kz=kin.GetWlab()*fBeamEnergy/kin.GetKlab() + (kin.GetQsquared())/(2.*kin.GetKlab());
  double kzprime=kin.GetWlab()*Eout/kin.GetKlab() + (-kin.GetQsquared())/(2.*kin.GetKlab());
  k_in[0]=fBeamEnergy; k_in[1]=k_out[1]=kx; k_in[2]=k_out[2]=0.;k_in[3]=kz;
  k_out[0]=Eout; k_out[3]=kzprime;
  
}

}  // namespace nucqcd


