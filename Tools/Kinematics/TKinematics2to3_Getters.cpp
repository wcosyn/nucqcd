/*
 * TKinematics2to3_Getters.cpp
 * TKinematics2to3: getters of vectors, angles and CM quantities, operator==, GetVarByGetter
 * Part of the implementation of TKinematics2to3 (split from TKinematics2to3.cpp, which holds the class
 * description, constructors and assignment); author of the original: Pieter Vancraeyveld.
 */

#include <Tools/Kinematics/TKinematics2to3.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/TLorentzQuaternion.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/MathUtil.hpp>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <Tools/Core/StringUtil.hpp>
#include <algorithm>
#include <optional>
#include <vector>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/Vector3.hpp>

namespace nucqcd {

// #include <numtoa.h>
using std::cout; using std::cerr;
using std::vector;
using std::strcmp;



//_____________________________________________________________________
int TKinematics2to3::GetNumberOfSteps() const
{
  return fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2]
    * fNumberOfSteps[3] * fNumberOfSteps[4] * fNumberOfSteps[5];
}

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfPhysicalSteps() const
{
  if( fNrOfPhysicals==-1 ) {
    TKinematics2to3 tk(*this); // temporary copy
    fNrOfPhysicals = 0;        // final result
    
    for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
      tk.GoTo(step);
      if(tk.IsPhysical()) ++fNrOfPhysicals;
    }
  }

  return fNrOfPhysicals;
}    

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfSteps(int varNumber) const
{
  
  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetNumberOfSteps(int):\n"
	      << "Index out-of-bounds!\n";
    return 0;
  }

  return fNumberOfSteps[varNumber-1];
}

//_____________________________________________________________________
int TKinematics2to3::GetNumberOfVariables() const
{
  // Determine how many independent variables are not fixed

  int nrOfVar = 0;
  for(int var=1; var<=6; ++var)
    if( !IsFixed(var) ) ++nrOfVar;

  return nrOfVar;
}

//_____________________________________________________________________
double TKinematics2to3::GetStepSize(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetStepSize(int):\n"
	      << "Index out-of-bounds!\n";
    return 0;
  }

  return fStepSize[varNumber-1];
}

//_____________________________________________________________________
int TKinematics2to3::GetStep(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetStep(int):\n"
	      << "Index out-of-bounds!\n";
    return -1;
  }

  // Determine step per variable
  int stepVar[6];
  int remainder=fStep;
  int denum=1;
  
  for(int i=6; i>0; --i)
    {
      denum = 1;
      for(int j=0; j<i-1; ++j)
	denum *= fNumberOfSteps[j];
      
      stepVar[i-1] = remainder / denum;
      remainder = remainder % denum;
    }
    
  return stepVar[varNumber-1];  
}

//_____________________________________________________________________
bool TKinematics2to3::IsFixed() const
{
  return !fIsVar[0] && !fIsVar[1] && !fIsVar[2]
    && !fIsVar[3] && !fIsVar[4] && !fIsVar[5];
}

//_____________________________________________________________________
bool TKinematics2to3::IsFixed(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::IsFixed(int):\n"
	      << "Index out-of-bounds!\n";
    return true;
  }

  return !fIsVar[varNumber-1];
}

//______________________________________________________________________________

// }

//_____________________________________________________________________
nucqcd::Vector3 TKinematics2to3::GetK3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the kaon in the LAB.

  if( axisConvention == kS )
    return GetK3(fAC);

  else if( axisConvention == kY ) {
    nucqcd::Vector3 rotatedKvector = fKvector->Vect();
    rotatedKvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedKvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::Vector3 rotatedKvector = fKvector->Vect();
    rotatedKvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedKvector;
  }

  return fKvector->Vect();
}

//_____________________________________________________________________
nucqcd::Vector3 TKinematics2to3::GetY3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the hyperon in the LAB.

  if( axisConvention == kS )
    return GetY3(fAC);

  else if( axisConvention == kY ) {
    nucqcd::Vector3 rotatedYvector = fYvector->Vect();
    rotatedYvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedYvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::Vector3 rotatedYvector = fYvector->Vect();
    rotatedYvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedYvector;
  }

  return fYvector->Vect();
}

//_____________________________________________________________________
nucqcd::Vector3 TKinematics2to3::GetN3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the nucleon in the LAB.

  if( axisConvention == kS )
    return GetN3(fAC);

  else if( axisConvention == kY ) {
    nucqcd::Vector3 rotatedNvector = fNvector->Vect();
    rotatedNvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedNvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::Vector3 rotatedNvector = fNvector->Vect();
    rotatedNvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedNvector;
  }

  return fNvector->Vect();
}

//_____________________________________________________________________
nucqcd::Vector3 TKinematics2to3::GetG3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the photon in the LAB.
  return fGvector->Vect();
}

//_____________________________________________________________________
nucqcd::Vector3 TKinematics2to3::GetD3(EAxisConvention axisConvention) const
{
  // Return the momentum 3-vector of the deuteron in the LAB.
  return fDvector->Vect();
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::GetK4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the kaon in the LAB.

  if( axisConvention == kS )
    return GetK4(fAC);

  else if( axisConvention == kY ) {
    nucqcd::LorentzVector rotatedKvector = *fKvector;
    rotatedKvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedKvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::LorentzVector rotatedKvector = *fKvector;
    rotatedKvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedKvector;
  }

  return *fKvector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::GetY4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the hyperon in the LAB.

  if( axisConvention == kS )
    return GetY4(fAC);

  else if( axisConvention == kY ) {
    nucqcd::LorentzVector rotatedYvector = *fYvector;
    rotatedYvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedYvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::LorentzVector rotatedYvector = *fYvector;
    rotatedYvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedYvector;
  }

  return *fYvector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::GetN4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the nucleon in the LAB.

  if( axisConvention == kS )
    return GetN4(fAC);

  else if( axisConvention == kY ) {
    nucqcd::LorentzVector rotatedNvector = *fNvector;
    rotatedNvector.RotateZ(-1.0*GetPhiy(kN));
    return rotatedNvector;
  }
  
  else if( axisConvention == kK ) {
    nucqcd::LorentzVector rotatedNvector = *fNvector;
    rotatedNvector.RotateZ(-1.0*GetPhik(kN));
    return rotatedNvector;
  }

  return *fNvector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::GetG4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the photon in the LAB.
  return *fGvector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::GetD4(EAxisConvention axisConvention) const
{
  // Return the energy-momentum 4-vector of the deuteron in the LAB
  return *fDvector;
}

//_____________________________________________________________________
double TKinematics2to3::GetStot() const
{
  return fMd*fMd + 2*GetWlab()*fMd - GetQsquared();
}

//_____________________________________________________________________
double TKinematics2to3::GetThk() const 
{ 
  return fKvector->Angle( fGvector->Vect() ); 
}
  
//_____________________________________________________________________
double TKinematics2to3::GetThy() const 
{ 
  return fYvector->Angle( fGvector->Vect() ); 
}
 
//_____________________________________________________________________
double TKinematics2to3::GetThn() const 
{ 
  return fNvector->Angle( fGvector->Vect() ); 
}

//_____________________________________________________________________
double TKinematics2to3::GetThk_deg() const 
{ 
  return GetThk() * nucqcd::rad_to_deg; 
}
  
//_____________________________________________________________________
double TKinematics2to3::GetThy_deg() const 
{ 
  return GetThy() * nucqcd::rad_to_deg; 
}
 
//_____________________________________________________________________
double TKinematics2to3::GetThn_deg() const 
{ 
  return GetThn() * nucqcd::rad_to_deg; 
}

//_____________________________________________________________________
double TKinematics2to3::GetPhik(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    return GetPhik(fAC); // the return was missing until 2026-11

  else if( axisConvention == kN )
    return fKvector->Vect().XYvector().Phi();

  else if( axisConvention == kK )
    return 0.0;

  // axisConvention == kY
  return nucqcd::Vector2::Phi_0_2pi( GetPhik(kN) - GetPhiy(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhiy(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    return GetPhiy(fAC);

  else if( axisConvention == kN )
    return fYvector->Vect().XYvector().Phi();

  else if( axisConvention == kY )
    return 0.0;

  // axisConvention == kK
  return nucqcd::Vector2::Phi_0_2pi( GetPhiy(kN) - GetPhik(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhin(EAxisConvention axisConvention) const
{
  if( axisConvention == kS )
    return GetPhin(fAC);

  else if( axisConvention == kN )
    return 0.0;

  else if( axisConvention == kK )
    return nucqcd::Vector2::Phi_0_2pi( -1.0*GetPhik(kN) );

  // axisConvention == kY
  return nucqcd::Vector2::Phi_0_2pi( -1.0*GetPhiy(kN) );
}

//_____________________________________________________________________
double TKinematics2to3::GetPhik_deg(EAxisConvention axisConvention) const
{
  return GetPhik(axisConvention) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetPhiy_deg(EAxisConvention axisConvention) const
{
  return GetPhiy(axisConvention) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetPhin_deg(EAxisConvention axisConvention) const
{
  return GetPhin(axisConvention) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetThky() const
{
  return fGvector->Angle( fKvector->Vect() + fYvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThkn() const
{
  return fGvector->Angle( fKvector->Vect() + fNvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThyn() const
{
  return fGvector->Angle( fYvector->Vect() + fNvector->Vect() );
}

//_____________________________________________________________________
double TKinematics2to3::GetThky_deg() const
{
  return GetThky() * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetThkn_deg() const
{
  return GetThkn() * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetThyn_deg() const
{
  return GetThyn() * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetTkg() const
{
  return ( *fKvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTyg() const
{
  return ( *fYvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTng() const
{
  return ( *fNvector - *fGvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTkd() const
{
  return ( *fKvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTyd() const
{
  return ( *fYvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetTnd() const
{
  return ( *fNvector - *fDvector ).Mag2();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_Ek() const
{
  return std::sqrt( fMk*fMk + fKYsystem->GetPk()*fKYsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_Ey() const
{
  return std::sqrt( fMy*fMy + fKYsystem->GetPk()*fKYsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_Ek() const
{
  return std::sqrt( fMk*fMk + fKNsystem->GetPk()*fKNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_En() const
{
  return std::sqrt( fMn*fMn + fKNsystem->GetPk()*fKNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_Ey() const
{
  return std::sqrt( fMy*fMy + fYNsystem->GetPk()*fYNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_En() const
{
  return std::sqrt( fMn*fMn + fYNsystem->GetPk()*fYNsystem->GetPk() );
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_costhkcm() const
{
  return fKYsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_thkcm() const
{
  return std::acos(fKYsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_costhycm() const
{
  return -1.0 * fKYsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_thycm() const
{
  return std::acos(-1.0 * fKYsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_costhkcm() const
{
  return fKNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_thkcm() const
{
  return std::acos(fKNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_costhncm() const
{
  return -1.0 * fKNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_thncm() const
{
  return std::acos(-1.0 * fKNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_costhycm() const
{
  return -1.0 * fYNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_thycm() const
{
  return std::acos(-1.0 * fYNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_costhncm() const
{
  return fYNsystem->GetCosthkcm();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_thncm() const
{
  return std::acos(fYNsystem->GetCosthkcm());
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phikcm(EAxisConvention ac) const
{
  // Transform Kvector to the KY cm 
  nucqcd::LorentzVector vectorKcm = BoostToKYcm(GetK4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorKcm.X())<nucqcd::strange_underflow && fabs(vectorKcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorKcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phiycm(EAxisConvention ac) const
{
  // Transform Yvector to the KY cm 
  nucqcd::LorentzVector vectorYcm = BoostToKYcm(GetY4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorYcm.X())<nucqcd::strange_underflow && fabs(vectorYcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorYcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phikcm(EAxisConvention ac) const
{
  // Transform Kvector to the KN cm 
  nucqcd::LorentzVector vectorKcm = BoostToKNcm(GetK4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorKcm.X())<nucqcd::strange_underflow && fabs(vectorKcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorKcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phincm(EAxisConvention ac) const
{
  // Transform Nvector to the KN cm 
  nucqcd::LorentzVector vectorNcm = BoostToKNcm(GetN4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorNcm.X())<nucqcd::strange_underflow && fabs(vectorNcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorNcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phiycm(EAxisConvention ac) const
{
  // Transform Yvector to the YN cm 
  nucqcd::LorentzVector vectorYcm = BoostToYNcm(GetY4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorYcm.X())<nucqcd::strange_underflow && fabs(vectorYcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorYcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phincm(EAxisConvention ac) const
{
  // Transform Nvector to the YN cm 
  nucqcd::LorentzVector vectorNcm = BoostToYNcm(GetN4(ac),ac);

  // and determine azimutal angle
  if(fabs(vectorNcm.X())<nucqcd::strange_underflow && fabs(vectorNcm.Y())<nucqcd::strange_underflow)
    return 0.;
  
  return (vectorNcm).Vect().XYvector().Phi();
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phikcm_deg(EAxisConvention ac) const
{
  return GetKYsystem_phikcm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetKYsystem_phiycm_deg(EAxisConvention ac) const
{
  return GetKYsystem_phiycm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phikcm_deg(EAxisConvention ac) const
{
  return GetKNsystem_phikcm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetKNsystem_phincm_deg(EAxisConvention ac) const
{
  return GetKNsystem_phincm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phiycm_deg(EAxisConvention ac) const
{
  return GetYNsystem_phiycm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
double TKinematics2to3::GetYNsystem_phincm_deg(EAxisConvention ac) const
{
  return GetYNsystem_phincm(ac) * nucqcd::rad_to_deg;
}

//_____________________________________________________________________
bool TKinematics2to3::operator==(const TKinematics2to3& rhs) const
{
  // Evaluate whether both objects currently refer to the same point in the
  // kinematical plane.

  // when the objects are the same, the answer is trivial
  if( this==&rhs ) return true;

  double areEqual = true;

  if( GetIsospin()!=rhs.GetIsospin() ) areEqual=false;

  if( fabs(GetStot())>1. ) {
    if( fabs((GetStot()-rhs.GetStot()))>1.e-6*fabs(GetStot()) )
      areEqual = false;
  } else {
    if( fabs(GetStot()-rhs.GetStot())>1.e-6 ) 
      areEqual = false;
  }

  if( fabs(GetSkn())>1. ) {
    if( fabs((GetSkn()-rhs.GetSkn()))>1.e-6*fabs(GetSkn()) )
      areEqual = false;
  } else {
    if( fabs((GetSkn()-rhs.GetSkn()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetSyn())>1. ) {
    if( fabs((GetSyn()-rhs.GetSyn()))>1.e-6*fabs(GetSyn()) )
      areEqual = false;
  } else {
    if( fabs((GetSyn()-rhs.GetSyn()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetTkg())>1. ) {
    if( fabs((GetTkg()-rhs.GetTkg()))>1.e-6*fabs(GetTkg()) )
      areEqual = false;
  } else {
    if( fabs((GetTkg()-rhs.GetTkg()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetTyd())>1. ) {
    if( fabs((GetTyd()-rhs.GetTyd()))>1.e-6*fabs(GetTyd()) )
      areEqual = false;
  } else {
    if( fabs((GetTyd()-rhs.GetTyd()))>1.e-6 )
      areEqual = false;
  }

  if( fabs(GetQsquared())>1. ) {
    if( fabs((GetQsquared()-rhs.GetQsquared()))>1.e-6*fabs(GetQsquared()) )
      areEqual = false;
  } else {
    if( fabs((GetQsquared()-rhs.GetQsquared()))>1.e-6 )
      areEqual = false;
  }

  return areEqual;  
}

//_____________________________________________________________________
double TKinematics2to3::GetVarByGetter(const std::string& variable) const
{
  // Return the value of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // This function is used by TKinematics2to3::UpdateKinematics() to check the solution.

  double(TKinematics2to3::*ptrToVarGetter)() const; // pointer to TKinematics getters
  
  const std::string varName = nucqcd::to_lower(variable);
  
  // Assign the correct Getter
  if( (varName=="pk") ) {
    ptrToVarGetter = &TKinematics2to3::GetPk;
  }
  else if( (varName=="py") ) {
    ptrToVarGetter = &TKinematics2to3::GetPy;
  }
  else if( (varName=="pn") ) {
    ptrToVarGetter = &TKinematics2to3::GetPn;
  }
  else if( (varName=="klab") ) {
    ptrToVarGetter = &TKinematics2to3::GetKlab;
  }
  else if( (varName=="ek") ) {
    ptrToVarGetter = &TKinematics2to3::GetEk;
  }
  else if( (varName=="ey") ) {
    ptrToVarGetter = &TKinematics2to3::GetEy;
  }
  else if( (varName=="en") ) {
    ptrToVarGetter = &TKinematics2to3::GetEn;
  }
  else if( (varName=="wlab") ) {
    ptrToVarGetter = &TKinematics2to3::GetWlab;
  }
  else if( (varName=="thk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk_deg;
  }
  else if( (varName=="thy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy_deg;
  }
  else if( (varName=="thn_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn_deg;
  }
  else if( (varName=="thk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThk;
  }
  else if( (varName=="thy") ) {
    ptrToVarGetter = &TKinematics2to3::GetThy;
  }
  else if( (varName=="thn") ) {
    ptrToVarGetter = &TKinematics2to3::GetThn;
  }
  else if( (varName=="costhk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthk;
  }
  else if( (varName=="costhy") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthy;
  }
  else if( (varName=="costhn") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthn;
  }
  else if( (varName=="phik_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_deg_kS;
  }
  else if( (varName=="phiy_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_deg_kS;
  }
  else if( (varName=="phin_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_deg_kS;
  }
  else if( (varName=="phik") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhik_kS;
  }
  else if( (varName=="phiy") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhiy_kS;
  }
  else if( (varName=="phin") ) {
    ptrToVarGetter = &TKinematics2to3::GetPhin_kS;
  }
  else if( (varName=="thky_deg") ||
	   (varName=="thyk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThky_deg;
  }
  else if( (varName=="thkn_deg") ||
	   (varName=="thnk_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn_deg;
  }
  else if( (varName=="thyn_deg") ||
	   (varName=="thny_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn_deg;
  }
  else if( (varName=="thky") ||
	   (varName=="thyk")) {
    ptrToVarGetter = &TKinematics2to3::GetThky;
  }
  else if( (varName=="thkn") ||
	   (varName=="thnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetThkn;
  }
  else if( (varName=="thyn") ||
	   (varName=="thny") ) {
    ptrToVarGetter = &TKinematics2to3::GetThyn;
  }
  else if( (varName=="costhky") ||
	   (varName=="costhyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthky;
  }
  else if( (varName=="costhkn") ||
	   (varName=="costhnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthkn;
  }
  else if( (varName=="costhyn") ||
	   (varName=="costhny") ) {
    ptrToVarGetter = &TKinematics2to3::GetCosthyn;
  }
  else if( (varName=="qsquared") ) {
    ptrToVarGetter = &TKinematics2to3::GetQsquared;
  }
  else if( (varName=="wtot") ) {
    ptrToVarGetter = &TKinematics2to3::GetWtot;
  }
  else if( (varName=="wky") ||
	   (varName=="wyk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWky;
  }
  else if( (varName=="wkn") ||
	   (varName=="wnk") ) {
    ptrToVarGetter = &TKinematics2to3::GetWkn;
  }
  else if( (varName=="wyn") ||
	   (varName=="wny") ) {
    ptrToVarGetter = &TKinematics2to3::GetWyn;
  }
  else if( (varName=="stot") ) {
    ptrToVarGetter = &TKinematics2to3::GetStot;
  }
  else if( (varName=="sky") ||
	   (varName=="syk") ) {
    ptrToVarGetter = &TKinematics2to3::GetSky;
  }
  else if( (varName=="skn") ||
	   (varName=="snk")) {
    ptrToVarGetter = &TKinematics2to3::GetSkn;
  }
  else if( (varName=="syn") ||
	   (varName=="sny")) {
    ptrToVarGetter = &TKinematics2to3::GetSyn;
  }
  else if( (varName=="tkg") ||
	   (varName=="tgk")) {
    ptrToVarGetter = &TKinematics2to3::GetTkg;
  }
  else if( (varName=="tyg") ||
	   (varName=="tgy")) {
    ptrToVarGetter = &TKinematics2to3::GetTyg;
  }
  else if( (varName=="tng") ||
	   (varName=="tgn")) {
    ptrToVarGetter = &TKinematics2to3::GetTng;
  }
  else if( (varName=="tkd") ||
	   (varName=="tdk")) {
    ptrToVarGetter = &TKinematics2to3::GetTkd;
  }
  else if( (varName=="tyd") ||
	   (varName=="tdy")) {
    ptrToVarGetter = &TKinematics2to3::GetTyd;
  }
  else if( (varName=="tnd") ||
	   (varName=="tdn")) {
    ptrToVarGetter = &TKinematics2to3::GetTnd;
  }
  else if( (varName=="kypk") ||
	   (varName=="ykpk")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Pk;
  }
  else if( (varName=="kypy") ||
	   (varName=="ykpy")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Py;
  }
  else if( (varName=="kyek") ||
	   (varName=="ykek")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ek;
  }
  else if( (varName=="kyey") ||
	   (varName=="ykey")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_Ey;
  }
  else if( (varName=="kys") ||
	   (varName=="yks")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_S;
  }
  else if( (varName=="kyw") ||
	   (varName=="ykw")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_W;
  }
  else if( (varName=="kycosthkcm") ||
	   (varName=="ykcosthkcm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhkcm;
  }
  else if( (varName=="kycosthycm") ||
	   (varName=="ykcosthycm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_costhycm;
  }
  else if( (varName=="kythkcm") ||
	   (varName=="ykthkcm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thkcm;
  }
  else if( (varName=="kythycm") ||
	   (varName=="ykthycm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_thycm;
  }
  else if( (varName=="kyphikcm") ||
	   (varName=="ykphikcm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_kS;
  }
  else if( (varName=="kyphiycm") ||
	   (varName=="ykphiycm")) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_kS;
  }
  else if( (varName=="kyhpikcm_deg") ||
	   (varName=="ykphikcm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phikcm_deg_kS;
  }
  else if( (varName=="kyphiycm_deg") ||
	   (varName=="ykphiycm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetKYsystem_phiycm_deg_kS;
  }
  else if( (varName=="knpk") ||
	   (varName=="nkpk") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pk;
  }
  else if( (varName=="knpn") ||
	   (varName=="nkpn") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Pn;
  }
  else if( (varName=="knek") ||
	   (varName=="nkek") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_Ek;
  }
  else if( (varName=="knen") ||
	   (varName=="nken") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_En;
  }
  else if( (varName=="kns") ||
	   (varName=="nks") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_S;
  }
  else if( (varName=="knw") ||
	   (varName=="nkw") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_W;
  }
  else if( (varName=="kncosthkcm") ||
	   (varName=="nkcosthkcm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhkcm;
  }
  else if( (varName=="kncosthncm") ||
	   (varName=="nkcosthncm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_costhncm;
  }
  else if( (varName=="knthkcm") ||
	   (varName=="nkthkcm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thkcm;
  }
  else if( (varName=="knthncm") ||
	   (varName=="nkthncm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_thncm;
  }
  else if( (varName=="knphikcm") ||
	   (varName=="nkphikcm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_kS;
  }
  else if( (varName=="knphincm") ||
	   (varName=="nkphincm") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_kS;
  }
  else if( (varName=="knphikcm_deg") ||
	   (varName=="nkphikcm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phikcm_deg_kS;
  }
  else if( (varName=="knphincm_deg") ||
	   (varName=="nkphincm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetKNsystem_phincm_deg_kS;
  }
  else if( (varName=="ynpy") ||
	   (varName=="nypy") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Py;
  }
  else if( (varName=="ynpn") ||
	   (varName=="nypn") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Pn;
  }
  else if( (varName=="yney") ||
	   (varName=="nyey") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_Ey;
  }
  else if( (varName=="ynen") ||
	   (varName=="nyen") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_En;
  }
  else if( (varName=="yns") ||
	   (varName=="nys") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_S;
  }
  else if( (varName=="ynw") ||
	   (varName=="nyw") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_W;
  }
  else if( (varName=="yncosthycm") ||
	   (varName=="nycosthycm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhycm;
  }
  else if( (varName=="yncosthncm") ||
	   (varName=="nycosthncm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_costhncm;
  }
  else if( (varName=="ynthycm") ||
	   (varName=="nythycm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thycm;
  }
  else if( (varName=="ynthncm") ||
	   (varName=="nythncm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_thncm;
  }
  else if( (varName=="ynphiycm") ||
	   (varName=="nyphiycm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_kS;
  }
  else if( (varName=="ynphincm") ||
	   (varName=="nyphincm") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_kS;
  }
  else if( (varName=="ynphiycm_deg") ||
	   (varName=="nyphiycm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phiycm_deg_kS;
  }
  else if( (varName=="ynphincm_deg") ||
	   (varName=="nyphincm_deg") ) {
    ptrToVarGetter = &TKinematics2to3::GetYNsystem_phincm_deg_kS;
  }
  else {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in  TKinematics2to3::GetVarByGetter(std::string): "
	      << "Kinematic variable \"" << variable << "\" unknown!");
  }

  return (this->*ptrToVarGetter)();
}

}  // namespace nucqcd
