/*
 * TKinematics2to3_Frames.cpp
 * TKinematics2to3: boosts to and from the CM frames of the subsystems, rotations to the axis conventions, bounds
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
nucqcd::LorentzVector TKinematics2to3::BoostToKYcm(const nucqcd::LorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  return BoostToKYcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostToKNcm(const nucqcd::LorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  return BoostToKNcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostToYNcm(const nucqcd::LorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  return BoostToYNcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostToGNcm(const nucqcd::LorentzVector& vector,
					    EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostToGNcm(ac) * vector;
}
//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostFromKYcm(const nucqcd::LorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromKYcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostFromKNcm(const nucqcd::LorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromKNcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostFromYNcm(const nucqcd::LorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostFromYNcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::BoostFromGNcm(const nucqcd::LorentzVector& vector,
					      EAxisConvention ac) const
{
  // Boost the 4vector 'vector' defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostFromGNcm(ac) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostToKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  nucqcd::LorentzRotation boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostToKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  nucqcd::LorentzRotation boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostToYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  nucqcd::LorentzRotation boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostToGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  nucqcd::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostToKYcm(ac);
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostFromKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  nucqcd::LorentzRotation boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostFromKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  nucqcd::LorentzRotation boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostFromYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  nucqcd::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  nucqcd::LorentzRotation boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::BoostFromGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  nucqcd::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToKYcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  return BoostSpinorToKYcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToKNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  return BoostSpinorToKNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToYNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  return BoostSpinorToYNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostToGNcm(const TSpinor& spinor,
				     EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostSpinorToGNcm(ac) * spinor;
}
//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromKYcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromKYcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromKNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromKNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromYNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame .

  return BoostSpinorFromYNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::BoostFromGNcm(const TSpinor& spinor,
				       EAxisConvention ac) const
{
  // Boost the Dirac spinor 'spinor' defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  return BoostSpinorFromGNcm(ac) * spinor;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToKYcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*KYvector).BoostVector();
  TSpinor::LorentzRotation boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToKNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*KNvector).BoostVector();
  TSpinor::LorentzRotation boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToYNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*YNvector).BoostVector();
  TSpinor::LorentzRotation boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorToGNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TSpinor::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM 
    = nucqcd::LorentzRotation().RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X())) 
    * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostSpinorToKYcm(ac);
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromKYcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*KYvector).BoostVector();
  TSpinor::LorentzRotation boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromKNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*KNvector).BoostVector();
  TSpinor::LorentzRotation boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromYNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TSpinor::LorentzRotation rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector 
    = -1.0 * (nucqcd::LorentzRotation().RotateY(rotationAngle)*YNvector).BoostVector();
  TSpinor::LorentzRotation boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::BoostSpinorFromGNcm(EAxisConvention ac) const
{
  // Object to boost a Dirac spinor defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TSpinor::LorentzRotation rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM 
    =  nucqcd::LorentzRotation().RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()))
    * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostSpinorFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-hyperon subsystem.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzQuaternion boostCMKY = boostVector;

  return boostCMKY * rotateToKY;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the kaon-nucleon subsystem.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzQuaternion boostCMKN = boostVector;

  return boostCMKN * rotateToKN;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the hyperon-nucleon subsystem.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzQuaternion boostCMYN = boostVector;

  return boostCMYN * rotateToYN;
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionToGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the LAB frame with the 'ac'
  // axis convention to the CM frame of the photon and the struck nucleon
  // inside the deuteron.
  // The z-axis lies along the photon and the y-axis is defined as being 
  // perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzQuaternion rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return rotateAxes * BoostQuaternionToKYcm(ac);
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromKYcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-hyperon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KYvector = GetK4(ac) + GetY4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKY,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKY; 
  double rotationAngle=0.; 
  if( 1.-fabs(KYvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KYvector.X())
      *KYvector.Vect().Angle(fGvector->Vect());
  rotateToKY.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KY CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKY*KYvector).BoostVector();
  TLorentzQuaternion boostCMKY = boostVector;

  return  rotateToKY.Inverse() * boostCMKY.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromKNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // kaon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector KNvector = GetK4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pKN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToKN; 
  double rotationAngle=0.; 
  if( 1.-fabs(KNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,KNvector.X())
      *KNvector.Vect().Angle(fGvector->Vect());
  rotateToKN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to KN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToKN*KNvector).BoostVector();
  TLorentzQuaternion boostCMKN = boostVector;

  return rotateToKN.Inverse() * boostCMKN.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromYNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the 
  // hyperon-nucleon subsystem with the 'ac' axis convention to the LAB frame.

  nucqcd::LorentzVector YNvector = GetY4(ac) + GetN4(ac);

  // Define rotation matrix that changes reference frame to
  // z' ~ pYN,  y' ~ y and x' ~ y' x z'
  TLorentzQuaternion rotateToYN; 
  double rotationAngle=0.; 
  if( 1.-fabs(YNvector.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,YNvector.X())
      *YNvector.Vect().Angle(fGvector->Vect());
  rotateToYN.RotateY( rotationAngle );
  
  // Define Lorentzboost that boosts in reference frame defined above
  // to YN CM system.
  nucqcd::Vector3 boostVector = -1.0 * (rotateToYN*YNvector).BoostVector();
  TLorentzQuaternion boostCMYN = boostVector;

  return rotateToYN.Inverse() * boostCMYN.Inverse();
}

//_____________________________________________________________________
TLorentzQuaternion TKinematics2to3::BoostQuaternionFromGNcm(EAxisConvention ac) const
{
  // Object to boost a 4vector defined in the CM frame of the photon 
  // and the struck nucleon inside the deuteron to the LAB frame
  // with the 'ac' axis convention.
  // 
  // In the CM frame, the z-axis lies along the photon and the y-axis is 
  // defined as being perpendicular to the gN and KY plane.
  // This is the reference frame where the matrix elements get evaluated
  // in strangecalc. 

  // Find the photon vector in the KY CM
  nucqcd::LorentzVector GvectorCM = BoostToKYcm(GetG4(ac),ac);

  // Define a rotation that puts z along the photon momentum
  TLorentzQuaternion rotateAxes;
  rotateAxes.RotateY(-std::copysign(GvectorCM.Theta(),GvectorCM.X()));

  // Find the kaon vector in the KY CM with this new reference frame
  nucqcd::LorentzVector KvectorCM = rotateAxes * BoostToKYcm(GetK4(ac),ac);

  // Define an extra rotation to a reference frame where the azimuthal
  // angle of the kaon is zero and as a consequence the y-axis is 
  // perpendicular to the KY plane.
  rotateAxes.RotateZ(-KvectorCM.Phi());

  return BoostQuaternionFromKYcm(ac) * rotateAxes.Inverse();
}

//_____________________________________________________________________
void TKinematics2to3::Help()
{
  std::cout << "I can't help you! "
	       << "In fact, can anyone?\n";
}

//_____________________________________________________________________
// FourVector<double> ToFourVector(const nucqcd::LorentzVector& v)
// {
//   return FourVector<double>(v.E(),v.X(),v.Y(),v.Z());
// }
// 
// //_____________________________________________________________________
// nucqcd::LorentzVector ToLorentzVector(const FourVector<double>& v)
// {
//   return nucqcd::LorentzVector(v[1],v[2],v[3],v[0]);
// }

//_____________________________________________________________________
double TKinematics2to3::GetPhotonThreshold() const
{
  // Calculate the threshold photon lab energy for the current isospin
  // channel.
  return (( std::pow(fMk+fMy+fMn,2.) + fQsquared - fMd*fMd ) / 2. / fMd) + kUflow;
}

//_____________________________________________________________________
void TKinematics2to3::GetBoundsMomentumParticle1(double th1,
						 double& lower,
						 double& upper) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has an upper and lower bound.
  //
  // These bounds are calculated by imposing the constraint W_23>=m_2+m_3.
  // This leads to a quadratic equation for |p_1|, which we solve.

  double eg = fGvector->E();
  double stot = fMd*fMd - fQsquared + 2.*fMd*fGvector->E();
  
  double costh1 = std::cos(th1);
  double m1=0., m23=0.;
  if( fAC==kK ) {
    m1 = fMk;
    m23 = fMy + fMn;
  } // particle1 = kaon
  else if( fAC==kY ) {
    m1 = fMy;
    m23 = fMk + fMn;
  } // particle1 = hyperon
  else {
    m1 = fMn;
    m23 = fMk + fMy;
  } // particle1 = nucleon

  // Define the quadratic equation...
  std::vector<double> roots = nucqcd::quadratic_real_roots(// |p_1|^2 coefficient
		(fQsquared+eg*eg)*costh1*costh1 - (eg+fMd)*(eg+fMd),
		// |p_1| coefficient
		std::sqrt(fQsquared+eg*eg)*costh1*(stot+m1*m1-(m23)*(m23)),
		// constant coefficient
		(std::pow(stot+m1*m1-(m23)*(m23),2.)/4.
		 -m1*m1*(eg+fMd)*(eg+fMd)));

  // ...and solve
  
  if( roots[0]>roots[1] ) {
    upper = std::max(roots[0],0.);
    lower = std::max(roots[1],0.);
  } else {
    upper = std::max(roots[1],0.);
    lower = std::max(roots[0],0.);
  }
}

//_____________________________________________________________________
double TKinematics2to3::GetLowerBoundMomentumParticle1(double th1) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has a lower bound.

  double lower,upper;
  GetBoundsMomentumParticle1(th1,lower,upper);
  
  return lower + kUflow;
}

//_____________________________________________________________________
double TKinematics2to3::GetUpperBoundMomentumParticle1(double th1) const
{
  // For the current invariant mass 'w' and a given scattering angle 
  // 'th1' of the particle specified by the axis convention, the lab 
  // momentum |p_1| has an upper bound.

  double lower,upper;
  GetBoundsMomentumParticle1(th1,lower,upper);
  
  return upper - kUflow;
}

//_____________________________________________________________________
double TKinematics2to3::GetMaximumScatterAngleParticle1() const
{
  // For the current invariant mass 'w' and lab momentum 'p_1' of the
  // particle specified by the axis convention, the lab scattering angle 
  // 'theta_1' has a maximum.
  //
  // This maximum is found by imposing the constraint W_23 >= m_2+m_3.
  
  // Find the 4-vector of particle 1 and sum masses of particles 2 and 3.
  nucqcd::LorentzVector *p1;
  double m23;
  if( fAC==kK ) {
    p1 = fKvector;
    m23 = fMy + fMn;
  } // particle1 = kaon
  else if( fAC==kY ) {
    p1 = fYvector;
    m23 = fMk + fMn;
  } // particle1 = hyperon
  else {
    p1 = fNvector;
    m23 = fMk + fMy;
  } // particle1 = nucleon

  double stot = fMd*fMd - fQsquared + 2.*fMd*fGvector->E();

  double costh1Min = 
    ( m23*m23 + 2.*p1->E() * (fGvector->E() + fMd)
      - stot - p1->M2() )
    / 2. / p1->P() / fGvector->P();

  if( costh1Min > 1.-kUflow ) return 0.;

  return std::acos(costh1Min) - kUflow;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::RotateToKac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.
  nucqcd::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToKac(fAC);

  else if( axisConvention == kY ) 
    rotation.RotateZ(-1.0*GetPhik(kY));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhik(kN));

  // else( axisConvention == kK )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::RotateToYac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.
  nucqcd::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToYac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhiy(kK));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhiy(kN));

  // else( axisConvention == kY )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
nucqcd::LorentzRotation TKinematics2to3::RotateToNac(EAxisConvention axisConvention) const
{
  // Object to rotate a 4vector defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.
  nucqcd::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateToNac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhin(kK));

  else if( axisConvention == kY )
    rotation.RotateZ(-1.0*GetPhin(kY));

  // else( axisConvention == kN )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToKac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToKac(fAC);

  else if( axisConvention == kY ) 
    rotation.RotateZ(-1.0*GetPhik(kY));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhik(kN));

  // else( axisConvention == kK )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToYac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToYac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhiy(kK));

  else if( axisConvention == kN )
    rotation.RotateZ(-1.0*GetPhiy(kN));

  // else( axisConvention == kY )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
TSpinor::LorentzRotation TKinematics2to3::RotateSpinorToNac(EAxisConvention axisConvention) const
{
  // Object to rotate a Dirac spinor defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.
  TSpinor::LorentzRotation rotation;
  
  if( axisConvention == kS )
    return RotateSpinorToNac(fAC);

  else if( axisConvention == kK ) 
    rotation.RotateZ(-1.0*GetPhin(kK));

  else if( axisConvention == kY )
    rotation.RotateZ(-1.0*GetPhin(kY));

  // else( axisConvention == kN )
  // unit rotation

  return rotation;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::RotateToKac(const nucqcd::LorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.  
  return RotateToKac(axisConvention) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::RotateToYac(const nucqcd::LorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.  
  return RotateToYac(axisConvention) * vector;
}

//_____________________________________________________________________
nucqcd::LorentzVector TKinematics2to3::RotateToNac(const nucqcd::LorentzVector& vector,
					    EAxisConvention axisConvention) const
{
  // Rotate the 4vector 'vector' defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.  
  return RotateToNac(axisConvention) * vector;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToKac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kK' axis convention.  
  return RotateSpinorToKac(axisConvention) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToYac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kY' axis convention.  
  return RotateSpinorToYac(axisConvention) * spinor;
}

//_____________________________________________________________________
TSpinor TKinematics2to3::RotateToNac(const TSpinor& spinor,
				     EAxisConvention axisConvention) const
{
  // Rotate the Dirac spinor 'spinor' defined in the LAB frame from the 'ac'
  // axis convention to the 'kN' axis convention.  
  return RotateSpinorToNac(axisConvention) * spinor;
}


}  // namespace nucqcd
