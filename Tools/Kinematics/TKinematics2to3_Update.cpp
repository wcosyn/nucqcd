/*
 * TKinematics2to3_Update.cpp
 * TKinematics2to3: solving the kinematics (UpdateKinematics and the subsystems)
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
void TKinematics2to3::InitializeSubsystems()
{
  // Create and set the KY-, KN- and YN-subsystems.

  if(!fKYsystem) fKYsystem = new TKinematics2to2("fKYsystem","KY CM subsystem",
						 fMd,fMk,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fKNsystem) fKNsystem = new TKinematics2to2("fKNsystem","KN CM subsystem",
						 fMd,fMk,fMn,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fYNsystem) fYNsystem = new TKinematics2to2("fYNsystem","YN CM subsystem",
						 fMd,fMn,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
}

//_____________________________________________________________________
void TKinematics2to3::UpdateKinematics()
{
  // Solve the kinematics!

  fIsPhysical = true;

  // declare 4vectors in the fAC axis convention.
  nucqcd::LorentzVector vector1,vector2,vector3;
  nucqcd::LorentzVector vector23;
  
  // pointer to subsystem
  TKinematics2to2 *subSystem;
  if( fAC==kK ) subSystem = fYNsystem;
  else if( fAC==kY ) subSystem = fKNsystem;
  else subSystem = fKYsystem;

  // declare final state particle masses
  double m1,m2,m3;
  if( fAC==kK ) { m1 = fMk;
    if( fKinSolution[1]=='Y' ) {
      m2 = fMy;
      m3 = fMn;
    } else {
      m2 = fMn;
      m3 = fMy;
    }
  } else if( fAC==kY ) { m1 = fMy;
    if( fKinSolution[1]=='K' ) {
      m2 = fMk;
      m3 = fMn;
    } else {
      m2 = fMn;
      m3 = fMk;
    }
  } else { m1 = fMn;
    if( fKinSolution[1]=='K' ) {
      m2 = fMk;
      m3 = fMy;
    } else {
      m2 = fMy;
      m3 = fMk;
    }
  }

  fDvector->SetXYZM(0.,0.,0.,fMd); // trivial

  // 1) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    One variable should either specify the photon energy
  //    or the invariant mass of the total system or the
  //    (2+3) subsystem (fEVar).
  //    In addition angular info needs to be given about a
  //    second final particle. Both theta (fTh2Var) and phi (fPhi2Var)
  //    should be either in the cm of the (2+3) subsystem or in the
  //    lab frame. This info will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "1YL": scenario 1, 2nd particle = hyperon with angles in LAB
  //    "1K*": scenario 1, 2nd particle = kaon with angles in CM
  //
  //    The first option with the angles in the LAB system can give
  //    multiple solutions. We advise the user not to use this option
  //    in ReadFormatString(const char*).
  //    TKinematics2to3WitLabAngles can be used for these cases.
  
  if( fKinSolution[0] == '1' ) {

    // Determine 4vector of photon & 1st particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 1st particle
    double e1,p1;
    if( fE1Var==&fPk || fE1Var==&fPy || fE1Var==&fPn ) {
      p1 = *fE1Var;
      e1 = std::sqrt( p1*p1 + m1*m1 );
    } else {
      e1 = *fE1Var;
      p1 = std::sqrt( e1*e1 - m1*m1 );
    }
    if(e1<0. || p1 <0.) fIsPhysical = false;

    // Determine photon energy
    //
    // Special attention is needed when energy is given as s23 or w23
    double eg=-1, s23=-1;
    bool egCalculated = true;
    if( fEVar==&fWlab )
      eg = fWlab;
    else if( fEVar==&fKlab )
      eg = std::sqrt( fKlab*fKlab - fQsquared );
    else if( fEVar==&fWtot )
      eg = ( fWtot*fWtot + fQsquared - fMd*fMd ) / 2. / fMd;
    else if( fEVar==&fStot )
      eg = ( fStot + fQsquared - fMd*fMd ) / 2. / fMd;
    else if( (fEVar==&fWky || fEVar==&fWkn || fEVar==&fWyn) &&
	     (fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn) )
      eg = ( *fEVar**fEVar - *fTh1Var -fMd*fMd + 2.*fMd*e1 ) / 2. / fMd;
    else if( (fEVar==&fSky || fEVar==&fSkn || fEVar==&fSyn) &&
	     (fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn) )
      eg = ( *fEVar - *fTh1Var -fMd*fMd + 2.*fMd*e1 ) / 2. / fMd;
    else {
      // if energy is given by s_{23} or W_{23}, it is less straightforward
      // to find the photon's energy.
      // First, we set s23
      if( fEVar==&fWky || fEVar==&fWkn || fEVar==&fWyn )
	s23 = *fEVar * *fEVar;
      else
	s23 = *fEVar;

      if( s23<(pow(m2+m3,2.)-nucqcd::strange_underflow) ) fIsPhysical = false;

      // we can calculate eg now, if t_{g1} is the angular variable
      if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
	eg = (s23 - *fTh1Var - fMd*fMd + 2.*e1*fMd) /2. /fMd;
      else
	egCalculated = false;
    }

    // Polar angle 1st particle
    double costh1;
    if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
      costh1 = ( *fTh1Var + fQsquared - m1*m1 + 2.*eg*e1 ) 
	/ 2. / p1 / std::sqrt(eg*eg+fQsquared);
    else if( fTh1Var==&fThk || fTh1Var==&fThy || fTh1Var==&fThn )
      costh1 = std::cos( *fTh1Var );
    else
      costh1 = *fTh1Var;
    
    if( std::abs(costh1) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::abs(costh1)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth1=-2.;
    if( 1.-std::abs(costh1) <= kUflow) sinth1=0.;
    else sinth1 = std::sqrt(1. - costh1*costh1); 

    // Photon energy for the remaining case,
    // that is when energy is given by s_{23} or W_{23}
    // and t_{g1} isn't the angular variable.
    //
    // We do this by solving s23 = (p_g + p_D - p_1)^2
    if( !egCalculated ) {
      if( fabs(costh1)<kUflow ) {
	eg = ( s23 + fQsquared - fMd*fMd - m1*m1 + 2.*fMd*e1 ) / 2. / ( fMd - e1 );
      } 
      else {
	std::vector<double> roots = nucqcd::quadratic_real_roots(pow((fMd-e1)/p1/costh1,2)-1.,
		      2.*(fMd-e1)/p1/costh1*(fMd*fMd+m1*m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,
		      pow((fMd*fMd+m1*m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,2)-fQsquared);

	if(roots.size() == 0 ) {
	  fIsPhysical = false;
	} 
	else { 
	  // test the solutions:
	  // * should be positive
	  // * sign of sqrt(eg^2+Q^2) should be positive
	  bool sol1=false, sol2=false;
	  
	  if( roots[0]>0. && (s23 + fQsquared -fMd*fMd -m1*m1 
			      +2.*e1*fMd +2.*roots[0]*(e1-fMd))/costh1>0. ) sol1=true;
	  
	  if( roots.size()==2 ) {
	    if( roots[1]>0. && (s23 + fQsquared -fMd*fMd -m1*m1 
				+2.*e1*fMd +2.*roots[1]*(e1-fMd))/costh1>0. ) sol2=true;
	  }
	  
	  if( roots.size()==1 ) {
	    eg = sol1 ? roots[0] : -1;
	  }
	  else {
	    if( sol1 && sol2 ) {
	      if( fabs(roots[0]-roots[1]) > nucqcd::strange_underflow*fabs(roots[0]))
		std::cout << "WARNING in TKinematics2to3::UpdateKinematics(): "
			     << "found multiple solutions. Choosing one arbitrarily.\n";
	      eg = roots[0];
	    }
	    else {
	      eg = sol1 ? roots[0] : roots[1];
	    }
	  }
	} // solution found
      } // non-zero costh1
    } // if( !egCalculated )

    double threshold = ((m1+m2+m3)*(m1+m2+m3)-fMd*fMd+fQsquared) / 2. / fMd;
    if( eg<threshold-nucqcd::strange_underflow ) fIsPhysical = false;

    // Set 4vectors
    fGvector->SetXYZT(0.,0.,std::sqrt(eg*eg+fQsquared),eg);
    vector1.SetXYZT(p1*sinth1,0.,p1*costh1,e1);
    
    // Determine (23) subsystem
    //~~~~~~~~~~~~~~~~~~~~~~~~~

    vector23 = *fGvector + *fDvector - vector1;
    subSystem->SetVar(1,vector23.Mag2()); // set S_23
    subSystem->SetVar(3,-1.*(*fGvector-vector1).Mag2()); // set -t_g1
       

    // Determine lorentz rotations
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define rotation matrix that changes reference frame to
    // z' ~ p23,  y' ~ y and x' ~ y' x z'
    nucqcd::LorentzRotation rotateTo23; 
    double rotationAngle=0.; 
    if( 1.-fabs(vector23.CosTheta())>kUflow )
      rotationAngle = -1.*std::copysign(1.,vector23.X())
	*vector23.Vect().Angle(fGvector->Vect());
    rotateTo23.RotateY(rotationAngle);
    // Define Lorentzboost that boosts in reference frame defined above
    // to 23 CM system.
    nucqcd::Vector3 boostVector = -1.0 * (rotateTo23*vector23).BoostVector();
    nucqcd::LorentzRotation boostCM23 = boostVector;

 
    // Determine 4vector of 2nd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // CM energy particle 2
    double p2cm = subSystem->GetPk();
    double e2cm = std::sqrt( p2cm*p2cm + m2*m2 );

    if( fKinSolution[2] == '*' ) {
      // CM angle particle 2
      double costh2cm;
      if( fTh2Var==&fKYthkcm || fTh2Var==&fKYthycm || fTh2Var==&fKNthncm ||
	  fTh2Var==&fKNthkcm || fTh2Var==&fYNthycm || fTh2Var==&fYNthncm )
	costh2cm = std::cos(*fTh2Var);
      else if( fTh2Var==&fKYcosthkcm || fTh2Var==&fKYcosthycm || fTh2Var==&fKNcosthncm ||
	       fTh2Var==&fKNcosthkcm || fTh2Var==&fYNcosthycm || fTh2Var==&fYNcosthncm )
	costh2cm = *fTh2Var;
      else {
	double wcm = ( boostCM23*rotateTo23 * *fGvector ).E();
	double kcm = ( boostCM23*rotateTo23 * *fGvector ).P();
	costh2cm = ( *fTh2Var + fQsquared - m2*m2 + 2.*wcm*e2cm ) / 2. / kcm / p2cm;
      }
      
      if( std::abs(costh2cm) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( std::abs(costh2cm)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	costh2cm = ( costh2cm>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // sin polar angle always positive
      // cos -> sin calculation is very unstable
      // when cos close to +-1.
      double sinth2cm=-2.;
      if( 1.-std::abs(costh2cm) <= kUflow ) sinth2cm=0.;
      else sinth2cm = std::sqrt(1. - costh2cm*costh2cm); 

      if( std::abs(sinth2cm) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( std::abs(sinth2cm)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	sinth2cm = ( sinth2cm>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // set 4vector for particle 2 in (23) CM system
      vector2.SetXYZM(p2cm*sinth2cm*std::cos(*fPhi2Var),
		      p2cm*sinth2cm*std::sin(*fPhi2Var),
		      p2cm*costh2cm, m2);

      // calculate 4vector particle 3
      // Remark:
      // one would logically now boost/rotate vector2 back to
      // the LAB frame and then calculate vector3 by taking the
      // difference with vector23. This however (the god of the
      // floating point knows why), introduces a small error (order 1e-9)
      // that keeps propagating. The procedure used below is more
      // accurate. If the reader of this comment understands why,
      // please let me (=Pieter) know.
      vector3.SetXYZM(-vector2.X(),-vector2.Y(),-vector2.Z(),m3);
      vector3 = rotateTo23.Inverse() * (boostCM23.Inverse() * vector3);

      // calculate 4vector particle 2 in the LAB system
      vector2 = vector23 - vector3;

    } // angles 2nd particle in CM system
    else {
      // LAB angle particle 2
      double costh2;
      if( fTh2Var==&fThk || fTh2Var==&fThy || fTh2Var==&fThn )
	costh2 = std::cos( *fTh2Var );
      else
	costh2 = *fTh2Var;
      
      if( std::abs(costh2) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( std::abs(costh2)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;
      
      // sin polar angle always positive
      // cos -> sin calculation is very unstable
      // when cos close to +-1.
      double sinth2=-2.;
      if( 1.-std::abs(costh2) <= kUflow ) sinth2=0.;
      else sinth2 = std::sqrt(1. - costh2*costh2); 
      
      if( std::abs(sinth2) <= 1.0 )
	fIsPhysical = true && fIsPhysical;
      else if( std::abs(sinth2)-1.0 < kUflow ) {
	fIsPhysical = true && fIsPhysical;
	sinth2 = ( sinth2>0.0 ? 1.0 : -1.0 );
      }
      else
	fIsPhysical = false;

      // We wish to determine the 3-momentum of the 2nd particle in the LAB frame
      // Because we have the LAB angles we know the LAB 4-vector except for the
      // 3-momentum. This 4-vector can be boosted to the CM frame.
      // Equating the energy component of this boosted 4-vector with the energy of
      // the 2nd particle in the CM frame (which we calculated supra), leads to a
      // quadratic equation. 
      // The equation can have 2 physical solutions!
      double e23 = vector23.E();
      double w23 = vector23.M();
      double p23 = (*fGvector-vector1).P();
      double costh23 = vector23.CosTheta();
      if( 1.-fabs(costh23)<kUflow ) costh23 = std::copysign(1.,costh23);
      double sinth23 = -2.;
      if( 1.-std::abs(costh23) <= kUflow ) sinth23=0.;
      else sinth23 = std::sqrt(1. - costh23*costh23);

      std::vector<double> roots = nucqcd::quadratic_real_roots(// coefficient p_2^2
		    p23*p23*sinth23*sinth23*sinth2*sinth2*std::cos(*fPhi2Var)*std::cos(*fPhi2Var)
		    +p23*p23*costh23*costh23*costh2*costh2
		    -2*p23*p23*sinth23*costh23*costh2*sinth2*std::cos(*fPhi2Var)
		    -e23*e23,
		    // coefficient p_2
		    2*e2cm*w23*p23*costh23*costh2
		    -2*e2cm*w23*p23*sinth23*sinth2*std::cos(*fPhi2Var),
		    // constant term
		    e2cm*e2cm*w23*w23
		    -e23*e23*m2*m2);
      
      double p2 = -1.;
      if(roots.size() == 0 ) {
	fIsPhysical = false;
      } 
      else { 
	// test the solutions:
	// * should be positive
	// * sign of sqrt(p2^2+m2^2) should be positive
	bool sol1=false, sol2=false;
	  
	if( roots[0]>0. && (e2cm*w23-roots[0]
			    *(p23*sinth23*sinth2*std::cos(*fPhi2Var)
			      -p23*costh23*costh2))>0. ) sol1=true;
	  
	if( roots.size()==2 ) {
	  if( roots[1]>0. && (e2cm*w23-roots[1]
			      *(p23*sinth23*sinth2*std::cos(*fPhi2Var)
				-p23*costh23*costh2))>0. ) sol2=true;
	}
	  
	if( roots.size()==1 ) {
	  p2 = sol1 ? roots[0] : -1;
	}
	else {
	  if( sol1 && sol2 ) {
	    if( fabs(roots[0]-roots[1]) > nucqcd::strange_underflow*fabs(roots[0]))
	      std::cout << "WARNING in TKinematics2to3::UpdateKinematics(): "
			   << "found multiple solutions. Choosing one arbitrarily.\n";
	    p2 = roots[0];
	    std::cout << "sol1= " << roots[0] << std::endl
			 << "sol2= " << roots[1] << std::endl;
	  }
	  else {
	    p2 = sol1 ? roots[0] : roots[1];
	  }
	}
      } // solution found
      
      if(p2<0.) fIsPhysical = false;

      // set 4vector for particle 2
      vector2.SetXYZM(sinth2*std::cos(*fPhi2Var)*p2,
		      sinth2*std::sin(*fPhi2Var)*p2,
		      costh2*p2, m2);

      vector3 = vector23 - vector2;
   
    } // angles 2nd particle in LAB system
   
  } // end 1 + (2+3) case

  // 2) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    In addition a second final particle is fully specified,
  //    energy, theta and phi (fE2Var,fTh2Var,fPhi2Var) in the lab frame.
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "2N-": scenario 2, 2nd particle = nucleon

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                            !!
  // !! No longer considered a valid way of solving the kinematics !!
  // !!                                                            !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else if( fKinSolution[0] == '2' ) {
   
    //     // Determine 4vector of 1st particle
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // Energy 1st particle
    //     double e1,p1;
    //     if( fE1Var==&fPk || fE1Var==&fPy || fE1Var==&fPn ) {
    //       p1 = *fE1Var;
    //       e1 = std::sqrt( p1*p1 + m1*m1 );
    //     } else {
    //       e1 = *fE1Var;
    //       p1 = std::sqrt( e1*e1 - m1*m1 );
    //     }
    //     if(e1<0. || p1 <0.) fIsPhysical = false;

    //     // Polar angle 1st particle
    //     double costh1;
    //     if( fTh1Var==&fThk || fTh1Var==&fThy || fTh1Var==&fThn )
    //       costh1 = std::cos( *fTh1Var );
    //     else
    //       costh1 = *fTh1Var;
    
    //     if( std::abs(costh1) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( std::abs(costh1)-1.0 < kUflow ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;

    //     // set 4vector
    //     vector1.SetXYZT(p1*std::sqrt(1.-costh1*costh1),0.,p1*costh1,e1);

 
    //     // Determine 4vector of 2nd particle
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // Energy 2nd particle
    //     double e2,p2;
    //     if( fE2Var==&fPk || fE2Var==&fPy || fE2Var==&fPn ) {
    //       p2 = *fE2Var;
    //       e2 = std::sqrt( p2*p2 + m2*m2 );
    //     } else {
    //       e2 = *fE2Var;
    //       p2 = std::sqrt( e2*e2 - m2*m2 );
    //     }
    //     if(e2<0. || p2 <0.) fIsPhysical = false;

    //     // Polar angle 2nd particle
    //     double costh2;
    //     if( fTh2Var==&fThk || fTh2Var==&fThy || fTh2Var==&fThn )
    //       costh2 = std::cos( *fTh2Var );
    //     else
    //       costh2 = *fTh2Var;
    //     double sinth2 = std::sqrt(1. -costh2*costh2);
    
    //     if( std::abs(costh2) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( std::abs(costh2)-1.0 < kUflow ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;

    //     // set 4vector
    //     vector2.SetXYZT(sinth2*std::cos(*fPhi2Var)*p2,
    // 		    sinth2*std::sin(*fPhi2Var)*p2,
    // 		    costh2*p2, e2);
    
    //     // Determine 4vector of 3rd particle and photon
    //     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // We do this by writing conservation of 4momentum: p_g+p_D=p_1+p_2+p_3
    //     // in components:
    //     //        wlab + fMd  = e1 + e2 + e3                     (1)
    //     //                 0  = p1x + p2x + p3*sinth3*cos(phi3)  (2)
    //     //                 0  = p2y + p3*sinth3*sin(phi3)        (3)
    //     // sqrt(wlab^2 + Q^2) = p1z + p2z + p3*costh3            (4)
    //     //
    //     // costh3 = +-sqrt(1 - ((p1x+p2x)^2+p2y^2)/p3^2 )        (5) = (2)^2+(3)^2

    //     // Azimuthal angle: tan(phi3) = (3)/(2) = -p2y / (-p1x-p2x)
    //     double phi3 = 
    //       nucqcd::Vector2::Phi_0_2pi(atan2(-1.*vector2.Y(),-1.*(vector1.X()+vector2.X())));

    //     // momentum 3rd particle. Iteratively solve (4) after inserting (1) and (5)
    //     double p3low2=pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y(),
    //       p3high2=p3low2;

    //     // first we look for a lower and upper limit
    //     while( fabs( -1.*(vector1.Z()+vector2.Z()) 
    // 		 + std::sqrt(pow(e1+e2+std::sqrt(p3high2+m3*m3)-fMd,2)+fQsquared) )
    // 	 - std::sqrt(p3high2-(pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y()))
    // 	 > 0. ) {
    //       p3low2 = p3high2;
    //       p3high2 *= 2.;
    //     }
    //     if(p3high2 <= p3low2) fIsPhysical = false;

    //     // solve for p3^2
    //     double delta = p3high2 - p3low2;
    
    //     while( delta/p3low2>1e-9 && p3low2<p3high2 ) {
     
    //       if( fabs( -1.*(vector1.Z()+vector2.Z()) 
    // 		+ std::sqrt(pow(e1+e2+std::sqrt(p3low2+delta/2.+m3*m3)-fMd,2)
    // 		       +fQsquared) )
    // 	  - std::sqrt(p3low2+delta/2.
    // 		 -(pow(vector1.X()+vector2.X(),2)+vector2.Y()*vector2.Y()))
    // 	  < 0. )
    // 	p3high2 -= delta/2.;
    //       else
    // 	p3low2 += delta/2.;
      
    //       delta = p3high2 - p3low2;
    //     }
    
    //     double p3=0.;
    //     if( p3low2>=0. ) p3 = std::sqrt(p3low2);
    //     else fIsPhysical = false;


    //     // energy 3rd particle
    //     double e3 = std::sqrt( p3*p3 + m3*m3 );

    //     // Determine photon 4vector
    //     //~~~~~~~~~~~~~~~~~~~~~~~~~
    //     // photon energy from (1)
    //     double eg = e1 + e2 + e3 - fMd;

    //     // set photon 4vector
    //     fGvector->SetXYZT(0.,0.,std::sqrt(eg*eg+fQsquared),eg);

    //     // Polar angle 3rd particle from (4)
    //     double costh3 = ( std::sqrt(eg*eg+fQsquared)-vector1.Z()-vector2.Z() ) / p3;
    //     if( std::abs(costh3) <= 1.0 )
    //       fIsPhysical = true && fIsPhysical;
    //     else if( std::abs(costh3)-1.0 < 1e-14 ) {
    //       fIsPhysical = true && fIsPhysical;
    //       costh3 = ( costh3>0.0 ? 1.0 : -1.0 );
    //     }
    //     else
    //       fIsPhysical = false;
    //     double sinth3 = std::sqrt(1. - costh3*costh3);

    //     // set 4vector 3rd particle
    //     vector3.SetXYZT(p3*sinth3*std::cos(phi3),
    // 		    p3*sinth3*std::sin(phi3),
    // 		    p3*costh3, e3);

  } // end 1 + 2 (+3) case

  // 3) All variables are mandelstam variables, except for Q^2
  //    obviously. 
  //    S_tot = (P_gamma + P_deut)^2 should always be given.
  //    fAC fixes the 1st particle. This leaves 7 combinations:
  //    (3a) s_12 + s_23 + t_g1 + t_d3
  //    (3b) s_12 + s_23 + t_g1 + t_d2
  //    (3c) s_12 + s_23 + t_d2 + t_d3
  //    (3d) s_12 + t_d2 + t_g1 + t_d3
  //    (3e) s_12 + s_23 + t_g1 + t_g2
  //    (3f) s_12 + s_23 + t_g2 + t_d3
  //    (3g) t_g2 + s_23 + t_g1 + t_d3
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "3Kb": scenario 3b with 2nd particle = kaon

  else if( fKinSolution[0] == '3' ) {

    double *s12,*s23,*tg1,*tg2,*td2,*td3;

    if( fAC==kK ) {
      tg1 = &fTgk;
      s23 = &fSyn;
      if( fKinSolution[1] == 'Y' ) {
	s12 = &fSky;
	tg2 = &fTgy;
	td2 = &fTdy;
	td3 = &fTdn;
      } else {
	s12 = &fSkn;
	tg2 = &fTgn;
	td2 = &fTdn;
	td3 = &fTdy;
      }
    } else if( fAC==kY ) {
      tg1 = &fTgy;
      s23 = &fSkn;
      if( fKinSolution[1] == 'K' ) {
	s12 = &fSky;
	tg2 = &fTgk;
	td2 = &fTdk;
	td3 = &fTdn;
      } else {
	s12 = &fSyn;
	tg2 = &fTgn;
	td2 = &fTdn;
	td3 = &fTdk;
      }
    } else {
      tg1 = &fTgn;
      s23 = &fSky;
      if( fKinSolution[1] == 'K' ) {
	s12 = &fSkn;
	tg2 = &fTgk;
	td2 = &fTdk;
	td3 = &fTdy;
      } else {
	s12 = &fSyn;
	tg2 = &fTgy;
	td2 = &fTdy;
	td3 = &fTdk;
      }
    }

    // determine the remaining mandelstam variables
    switch(fKinSolution[2]) {
    case 'a':
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'b':
      *td3 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td2 - *s23;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'c':
      *tg1 = *td2 + *td3 + *s23 - fMd*fMd - m2*m2 - m3*m3;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'd':
      *s23 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td2 - *td3;
      *tg2 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg1 - *s12;
      break;
    case 'e':
      *td3 = *tg2 + *tg1 + *s12 + fQsquared - m1*m1 - m2*m2;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'f':
      *tg1 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg2 - *s12;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    case 'g':
      *s12 = m1*m1 + m2*m2 - fQsquared + *td3 - *tg2 - *tg1;
      *td2 = m2*m2 + m3*m3 + fMd*fMd + *tg1 - *td3 - *s23;
      break;
    default:
      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::UpdateKinematics(): "
	   << "Invalid option.");
    }

    double threshold = (m1+m2+m3)*(m1+m2+m3);
    if( fStot<threshold-nucqcd::strange_underflow ) fIsPhysical = false;

    // Determine photon 4vector
    //~~~~~~~~~~~~~~~~~~~~~~~~~
    // Photon energy
    double eg = ( fStot + fQsquared - fMd*fMd ) / 2. / fMd;
    if(eg<=0.) fIsPhysical = false;
    double klab = std::sqrt(eg*eg+fQsquared);

    // set photon 4vector
    fGvector->SetXYZT(0.,0.,klab,eg);

    
    // Determine 4vector 1st particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 1st particle
    double e1 = (2.*fMd*eg - *s23 + *tg1 + fMd*fMd) / 2. / fMd;
    if(e1<m1) fIsPhysical = false;
    double p1 = std::sqrt(e1*e1 - m1*m1);

    // Polar angle 1st particle
    double costh1 = ( *tg1 + fQsquared - m1*m1 + 2.*eg*e1 ) / 2. / klab / p1;
    
    if( std::abs(costh1) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::abs(costh1)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh1 = ( costh1>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;
    
    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth1=-2.;
    if( 1.-std::abs(costh1) <= kUflow ) sinth1=0.;
    else sinth1 = std::sqrt(1. - costh1*costh1); 

    // set 1st particle 4vector
    vector1.SetXYZT(p1*sinth1,0.,p1*costh1,e1);


    // Determine (23) subsystem energy wise
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector23 = *fGvector + *fDvector - vector1;
    subSystem->SetVar(1,*s23);
    subSystem->SetVar(3,-1.* *tg1);


    // Determine lorentz rotations
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define rotation matrix that changes reference frame to
    // z' ~ p23,  y' ~ y and x' ~ y' x z'
    nucqcd::LorentzRotation rotateTo23; 
    double rotationAngle=0.; 
    if( 1.-fabs(vector23.CosTheta())>kUflow )
      rotationAngle = -1.*std::copysign(1.,vector23.X())
	*vector23.Vect().Angle(fGvector->Vect());
    rotateTo23.RotateY(rotationAngle);
    // Define Lorentzboost that boosts in reference frame defined above
    // to 23 CM system.
    nucqcd::Vector3 boostVector = -1.0 * (rotateTo23*vector23).BoostVector();
    nucqcd::LorentzRotation boostCM23 = boostVector;


    // Determine 4vector 2nd and 3rd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Energy 2nd particle
    double e2 = ( fMd*fMd + m2*m2 - *td2 ) / 2. / fMd;
    if(e2<m2) fIsPhysical = false;
    double p2 = std::sqrt(e2*e2 - m2*m2);

    // Energy 3rd particle
    double e3 = eg + fMd - e1 - e2;
    if(e3<m3) fIsPhysical = false;
    double p3 = std::sqrt(e3*e3 - m3*m3);

    // Polar angle 2nd particle
    double costh2 = ( *tg2 + fQsquared - m2*m2 + 2.*eg*e2 ) / 2. / klab / p2;
    
    if( std::abs(costh2) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::abs(costh2)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh2 = ( costh2>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth2=-2.;
    if( 1.-std::abs(costh2) <= kUflow ) sinth2=0.;
    else sinth2 = std::sqrt(1. - costh2*costh2);
    
    //express conservation of 4momentum: p_g+p_D=p_1+p_2+p_3 in components:
    //
    //    p3*sinth3*cos(phi3) = -p1x - p2*sinth2*cos(phi2)  (1)
    //    p3*sinth3*sin(phi3) = -p2*sinth2*cos(phi2)        (2)
    
    // Polar angle 3rd particle
    double costh3 = ( vector23.Z() - p2*costh2 ) / p3;
    
    if( std::abs(costh3) <= 1.0 )
      fIsPhysical = true && fIsPhysical;
    else if( std::abs(costh3)-1.0 < kUflow ) {
      fIsPhysical = true && fIsPhysical;
      costh3 = ( costh3>0.0 ? 1.0 : -1.0 );
    }
    else
      fIsPhysical = false;

    // sin polar angle always positive
    // cos -> sin calculation is very unstable
    // when cos close to +-1.
    double sinth3=-2.;
    if( 1.-std::abs(costh3) <= kUflow ) sinth3=0.;
    else sinth3 = std::sqrt(1. - costh3*costh3);

    // Find azimuthal angle particle 2
    double phi2;
    if( vector1.X() != 0. ) {
      // Azimuthal angle particle 2 from (1)^2 + (2)^2
      phi2 =std::acos( ( p3*p3*sinth3*sinth3 - vector1.X()*vector1.X()
		    - p2*p2*sinth2*sinth2 )
		  / 2. / vector1.X() / p2 / sinth2 ) ;
      
      // with the formula above we can only determine Phi_2 upto its sign.
      // This might seem wrong but it isn't. There's a point symmetry.
    }

    else {
      // if particle 1 has no x-component, the x- and y-axis are not properly defined
      phi2=0;      
    }

    // set 4vector for particle 2 again
    vector2.SetXYZT(sinth2*std::cos(phi2)*p2,
		    sinth2*std::sin(phi2)*p2,
		    costh2*p2, e2);
    

    // Determine 4vector of 3rd particle
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vector3 = vector23 - vector2;

  } // end of mandelstam case


  // Put the final state 4vectors in the correct members
  if( fAC==kK ) {
    *fKvector = vector1;
    if( fKinSolution[1] == 'Y' ) {
      *fYvector = vector2;
      *fNvector = vector3;
    } else {
      *fNvector = vector2;
      *fYvector = vector3;
    }
  } else if( fAC==kY ) {
    *fYvector = vector1;
    if( fKinSolution[1] == 'K' ) {
      *fKvector = vector2;
      *fNvector = vector3;
    } else {
      *fNvector = vector2;
      *fKvector = vector3;
    }
  } else {
    *fNvector = vector1;
    if( fKinSolution[1] == 'K' ) {
      *fKvector = vector2;
      *fYvector = vector3;
    } else {
      *fYvector = vector2;
      *fKvector = vector3;
    }
  }

  // The final state 4vectors have been calculated in the user's frame of choice.
  // When necessary we need to rotate them to our choice of frame (=kN).
  if( fAC!=kN ) {
    double angle = -1. * fNvector->Vect().XYvector().Phi();
    fKvector->RotateZ(angle);
    fYvector->RotateZ(angle);
    fNvector->RotateZ(angle);
  }
  
  UpdateSubsystems();

  // Do a final check whether we are in a physical point
  // by checking conservation of energy and momentum
  // nucqcd::LorentzVector sum = *fGvector + *fDvector - *fKvector - *fYvector - *fNvector;
  // if( fabs(sum.E()) > nucqcd::strange_underflow || fabs(sum.X()) > nucqcd::strange_underflow ||
  //     fabs(sum.Y()) > nucqcd::strange_underflow || fabs(sum.Z()) > nucqcd::strange_underflow )
  //   fIsPhysical = false;
  // cout << "after " << fIsPhysical << " " << sum.E() << " " << sum.X() << " " << sum.Y() << " " << sum.Z() << std::endl;

  fIsPhysical *= fKYsystem->IsPhysical() 
    * fKNsystem->IsPhysical() 
    * fYNsystem->IsPhysical();
}

//_____________________________________________________________________
void TKinematics2to3::UpdateSubsystems()
{
  // Recalculate the KY-, KN- and YN-subsystems.

  fKYsystem->SetVar(1,(*fKvector+*fYvector).Mag2()); // set s_ky
  fKYsystem->SetVar(2,(*fGvector-*fNvector-*fKvector).Mag2()); // set (p_g-p_n-p_k)^2
  fKYsystem->SetVar(3,-1.*(*fGvector-*fNvector).Mag2()); // set -(p_g-p_n)^2

  fKNsystem->SetVar(1,(*fKvector+*fNvector).Mag2()); // set s_kn
  fKNsystem->SetVar(2,(*fGvector-*fYvector-*fKvector).Mag2()); // set (p_g-p_y-p_k)^2
  fKNsystem->SetVar(3,-1.*(*fGvector-*fYvector).Mag2()); // set -(p_g-p_y)^2

  fYNsystem->SetVar(1,(*fYvector+*fNvector).Mag2()); // set s_yn
  fYNsystem->SetVar(2,(*fGvector-*fKvector-*fNvector).Mag2()); // set (p_g-p_k-p_N)^2
  fYNsystem->SetVar(3,-1.*(*fGvector-*fKvector).Mag2()); // set -(p_g-p_k)^2
}

//_____________________________________________________________________
void TKinematics2to3::UpdateVariables()
{
  // After UpdateKinematics() not all data member are up-to-date.
  // This function makes sure they all are.

  fWlab = GetWlab();
  fKlab = GetKlab();
  fEk = GetEk();
  fEy = GetEy();
  fEn = GetEn();
  fPk = GetPk();
  fPy = GetPy();
  fPn = GetPn();
  fThk = GetThk();
  fThy = GetThy();
  fThn = GetThn();
  fCosthk = GetCosthk();
  fCosthy = GetCosthy();
  fCosthn = GetCosthn();
  fPhik = GetPhik(kN);
  fPhiy = GetPhiy(kN);
  fPhin = GetPhin(kN);
  fWtot = GetWtot();
  fWky = GetWky();
  fWkn = GetWkn();
  fWyn = GetWyn();
  fStot = GetStot();
  fSky = GetSky();
  fSkn = GetSkn();
  fSyn = GetSyn();
  fTgk = GetTkg();
  fTgy = GetTyg();
  fTgn = GetTng();
  fTdk = (*fDvector-*fKvector).Mag2();
  fTdy = (*fDvector-*fYvector).Mag2();
  fTdn = (*fDvector-*fNvector).Mag2();
  fKYthkcm = std::acos(GetKYsystem_costhkcm());
  fKYthycm = std::acos(GetKYsystem_costhycm());
  fKYcosthkcm = GetKYsystem_costhkcm();
  fKYcosthycm = GetKYsystem_costhycm();
  fKYphikcm = GetKYsystem_phikcm();
  fKYphiycm = GetKYsystem_phiycm();
  fKNthkcm = std::acos(GetKNsystem_costhkcm());
  fKNthncm = std::acos(GetKNsystem_costhncm());
  fKNcosthkcm = GetKNsystem_costhkcm();
  fKNcosthncm = GetKNsystem_costhncm();
  fKNphikcm = GetKNsystem_phikcm();
  fKNphincm = GetKNsystem_phincm();
  fYNthycm = std::acos(GetYNsystem_costhycm());
  fYNthncm = std::acos(GetYNsystem_costhncm());
  fYNcosthycm = GetYNsystem_costhycm();
  fYNcosthncm = GetYNsystem_costhncm();
  fYNphiycm = GetYNsystem_phiycm();
  fYNphincm = GetYNsystem_phincm();
}

}  // namespace nucqcd
