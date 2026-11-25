/*
 * Kinematics2to3_Update.cpp
 * Kinematics2to3: solving the kinematics (UpdateKinematics and the subsystems)
 * Part of the implementation of Kinematics2to3 (split from Kinematics2to3.cpp, which holds the class
 * description, constructors and assignment); author of the original: Pieter Vancraeyveld.
 */

#include <Tools/Kinematics/Kinematics2to3.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/LorentzQuaternion.hpp>
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
void Kinematics2to3::InitializeSubsystems()
{
  // Create and set the KY-, KN- and YN-subsystems.

  if(!fKYsystem) fKYsystem = new Kinematics2to2("fKYsystem","KY CM subsystem",
						 fMd,fMk,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fKNsystem) fKNsystem = new Kinematics2to2("fKNsystem","KN CM subsystem",
						 fMd,fMk,fMn,"s:t:qsquared",
						 0.0,0.0,0.0);
  if(!fYNsystem) fYNsystem = new Kinematics2to2("fYNsystem","YN CM subsystem",
						 fMd,fMn,fMy,"s:t:qsquared",
						 0.0,0.0,0.0);
}


//_____________________________________________________________________
/*! \brief Intermediate results of UpdateKinematics, passed between the steps of the solution */
struct Kinematics2to3::UpdateWork {
  nucqcd::LorentzVector vector1; ///< the final state particle that fixes the x axis (fAC axis convention)
  nucqcd::LorentzVector vector2; ///< the second final state particle
  nucqcd::LorentzVector vector3; ///< the third final state particle
  nucqcd::LorentzVector vector23;                  ///< (2+3) subsystem
  double m1 = 0.; ///< mass of particle 1
  double m2 = 0.; ///< mass of particle 2
  double m3 = 0.; ///< mass of particle 3
  Kinematics2to2 *subSystem = nullptr;            ///< the CM subsystem of particles 2 and 3
  nucqcd::LorentzRotation rotateTo23;              ///< rotation to the frame with z' along p23
  nucqcd::LorentzRotation boostCM23;               ///< boost from there to the (2+3) CM
  double p2cm = 0.; ///< momentum of particle 2 in the (2+3) CM
  double e2cm = 0.; ///< energy of particle 2 in the (2+3) CM
  double *s12 = nullptr; ///< scenario 3: mandelstam variable s_12
  double *s23 = nullptr; ///< scenario 3: mandelstam variable s_23
  double *tg1 = nullptr; ///< scenario 3: mandelstam variable t_g1
  double *tg2 = nullptr; ///< scenario 3: mandelstam variable t_g2
  double *td2 = nullptr; ///< scenario 3: mandelstam variable t_d2
  double *td3 = nullptr; ///< scenario 3: mandelstam variable t_d3
  double eg = 0.; ///< scenario 3: photon energy
  double klab = 0.; ///< scenario 3: photon momentum
  double e1 = 0.; ///< scenario 3: energy of particle 1
};

//_____________________________________________________________________
void Kinematics2to3::ClampCosine(double &c)
{
  // A cosine just outside [-1,1] (by less than kUflow) is rounding: put it at the edge; further out the
  // point is not physical.
  // (The point stays physical in the first two cases: the former `fIsPhysical = true && fIsPhysical`.)
  if( std::abs(c) <= 1.0 ) {}
  else if( std::abs(c)-1.0 < kUflow ) c = ( c>0.0 ? 1.0 : -1.0 );
  else fIsPhysical = false;
}

//_____________________________________________________________________
double Kinematics2to3::SineFromCosine(double c)
{
  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  if( 1.-std::abs(c) <= kUflow ) return 0.;
  return std::sqrt(1. - c*c);
}

//_____________________________________________________________________
void Kinematics2to3::RotateTo23Frame(UpdateWork &w)
{
  // Determine lorentz rotations
  // Define rotation matrix that changes reference frame to
  // z' ~ p23,  y' ~ y and x' ~ y' x z'
  double rotationAngle=0.; 
  if( 1.-fabs(w.vector23.CosTheta())>kUflow )
    rotationAngle = -1.*std::copysign(1.,w.vector23.X())
      *w.vector23.Vect().Angle(fGvector->Vect());
  w.rotateTo23.RotateY(rotationAngle);
  // Define Lorentzboost that boosts in reference frame defined above
  // to 23 CM system.
  nucqcd::Vector3 boostVector = -1.0 * (w.rotateTo23*w.vector23).BoostVector();
  w.boostCM23 = boostVector;
}

//_____________________________________________________________________
void Kinematics2to3::UpdateKinematics()
{
  // Solve the kinematics!

  fIsPhysical = true;

  UpdateWork w;

  
  // the CM subsystem of particles 2 and 3
  if( fAC==kK ) w.subSystem = fYNsystem;
  else if( fAC==kY ) w.subSystem = fKNsystem;
  else w.subSystem = fKYsystem;

  // the masses of the final state particles
  if( fAC==kK ) { w.m1 = fMk;
    if( fKinSolution[1]=='Y' ) {
      w.m2 = fMy;
      w.m3 = fMn;
    } else {
      w.m2 = fMn;
      w.m3 = fMy;
    }
  } else if( fAC==kY ) { w.m1 = fMy;
    if( fKinSolution[1]=='K' ) {
      w.m2 = fMk;
      w.m3 = fMn;
    } else {
      w.m2 = fMn;
      w.m3 = fMk;
    }
  } else { w.m1 = fMn;
    if( fKinSolution[1]=='K' ) {
      w.m2 = fMk;
      w.m3 = fMy;
    } else {
      w.m2 = fMy;
      w.m3 = fMk;
    }
  }

  fDvector->SetXYZM(0.,0.,0.,fMd); // trivial

  // The solution depends on the variables the user gave (fKinSolution[0], see ReadFormatString):
  // 1) an energy, the angle of the particle on the x axis and angles of a second particle
  // 3) all variables are mandelstam variables
  // (2, a second particle fully specified in the lab frame, is refused by ReadFormatString.)
  if( fKinSolution[0] == '1' ) SolveEnergyScenario(w);
  else if( fKinSolution[0] == '3' ) SolveMandelstamScenario(w);

// Put the final state 4vectors in the correct members
if( fAC==kK ) {
  *fKvector = w.vector1;
  if( fKinSolution[1] == 'Y' ) {
    *fYvector = w.vector2;
    *fNvector = w.vector3;
  } else {
    *fNvector = w.vector2;
    *fYvector = w.vector3;
  }
} else if( fAC==kY ) {
  *fYvector = w.vector1;
  if( fKinSolution[1] == 'K' ) {
    *fKvector = w.vector2;
    *fNvector = w.vector3;
  } else {
    *fNvector = w.vector2;
    *fKvector = w.vector3;
  }
} else {
  *fNvector = w.vector1;
  if( fKinSolution[1] == 'K' ) {
    *fKvector = w.vector2;
    *fYvector = w.vector3;
  } else {
    *fYvector = w.vector2;
    *fKvector = w.vector3;
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
void Kinematics2to3::SolveEnergyScenario(UpdateWork &w)
{
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
  //    Kinematics2to3WitLabAngles can be used for these cases.

  SolvePhotonAndParticle1(w);
  RotateTo23Frame(w);

  // Determine 4vector of 2nd particle
  // CM energy particle 2
  w.p2cm = w.subSystem->GetPk();
  w.e2cm = std::sqrt( w.p2cm*w.p2cm + w.m2*w.m2 );

  if( fKinSolution[2] == '*' ) SolveParticle2InCM(w);
  else SolveParticle2InLab(w);
}

//_____________________________________________________________________
void Kinematics2to3::SolvePhotonAndParticle1(UpdateWork &w)
{
  // Determine 4vector of photon & 1st particle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Energy 1st particle
  double e1,p1;
  if( fE1Var==&fPk || fE1Var==&fPy || fE1Var==&fPn ) {
    p1 = *fE1Var;
    e1 = std::sqrt( p1*p1 + w.m1*w.m1 );
  } else {
    e1 = *fE1Var;
    p1 = std::sqrt( e1*e1 - w.m1*w.m1 );
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

    if( s23<(pow(w.m2+w.m3,2.)-nucqcd::strange_underflow) ) fIsPhysical = false;

    // we can calculate eg now, if t_{g1} is the angular variable
    if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
	eg = (s23 - *fTh1Var - fMd*fMd + 2.*e1*fMd) /2. /fMd;
    else
	egCalculated = false;
  }

  // Polar angle 1st particle
  double costh1;
  if( fTh1Var==&fTgk || fTh1Var==&fTgy || fTh1Var==&fTgn )
    costh1 = ( *fTh1Var + fQsquared - w.m1*w.m1 + 2.*eg*e1 ) 
	/ 2. / p1 / std::sqrt(eg*eg+fQsquared);
  else if( fTh1Var==&fThk || fTh1Var==&fThy || fTh1Var==&fThn )
    costh1 = std::cos( *fTh1Var );
  else
    costh1 = *fTh1Var;

  ClampCosine(costh1);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth1 = SineFromCosine(costh1);

  // Photon energy for the remaining case,
  // that is when energy is given by s_{23} or W_{23}
  // and t_{g1} isn't the angular variable.
  //
  // We do this by solving s23 = (p_g + p_D - p_1)^2
  if( !egCalculated ) {
    if( fabs(costh1)<kUflow ) {
	eg = ( s23 + fQsquared - fMd*fMd - w.m1*w.m1 + 2.*fMd*e1 ) / 2. / ( fMd - e1 );
    } 
    else {
	std::vector<double> roots = nucqcd::quadratic_real_roots(pow((fMd-e1)/p1/costh1,2)-1.,
		      2.*(fMd-e1)/p1/costh1*(fMd*fMd+w.m1*w.m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,
		      pow((fMd*fMd+w.m1*w.m1-2.*fMd*e1-s23-fQsquared)/2./p1/costh1,2)-fQsquared);

	if(roots.size() == 0 ) {
	  fIsPhysical = false;
	} 
	else { 
	  // test the solutions:
	  // * should be positive
	  // * sign of sqrt(eg^2+Q^2) should be positive
	  bool sol1=false, sol2=false;

	  if( roots[0]>0. && (s23 + fQsquared -fMd*fMd -w.m1*w.m1 
			      +2.*e1*fMd +2.*roots[0]*(e1-fMd))/costh1>0. ) sol1=true;

	  if( roots.size()==2 ) {
	    if( roots[1]>0. && (s23 + fQsquared -fMd*fMd -w.m1*w.m1 
				+2.*e1*fMd +2.*roots[1]*(e1-fMd))/costh1>0. ) sol2=true;
	  }

	  if( roots.size()==1 ) {
	    eg = sol1 ? roots[0] : -1;
	  }
	  else {
	    if( sol1 && sol2 ) {
	      if( fabs(roots[0]-roots[1]) > nucqcd::strange_underflow*fabs(roots[0]))
		std::cout << "WARNING in Kinematics2to3::UpdateKinematics(): "
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

  double threshold = ((w.m1+w.m2+w.m3)*(w.m1+w.m2+w.m3)-fMd*fMd+fQsquared) / 2. / fMd;
  if( eg<threshold-nucqcd::strange_underflow ) fIsPhysical = false;

  // Set 4vectors
  fGvector->SetXYZT(0.,0.,std::sqrt(eg*eg+fQsquared),eg);
  w.vector1.SetXYZT(p1*sinth1,0.,p1*costh1,e1);

  // Determine (23) subsystem
  //~~~~~~~~~~~~~~~~~~~~~~~~~

  w.vector23 = *fGvector + *fDvector - w.vector1;
  w.subSystem->SetVar(1,w.vector23.Mag2()); // set S_23
  w.subSystem->SetVar(3,-1.*(*fGvector-w.vector1).Mag2()); // set -t_g1
}

//_____________________________________________________________________
void Kinematics2to3::SolveParticle2InCM(UpdateWork &w)
{
  // CM angle particle 2
  double costh2cm;
  if( fTh2Var==&fKYthkcm || fTh2Var==&fKYthycm || fTh2Var==&fKNthncm ||
	  fTh2Var==&fKNthkcm || fTh2Var==&fYNthycm || fTh2Var==&fYNthncm )
	costh2cm = std::cos(*fTh2Var);
  else if( fTh2Var==&fKYcosthkcm || fTh2Var==&fKYcosthycm || fTh2Var==&fKNcosthncm ||
	       fTh2Var==&fKNcosthkcm || fTh2Var==&fYNcosthycm || fTh2Var==&fYNcosthncm )
	costh2cm = *fTh2Var;
  else {
	double wcm = ( w.boostCM23*w.rotateTo23 * *fGvector ).E();
	double kcm = ( w.boostCM23*w.rotateTo23 * *fGvector ).P();
	costh2cm = ( *fTh2Var + fQsquared - w.m2*w.m2 + 2.*wcm*w.e2cm ) / 2. / kcm / w.p2cm;
  }

  ClampCosine(costh2cm);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth2cm = SineFromCosine(costh2cm);

  ClampCosine(sinth2cm);

  // set 4vector for particle 2 in (23) CM system
  w.vector2.SetXYZM(w.p2cm*sinth2cm*std::cos(*fPhi2Var),
		      w.p2cm*sinth2cm*std::sin(*fPhi2Var),
		      w.p2cm*costh2cm, w.m2);

  // calculate 4vector particle 3
  // Remark:
  // one would logically now boost/rotate w.vector2 back to
  // the LAB frame and then calculate w.vector3 by taking the
  // difference with w.vector23. This however (the god of the
  // floating point knows why), introduces a small error (order 1e-9)
  // that keeps propagating. The procedure used below is more
  // accurate. If the reader of this comment understands why,
  // please let me (=Pieter) know.
  w.vector3.SetXYZM(-w.vector2.X(),-w.vector2.Y(),-w.vector2.Z(),w.m3);
  w.vector3 = w.rotateTo23.Inverse() * (w.boostCM23.Inverse() * w.vector3);

  // calculate 4vector particle 2 in the LAB system
  w.vector2 = w.vector23 - w.vector3;
}

//_____________________________________________________________________
void Kinematics2to3::SolveParticle2InLab(UpdateWork &w)
{
  // LAB angle particle 2
  double costh2;
  if( fTh2Var==&fThk || fTh2Var==&fThy || fTh2Var==&fThn )
	costh2 = std::cos( *fTh2Var );
  else
	costh2 = *fTh2Var;

  ClampCosine(costh2);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth2 = SineFromCosine(costh2);

  ClampCosine(sinth2);

  // We wish to determine the 3-momentum of the 2nd particle in the LAB frame
  // Because we have the LAB angles we know the LAB 4-vector except for the
  // 3-momentum. This 4-vector can be boosted to the CM frame.
  // Equating the energy component of this boosted 4-vector with the energy of
  // the 2nd particle in the CM frame (which we calculated supra), leads to a
  // quadratic equation. 
  // The equation can have 2 physical solutions!
  double e23 = w.vector23.E();
  double w23 = w.vector23.M();
  double p23 = (*fGvector-w.vector1).P();
  double costh23 = w.vector23.CosTheta();
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
		    2*w.e2cm*w23*p23*costh23*costh2
		    -2*w.e2cm*w23*p23*sinth23*sinth2*std::cos(*fPhi2Var),
		    // constant term
		    w.e2cm*w.e2cm*w23*w23
		    -e23*e23*w.m2*w.m2);

  double p2 = -1.;
  if(roots.size() == 0 ) {
	fIsPhysical = false;
  } 
  else { 
	// test the solutions:
	// * should be positive
	// * sign of sqrt(p2^2+w.m2^2) should be positive
	bool sol1=false, sol2=false;

	if( roots[0]>0. && (w.e2cm*w23-roots[0]
			    *(p23*sinth23*sinth2*std::cos(*fPhi2Var)
			      -p23*costh23*costh2))>0. ) sol1=true;

	if( roots.size()==2 ) {
	  if( roots[1]>0. && (w.e2cm*w23-roots[1]
			      *(p23*sinth23*sinth2*std::cos(*fPhi2Var)
				-p23*costh23*costh2))>0. ) sol2=true;
	}

	if( roots.size()==1 ) {
	  p2 = sol1 ? roots[0] : -1;
	}
	else {
	  if( sol1 && sol2 ) {
	    if( fabs(roots[0]-roots[1]) > nucqcd::strange_underflow*fabs(roots[0]))
	      std::cout << "WARNING in Kinematics2to3::UpdateKinematics(): "
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
  w.vector2.SetXYZM(sinth2*std::cos(*fPhi2Var)*p2,
		      sinth2*std::sin(*fPhi2Var)*p2,
		      costh2*p2, w.m2);

  w.vector3 = w.vector23 - w.vector2;

}

//_____________________________________________________________________
void Kinematics2to3::SolveMandelstamScenario(UpdateWork &w)
{
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

  SelectMandelstamVariables(w);
  SolveMandelstamParticle1(w);
  RotateTo23Frame(w);
  SolveMandelstamParticles23(w);
}

//_____________________________________________________________________
void Kinematics2to3::SelectMandelstamVariables(UpdateWork &w)
{
  // the mandelstam variables for the axis convention and second particle, and the remaining ones from the given four
  if( fAC==kK ) {
    w.tg1 = &fTgk;
    w.s23 = &fSyn;
    if( fKinSolution[1] == 'Y' ) {
	w.s12 = &fSky;
	w.tg2 = &fTgy;
	w.td2 = &fTdy;
	w.td3 = &fTdn;
    } else {
	w.s12 = &fSkn;
	w.tg2 = &fTgn;
	w.td2 = &fTdn;
	w.td3 = &fTdy;
    }
  } else if( fAC==kY ) {
    w.tg1 = &fTgy;
    w.s23 = &fSkn;
    if( fKinSolution[1] == 'K' ) {
	w.s12 = &fSky;
	w.tg2 = &fTgk;
	w.td2 = &fTdk;
	w.td3 = &fTdn;
    } else {
	w.s12 = &fSyn;
	w.tg2 = &fTgn;
	w.td2 = &fTdn;
	w.td3 = &fTdk;
    }
  } else {
    w.tg1 = &fTgn;
    w.s23 = &fSky;
    if( fKinSolution[1] == 'K' ) {
	w.s12 = &fSkn;
	w.tg2 = &fTgk;
	w.td2 = &fTdk;
	w.td3 = &fTdy;
    } else {
	w.s12 = &fSyn;
	w.tg2 = &fTgy;
	w.td2 = &fTdy;
	w.td3 = &fTdk;
    }
  }

  // determine the remaining mandelstam variables
  switch(fKinSolution[2]) {
  case 'a':
    *w.tg2 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg1 - *w.s12;
    *w.td2 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td3 - *w.s23;
    break;
  case 'b':
    *w.td3 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td2 - *w.s23;
    *w.tg2 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg1 - *w.s12;
    break;
  case 'c':
    *w.tg1 = *w.td2 + *w.td3 + *w.s23 - fMd*fMd - w.m2*w.m2 - w.m3*w.m3;
    *w.tg2 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg1 - *w.s12;
    break;
  case 'd':
    *w.s23 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td2 - *w.td3;
    *w.tg2 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg1 - *w.s12;
    break;
  case 'e':
    *w.td3 = *w.tg2 + *w.tg1 + *w.s12 + fQsquared - w.m1*w.m1 - w.m2*w.m2;
    *w.td2 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td3 - *w.s23;
    break;
  case 'f':
    *w.tg1 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg2 - *w.s12;
    *w.td2 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td3 - *w.s23;
    break;
  case 'g':
    *w.s12 = w.m1*w.m1 + w.m2*w.m2 - fQsquared + *w.td3 - *w.tg2 - *w.tg1;
    *w.td2 = w.m2*w.m2 + w.m3*w.m3 + fMd*fMd + *w.tg1 - *w.td3 - *w.s23;
    break;
  default:
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in Kinematics2to3::UpdateKinematics(): "
	   << "Invalid option.");
  }
}

//_____________________________________________________________________
void Kinematics2to3::SolveMandelstamParticle1(UpdateWork &w)
{
  double threshold = (w.m1+w.m2+w.m3)*(w.m1+w.m2+w.m3);
  if( fStot<threshold-nucqcd::strange_underflow ) fIsPhysical = false;

  // Determine photon 4vector
  //~~~~~~~~~~~~~~~~~~~~~~~~~
  // Photon energy
  w.eg = ( fStot + fQsquared - fMd*fMd ) / 2. / fMd;
  if(w.eg<=0.) fIsPhysical = false;
  w.klab = std::sqrt(w.eg*w.eg+fQsquared);

  // set photon 4vector
  fGvector->SetXYZT(0.,0.,w.klab,w.eg);

  // Determine 4vector 1st particle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Energy 1st particle
  w.e1 = (2.*fMd*w.eg - *w.s23 + *w.tg1 + fMd*fMd) / 2. / fMd;
  if(w.e1<w.m1) fIsPhysical = false;
  double p1 = std::sqrt(w.e1*w.e1 - w.m1*w.m1);

  // Polar angle 1st particle
  double costh1 = ( *w.tg1 + fQsquared - w.m1*w.m1 + 2.*w.eg*w.e1 ) / 2. / w.klab / p1;

  ClampCosine(costh1);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth1 = SineFromCosine(costh1);

  // set 1st particle 4vector
  w.vector1.SetXYZT(p1*sinth1,0.,p1*costh1,w.e1);

  // Determine (23) subsystem energy wise
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  w.vector23 = *fGvector + *fDvector - w.vector1;
  w.subSystem->SetVar(1,*w.s23);
  w.subSystem->SetVar(3,-1.* *w.tg1);
}

//_____________________________________________________________________
void Kinematics2to3::SolveMandelstamParticles23(UpdateWork &w)
{
  // Determine 4vector 2nd and 3rd particle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Energy 2nd particle
  double e2 = ( fMd*fMd + w.m2*w.m2 - *w.td2 ) / 2. / fMd;
  if(e2<w.m2) fIsPhysical = false;
  double p2 = std::sqrt(e2*e2 - w.m2*w.m2);

  // Energy 3rd particle
  double e3 = w.eg + fMd - w.e1 - e2;
  if(e3<w.m3) fIsPhysical = false;
  double p3 = std::sqrt(e3*e3 - w.m3*w.m3);

  // Polar angle 2nd particle
  double costh2 = ( *w.tg2 + fQsquared - w.m2*w.m2 + 2.*w.eg*e2 ) / 2. / w.klab / p2;

  ClampCosine(costh2);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth2 = SineFromCosine(costh2);

  //express conservation of 4momentum: p_g+p_D=p_1+p_2+p_3 in components:
  //
  //    p3*sinth3*cos(phi3) = -p1x - p2*sinth2*cos(phi2)  (1)
  //    p3*sinth3*sin(phi3) = -p2*sinth2*cos(phi2)        (2)

  // Polar angle 3rd particle
  double costh3 = ( w.vector23.Z() - p2*costh2 ) / p3;

  ClampCosine(costh3);

  // sin polar angle always positive
  // cos -> sin calculation is very unstable
  // when cos close to +-1.
  double sinth3 = SineFromCosine(costh3);

  // Find azimuthal angle particle 2
  double phi2;
  if( w.vector1.X() != 0. ) {
    // Azimuthal angle particle 2 from (1)^2 + (2)^2
    phi2 =std::acos( ( p3*p3*sinth3*sinth3 - w.vector1.X()*w.vector1.X()
		    - p2*p2*sinth2*sinth2 )
		  / 2. / w.vector1.X() / p2 / sinth2 ) ;

    // with the formula above we can only determine Phi_2 upto its sign.
    // This might seem wrong but it isn't. There's a point symmetry.
  }

  else {
    // if particle 1 has no x-component, the x- and y-axis are not properly defined
    phi2=0;      
  }

  // set 4vector for particle 2 again
  w.vector2.SetXYZT(sinth2*std::cos(phi2)*p2,
		    sinth2*std::sin(phi2)*p2,
		    costh2*p2, e2);

  // Determine 4vector of 3rd particle
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  w.vector3 = w.vector23 - w.vector2;
}



//_____________________________________________________________________
void Kinematics2to3::UpdateSubsystems()
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
void Kinematics2to3::UpdateVariables()
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
