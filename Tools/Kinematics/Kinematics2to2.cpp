/*
 * Kinematics2to2.cpp
 *
 * Author:        Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Kinematics2to2                                                       //
//                                                                       //
// Kinematics2to2 allows to calculate kinematics for 2 -> 2 scattering. //
//                                                                       //
// All quantities are in MeV or rad.                                     //
//                                                                       //
// Kinematics2to2 is an extension of the Kinematics class defined      //
// in the strangecalc-wrapper project.                                   //
// The user can go beyond the 6 standard 'isospin' channels by setting   //
// the nucleon, kaon and hyperon mass and choosing Q^2 intelligently.    //
//                                                                       // 
// Because Kinematics was not intended for general 2->2 scattering, the //
// getters and setters can have awkward names. This should not be an     //
// insurmountable problem however.                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <Tools/Kinematics/Kinematics2to2.hpp>

namespace nucqcd {



//_____________________________________________________________________
Kinematics2to2::Kinematics2to2(const char *name,const char *title,
				 double mN, double mK, double mY,
				 const char *formatString,
				 double var1,double var2,double var3)
: Kinematics(name,title,1,formatString,var1,var2,var3)
{
  // Constructor
  //
  // Q^2,mN = masses of initial state particles
  // mK,mY  = masses of final state particles
  //
  // Remember all quantities should be in MeV.
  //
  // See Kinematics::Kinematics for more info

  SetMasses(mN,mK,mY);
}

//_____________________________________________________________________

//_____________________________________________________________________
void Kinematics2to2::SetMasses(double mN, double mK, double mY)
{
  // mN = mass of initial state particle
  // mK,mY  = masses of final state particles 
  //
  // Masses should be in MeV.

  SetMn(mN);
  SetMk(mK);
  SetMy(mY);
  UpdateKinematics();
}

}  // namespace nucqcd

