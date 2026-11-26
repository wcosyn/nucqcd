/*! \file Kinematics2to2.hpp
 * \brief 2 -> 2 kinematics with masses set by the user
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef KINEMATICS2TO2_H
#define KINEMATICS2TO2_H

#include <Tools/Kinematics/Kinematics.hpp>

namespace nucqcd {


/*! \brief Kinematics with the masses of N, K and Y set by the user (e.g. A(e,e'p) with N the target
 * nucleus, K the residual nucleus and Y the knocked-out nucleon) and three independent variables */
class Kinematics2to2 : public Kinematics
{
 public:
  /*! constructor
   * \param mN,mK,mY [MeV] masses of the target N, the meson K and the baryon Y
   * \param formatString format string with the three independent variables, see Kinematics::SetFormat
   * \param var1,var2,var3 initial values of the independent variables */
  Kinematics2to2(double mN, double mK, double mY,
		  const char *formatString,
		  double var1,double var2,double var3);
  
  virtual void SetMasses(double mN, double mK, double mY); ///< sets the masses [MeV] and updates the kinematics

 protected:
  virtual double Underflow() const { return 1e-6; } ///< underflow parameter (1e-6)
  
    
};



}  // namespace nucqcd

/** @} */
#endif
