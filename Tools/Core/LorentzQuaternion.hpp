/*! \file LorentzQuaternion.hpp
 * \brief Lorentz transformations as biquaternions
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on November 13, 2010
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef LORENTZQUATERNION_H
#define LORENTZQUATERNION_H

namespace nucqcd { class LorentzVector; }
namespace nucqcd { class Vector3; }

#include <complex>

namespace nucqcd {


// Forward declarations
class Rotation;

/*! \brief Lorentz transformation (boosts and rotations) as a biquaternion Q = [Q0, Q1, Q2, Q3]
 *
 * Products of biquaternions chain transformations; a four-vector is transformed through its quaternion
 * equivalent. See J. Lambek, The Mathematical Intelligencer 17 (4), 7-15 (1995), for the decomposition
 * into a rotation and a boost.
 */
class LorentzQuaternion
{
 public:
  LorentzQuaternion(); ///< identity transformation
  LorentzQuaternion(const double, const nucqcd::Vector3&); ///< rotation over an angle [rad] around an axis
  LorentzQuaternion(const nucqcd::Rotation&); ///< rotation
  LorentzQuaternion(const nucqcd::Vector3&); ///< pure boost with velocity beta
  LorentzQuaternion(const double betaX, const double betaY, const double betaZ); ///< pure boost with velocity (betaX, betaY, betaZ)
  
  LorentzQuaternion& Boost(const nucqcd::Vector3&); ///< chains a pure boost
  LorentzQuaternion& Boost(const double,const double,const double); ///< chains a pure boost with velocity (betaX, betaY, betaZ)
  LorentzQuaternion& Rotate(const double, const nucqcd::Vector3&); ///< chains a rotation over an angle [rad] around an axis
  LorentzQuaternion& Rotate(const double, const nucqcd::Vector3*); ///< chains a rotation over an angle [rad] around an axis
  LorentzQuaternion& RotateX(const double); ///< chains a rotation around the x axis [rad]
  LorentzQuaternion& RotateY(const double); ///< chains a rotation around the y axis [rad]
  LorentzQuaternion& RotateZ(const double); ///< chains a rotation around the z axis [rad]

  LorentzQuaternion& Invert(); ///< inverts the transformation
  LorentzQuaternion Inverse() const; ///< inverse transformation
  LorentzQuaternion& Conjugate(); ///< conjugates: [Q0, -Q1, -Q2, -Q3]
  LorentzQuaternion& ComplexConjugate(); ///< complex conjugates all components
  LorentzQuaternion& Hermitian(); ///< Hermitian conjugate: conjugate and complex conjugate
  static LorentzQuaternion Conjugate(const LorentzQuaternion&); ///< conjugate [Q0, -Q1, -Q2, -Q3]
  static LorentzQuaternion ComplexConjugate(const LorentzQuaternion&); ///< complex conjugate
  static LorentzQuaternion Hermitian(const LorentzQuaternion&); ///< Hermitian conjugate

  LorentzQuaternion operator*(const LorentzQuaternion&) const; ///< quaternion product (chained transformation)
  nucqcd::LorentzVector operator*(const nucqcd::LorentzVector&) const; ///< transforms a four-vector
  nucqcd::Vector3 operator*(const nucqcd::Vector3&) const; ///< transforms a three-vector
  LorentzQuaternion& operator*=(const LorentzQuaternion&); ///< a *= b sets a = a*b
  LorentzQuaternion& Transform(const LorentzQuaternion&); ///< a.Transform(b) sets a = b*a

  bool IsIdentity() const; ///< is this the identity?
  void Decompose(LorentzQuaternion& rotation, LorentzQuaternion& boost) const; ///< decomposes into a rotation and a boost (Lambek 1995)
void Decompose(double& theta, nucqcd::Vector3& axis, nucqcd::Vector3& boost) const; ///< decomposes into a rotation (angle [rad] and axis) and a boost (velocity)

  friend std::ostream& operator<<(std::ostream&, const LorentzQuaternion&); ///< prints the components

 protected:
  void SetBoost(const double betaX, const double betaY, const double betaZ); ///< sets a pure boost
  void SetRotation(const double, const nucqcd::Vector3&); ///< sets a rotation over an angle [rad] around an axis

 private:
  LorentzQuaternion(const std::complex<double>& q0,const std::complex<double>& q1,const std::complex<double>& q2,const std::complex<double>& q3); ///< constructor from the four components
 
 private:
  std::complex<double> fQ0; ///< component Q0
  std::complex<double> fQ1; ///< component Q1
  std::complex<double> fQ2; ///< component Q2
  std::complex<double> fQ3; ///< component Q3
 

}; // class LorentzQuaternion


}  // namespace nucqcd

/** @} */
#endif // LORENTZQUATERNION_H
