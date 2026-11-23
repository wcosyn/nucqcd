/*! \file TLorentzQuaternion.hpp
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

// Forward declarations
namespace nucqcd { class Rotation; }

/*! \brief Lorentz transformation (boosts and rotations) as a biquaternion Q = [Q0, Q1, Q2, Q3]
 *
 * Products of biquaternions chain transformations; a four-vector is transformed through its quaternion
 * equivalent. See J. Lambek, The Mathematical Intelligencer 17 (4), 7-15 (1995), for the decomposition
 * into a rotation and a boost.
 */
class TLorentzQuaternion
{
 public:
  TLorentzQuaternion(); ///< identity transformation
  TLorentzQuaternion(const double, const nucqcd::Vector3&); ///< rotation over an angle [rad] around an axis
  TLorentzQuaternion(const nucqcd::Rotation&); ///< rotation
  TLorentzQuaternion(const nucqcd::Vector3&); ///< pure boost with velocity beta
  TLorentzQuaternion(const double betaX, const double betaY, const double betaZ); ///< pure boost with velocity (betaX, betaY, betaZ)
  
  TLorentzQuaternion& Boost(const nucqcd::Vector3&); ///< chains a pure boost
  TLorentzQuaternion& Boost(const double,const double,const double); ///< chains a pure boost with velocity (betaX, betaY, betaZ)
  TLorentzQuaternion& Rotate(const double, const nucqcd::Vector3&); ///< chains a rotation over an angle [rad] around an axis
  TLorentzQuaternion& Rotate(const double, const nucqcd::Vector3*); ///< chains a rotation over an angle [rad] around an axis
  TLorentzQuaternion& RotateX(const double); ///< chains a rotation around the x axis [rad]
  TLorentzQuaternion& RotateY(const double); ///< chains a rotation around the y axis [rad]
  TLorentzQuaternion& RotateZ(const double); ///< chains a rotation around the z axis [rad]

  TLorentzQuaternion& Invert(); ///< inverts the transformation
  TLorentzQuaternion Inverse() const; ///< inverse transformation
  TLorentzQuaternion& Conjugate(); ///< conjugates: [Q0, -Q1, -Q2, -Q3]
  TLorentzQuaternion& ComplexConjugate(); ///< complex conjugates all components
  TLorentzQuaternion& Hermitian(); ///< Hermitian conjugate: conjugate and complex conjugate
  static TLorentzQuaternion Conjugate(const TLorentzQuaternion&); ///< conjugate [Q0, -Q1, -Q2, -Q3]
  static TLorentzQuaternion ComplexConjugate(const TLorentzQuaternion&); ///< complex conjugate
  static TLorentzQuaternion Hermitian(const TLorentzQuaternion&); ///< Hermitian conjugate

  TLorentzQuaternion operator*(const TLorentzQuaternion&) const; ///< quaternion product (chained transformation)
  nucqcd::LorentzVector operator*(const nucqcd::LorentzVector&) const; ///< transforms a four-vector
  nucqcd::Vector3 operator*(const nucqcd::Vector3&) const; ///< transforms a three-vector
  TLorentzQuaternion& operator*=(const TLorentzQuaternion&); ///< a *= b sets a = a*b
  TLorentzQuaternion& Transform(const TLorentzQuaternion&); ///< a.Transform(b) sets a = b*a

  bool IsIdentity() const; ///< is this the identity?
  void Decompose(TLorentzQuaternion& rotation, TLorentzQuaternion& boost) const; ///< decomposes into a rotation and a boost (Lambek 1995)
void Decompose(double& theta, nucqcd::Vector3& axis, nucqcd::Vector3& boost) const; ///< decomposes into a rotation (angle [rad] and axis) and a boost (velocity)

  friend std::ostream& operator<<(std::ostream&, const TLorentzQuaternion&); ///< prints the components

 protected:
  void SetBoost(const double betaX, const double betaY, const double betaZ); ///< sets a pure boost
  void SetRotation(const double, const nucqcd::Vector3&); ///< sets a rotation over an angle [rad] around an axis

 private:
  TLorentzQuaternion(const std::complex<double>& q0,const std::complex<double>& q1,const std::complex<double>& q2,const std::complex<double>& q3); ///< constructor from the four components
 
 private:
  std::complex<double> fQ0; ///< component Q0
  std::complex<double> fQ1; ///< component Q1
  std::complex<double> fQ2; ///< component Q2
  std::complex<double> fQ3; ///< component Q3
 

}; // class TLorentzQuaternion

/** @} */
#endif // LORENTZQUATERNION_H
