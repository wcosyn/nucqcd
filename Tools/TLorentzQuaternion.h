/*! \file TLorentzQuaternion.h
 * \brief Lorentz transformations as biquaternions
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on November 13, 2010
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef LORENTZQUATERNION_H
#define LORENTZQUATERNION_H

// Include statements for classes that could be merely forward declared
// but need a full declaration to allow the Dict.cpp file to be compiled. 
#ifdef ROOT_Rtypes
#include <TLorentzVector.h>
#include <TVector3.h>
#else
class TLorentzVector;
class TVector3;
#endif

#include <TObject.h>
#include <complex>

// Forward declarations
class TRotation;

/*! \brief Lorentz transformation (boosts and rotations) as a biquaternion Q = [Q0, Q1, Q2, Q3]
 *
 * Products of biquaternions chain transformations; a four-vector is transformed through its quaternion
 * equivalent. See J. Lambek, The Mathematical Intelligencer 17 (4), 7-15 (1995), for the decomposition
 * into a rotation and a boost.
 */
class TLorentzQuaternion : public TObject
{
 public:
  TLorentzQuaternion(); ///< identity transformation
  TLorentzQuaternion(TRootIOCtor*); ///< ROOT I/O constructor
  TLorentzQuaternion(const double, const TVector3&); ///< rotation over an angle [rad] around an axis
  TLorentzQuaternion(const TRotation&); ///< rotation
  TLorentzQuaternion(const TVector3&); ///< pure boost with velocity beta
  TLorentzQuaternion(const double betaX, const double betaY, const double betaZ); ///< pure boost with velocity (betaX, betaY, betaZ)
  
  TLorentzQuaternion& Boost(const TVector3&); ///< chains a pure boost
  TLorentzQuaternion& Boost(const double,const double,const double); ///< chains a pure boost with velocity (betaX, betaY, betaZ)
  TLorentzQuaternion& Rotate(const double, const TVector3&); ///< chains a rotation over an angle [rad] around an axis
  TLorentzQuaternion& Rotate(const double, const TVector3*); ///< chains a rotation over an angle [rad] around an axis
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
  TLorentzVector operator*(const TLorentzVector&) const; ///< transforms a four-vector
  TVector3 operator*(const TVector3&) const; ///< transforms a three-vector
  TLorentzQuaternion& operator*=(const TLorentzQuaternion&); ///< a *= b sets a = a*b
  TLorentzQuaternion& Transform(const TLorentzQuaternion&); ///< a.Transform(b) sets a = b*a

  bool IsIdentity() const; ///< is this the identity?
  void Decompose(TLorentzQuaternion& rotation, TLorentzQuaternion& boost) const; ///< decomposes into a rotation and a boost (Lambek 1995)
void Decompose(double& theta, TVector3& axis, TVector3& boost) const; ///< decomposes into a rotation (angle [rad] and axis) and a boost (velocity)

  friend std::ostream& operator<<(std::ostream&, const TLorentzQuaternion&); ///< prints the components

 protected:
  void SetBoost(const double betaX, const double betaY, const double betaZ); ///< sets a pure boost
  void SetRotation(const double, const TVector3&); ///< sets a rotation over an angle [rad] around an axis

 private:
  TLorentzQuaternion(const std::complex<double>& q0,const std::complex<double>& q1,const std::complex<double>& q2,const std::complex<double>& q3); ///< constructor from the four components
 
 private:
  std::complex<double> fQ0; ///< component Q0
  std::complex<double> fQ1; ///< component Q1
  std::complex<double> fQ2; ///< component Q2
  std::complex<double> fQ3; ///< component Q3
 
  ClassDef(TLorentzQuaternion,1); // Lorentz transformations (boosts and rotations) represented by a biquaternion

}; // class TLorentzQuaternion

/** @} */
#endif // LORENTZQUATERNION_H
