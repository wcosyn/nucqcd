/*
 * LorentzQuaternion.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * Work started on November 13, 2010
 * 
 */

#include <Tools/Core/LorentzQuaternion.hpp>
 #include <Tools/Core/constants.hpp>
#include <complex>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/LorentzVector.hpp>
#include <Tools/Core/Vector3.hpp>
#include <cmath>
#include <cassert>
#include <iostream>

namespace nucqcd {

using std::complex;
using std::sqrt;
using std::cosh;
using std::sinh;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::cosh;
using std::sinh;
using std::atan2;

//////////////////////////////////////////////////////////////////////////
//
// LorentzQuaternion
//
// Representation of Lorentz transformationS using biquaternions.
//
// |------------|
// | Definition |
// |------------|
// A biquaternion Q is a set of four complex numbers:
//BEGIN_LATEX Q = [Q_{0},Q_{1},Q_{2},Q_{3}] = [Q_{0}, #vec{Q} ]END_LATEX
// for which we define addition and multiplication:
//BEGIN_LATEX
//Q + Q' = [ Q_{0}+Q_{0}, #vec{Q}+#vec{Q}' ]
//QQ' = [ Q_{0}Q'_{0} - #vec{Q}#upoint#vec{Q}', Q_{0}#vec{Q}' + Q'_{0}#vec{Q} + #vec{Q}#times#vec{Q}']
//END_LATEX
// In addition we define conjugation (Conjugate()), complex conjugation (ComplexConjugate())
// and hermitian conjugation (Hermitian()):
//BEGIN_LATEX
//#bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]
//Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]
//Q^{+} = #bar{Q}*
//END_LATEX
//
// |------------------------|
// | Lorentz transformation |
// |------------------------|
// A general Lorentz transformation can be represented by a biquaterion of
// unit norm, i.e.
//BEGIN_LATEX Q#bar{Q}=1 END_LATEX
// Consecutive transformations can be chained by multiplying from the left.
// Representing a general 4-vector v as a quaternion V
//BEGIN_LATEX v#equiv(t,#vec{x}) #Rightarrow V#equiv[t,i#vec{x}]END_LATEX
// the 4-vector can be Lorentz transformed applying
//BEGIN_LATEX V' =  Q V Q^{+}END_LATEX
// with Q the transformation.
//
// A pure rotation over angle theta around a unit axis n is represented by
//BEGIN_LATEX
//Q_{rotation}(#theta,#vec{n}) = [ cos#frac{#theta}{2}, #vec{n}sin#frac{#theta}{2}]
//END_LATEX
// A pure boost with boost vector beta and rapidity ksi is represented by
//BEGIN_LATEX
//Q_{boost}(#vec{#beta}) = [ cosh#frac{#xi}{2}, #frac{#vec{#beta}}{#beta}sinh#frac{#xi}{2}]
//END_LATEX
//
// |---------------|
// | Decomposition |
// |---------------|
// Any Lorentz transformation can be decomposed in a pure rotation r followed by
// a pure boost b. When the transformation is represented by a biquaternion Q
// this decomposition can easily be performed:
//BEGIN_LATEX
//Q #equiv U + i V  (U,V ordinary quaternions)
//#mu = #sqrt{UU^{+}}
//r := #frac{U}{#mu}
//s := [ #mu, -i #frac{UV^{+}}{#mu} ] 
//END_LATEX
// Decomposition of a LorentzQuaternion is done with the methods Decompose(..).
//
//////////////////////////////////////////////////////////////////////////


//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion(const std::complex<double>& q0,
				       const std::complex<double>& q1,
				       const std::complex<double>& q2,
				       const std::complex<double>& q3)
: fQ0(q0), fQ1(q1), fQ2(q2), fQ3(q3)
{
  // Private constructor
}

//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion()
  : fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize identity transformation
}

//_____________________________________________________________________

//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion(const double theta, const nucqcd::Vector3& axis)
  : fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize rotation over angle 'theta' around axis defined by 'axis'
  SetRotation(theta,axis);
}

//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion(const nucqcd::Rotation& rotation)
  : fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize rotation
  double theta;
  nucqcd::Vector3 axis;
  rotation.AngleAxis(theta,axis);

  SetRotation(theta,axis);
}

//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion(const nucqcd::Vector3& beta)
  : fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize pure Lorentz boost defined by boost vector 'beta'
  SetBoost(beta.X(),beta.Y(),beta.Z());
}

//_____________________________________________________________________
LorentzQuaternion::LorentzQuaternion(const double betaX, 
				       const double betaY, 
				       const double betaZ)
  : fQ0(complex<double>(1.)), fQ1(complex<double>(0.)),
    fQ2(complex<double>(0.)), fQ3(complex<double>(0.))
{
  // Initialize pure Lorentz boost defined by components of boost vector
  SetBoost(betaX,betaY,betaZ);
}

//_____________________________________________________________________
void LorentzQuaternion::SetBoost(const double betaX, 
				  const double betaY, 
				  const double betaZ)
{
  // Create a pure boost (defined by the components of the boost vector).
  const double beta = sqrt(betaX*betaX
			   +betaY*betaY
			   +betaZ*betaZ);
  const double coshksi2 = sqrt( (1./sqrt(1.-beta*beta) + 1.)/2. );
  const double sinhksi2 = sqrt( (1./sqrt(1.-beta*beta) - 1.)/2. );

  fQ0 = complex<double>(coshksi2);
  if( fabs(beta)>nucqcd::strange_underflow ) {
    fQ1 = complex<double>(0.,betaX/beta*sinhksi2);
    fQ2 = complex<double>(0.,betaY/beta*sinhksi2);
    fQ3 = complex<double>(0.,betaZ/beta*sinhksi2);
  }
}

//_____________________________________________________________________
void LorentzQuaternion::SetRotation(const double theta, const nucqcd::Vector3& axis)
{
  // Create a rotation over angle 'theta' around axis defined by 'axis'
  const double norm = axis.Mag();
  if( fabs(norm)<nucqcd::strange_underflow ) {
    cerr << "ERROR in LorentzQuaternion::SetRotation(const double, "
	 << "const nucqcd::Vector3&): "
	 << "invalid rotation axis." << endl;
    assert(1==0);
  }
  const double sinThetaOver2 = sin(theta/2.);

  fQ0 = complex<double>(cos(theta/2.));
  fQ1 = complex<double>(axis.X()/norm*sinThetaOver2);
  fQ2 = complex<double>(axis.Y()/norm*sinThetaOver2);
  fQ3 = complex<double>(axis.Z()/norm*sinThetaOver2);
}

//_____________________________________________________________________
ostream& operator<<(ostream& stream, const LorentzQuaternion& q)
{
  stream << "[" << q.fQ0 << "," << q.fQ1 << "," << q.fQ2 << "," << q.fQ3 << "]";
  return stream;
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Conjugate()
{
  // Conjugate the biquaternion
  //BEGIN_LATEX #bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]END_LATEX
  static const complex<double> minus(-1.,0.);
  
  fQ1 *= minus;
  fQ2 *= minus;
  fQ3 *= minus;

  return *this;
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::ComplexConjugate()
{
  // Complex conjugate the biquaternion
  //BEGIN_LATEX Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]END_LATEX
  fQ0 = conj(fQ0);
  fQ1 = conj(fQ1);
  fQ2 = conj(fQ2);
  fQ3 = conj(fQ3);

  return *this;
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Hermitian()
{
  // Hermitian conjugate the biquaternion
  //BEGIN_LATEX Q^{+} = #bar{Q}*END_LATEX
  static const complex<double> minus(-1.,0.);

  fQ0 = conj(fQ0);
  fQ1 = conj(fQ1) * minus;
  fQ2 = conj(fQ2) * minus;
  fQ3 = conj(fQ3) * minus;

  return *this;
}

//_____________________________________________________________________
LorentzQuaternion LorentzQuaternion::Conjugate(const LorentzQuaternion& rhs)
{
  // Conjugate the biquaternion
  //BEGIN_LATEX #bar{Q} = [Q_{0},-Q_{1},-Q_{2},-Q_{3}]END_LATEX
  LorentzQuaternion lhs = rhs;
  return lhs.Conjugate();
}

//_____________________________________________________________________
LorentzQuaternion LorentzQuaternion::ComplexConjugate(const LorentzQuaternion& rhs)
{
  // Complex conjugate the biquaternion
  //BEGIN_LATEX Q* = [Q*_{0},Q*_{1},Q*_{2},Q*_{3}]END_LATEX
  LorentzQuaternion lhs = rhs;
  return lhs.ComplexConjugate();
}

//_____________________________________________________________________
LorentzQuaternion LorentzQuaternion::Hermitian(const LorentzQuaternion& rhs)
{
  // Hermitian conjugate the biquaternion
  //BEGIN_LATEX Q^{+} = #bar{Q}*END_LATEX
  LorentzQuaternion lhs = rhs;
  return lhs.Hermitian();
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Invert()
{
  // Invert the Lorentz transformation
  return Conjugate();
}

//_____________________________________________________________________
LorentzQuaternion LorentzQuaternion::Inverse() const
{
  // Returns the inverse Lorentz transformation
  return Conjugate(*this);
}

//_____________________________________________________________________
LorentzQuaternion LorentzQuaternion::operator*(const LorentzQuaternion& rhs) const
{
  // Multiply with a biquaternion from the right:
  //BEGIN_LATEX QQ' = [ Q_{0}Q'_{0} - #vec{Q}#upoint#vec{Q}', Q_{0}#vec{Q}' + Q'_{0}#vec{Q} + #vec{Q}#times#vec{Q}'] END_LATEX
  LorentzQuaternion result;
  
  result.fQ0 = fQ0*rhs.fQ0 - fQ1*rhs.fQ1 - fQ2*rhs.fQ2 - fQ3*rhs.fQ3;
  result.fQ1 = fQ0*rhs.fQ1 + rhs.fQ0*fQ1 + fQ2*rhs.fQ3 - fQ3*rhs.fQ2;
  result.fQ2 = fQ0*rhs.fQ2 + rhs.fQ0*fQ2 + fQ3*rhs.fQ1 - fQ1*rhs.fQ3;
  result.fQ3 = fQ0*rhs.fQ3 + rhs.fQ0*fQ3 + fQ1*rhs.fQ2 - fQ2*rhs.fQ1;

  return result;
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::operator*=(const LorentzQuaternion& rhs)
{
  // Matrix multiplication
  // Note a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;
  return operator=( (*this) * rhs );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Transform(const LorentzQuaternion& rhs)
{
  // Chain LorentzTransformations.
  // Note a.transform(b); <=> a = b * a; while a *= b; <=> a = a * b;
  return operator=( rhs * (*this) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Boost(const nucqcd::Vector3& boostVector)
{
  // Chain pure Lorentz boost to transformation.
  return Transform( LorentzQuaternion(boostVector) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Boost(const double betaX, 
					      const double betaY, 
					      const double betaZ)
{
  // Chain pure Lorentz boost to transformation.
  return Transform( LorentzQuaternion(betaX,betaY,betaZ) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Rotate(const double theta,
					       const nucqcd::Vector3& axis)
{
  // Chain a rotation over 'theta' around 'axis' to transformation
  return Transform( LorentzQuaternion(theta,axis) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::Rotate(const double theta,
					       const nucqcd::Vector3* axis)
{
  // Chain a rotation over 'theta' around 'axis' to transformation
  return Transform( LorentzQuaternion(theta,*axis) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::RotateX(const double theta)
{
  // Chain a rotation over 'theta' around x-axis to transformation
  return Transform( LorentzQuaternion(theta,nucqcd::Vector3(1.,0.,0.)) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::RotateY(const double theta)
{
  // Chain a rotation over 'theta' around y-axis to transformation
  return Transform( LorentzQuaternion(theta,nucqcd::Vector3(0.,1.,0.)) );
}

//_____________________________________________________________________
LorentzQuaternion& LorentzQuaternion::RotateZ(const double theta)
{
  // Chain a rotation over 'theta' around z-axis to transformation
  return Transform( LorentzQuaternion(theta,nucqcd::Vector3(0.,0.,1.)) );
}

//_____________________________________________________________________
nucqcd::LorentzVector LorentzQuaternion::operator*(const nucqcd::LorentzVector& v) const
{
  // Transform the 4-vector

  // Define quaternion equivalent of 4-vector
  LorentzQuaternion q(complex<double>(v.E()),
		       complex<double>(0.,v.X()),
		       complex<double>(0.,v.Y()),
		       complex<double>(0.,v.Z()));

  // Transform the quaternion
  LorentzQuaternion qp = (*this) * q * Hermitian(*this);

  // Go back from quaternion to 4-vector
  return nucqcd::LorentzVector(qp.fQ1.imag(),
			qp.fQ2.imag(),
			qp.fQ3.imag(),
			qp.fQ0.real());
}

//_____________________________________________________________________
nucqcd::Vector3 LorentzQuaternion::operator*(const nucqcd::Vector3& v) const
{
  // Transform the 3-vector
  nucqcd::LorentzVector vp = (*this) * nucqcd::LorentzVector(v.X(),v.Y(),v.Z(),0.);
  return nucqcd::Vector3(vp.X(),vp.Y(),vp.Z());
}

//_____________________________________________________________________
bool LorentzQuaternion::IsIdentity() const
{
  return (fabs(fQ1.real())<nucqcd::strange_underflow) && (fabs(fQ1.imag())<nucqcd::strange_underflow)
    && (fabs(fQ2.real())<nucqcd::strange_underflow) && (fabs(fQ2.imag())<nucqcd::strange_underflow)
    && (fabs(fQ3.real())<nucqcd::strange_underflow) && (fabs(fQ3.imag())<nucqcd::strange_underflow)
    && (fabs(fQ0.imag())<nucqcd::strange_underflow) && (fabs(fQ0.real()-1.)<nucqcd::strange_underflow);
}

//_____________________________________________________________________
void LorentzQuaternion::Decompose(LorentzQuaternion& rotation, 
				   LorentzQuaternion& boost) const
{
  // Decompose the current Lorentz transformation into a rotation + boost.
  // Details can be found in J.Lambek, The mathematical intelligencer, 
  // Volume 17, Number 4, 7-15 (1995) (http://dx.doi.org/10.1007/BF03024783)

  // Decompose the quaternion in a real and imaginary part
  LorentzQuaternion u(fQ0.real(),
		       fQ1.real(),
		       fQ2.real(),
		       fQ3.real());
  LorentzQuaternion v(fQ0.imag(),
		       fQ1.imag(),
		       fQ2.imag(),
		       fQ3.imag());
  
  // Find the norm of the real part
  const double mu = sqrt( (u*Conjugate(u)).fQ0.real() );

  // Define the rotation
  rotation= u;
  rotation.fQ0 /= mu;
  rotation.fQ1 /= mu;
  rotation.fQ2 /= mu;
  rotation.fQ3 /= mu;

  // Define the boost
  boost = u*Conjugate(v);
  const complex<double> factor(0.,-1./mu);
  boost.fQ0 = mu;
  boost.fQ1 *= factor;
  boost.fQ2 *= factor;
  boost.fQ3 *= factor;
}

//_____________________________________________________________________
void LorentzQuaternion::Decompose(double& theta, nucqcd::Vector3& axis,
				   nucqcd::Vector3& boost) const
{
  // Decompose the current Lorentz transformation into a rotation + boost.
  // The rotation is given as an angle 'theta' with its rotation axis.
  // The boost is given by its boost vector.

  // Decompose the transformation
  LorentzQuaternion r,b;
  Decompose(r,b);

  // Find the angle+axis of the rotation
  const double sintheta2 = sqrt(r.fQ1.real()*r.fQ1.real()
				+r.fQ2.real()*r.fQ2.real()
				+r.fQ3.real()*r.fQ3.real());
  theta = 2.*atan2(sintheta2,r.fQ0.real());
  if( fabs(sintheta2)>nucqcd::strange_underflow )
    axis = nucqcd::Vector3(r.fQ1.real()/sintheta2,
		    r.fQ2.real()/sintheta2,
		    r.fQ3.real()/sintheta2);
  else 
    axis = nucqcd::Vector3(0.,0.,1.);
    
  // Find the boost vector
  const double sinhksi2 = sqrt(b.fQ1.imag()*b.fQ1.imag()
			       +b.fQ2.imag()*b.fQ2.imag()
			       +b.fQ3.imag()*b.fQ3.imag());
  const double beta = tanh(2.*log(b.fQ0.real()+sinhksi2));
  if( fabs(sinhksi2)>nucqcd::strange_underflow)
    boost = nucqcd::Vector3(beta*b.fQ1.imag()/sinhksi2,
		     beta*b.fQ2.imag()/sinhksi2,
		     beta*b.fQ3.imag()/sinhksi2 );
  else 
    boost = nucqcd::Vector3(0.,0.,0.);
}

}  // namespace nucqcd

