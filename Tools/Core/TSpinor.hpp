/*! \file TSpinor.hpp
 * \brief Dirac spinors (TSpinor, TVSpinor), helicity spinors (THelicitySpinor, THelicityVSpinor),
 * spin polarization states and Lorentz transformations in spinor space
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on August 10, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TSPINOR_H
#define TSPINOR_H

// Include statements for classes that could be merely forward declared
// but need a full declaration to allow the Dict.cpp file to be compiled. 
#ifdef ROOT_Rtypes
#include <Tools/Core/TLorentzQuaternion.hpp>
#include <TLorentzRotation.h>
#else
class TLorentzQuaternion;
class TLorentzRotation;
#endif

#include <Tools/Core/Matrix.hpp>
#include <Tools/Core/FourVector.hpp>
#include <complex>
#include <Tools/Core/GammaStructure.hpp>
#include <TVector3.h>
#include <iostream>

// Forward declarations
class TLorentzVector;
class TRotation;

/*! \brief Dirac spinor u(p,s) of a positive-energy spin-1/2 particle
 *
 * Conventions of Bjorken & Drell and Gross (eq. (A.15) of Gross, Relativistic Quantum Mechanics and
 * Field Theory), in the Dirac representation:
 * \f[ u(\vec p,s) = \sqrt{E+m} \begin{pmatrix} \chi^s \\ \frac{\vec\sigma\cdot\vec p}{E+m}\chi^s \end{pmatrix},
 *     \qquad E = \sqrt{\vec p^{\,2}+m^2}, \f]
 * normalized as \f$ \bar u u = 2m \f$ (kDoubleMass, the default) or \f$ \bar u u = 1 \f$ (kUnity). Only
 * the three-momentum is used: the energy is always the on-shell one. The two-component spinor
 * \f$ \chi^s \f$ has spin projection s along the quantization axis of the Polarization, see
 * Polarization::GetPauliSpinor. The negative-energy spinor v is TVSpinor; helicity spinors are
 * THelicitySpinor and THelicityVSpinor.
 */
class TSpinor
{
 public:
  // forward declarations
  class Polarization;

  /*! normalization of the spinors */
  enum Normalization { kUnity=1,     ///< Bar(U)U = 1
		       kDoubleMass=2 ///< Bar(U)U = 2M
  };

 protected:
  TSpinor(double px,double py,double pz,double mass, const Polarization&, Normalization,int energyState); ///< constructor for positive (energyState > 0) or negative (< 0) energy solutions, used by TVSpinor; momentum (px,py,pz) [MeV], mass [MeV]
 public:
  TSpinor(const TLorentzVector&, double mass, const Polarization&, Normalization=kDoubleMass); ///< constructor: four-momentum [MeV] (only its spatial part is used), mass [MeV], polarization, normalization
  TSpinor(const FourVector<double>&, double mass, const Polarization&, Normalization=kDoubleMass); ///< as the TLorentzVector constructor
  TSpinor(const TSpinor&); ///< copy constructor
  virtual ~TSpinor(); ///< destructor
  TSpinor& operator=(const TSpinor&); ///< assignment
  virtual TSpinor* Clone(const char* ="") const; ///< virtual copy constructor

  operator const Matrix<4,1>&() const { return *fComponent; } ///< the four components as a 4x1 matrix
  /// outer product of the spinor (4x1) with a 1xK matrix
  template<int K>
    Matrix<4,K> operator*(const Matrix<1,K>& rhs) { return (Matrix<4,1>)(*this) * rhs; }
  
  /// prints the components
  friend 
    std::ostream& operator<<(std::ostream&, const TSpinor&);
  
  static Matrix<1,4> Bar(const TSpinor&); ///< Dirac adjoint \f$ \bar u = u^\dagger \gamma^0 \f$
  static const FourVector< Matrix<2,2> > kSigmaPauli; ///< four-vector of Pauli matrices (1, sigma_x, sigma_y, sigma_z)
  
 protected:
  void InitializeSpinor(double px,double py,double pz,double mass,const Polarization&, Normalization,int energyState); ///< fills the components, see the class description
  
 protected:
  /// the four components of the spinor
  Matrix<4,1> *fComponent; //! a spinor is a 4-component vector
 private:

  
  /*! \brief Spin polarization state of a spin-1/2 particle: spin up or down along the quantization
   * axis with polar angle theta and azimuthal angle phi (the z axis by default in most of the code) */
 public:
  class Polarization
  {
  public:
    /*! spin projection +1/2 (kUp) or -1/2 (kDown) along the quantization axis */
    enum State { kUp=1,
		 kDown=-1 };

  public:   
    Polarization(double,double,State); ///< constructor: theta [rad], phi [rad] of the quantization axis, and the state
    Polarization(const Polarization&); ///< copy constructor
    virtual ~Polarization(); ///< destructor
    Polarization& operator=(const Polarization&); ///< assignment

    bool operator==(const Polarization&); ///< same axis and state
    bool operator!=(const Polarization& pol) { return !((*this)==pol); } ///< not the same axis and state
    bool operator==(State polarizationState); ///< same state
    bool operator!=(State polarizationState) { return !((*this)==polarizationState); } ///< different state
    operator int() const; ///< twice the spin projection (+1 or -1)
    operator State() const; ///< the state

    void SetTheta(double theta); ///< sets the polar angle of the quantization axis [rad]
    void SetPhi(double phi); ///< sets the azimuthal angle of the quantization axis [rad]
    void SetState(State polState) { fPolarizationState=polState; } ///< sets the state

    Polarization& operator++(); ///< raises the state to kUp
    Polarization& operator--(); ///< lowers the state to kDown

    Polarization Inverse() const; ///< a copy with the opposite state
    Polarization& Invert(); ///< flips the state

    Matrix<2,1> GetPauliSpinor(int eSolution=1) const; ///< two-component spinor: for particles (eSolution > 0) \f$ \chi^{+1/2} = (\cos\frac\theta2, e^{i\phi}\sin\frac\theta2) \f$, \f$ \chi^{-1/2} = (-e^{-i\phi}\sin\frac\theta2, \cos\frac\theta2) \f$; for antiparticles (eSolution < 0) \f$ \eta^s = \chi^{-s} \f$

  private:
    double            fTheta; ///< [rad] polar angle of the quantization axis
    double            fPhi; ///< [rad] azimuthal angle of the quantization axis
    State fPolarizationState; ///< spin state


  }; // class TSpinor::Polarization

  /*! \brief Lorentz transformation (boosts and rotations) in spinor space, S(Lambda) as a 4x4 GammaStructure */
 public:
  class LorentzRotation
  {
  public:
    LorentzRotation(); ///< unity
    LorentzRotation(const TRotation&); ///< pure rotation
    LorentzRotation(const TVector3&); ///< pure boost with velocity beta
    LorentzRotation(double betaX, double betaY, double betaZ); ///< pure boost with velocity (betaX, betaY, betaZ)

    // Const Cast to GammaStructure
    operator const GammaStructure&() const { return fTransformation; } ///< the transformation matrix

    // Cast to GammaStructure
    operator GammaStructure&() { return fTransformation; } ///< the transformation matrix

        /// prints the transformation matrix
    friend 
      std::ostream& operator<<(std::ostream&, const LorentzRotation&);

    LorentzRotation& Boost(const TVector3& boostVector) { return Transform(LorentzRotation(boostVector)); } ///< followed by a boost with velocity boostVector
    LorentzRotation& Boost(double betaX, double betaY, double betaZ) { return Transform(LorentzRotation(betaX,betaY,betaZ)); } ///< followed by a boost with velocity (betaX, betaY, betaZ)
    LorentzRotation  Inverse() const; ///< inverse transformation
    LorentzRotation& Invert(); ///< inverts the transformation
    TSpinor          operator*(const TSpinor&) const; ///< transforms a spinor
    LorentzRotation  operator*(const LorentzRotation&) const; ///< compound transformation
    LorentzRotation& operator*=(const LorentzRotation&); ///< compound transformation
    GammaStructure   operator*(const GammaStructure&) const; ///< product with a 4x4 matrix
    LorentzRotation& Rotate(double angle, const TVector3& axis); ///< followed by a rotation over angle [rad] around axis
    LorentzRotation& Rotate(double angle, const TVector3* axis) { return Rotate(angle,*axis); } ///< followed by a rotation over angle [rad] around axis
    LorentzRotation& RotateX(double angle) { return Rotate(angle,TVector3(1.,0.,0.)); } ///< followed by a rotation around the x axis [rad]
    LorentzRotation& RotateY(double angle) { return Rotate(angle,TVector3(0.,1.,0.)); } ///< followed by a rotation around the y axis [rad]
    LorentzRotation& RotateZ(double angle) { return Rotate(angle,TVector3(0.,0.,1.)); } ///< followed by a rotation around the z axis [rad]
    LorentzRotation& Transform(const LorentzRotation&); ///< a.Transform(b) sets a = b*a

  protected:
    void SetBoost(double betaX, double betaY, double betaZ); ///< sets a pure boost, \f$ S = \cosh\frac\xi2 + \sinh\frac\xi2\, \hat\beta_i\gamma^0\gamma^i \f$ with rapidity xi
    void SetRotation(double thetaX, double thetaY, double thetaZ); ///< sets a pure rotation over the angle vector (thetaX, thetaY, thetaZ) [rad]
		     
  private:
    GammaStructure fTransformation; ///< transformation matrix


  }; // class TSpinor::LorentzRotation  

}; // class TSpinor

/*! \brief Dirac spinor v(p,s) of a negative-energy solution (antiparticle), with
 * \f$ v = \sqrt{E+m}\,(\frac{\vec\sigma\cdot\vec p}{E+m}\eta^s, \eta^s) \f$ and \f$ \eta^s = \chi^{-s} \f$, see TSpinor */
class TVSpinor : public TSpinor
{
 public:
  TVSpinor(const TLorentzVector&, double mass, const Polarization&, Normalization=kDoubleMass); ///< constructor, see TSpinor
  TVSpinor(const FourVector<double>&, double mass, const Polarization&, Normalization=kDoubleMass); ///< constructor, see TSpinor
  TVSpinor(const TVSpinor&); ///< copy constructor
  virtual ~TVSpinor(); ///< destructor
  TVSpinor& operator=(const TVSpinor&); ///< assignment
  virtual TVSpinor* Clone(const char* ="") const; ///< virtual copy constructor


}; // class TVSpinor

/*! \brief Helicity spinor of a positive-energy spin-1/2 particle (spin along its momentum), with the
 * helicity transformation and Wick rotations of E. Leader, Spin in Particle Physics (sections 1.2.2 and 2.1) */
class THelicitySpinor : public TSpinor
{
 public:
  THelicitySpinor(const TLorentzVector&, double mass, const Polarization::State&, Normalization=kDoubleMass); ///< constructor: four-momentum [MeV], mass [MeV], helicity state, normalization
  THelicitySpinor(const FourVector<double>&, double mass, const Polarization::State&, Normalization=kDoubleMass); ///< as the TLorentzVector constructor
  THelicitySpinor(const THelicitySpinor&); ///< copy constructor
  virtual ~THelicitySpinor(); ///< destructor
  THelicitySpinor& operator=(const THelicitySpinor&); ///< assignment
  virtual THelicitySpinor* Clone(const char* ="") const; ///< virtual copy constructor

  static TLorentzRotation HelicityTransformation(const TLorentzVector&); ///< helicity transformation \f$ h(p) = r_z(\phi) r_y(\theta) l_z(|\vec p|/p^0) \f$ (Leader, section 1.2.2)
  static TSpinor::LorentzRotation HelicitySpinorTransformation(const TLorentzVector&); ///< helicity transformation h(p) in spinor space
  static TLorentzQuaternion HelicityQuaternionTransformation(const TLorentzVector&); ///< helicity transformation h(p) as a quaternion
  static TLorentzRotation WickRotation(const TLorentzRotation&, const TLorentzVector&); ///< Wick helicity rotation \f$ r(l,p) = h(lp)^{-1}\, l\, h(p) \f$ (Leader, eq. (2.1.7)); l acts on four-vectors, not on the frame
  static TLorentzQuaternion WickRotation(const TLorentzQuaternion&, const TLorentzVector&); ///< Wick helicity rotation, quaternion version
  static void WickRotationEulerAngles(const TLorentzRotation&, const TLorentzVector&,double& alpha, double& beta, double& gamma); ///< Euler angles (alpha, beta, gamma), Y convention of TRotation, of the Wick helicity rotation
  static Matrix<2,2> RotationMatrix(const TLorentzRotation&, const TLorentzVector&); ///< spin-1/2 Wigner matrix of the Wick helicity rotation: \f$ U(l)|p,\lambda\rangle = \sum_{\lambda'} M_{\lambda\lambda'}|lp,\lambda'\rangle \f$
  static Matrix<2,2> RotationMatrix(const TLorentzQuaternion&, const TLorentzVector&); ///< as above, quaternion version

 private:
  void ToHelicitySpinor(const Polarization::State&, const double phi); ///< rotates the Pauli spinors of the TSpinor constructor to the helicity basis


}; // class THelicitySpinor

/*! \brief Helicity spinor of a negative-energy solution (antiparticle), see THelicitySpinor and TVSpinor */
class THelicityVSpinor : public TSpinor
{
 public:
  THelicityVSpinor(const TLorentzVector&, double mass, const Polarization::State&, Normalization=kDoubleMass); ///< constructor, see THelicitySpinor
  THelicityVSpinor(const FourVector<double>&, double mass, const Polarization::State&, Normalization=kDoubleMass); ///< constructor, see THelicitySpinor
  THelicityVSpinor(const THelicityVSpinor&); ///< copy constructor
  virtual ~THelicityVSpinor(); ///< destructor
  THelicityVSpinor& operator=(const THelicityVSpinor&); ///< assignment
  virtual THelicityVSpinor* Clone(const char* ="") const; ///< virtual copy constructor

 private:
  void ToHelicityVSpinor(const Polarization::State&, const double phi); ///< rotates the Pauli spinors to the helicity basis


}; // class THelicityVSpinor

//______________________________________________________________________________
Matrix<4,1> operator*(const Matrix<4,4>& lhs, const TSpinor& rhs); ///< 4x4 matrix times spinor
std::complex<double> operator*(const Matrix<1,4>& lhs, const TSpinor& rhs); ///< row (e.g. a bar spinor) times spinor
Matrix<4,1> operator*(const GammaStructure& lhs, const TSpinor& rhs); ///< gamma-matrix structure times spinor

/** @} */
#endif // TSPINOR_H
