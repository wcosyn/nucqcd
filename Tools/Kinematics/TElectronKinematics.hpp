/*! \file TElectronKinematics.hpp
 * \brief Kinematics of the scattered electron in electroproduction
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on March 24, 2010
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TELECTRONKINEMATICS_H
#define TELECTRONKINEMATICS_H

#include <Tools/Core/FourVector.hpp>
// Forward declarations
class TKinematics2to2;
class TKinematics2to3;

/*! \brief Electron kinematics for electroproduction: fixed by one of beam energy, virtual-photon
 * polarization epsilon, or cosine of the electron scattering angle; the other two follow from Q^2 and
 * the energy transfer of the hadronic kinematics (TKinematics2to2 or TKinematics2to3, lab frame)
 *
 * Create it with one of the named constructors; the caller owns the object.
 */
class TElectronKinematics
{
 public:
  TElectronKinematics(const TElectronKinematics&); ///< copy constructor
  TElectronKinematics& operator=(const TElectronKinematics&); ///< assignment
  virtual ~TElectronKinematics() {} // Destructor

  virtual TElectronKinematics* Clone(const char * ="") const { return new TElectronKinematics(*this); } ///< virtual copy constructor

  static TElectronKinematics* CreateWithBeamEnergy(double); ///< fixes the beam energy [MeV]; the caller owns the object
  static TElectronKinematics* CreateWithEpsilon(double); ///< fixes the transverse linear polarization epsilon of the virtual photon; the caller owns the object
  static TElectronKinematics* CreateWithCosScatterAngle(double); ///< fixes the cosine of the electron scattering angle, in ]-1,1]; the caller owns the object

  void SetBeamEnergy(double); ///< fixes the beam energy [MeV]
  void SetEpsilon(double); ///< fixes epsilon
  void SetCosScatterAngle(double); ///< fixes the cosine of the electron scattering angle, in ]-1,1]

  double GetBeamEnergy(const TKinematics2to2&) const; ///< [MeV] beam energy for the hadronic kinematics
  double GetEpsilon(const TKinematics2to2&) const; ///< transverse linear polarization epsilon of the virtual photon for the hadronic kinematics
  double GetCosScatterAngle(const TKinematics2to2&) const; ///< cosine of the electron scattering angle for the hadronic kinematics
  double GetTan2HalfAngle(const TKinematics2to2&) const; ///< tan^2 of half the electron scattering angle for the hadronic kinematics
  double GetBeamEnergy(const double Q2, const double nu) const; ///< [MeV] beam energy for Q^2 [MeV^2] and energy transfer nu [MeV]
  double GetEpsilon(const double Q2, const double nu) const; ///< epsilon for Q^2 [MeV^2] and nu [MeV]
  double GetCosScatterAngle(const double Q2, const double nu) const; ///< cosine of the electron scattering angle for Q^2 [MeV^2] and nu [MeV]
  double GetTan2HalfAngle(const double Q2, const double nu) const; ///< tan^2 of half the electron scattering angle for Q^2 [MeV^2] and nu [MeV]
  double GetBeamEnergy(const TKinematics2to3&) const; ///< [MeV] beam energy for the hadronic kinematics
  double GetEpsilon(const TKinematics2to3&) const; ///< transverse linear polarization epsilon of the virtual photon for the hadronic kinematics
  double GetCosScatterAngle(const TKinematics2to3&) const; ///< cosine of the electron scattering angle for the hadronic kinematics
  double GetTan2HalfAngle(const TKinematics2to3&) const; ///< tan^2 of half the electron scattering angle for the hadronic kinematics

  /*! constructs the four vectors for the incoming and outgoing lepton
   * \param kin kinematics of the hadron piece
   * \param[out] k_in [MeV] incoming lepton fourvector
   * \param[out] k_out [MeV] outgoing lepton fourvector
   */
  void GetLeptonVectors(const TKinematics2to2& kin, FourVector<double> &k_in, FourVector<double> &k_out); 
  double GetLeptonMass() const{return 0.511;} /*!< [MeV] returns lepton mass */
  
  
 protected:
  TElectronKinematics(); ///< protected: use the named constructors
  bool SolveKinematics(const TKinematics2to2&) const; ///< solves the electron kinematics; false when unphysical
  bool SolveKinematics(const double Q2, const double nu) const; ///< solves the electron kinematics for Q^2 [MeV^2] and nu [MeV]; false when unphysical
  bool SolveKinematics(const TKinematics2to3&) const; ///< solves the electron kinematics; false when unphysical
  
 private:
  /// which variable fixes the electron kinematics
  enum InputVariable { kBeamEnergy = 0,  ///< beam energy
		       kEpsilon = 1,     ///< epsilon
		       kScatterAngle = 2 ///< scattering angle
  };

  /// the input variable
  InputVariable  fInput;        // specifies the input variable
  /// [MeV] beam energy
  mutable double fBeamEnergy;   // energy incoming electron beam in [MeV]
  /// transverse linear polarization epsilon of the virtual photon
  mutable double fEpsilon;      // virtual photon's transverse linear polarization
  /// cosine of the electron scattering angle, in ]-1,1]
  mutable double fScatterAngle; // cos of electron scattering angle in ]-1,1]
  /// tan^2 of half the electron scattering angle
  mutable double ftan2HalfAngle;


}; // class TElectronKinematics

/** @} */
#endif // TELECTRONKINEMATICS_H
