/*! \file LeptonKinematics.hpp 
 * \brief Contains LeptonKinematics class declaration
 * \author Wim Cosyn
 * \date 12/12/2013
 * \addtogroup MePhys
 * @{
 */

#ifndef LEPTONKINEMATICS_H
#define LEPTONKINEMATICS_H

#include <Tools/Core/FourVector.hpp>

namespace nucqcd {


// Forward declarations
class Kinematics2to2;
class Kinematics2to3;

/*! \brief A Class that stores the lepton kinematics for a certain charged weak scattering process with incoming neutrino
 * In electroproduction reactions we separate the electron scattering
 * vertex from the virtual photon - deuteron system.
 * This class deals with the kinematics of the former. In conjunction
 * with a Kinematics2to3 object the kinematics are completely fixed by
 * the incoming electron's beam energy
 * or the outgoing electron's scattering angle
 * Objects should be created with the functions CreateWithBeamEnergy or CreateWithCosScatterAngle, 
 * these return a pointer to a newly created instance.  The user has to delete the Object after he is done with it!
 * Kinematics are solved when you pass a Kinematics2to2 or Kinematics2to3 object to the Getters.
 */

class LeptonKinematics
{
 public:
  LeptonKinematics(const LeptonKinematics&); /*!< copy constructor */
  LeptonKinematics& operator=(const LeptonKinematics&); /*!< assignment operator */
  virtual ~LeptonKinematics() {} /*!< Destructor */

  virtual LeptonKinematics* Clone() const { return new LeptonKinematics(*this); } /*!< ROOT Cloner */

  /// lepton type
  enum Lepton { electron = 0, ///< electron
		       muon = 1,     ///< muon
		       tau = 2       ///< tau
  };

  /*! indirect constructorcreator 
   * \param lep lepton type
   * \param beam [MeV] incoming beam energy */
  static LeptonKinematics* CreateWithBeamEnergy(Lepton lep, double beam); 
//   static LeptonKinematics* CreateWithEpsilon(double);
  /*! indirect constructorcreator 
   * \param lep lepton type
   * \param costheta_l [-1,1] cosine of the lepton scattering angle */
  static LeptonKinematics* CreateWithCosScatterAngle(Lepton lep, double costheta_l);

  void SetBeamEnergy(double beam ); /*!< \param beam [MeV] beam energy */
//   void SetEpsilon(double); 
  void SetCosScatterAngle(double costheta_l); /*!< \param costheta_l [-1,1] cosine of the lepton scattering angle */

  double GetBeamEnergy(const Kinematics2to2&) const; /*!< \return [MeV] energy of incoming beam */
//   double GetEpsilon(const Kinematics2to2&) const; 
  double GetCosScatterAngle(const Kinematics2to2&) const; /*!< \return costheta_l [-1,1] cosine of the lepton scattering angle */
//   double GetTan2HalfAngle(const Kinematics2to2&) const;
  double GetBeamEnergy(const Kinematics2to3&) const; ///< [MeV] energy of the incoming beam
//   double GetEpsilon(const Kinematics2to3&) const; /*!< \return costheta_l [-1,1] cosine of the lepton scattering angle */
  double GetCosScatterAngle(const Kinematics2to3&) const; ///< cosine of the lepton scattering angle, in [-1,1]
  double GetLeptonMass() const{return mass;} /*!< [MeV] returns lepton mass */
//   double GetTan2HalfAngle(const Kinematics2to3&) const;
  static const double masse; /*!< [MeV] electronmass */
  static const double massmu; /*!< [MeV] muon mass */
  static const double masstau; /*!< [MeV] tau mass */
  
  /*! constructs the four vectors for the incoming and outgoing lepton
   * \param kin kinematics of the hadron piece
   * \param[out] k_in [MeV] incoming lepton fourvector
   * \param[out] k_out [MeV] outgoing lepton fourvector
   */
  void GetLeptonVectors(const Kinematics2to2& kin, FourVector<double> &k_in, FourVector<double> &k_out); 
  
 protected:
  /*! basic constructor
   * \param type outgoing lepton type
   */
  LeptonKinematics(Lepton type); 
//   LeptonKinematics(); /*!< default constructor */
  /*! Solves the kinematics 
   * \param kin contains all gauge boson + A -> B+C kinematics */
  bool SolveKinematics(const Kinematics2to2& kin) const; 
  /*! Solves the kinematics 
   * \param kin contains all gauge boson + A -> B+C+D kinematics */
  bool SolveKinematics(const Kinematics2to3& kin) const;
  
 private:
  /*! Controls input variable */
  enum InputVariable { kBeamEnergy = 0,
		       kEpsilon = 1,
		       kScatterAngle = 2 };

  InputVariable  fInput;        /*!< specifies the input variable */
  mutable double fBeamEnergy;   /*!< energy incoming electron beam in [MeV] */
//   mutable double fEpsilon;      /*!< virtual photon's transverse linear polarization */
  mutable double fScatterAngle; /*!< cos of electron scattering angle in ]-1,1] */
//   mutable double ftan2HalfAngle; /*!< tan squared of half of the lepton scattering angle */
  double mass; /*!< [MeV] mass of outgoing lepton */
  

}; 

}  // namespace nucqcd

/** @} */
#endif 
