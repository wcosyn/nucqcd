/*! \file InclusiveCrossRes.hpp 
 * \brief Contains declaration of class InclusiveCrossRes, computes deuteron inclusive cross sections based on resonance model
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup DIS
 * @{
 */
#ifndef INCLUSIVECROSSRES_HPP
#define INCLUSIVECROSSRES_HPP

#include <string>
#include <memory>
#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/StructureFunctions/DeuteronStructure.hpp>
#include <Tools/WaveFunctions/Deuteron.hpp>

namespace nucqcd {



//forward declaration
class Resonance;

/*! \brief A class that computes inclusive deuteron cross sections based on a resonance model
 *
 * The final-state interactions in inclusive DIS on the deuteron are modelled with a set of
 * resonances (Resonance) that rescatter off the spectator nucleon, see W. Cosyn, W. Melnitchouk
 * and M. Sargsian, Phys. Rev. C 89, 014612 (2014).
 */

class InclusiveCrossRes{
  
public:
    /*! Constructor
   * \param elec Electron kinematics, see ElectronKinematics
     * \param proton photon interacting with proton [1] or neutron [0]
   * \param wavename Deuteron wave function name, see Deuteron
   * \param strucname which structure function parametrization ["CB"=Christy-Bosted (SLAC above W^2 = 9 GeV^2, see NuclStructure), "SLAC", "Alekhin", or any other name of NuclStructure: "CTEQ", "MSTW", "LHAPDF:<set>", combinations such as "CB+LHAPDF:CT18LO@4"], see NuclStructure
   * \param symm parameter that controls symmetric or not rescatterings in the FSI
   * \param offshell which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression
   * \param fixprop fixed masses in the propagators [1] or not [0]
   * \param kin_choice affects parts of phase space of spectator that are included, too high spectator momenta are excluded
   * because of energy conservation, provides smooth transition between resonance tresholds
   */  
  InclusiveCrossRes(bool proton, std::string strucname, std::string wavename, ElectronKinematics &elec, 
		    int symm, int offshell, bool fixprop, int kin_choice);
    /*! Calculates the plane-wave inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   * \return [] inclusive cross section without prefactors (structure functions times momentum distribution integrated)
   */
  double calc_F2Dinc(double Q2,double x);
  /*! Calculates the fsi inclusive cross section without any prefactors, just the deuteron structure functions time
   * momentum distribution
   * \param fsi1 [] inclusive fsi cross section without prefactors (structure functions times momentum distribution integrated)
   * calculated in one of two ways (see notes)
   * \param fsi2 [] other fsi formula
   * \param Q2 [MeV^2] Q^2 of virtual photon
   * \param x [] bjorken x
   */
  void calc_F2DincFSI(double &fsi1, double &fsi2, double Q2,double x);
  //void calc_F2DincFSI2(double &fsi1, double &fsi2, double Q2,double x);
  /*! add a Resonance instance to the vector containing all the resonances 
   * \param res Resonance instance you want to add
   */
  void addResonance(Resonance &res) {resonance_vec.push_back(res);} 
  std::vector<Resonance> & getResonance_vec() {return resonance_vec;} /*!< returns the vector with all the resonances */
  
private:
  double massi; /*!< [MeV] mass of nucleon that gets hit by photon */
  double massr; /*!< [MeV] mass of spectator nucleon */
  int symm; /*!< parameter that controls symmetric or asymmetric rescatterings in the FSI amplitude */
  bool fixprop; /*!< fixed masses in propagator [1] or not [0] */
  int t_choice;  /*!< affects how t is computed <BR> [0] = (pr-pr')^2 <BR> [1] = (pX-pX')^2 */
  int kin_choice; /*!< affects parts of phase space of spectator that are included,
      too high spectator momenta are excluded because of energy conservation, provides smooth transition between resonance tresholds */
  
  double przprime; /*!< [MeV] longitudinal spectator momentum after rescattering */
  double prz; /*!< [MeV] lingitudinal spectator momentum before rescattering */
  double Wxprime2; /*!< [MeV^2] invariant mass of the X that is paired with the spectator that gets integrated first */
  double otherWx2; /*!< [MeV^2] invariant mass of the other X that depends on the first ps integration */
  /*! which offshell parametrization do you want to use? <BR>
   * - 0: based on off-shell mass suppression (See M. Sargsian PRC82, 014612)
   * - 1: based on dipole FF suppression with cutoff lambda (See S. Jesschonek PRC78, 014007)
   * - 2: suppression with a off-shell beta parameter
   * - 3: no off-shell amplitude, fully suppressed
   * - 4: full off-shell amplitude, no suppression 
   */
  int offshellset;

  

  std::unique_ptr<Deuteron::Wavefunction> wf; /*!< pointer to deuteron wave funtion, see Deuteron */
  ElectronKinematics electron; /*!< electron kinematis, see ElectronKinematics */
  DeuteronStructure structure;  /*!< deuteron structure functions object, see DeuteronStructure */

  std::vector<Resonance> resonance_vec; /*!< vector that contains Resonance instances, used in the FSI amplitude */
  
  /*! integrand over the spectator momentum for the plane-wave F2 (calc_F2Dinc)
   * \param pr [MeV] spectator momentum \param[out] result integrand
   * \param ap x, Q^2 [MeV^2], and a pointer to the estimate for the inner cos theta integration */
  void int_pr(double pr, double *result, va_list ap);
  /*! integrand over the cosine of the spectator angle with q for the plane-wave F2
   * \param costheta cos of the angle between spectator momentum and q \param[out] result integrand
   * \param ap spectator momentum [MeV], x, Q^2 [MeV^2] */
  void int_costheta_incl(double costheta, double *result, va_list ap);
  /*! integrand over the spectator momentum for the FSI part of F2 (calc_F2DincFSI)
   * \param pr [MeV] spectator momentum \param[out] result integrand for the two FSI formulas of calc_F2DincFSI
   * \param ap x, Q^2 [MeV^2], the indices of the two resonances (size_t), and a pointer to the
   * estimate for the inner cos theta integration */
  void int_pr_fsi(double pr, double *result, va_list ap);
  /*! integrand over the cosine of the spectator angle with q for the FSI part of F2
   * \param costheta cos of the angle between spectator momentum and q \param[out] result integrand for the two FSI formulas of calc_F2DincFSI
   * \param ap spectator momentum [MeV], x, Q^2 [MeV^2], the indices of the two resonances (size_t) */
  void int_costheta_incl_fsi(double costheta, double *result, va_list ap);
  /*! integrand over the transverse momentum of the rescattering loop
   * \param qt [MeV] transverse momentum transferred in the rescattering \param[out] results integrand for the two FSI formulas of calc_F2DincFSI
   * \param ap the two kinematics objects (Kinematics2to2*), the deuteron wave function components
   * (complex<double>*), the indices of the two resonances (size_t), and a pointer to the estimate for
   * the inner azimuthal integration */
  void int_qt(double qt, double *results, va_list ap);
  /*! integrand over the azimuthal angle of the transverse momentum of the rescattering loop
   * \param qphi [rad] azimuthal angle \param[out] results integrand for the two FSI formulas of calc_F2DincFSI
   * \param ap qt [MeV], the two kinematics objects (Kinematics2to2*), the deuteron wave function
   * components (complex<double>*), the indices of the two resonances (size_t) */
  void int_qphi(double qphi, double *results, va_list ap); 
  /*! recursive method to find the pole in the fsi integration, longitudinal part,
   * also determines intermediate mass
   * \param pt [MeV] final transverse spectator momentum
   * \param Er [MeV] final spectator on-shell energy
   * \param pkin kinematics object containing the gamma+D->X+N kinematics <BR>
   * in Kinematics2to2 language: deuteron is N, nucleon is Kaon, X is hyperon
   * \param first which resonance comes first in the amplitude
   * \param res1 index of the initial resonance
   * \param res2 index of the resonance after rescattering
   */  
  void get_prz(double pt, double Er, Kinematics2to2 *pkin, int first, size_t res1, size_t res2);
 /*! gives you the scatter amplitude of the final-state interaction
   * \param t [MeV^2] momentum transfer squared
   * \param Q2 [MeV^2] Q^2 of the virtual photon
   * \param res1 index of the initial resonance
   * \param res2 index of the resonance after rescattering
   * \return \f$ \sigma_{tot} (I+\epsilon) e^{\beta t/2} multiplied by the relevant coefficients of the resonances \f$
   */  
  std::complex<double> scatter(double t, double Q2, size_t res1, size_t res2);
};


}  // namespace nucqcd

/*! @} */
#endif