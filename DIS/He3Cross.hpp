/*! \file He3Cross.hpp 
 * \brief Contains declaration of class He3Cross, computes He3 DIS cross sections with two tagged spectator nucleons
 * \author Wim Cosyn
 * \date 26/09/2014
 * 
 * \addtogroup DIS
 * @{
 */

#ifndef HE3CROSS_HPP
#define HE3CROSS_HPP

#include <Tools/WaveFunctions/He3wf.hpp>
#include <Tools/StructureFunctions/DeuteronStructure.hpp>
#include <Tools/Kinematics/LightConeKin2to3.hpp>

#include <string>

namespace nucqcd {


/*! \brief A class that computes tagged spectator he3 cross sections
 *
 * Plane-wave (virtual nucleon) DIS on He3 with both spectator nucleons detected, following the LDRD note
 * "Tagged Deep Inelastic Scattering for general nucleus" (he3.tex):
 * \f[ \frac{d\sigma}{dx_A dQ^2 d\phi_{e'}} = \frac{y_A^2\alpha^2}{Q^6}\, m_N L_{\mu\nu}W^{\mu\nu}_N\, S(p_{s_1},p_{s_2})\,
 *     [(A-1)!\,\Theta]\, d^3p_{s_1} d^3p_{s_2}. \f]
 * The two spectators are not identified (for a struck proton: one proton and one neutron, either at
 * either momentum), so they are treated like the identical pp pair of a struck neutron:
 * (A-1)! Theta = 1, and the cross section is integrated over the full (ordered) spectator phase space.
 *
 * Normalization of S: the note normalizes the wave function to 1 for each struck nucleon (its
 * optical-theorem check sums sigma_hN over the nucleons N), whereas the He3 wave function of He3wf is
 * antisymmetric, including isospin, and normalized to 1 summed over all isospin arrangements, so that
 * each arrangement carries 1/A of it. Hence S = A sum |Psi_He3wf|^2, and the momentum distribution
 * integrates to Z = 2 (struck proton) and N = 1 (struck neutron). Until 2026-11 the factor A = 3 was
 * missing (the tagged cross section was 3 times too small), and before that S also missed part of the D
 * and P waves (see He3wf).
 */
class He3Cross{
public:
    /*! Constructor
   * \param wfname He3 wave function name, see He3wf class constructor for possibilities ["AV18", etc.]
   * \param inputdir share dir [full path please, else fortran code breaks for some reason]
   * \param strucname which structure function parametrization ["CB"=Christy-Bosted (SLAC above W^2 = 9 GeV^2), "SLAC", "Alekhin", or any other name of NuclStructure: "CTEQ", "MSTW", "LHAPDF:<set>", combinations such as "CB+LHAPDF:CT18LO@4"], see NuclStructure
   * \param proton photon interacts with proton [1] or neutron [0]
   */
  He3Cross(std::string wfname, std::string inputdir, std::string strucname, bool proton);
  He3Cross(const He3Cross&); /*!< Copy Constructor */
  He3Cross& operator=(const He3Cross&); /*!< assignment operator */
  ~He3Cross(); /*!<Destructor */
  
  /*! He3 momentum distribution of the two (unidentified) spectators in the He3 rest frame
   * \param sp1 [MeV] momentum of spectator 1 \param sp2 [MeV] momentum of spectator 2
   * \return [MeV^-6] S(sp1,sp2) = A sum over all spins (m_A averaged) and the spectator isospins of
   * |Psi|^2, with the struck nucleon (momentum -sp1-sp2) the proton or neutron of the constructor;
   * \f$ \int d^3p_1 d^3p_2\, S = Z = 2 \f$ (struck proton) or N = 1 (struck neutron), see the class description */
  double getMomentumDistribution(const nucqcd::Vector3 &sp1, const nucqcd::Vector3 &sp2);
  /*! Calculate He3 momentum distribution for a certain kinematics of the spectators
   * \param[in] kin LightconeKin2to3 object that has all the kinematics
   * \return [MeV-4] getMomentumDistribution at the spectator momenta (shifted by a third of the He3
   * momentum), times E_sp1_lab*E_sp2_lab
   */  
  double getDensity(LightConeKin2to3 &kin);
    /*! get the average cross section
   * \param kin kinematics object containing the gamma+He3->X+N1+N2 kinematics <BR>
   * \return [nb/GeV^6] tagged spectators cross section \f$ \frac{d\sigma}{dx_A dQ^2 \frac{d^3p_{s_1}}{E_{s_1}}\frac{d^3p_{s_2}}{E_{s_2}}} \f$ 
   * NOT averaged over phi, for unidentified spectators (see the class description and the he3 note)
   */
  double getCross (LightConeKin2to3 &kin);
  
private:
  He3wf wavefunction; /*!< He3 wave function, see He3wf */
  DeuteronStructure strucfunc; /*!< structure functions of the struck nucleon, see DeuteronStructure */
  bool proton; /*!< photon interacts with a proton [1] or a neutron [0] */

};

}  // namespace nucqcd

/*! @} */
#endif  
