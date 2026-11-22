/*! \defgroup MePhys libMePhys: Library that contains all classes for general medium physics purposes
 * \author Wim Cosyn, Pieter Vancraeyveld
 * \date 29/08/2012
 * \brief This code implements classes to do a lot of general medium physics stuff
 * 
 * \details 
 * 
 * I jacked and hacked most of the stuff included here from the strangecalc libraries.
 * Generally I did not doxygen-ize these classes.  The relevant API can be found at 
 * http://rprmodel.ugent.be/api/. <BR><BR>
 * 
 * I also added some new classes: <BR>
 * - A class that contains nucleon structure functions [NuclStructure] <BR>
 * - A class that contains deuteron structure functions [DeuteronStructure] <BR>
 * - A leptonkinematics class for neutrino scattering [TLeptonKinematics] <BR>
 * - A class for he3 wave functions [He3wf] <BR>
 * -A class that has lightcone kinematics for a two body scattering reaction  [LightConeKin2to2] <BR><BR>
 * -A class that has lightcone kinematics for a two to three body scattering reaction  [LightConeKin2to3] <BR><BR>
 * 
 * Utilfunctions.hpp contains a few helper functions (mostly math related)<BR>
 * --------------------------------------<BR>
 * 
 * Some changes I did in Pieter's files and classes: <BR>
  * -moved from kinematics2to3.h and.cpp to fourvector.h: <BR>
  * // Global functions for conversion <BR>
  * // from ROOT 4vectors to strangecalc 4vectors.<BR>
  * template<typename T> <BR>
  * FourVector<T> operator*(const TLorentzRotation&, const FourVector<T>&); <BR>
  * <BR>
  * FourVector<double> ToFourVector(const TLorentzVector&); <BR>
  * TLorentzVector ToLorentzVector(const FourVector<double>&); <BR>
  * <BR>
  * -moved definitions of ToFourVector and ToLorentzVector from kinematics2to3.cpp to a new file FourVector.cpp <BR>
  * <BR>
  * -commented out include "2to3kinematics.h" in TDeuteron.cpp <BR>
  * <BR>
  * -commented out all hash functions <BR>
  * -commented out numtoa.h includes <BR>
  * <BR>
  * -moved all the necessary wrapper classes in the dir with deuteron code, changed include <> to include "" <BR>
  * <BR>
  * -commented out all includes of Structures.h as only defines of constants were needed <BR>
  * <BR>
  * -included needed constants (nucqcd::strange_underflow, nucqcd::mass_p, nucqcd::mass_n) to my file constants.hpp (some  double definitions now of course) <BR>
  * <BR>
  * -replaced all includes of Structures.h with include "constants.hpp" (was not even needed in TYukawaPWF.cpp) <BR>
  * <BR>
  * -changed TYukawaPWF to also include parmetrization with non-equidistant Mass(i) (AV18 f.i.) <BR>
  *  -added double *fM member <BR>
  * -added fM(NULL) to initializer list of constructor<BR>
  * -added new constructor with array of m values<BR>
  * -check in constructor to delete fM when necessary<BR>
  * -modified Mass(int i) for this case<BR>
  * <BR>
  *-added AV18 and Av18b Deuteron wf to TDeuteron.h & TDeuteron.cpp <BR>
  *<BR>
  *-added coordinate state wavefunction functions to TDeuteron<BR>
  *-some more vector types to DeuteronRState and DeuteronPState<BR>
  *-added off-shell momentum wave function parts to TWaveFunctionImplementation and derived classes
  * (GetUpoff etc, DeuteronPStateOff functions...)<BR>
  *<BR>
  *-commented out hash functions in TKinematics.cpp<BR>

  *-added hypron kin vars in TKinematics.cpp: updated constructors, getters, updatekinematics<BR>
  *-pklab: ++energycount changed to ++anglecount<BR>
  *-l1250 added pklab kinematics<BR>
  *<BR>
  *-added tan2 to Electronkin<BR>
  *<BR>
  *-electronkinematics changed all 2to3 to 2to2<BR>
  *<BR>
  *-added TKinematics2to3 to Tools, commented out functions that I moved previously...<BR>
  *<BR>
  *-added masses to constructors to Tkinematics2to3(withlabangles)<BR>
  *changed setisospin function for masses functionality<BR>
 *
 */



/*! \file NuclStructure.hpp 
 * \brief Contains declaration of class NuclStructure: nucleon structure functions at fixed (possibly off-shell) kinematics
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef NUCLSTRUCTURE_HPP
#define NUCLSTRUCTURE_HPP

#include "NucleonStructure.hpp"

#include <cmath>
#include <string>

/*! \brief Nucleon structure functions at fixed kinematics (Q^2, x and W^2, which may be off-shell)
 *
 * A thin wrapper around NucleonStructure: it stores the kinematics and the parametrization, and
 * evaluates them with NucleonStructure::getF_xQW, so the parametrizations are implemented once.
 * x and W^2 can be independent (off-shell nucleons in DeuteronStructure): Christy-Bosted and SLAC
 * use W^2, the others x. Until 2026-11 this class had its own copy of the parametrizations; F1 for
 * SLAC and CTEQ was then infinite, it is now the value of NucleonStructure (from F2 with R1998).
 */
class NuclStructure{
  
public:
  /*! Constructor
   * \param proton proton [1] or neutron [0]
   * \param var1 first variable, see below
   * \param var2 second variable, see below
   * \param switchvar Which input variables? Possibilities: <BR>
   * [0]: var1 is Q2 [MeV^2] and var2 is Bjorken x [] <BR>
   * [1]: var1 is Q2 [MeV^2] and var2 is W^2 [MeV^2] <BR>
   * [2]: var1 is Bjorken x [] and var2 is W^2 [MeV^2] <BR>
   * (the third variable follows from the on-shell relation W^2 = Q^2(1/x-1) + M^2)
   * \param name Name of the parametrization, see NucleonStructure: "CB" (Christy & Bosted; above
   * W^2 = 9 GeV^2 it switches to "SLAC", also in getName), "SLAC", "Alekhin", "CTEQ", "HMRS", "MSTW"
   */
  NuclStructure(const bool proton, double var1, double var2, int switchvar, std::string name);
  /*! Constructor with independent kinematics (off-shell)
   * \param proton proton [1] or neutron [0]
   * \param Q2in [MeV^2] Q^2 of the virtual photon
   * \param xin [] Bjorken x
   * \param Wsqin [MeV^2] invariant mass squared W^2
   * \param name Name of the parametrization, see the other constructor
   */
  NuclStructure(const bool proton, double Q2in, double xin, double Wsqin, std::string name);

  /*! Returns the F1 and F2 structure functions */
  void getF(double &F1, double &F2);
  double getF1(); /*!< returns the F1 structure function */
  double getF2(); /*!< returns the F2 structure function */
  const std::string getName() const{return name;} /*!< returns the name of the chosen parametrization */
  double getW(){return sqrt(Wsq);} /*!< [MeV] returns the invariant mass of X*/
  bool getProton() const{return proton;} ///< structure functions of the proton (true) or the neutron (false)
  
  /*! g1 from the grsv2000 parametrization, see NucleonStructure::getG1_grsv2000
   * \param proton 1=proton 0=neutron \param x Bjorken x \param Q2 [MeV^2] \return g1 */
  static double getG1_grsv2000(const bool proton, double x, double Q2){ return NucleonStructure::getG1_grsv2000(proton,x,Q2); }
  /*! R = sigma_L/sigma_T, 1990 SLAC fit, see NucleonStructure::getr1990
   * \param x Bjorken x \param Q2 [MeV^2] \return R */
  static double getr1990(double x, double Q2){ return NucleonStructure::getr1990(x,Q2); }
  /*! R = sigma_L/sigma_T, 1998 SLAC fit, see NucleonStructure::getr1998
   * \param x Bjorken x \param Q2 [MeV^2] \return R */
  static double getr1998(double x, double Q2){ return NucleonStructure::getr1998(x,Q2); }
  /*! g1+g2 from grsv2000 with the Wandzura-Wilczek relation, see NucleonStructure::getG1plusg2_grsv2000
   * \param proton 1=proton 0=neutron \param x Bjorken x \param Q2 [MeV^2] \return g1+g2 */
  static double getG1plusg2_grsv2000(const bool proton, double x, double Q2){ return NucleonStructure::getG1plusg2_grsv2000(proton,x,Q2); }
  
private:
  std::string name; /*!< name of the parametrization */
  bool proton; /*!< proton [1] or neutron [0] structure functions */
  double mass; /*!< [MeV] nucleon mass */
  double x; /*!< Bjorken x */
  double Q2; /*!< [MeV^2] Q^2 of the virtual photon */
  double Wsq; /*!< [MeV^2] invariant mass squared W^2 */
  NucleonStructure structure; /*!< the parametrization */
};
/*! @} */
#endif
