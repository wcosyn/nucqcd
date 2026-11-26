/*! \mainpage Glauber ISI/FSI RMSGA C++ Project
 * \defgroup Glauber libGlauber: Library that contains all Glauber related classes
 * \author Wim Cosyn
 * \date 16/08/2011
 * \brief This module implements classes for the Glauber RMSGA formalism, including CT and SRC effects.
 * 
 * \details 
 * - It contains classes for a mean field nucleus, a mean field nucleus with densities (used in thickness calculations). <BR>
 * [MeanFieldNucleus, MeanFieldNucleusThick] <BR>
 * - A class for correlated FSI calculations containing all the needed functions and a grid for the gamma functions. <BR>
 * [FsiCorrelator]<BR>
 * - Two abstract classes: one for a general FSI grid (with the interface for all the variants of a grid: FSI, SRC, decay, CT, see GridVariant),
 * one that adds the color transparency grid. (Until 2026-11 there were seven, one layer for every variant.) <BR>
 * [AbstractFsiGrid, AbstractFsiCTGrid]<BR>
 * - Three glauber classes : one without thickness, two with (also adds SRC to the ISI/FSI, one adds decay properties for ejected particles), the two
 * with thickness sharing GlauberThickGridCommon. <BR>
 * [GlauberGrid, GlauberGridThick, GlauberDecayGridThick, GlauberThickGridCommon] <BR>
 * - A special class for the glauber grid of one particle, exploiting the symmetry along the z-axis (no phi dependence).  <BR>
 * [OneGlauberGrid] <BR>
 * - A class for spinors that represent the MF wave function of a nucleon in a certain nucleus at a coordinate. <BR>
 * [MeanFieldSpinor] <BR>
 * - A class that includes all properties of certain "ejected" particles (scattering parameters, kinetics, etc). <BR>
 * [FastParticle] <BR>
 * - A class that stores grids for distorted momentum distributions <BR>
 * [DistMomDistrGrid] <BR>
 * - A class that computes nucleon form factors, dirac couplings, with medium modifications (CQM or QSM) if desired. <BR>
 * [NucleonEMOperator] <BR>
 * - A class that computes weak nucleon form factors, dirac couplings, with medium modifications (CQM or QSM) if desired. <BR>
 * [NucleonWeakEMOperator] <BR>
 * -A class that implements the ROMEA model for FSI based on the optical potentials from Cooper & Hana (Phys. Rev. C 80, 034605) <BR>
 * [ROMEAGrid] <BR>
*/



/*! \file AbstractFsiGrid.hpp 
 * \brief Contains declaration of abstract class AbstractFsiGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * 
 * \addtogroup Glauber
 * @{
 */

#ifndef ABSTRACTFSIGRID_H
#define ABSTRACTFSIGRID_H

#include <string>
#include <complex>
#include <vector>
#include <Tools/Core/Vector3.hpp>

#include <Glauber/FastParticle.hpp>
#include <Tools/Core/Grid3D.hpp>
#include <Tools/Core/IntegratorFlag.hpp>
#include <Glauber/MeanFieldNucleus.hpp>

namespace nucqcd {


/*! The variants of a Glauber grid: which of the factors short-range correlations (SRC), decay of the ejected particle (decay)
 * and color transparency (CT) are included, besides the FSI factor (Fsi). The value of an enumerator is the grid number
 * of GlauberDecayGridThick, the class with all of them (getFsiGridN_interp): bit 0 SRC, bit 1 decay, bit 2 CT. See
 * AbstractFsiGrid::getGridNumber for which of them a class has. */
enum class GridVariant {
  Fsi = 0,         /*!< FSI */
  Src = 1,         /*!< FSI+SRC */
  Decay = 2,       /*!< FSI+decay */
  SrcDecay = 3,    /*!< FSI+SRC+decay */
  Ct = 4,          /*!< FSI+CT */
  SrcCt = 5,       /*!< FSI+SRC+CT */
  CtDecay = 6,     /*!< FSI+CT+decay */
  SrcCtDecay = 7   /*!< FSI+SRC+CT+decay */
};

/*! \brief An abstract class for a ISI/FSI grid
 * 
 * Contains all necessary functions to operate a general fsi grid. Has a pointer to a nucleus for which the fsi grid is.  
 * A vector contains all the particles that undergo ISI or FSI.  Has interpolation and print functions for the grid.
 * 
 * Typically an object that is inherited from this class is operated as follows.<BR>
 * 1. Initialize object with constructor AbstractFsiGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus with addKnockout() <BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever...
 *   (use getFsiGridFull_interp3() for instance) <BR>
 * 
 * 
 */
class AbstractFsiGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in cos theta
   * \param phi_grid gridsizein phi
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param prec the precision you want in the integrations
   * \param integrator which integrator: numint::IntegrationMethod::Romberg (numint::cube_romb) or numint::IntegrationMethod::Adaptive (numint::cube_adaptive); the integers 1 and 2 convert with nucqcd::integration_method()
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   * The grids are cached in dir/grids/, or in $NUCQCD_GRID_DIR if that is set (see GridCache.hpp).
   */
  AbstractFsiGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus *pnucl, 
		  double prec, IntegrationMethod integrator, std::string dir);
  virtual ~AbstractFsiGrid(); /*!< Destructor */
  // Not copyable. The Glauber grids own their arrays as Grid3D (Tools/Core/Grid3D.hpp) since 2026-11, so copying them would no
  // longer free anything twice, but nothing copies them and their other members (FsiCorrelator, ...) were not checked;
  // DistMomDistrGrid still has raw arrays (docs/ROADMAP.md, step 6, "Glauber grid hierarchy").
  AbstractFsiGrid(const AbstractFsiGrid&) = delete;
  AbstractFsiGrid& operator=(const AbstractFsiGrid&) = delete;
  
  const std::string getFsi_Filename() const{ return fsi_filename;} /*!< returns filename for regular fsi grid */
  const std::string getDir() const{ return dir;}/*!< returns the grid cache directory (with trailing /) */ 
  int getRgrid() const { return rgrid;}/*!<  returns the gridsize in r*/
  int getCthgrid() const { return cthgrid;}/*!<  returns the gridsize in theta*/
  int getPhigrid() const {return phigrid;}/*!<  returns the gridsize in phi*/
  double getInvRstep() const {return invrstep;}/*!<  returns the inverse of the stepsize of the grid in r*/
  double getInvCthstep() const {return invcthstep;}/*!<  returns the inverse of the stepsize of the grid in theta*/
  double getInvPhistep() const {return invphistep;}/*!<  returns the inverse of the stepsize of the grid in phi*/
  double getS_interp() const {return s_interp;}/*!<  returns s, measure where the point for interpolation is between two indices in r*/
  double getT_interp() const {return t_interp;}/*!<  returns t, measure where the point for interpolation is between two indices in theta*/
  double getU_interp() const {return u_interp;}/*!<  returns u, measure where the point for interpolation is between two indices in phi*/
  double getComp_s_interp() const {return comp_s_interp;}/*!<  returns 1-s*/
  double getComp_t_interp() const {return comp_t_interp;}/*!<  returns 1-t*/
  double getComp_u_interp() const {return comp_u_interp;}/*!<  returns 1-u*/
  int getRindex() const {return rindex;}/*!<  returns the index for r used in the 3d interpolation*/
  int getCthindex() const {return cthindex;}/*!<  returns the index for theta used in the 3d interpolation*/
  int getPhiindex() const {return phiindex;}/*!<  returns the index for phi used in the 3d interpolation */
  MeanFieldNucleus *getPnucleus() const {return pnucleus;}/*!<  returns the pointer to the MeanFieldNucleus instance */
  std::vector<FastParticle> & getParticles() {return particles;}/*!<  returns the vector with all particles subject to isi/fsi */
  std::vector<int> getKnockoutLevels() const {return knockoutlevels;}/*!< returns the vector with all the shell index values of knocked out nucleons */
  std::vector<int> getKnockoutM() const {return knockoutm;}/*!<  returns the vector with all the m values of knocked out nucleons */
  bool getAllinplane() const {return allinplane;}/*!<  returns 1 if all the isi/fsi particles lie in the x-z plane */
  double getPrec() const {return prec;}/*!< returns the integration precision */
  IntegrationMethod getIntegrator() const{return integrator;} /*!< returns the integrator type */
  
  /*! add particle subject to isi/fsi, checks to see that the vector contains max one incoming particle.
   * The grid stores and uses the particle moved to the centre of its grid-cache bin
   * (FastParticle::gridBinCentre), so getParticles() returns the binned particles.
   *\param[in] newparticle instance of FastParticle you want to add*/
  virtual void addParticle(FastParticle &newparticle); 
  void clearParticles(); /*!< Clear contents of isi/fsi particles vector*/
  void printParticles() const; /*!< Print contents of isi/fsi particles vector*/
  /*! add quantum numbers of particle that is knocked out from nucleus
   *\param level shellindex level
   *\param m mj quantum number (times two!!)  [-2j,2j]   
   */
  void addKnockout(const int level, const int m); 
  void clearKnockout(); /*!< clear the vectors of the knocked out particles*/
  void printKnockout() const; /*!< print info about the particles that are knocked out from the nucleus*/

  int getTotalProtonOut() const {return totalprotonout;} /*!< get total protons subject to fsi/isi (isi counts negative)*/
  int getTotalNeutronOut() const {return totalneutronout;} /*!< get total neutrons subject to fsi/isi (isi counts negative)*/
  int getProtonKnockout() const {return protonknockout;} /*!< get total protons that are knocked out from nucleus*/
  int getNeutronKnockout() const {return neutronknockout;} /*!< get total neutrons that are knocked out from nucleus*/

  int getNumber_of_grids() const{return number_of_grids;}/*!< get total number of grids the instance contains*/

  virtual void printFsi_grid()=0;  /*!< Prints the FSI grid for a certain situation, pure virtual function!! */
  virtual void print_grid(int gridindex)=0;  /*!< Prints the FSI grid for a certain situation, pure virtual function!! */

  void setRinterp(const double r); /*!< sets the rindex variable used in the 3d interpolation */
  void setCthinterp(const double costheta); /*!< sets the thindex variable used in the 3d interpolation */
  void setPhiinterp(const double phi); /*!< sets the phiindex variable used in the 3d interpolation */
  
  /*! sets the interpolation indices and weights in r, cos(theta) and phi
   * \param r [fm] radial coordinate \param theta cos of the polar angle (despite the name)
   * \param phi [rad] azimuthal angle; with all particles in the x-z plane, phi > pi is mirrored to 2pi-phi */
  void setinterp(const double r, const double theta, const double phi);
  /*! sets the interpolation indices and weights in cos(theta) and phi (r set before)
   * \param costheta cos of the polar angle \param phi [rad] azimuthal angle */
  void setinterp(const double costheta, const double phi);
  /*! sets the interpolation index and weight in phi (r and cos(theta) set before) \param phi [rad] azimuthal angle */
  void setinterp(const double phi);
  
  /*! The grid number (as in getFsiGridN_interp) of a variant of this class, -1 if it has no such grid. The default is a
   * class with only the plain FSI grid. Overridden by the classes with more grids.
   * \param variant the variant */
  virtual int getGridNumber(GridVariant variant) const { return variant == GridVariant::Fsi ? 0 : -1; }
  /*! \return true if this class has a grid for the variant, see getGridNumber */
  bool hasVariant(GridVariant variant) const { return getGridNumber(variant) >= 0; }

  // The named interpolation functions: one for every variant, with the point set before (no argument), with the point
  // as a vector, as (r,costheta,phi), as (costheta,phi) (r set before) and as phi (r and costheta set before). They are
  // the interpolation of getFsiGridN_interp with the grid number of the variant.
  /** @name The FSI grid */
  /** @{ */
  /*! FSI grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::Fsi)); }
  std::complex<double> getFsiGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::Fsi), rvec); } /*!< FSI grid at the point rvec */
  std::complex<double> getFsiGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::Fsi), r, costheta, phi); } /*!< FSI grid at (r,costheta,phi) */
  std::complex<double> getFsiGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::Fsi), costheta, phi); } /*!< FSI grid at (costheta,phi), r set before */
  std::complex<double> getFsiGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::Fsi), phi); } /*!< FSI grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+CT grid */
  /** @{ */
  /*! FSI+CT grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiCtGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::Ct)); }
  std::complex<double> getFsiCtGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::Ct), rvec); } /*!< FSI+CT grid at the point rvec */
  std::complex<double> getFsiCtGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::Ct), r, costheta, phi); } /*!< FSI+CT grid at (r,costheta,phi) */
  std::complex<double> getFsiCtGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::Ct), costheta, phi); } /*!< FSI+CT grid at (costheta,phi), r set before */
  std::complex<double> getFsiCtGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::Ct), phi); } /*!< FSI+CT grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+SRC grid */
  /** @{ */
  /*! FSI+SRC grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiSrcGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::Src)); }
  std::complex<double> getFsiSrcGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::Src), rvec); } /*!< FSI+SRC grid at the point rvec */
  std::complex<double> getFsiSrcGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::Src), r, costheta, phi); } /*!< FSI+SRC grid at (r,costheta,phi) */
  std::complex<double> getFsiSrcGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::Src), costheta, phi); } /*!< FSI+SRC grid at (costheta,phi), r set before */
  std::complex<double> getFsiSrcGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::Src), phi); } /*!< FSI+SRC grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+SRC+CT grid */
  /** @{ */
  /*! FSI+SRC+CT grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiSrcCtGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::SrcCt)); }
  std::complex<double> getFsiSrcCtGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::SrcCt), rvec); } /*!< FSI+SRC+CT grid at the point rvec */
  std::complex<double> getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::SrcCt), r, costheta, phi); } /*!< FSI+SRC+CT grid at (r,costheta,phi) */
  std::complex<double> getFsiSrcCtGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::SrcCt), costheta, phi); } /*!< FSI+SRC+CT grid at (costheta,phi), r set before */
  std::complex<double> getFsiSrcCtGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::SrcCt), phi); } /*!< FSI+SRC+CT grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+decay grid */
  /** @{ */
  /*! FSI+decay grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiDecayGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::Decay)); }
  std::complex<double> getFsiDecayGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::Decay), rvec); } /*!< FSI+decay grid at the point rvec */
  std::complex<double> getFsiDecayGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::Decay), r, costheta, phi); } /*!< FSI+decay grid at (r,costheta,phi) */
  std::complex<double> getFsiDecayGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::Decay), costheta, phi); } /*!< FSI+decay grid at (costheta,phi), r set before */
  std::complex<double> getFsiDecayGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::Decay), phi); } /*!< FSI+decay grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+CT+decay grid */
  /** @{ */
  /*! FSI+CT+decay grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiCtDecayGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::CtDecay)); }
  std::complex<double> getFsiCtDecayGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::CtDecay), rvec); } /*!< FSI+CT+decay grid at the point rvec */
  std::complex<double> getFsiCtDecayGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::CtDecay), r, costheta, phi); } /*!< FSI+CT+decay grid at (r,costheta,phi) */
  std::complex<double> getFsiCtDecayGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::CtDecay), costheta, phi); } /*!< FSI+CT+decay grid at (costheta,phi), r set before */
  std::complex<double> getFsiCtDecayGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::CtDecay), phi); } /*!< FSI+CT+decay grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+SRC+decay grid */
  /** @{ */
  /*! FSI+SRC+decay grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiSrcDecayGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::SrcDecay)); }
  std::complex<double> getFsiSrcDecayGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::SrcDecay), rvec); } /*!< FSI+SRC+decay grid at the point rvec */
  std::complex<double> getFsiSrcDecayGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::SrcDecay), r, costheta, phi); } /*!< FSI+SRC+decay grid at (r,costheta,phi) */
  std::complex<double> getFsiSrcDecayGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::SrcDecay), costheta, phi); } /*!< FSI+SRC+decay grid at (costheta,phi), r set before */
  std::complex<double> getFsiSrcDecayGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::SrcDecay), phi); } /*!< FSI+SRC+decay grid at phi, r and costheta set before */
  /** @} */
  /** @name The FSI+SRC+CT+decay grid */
  /** @{ */
  /*! FSI+SRC+CT+decay grid at the point set with setinterp
   * \throw std::logic_error if this class has no such grid (see getGridNumber) */
  std::complex<double> getFsiSrcCtDecayGridFull_interp() { return getFsiGridN_interp(variantNumber(GridVariant::SrcCtDecay)); }
  std::complex<double> getFsiSrcCtDecayGridFull_interpvec(const nucqcd::Vector3 &rvec) { return getFsiGridN_interpvec(variantNumber(GridVariant::SrcCtDecay), rvec); } /*!< FSI+SRC+CT+decay grid at the point rvec */
  std::complex<double> getFsiSrcCtDecayGridFull_interp3(const double r, const double costheta, const double phi) { return getFsiGridN_interp3(variantNumber(GridVariant::SrcCtDecay), r, costheta, phi); } /*!< FSI+SRC+CT+decay grid at (r,costheta,phi) */
  std::complex<double> getFsiSrcCtDecayGridFull_interp2(const double costheta, const double phi) { return getFsiGridN_interp2(variantNumber(GridVariant::SrcCtDecay), costheta, phi); } /*!< FSI+SRC+CT+decay grid at (costheta,phi), r set before */
  std::complex<double> getFsiSrcCtDecayGridFull_interp1(const double phi) { return getFsiGridN_interp1(variantNumber(GridVariant::SrcCtDecay), phi); } /*!< FSI+SRC+CT+decay grid at phi, r and costheta set before */
  /** @} */

  //interpolation functions, one without arguments is a pure virtual function
  /*!returns the value of a grid of choice for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiGridN_interpvec(const int grid, const nucqcd::Vector3 &rvec); 
  /*!returns the value of a grid of choice for a certain situation at coordinate (r,costheta,phi) */
  std::complex<double> getFsiGridN_interp3(const int grid, const double r, const double costheta, const double phi); 
  /*!returns the value of a grid of choice for a certain situation at coordinate (costheta,phi), r has been set previously */
  std::complex<double> getFsiGridN_interp2(const int grid, const double costheta, const double phi);
  /*!returns the value of a grid of choice for a certain situation at coordinate (phi), r&theta have been set previously */
  std::complex<double> getFsiGridN_interp1(const int grid, const double phi);
  /*!returns the value of a grid of choice for a certain situation at coordinate (r,theta,phi) that has been set previously
   ,pure virtual function!!*/
  virtual std::complex<double> getFsiGridN_interp(const int grid)=0;
  /*! 3d interpolation of a grid after a certain coordinate (r,costheta,phi) has been set! */
  virtual std::complex<double> getInterp(const Grid3D<std::complex<double> > &grid);  //3d interpolation of grid
  
  
  virtual void fillGrids();  /*!< fills the fsi grids that are used for interpolation */
  /*! updates the fsi grids that are used for interpolation, 
   different from fillGrids.  If new situation of knockoutshells and ejected particles is the same
   we don't have to fill the grids or compute them again and the situation remains as before.*/
  virtual void updateGrids();  
  
  double getR_hit() const{return r_hit;} /*!< [fm] r of the hard interaction point */
  double getCostheta_hit() const{return costheta_hit;} /*!< cos(theta) of the hard interaction point */
  double getSintheta_hit() const{return sintheta_hit;} /*!< sin(theta) of the hard interaction point */
  double getPhi_hit() const{return phi_hit;} /*!< [rad] phi of the hard interaction point */
  double getCosphi_hit() const{return cosphi_hit;} /*!< cos(phi) of the hard interaction point */
  double getSinphi_hit() const{return sinphi_hit;} /*!< sin(phi) of the hard interaction point */
  bool getFilledallgrid() const{return filledallgrid;} /*!< true when all grids (FSI and FSI+CT) are filled for the current particles */
  
protected:
  virtual void setFilenames(std::string dir); /*!< set filenames of the grids \param dir grid cache directory, plus the file name prefixes of the derived classes */ 
  bool filledgrid; /*!< denotes if the grid has been filled */
  bool filledallgrid; /*!< denotes if allthe (possible) grids have been filled */
  double r_hit;  /*!< r coordinate of the hard interaction point, what the grid depens on */
  //double theta_hit;/*!< theta coordinate of the hard interaction point, what the grid depens on */
  double costheta_hit;/*!< cos of theta coordinate of the hard interaction point, what the grid depens on */
  double sintheta_hit;/*!< sin of theta coordinate of the hard interaction point, what the grid depens on */
  double phi_hit;/*!< phi coordinate of the hard interaction point, what the grid depens on */
  double cosphi_hit;/*!< cos of phi coordinate of the hard interaction point, what the grid depens on */
  double sinphi_hit;/*!< sin of phi coordinate of the hard interaction point, what the grid depens on */
  std::string fsi_filename; /*!< filename for regular fsi grid */
  int number_of_grids; /*!<number of fsi grids */
  double prec; /*!< integration precision */
  IntegrationMethod integrator; /*!< which integrator: numint::IntegrationMethod::Romberg (numint::cube_romb) or numint::IntegrationMethod::Adaptive (numint::cube_adaptive); the integers 1 and 2 convert with nucqcd::integration_method() */
  
private:
  std::string dir;  /*!< grid cache directory (with trailing /) */ 
  
  int rgrid; /*!<  the gridsize in r*/
  int cthgrid; /*!<  the gridsize in theta*/
  int phigrid;/*!<  the gridsize in phi*/
  double invrstep; /*!<  the inverse of the stepsize of the grid in r*/
  double invcthstep;/*!<  the inverse of the stepsize of the grid in theta*/
  double invphistep;/*!<  the inverse of the stepsize of the grid in phi*/
  double s_interp; /*!<  s, measure where the point for interpolation is between two indices in r*/
  double t_interp; /*!<  t, measure where the point for interpolation is between two indices in theta*/
  double u_interp; /*!<  u, measure where the point for interpolation is between two indices in phi*/
  double comp_s_interp; /*!<  1-s*/
  double comp_t_interp; /*!<  1-t*/
  double comp_u_interp; /*!<  1-u*/
  int rindex; /*!<  the index for r used in the 3d interpolation*/
  int cthindex; /*!<  the index for cos(theta) used in the 3d interpolation*/
  int phiindex; /*!<  the index for phi used in the 3d interpolation */
  bool allinplane; /*!<  1 if all the isi/fsi particles lie in the x-z plane, can be used to exploit symmetry in phi */
  
  int totalprotonout; /*!<  total protons subject to fsi/isi (isi counts negative)*/
  int totalneutronout;  /*!<  total neutrons subject to fsi/isi (isi counts negative)*/
  int protonknockout;  /*!<  total protons that are knocked out from nucleus*/
  int neutronknockout; /*!<  total neutrons that are knocked out from nucleus*/
    
  std::vector<FastParticle> particles; /*!<  the vector with all particles subject to isi/fsi */
  std::vector<int> knockoutlevels; /*!< the vector with all the shell index values of knocked out nucleons */
  std::vector<int> knockoutm; /*!<  the vector with all the m values of knocked out nucleons */
  
  /*! the grid number of a variant \throw std::logic_error if this class has no such grid */
  int variantNumber(GridVariant variant) const;

  //pure virtual functions!!!
  virtual void constructAllGrids()=0; /*!< construct all grids */
  virtual void readinFsiGrid(std::ifstream &infile)=0; /*!< read in all grids */
  virtual void writeoutFsiGrid(std::ofstream &outfile)=0; /*!< write out all grids */
    
  MeanFieldNucleus *pnucleus;/*!<  the pointer to the MeanFieldNucleus instance */
  
  
};

}  // namespace nucqcd

/*! @} */
#endif