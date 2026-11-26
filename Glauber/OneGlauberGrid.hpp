/*! \file OneGlauberGrid.hpp 
 * \brief Contains declaration of class GlauberGrid
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */
#ifndef ONEGLAUBERGRID_H
#define ONEGLAUBERGRID_H

#include <cstdarg>

#include <Glauber/GlauberGrid.hpp>
#include <Integrator/numint/numint.hpp>

namespace nucqcd {


/*! \brief A class for a RMSGA ISI/FSI grid with one particle, exploiting the symmetry
 * 
 * The only ejected particle is oriented along the z-axis.  So care has to be taken when integrating back over 
 * an axis system where z is along q for instance!!! <BR>
 * 
 * The symmetry means we can seperate two integrations, greatly speeding up the calculation
 * The integration can be split up due to the fact that the density is phi independent, e.g.
 * \f$ \rho(\vec{b},z) = \rho(|\vec{b}|,z) \f$
 * hence the phi integration (angle between b' and b) can be integrated seperately from the 
 * integration over \f$|\vec{b}| \f$ and z
 *
 * Typically an object that is an instance from this class is operated as follows.<BR>
 * 1. Initialize object with constructor OneGlauberGrid()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus  (using addKnockout()
)<BR>
 * 5. Interpolate grid for a certain point or print the grid or whatever... <BR>
 * 
 * Every function that takes a \param grid argument: <BR>
 * 0: RMSGA <BR>
 * 1: RMSGA+CT <BR>
 */

class OneGlauberGrid : public GlauberGrid{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param pnucl pointer to an instance of MeanFieldNucleus
   * \param prec precision you want in the integrations
   * \param integrator which integrator: numint::IntegrationMethod::Romberg (numint::cube_romb) or numint::IntegrationMethod::Adaptive (numint::cube_adaptive); the integers 1 and 2 convert with nucqcd::integration_method()
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   */
  OneGlauberGrid(const int r_grid, const int cth_grid, MeanFieldNucleus *pnucl, 
		 double prec, IntegrationMethod integrator, std::string dir);
  virtual ~OneGlauberGrid();/*!< Destructor */
  /*! 2d interpolation of a grid after a certain coordinate (r,costheta) has been set! */
  virtual std::complex<double> getInterp(std::complex<double> ***grid);  
  /*! add particle subject to isi/fsi, checks to see the vector only contains max one particle 
   *\param newparticle pointer to instance of FastParticle you want to add*/
  virtual void addParticle(FastParticle& newparticle);
  
    
private:
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi, irrelevant here (always 0)
   */
  virtual void calcGlauberphasesBoth(const int i, const int j, const int k); 
  /*! calculates the glauberphases for one gridpoint (only FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi, irrelevant here (always 0)
   */
  virtual void calcGlauberphasesCt(const int i, const int j, const int k);
   /*! set filenames of the grids, includes "One" prefix 
    *\param dir dir where all input/output is located */ 
  virtual void setFilenames(std::string dir);
  
  
  
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.level,p.mm);
    }
    OneGlauberGrid *grid;/*!< pointer to the grid where the integration is performed */
    int level;/*!< knockout level */
    int mm; /*!< index in m_j */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param level knockout level
    * \param mm index in m_j
    */
    void (*f)(numint::vector_d & res, double x, double y, double z, OneGlauberGrid & grid, int level, int mm);
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_bound(numint::vector_d &, double b, double z, double phi, OneGlauberGrid & grid, int level, int mm);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_one_bound_ct(numint::vector_d &, double b, double z, double phi, OneGlauberGrid & grid, int level, int mm);

};

// 




}  // namespace nucqcd

/** @} */
#endif
