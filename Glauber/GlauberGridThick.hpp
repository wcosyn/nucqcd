/*! \file GlauberGridThick.hpp 
 * \brief Contains declaration of class GlauberGridThick
 * \author Wim Cosyn
 * \date 16/08/2011
 * \addtogroup Glauber
 * @{
 */

#ifndef GLAUBERGRIDTHICK_H
#define GLAUBERGRIDTHICK_H

#include <Integrator/numint/numint.hpp>


#include <cstdarg>
#include <Glauber/GlauberThickGridCommon.hpp>


namespace nucqcd {


/*! \brief An class for a Glauber ISI/FSI grid using thickness approximation, with the variants FSI, FSI+SRC, FSI+CT and FSI+SRC+CT (see GlauberThickGridCommon)
 * 
 * Four grids are actually calculated.  fsi_grid contains both the regular fsi grid as well as the grid that includes SRC.  
 * Same thing for the fsi_ct_grid
 * 
 * Typically an object that is an instance from this class is operated as follows.<BR>
 * 1. Initialize object with constructor GlauberGridThick()<BR>
 * 2. Add all particles subject to ISI/FSI with addParticle() <BR>
 * 3. Call fillGrids() or updateGrids() <BR>
 * 4. Add particles that are knocked out from nucleus  (using addKnockout()
)<BR>
 *  5. Interpolate grid for a certain point or print the grid or whatever... 
(use getFsiGridFull_interp3() for instance) <BR>
 * 
 * Every function that takes a \param grid argument: <BR>
 * 0: RMSGA <BR>
 * 1: RMSGA+SRC <BR>
 * 2: RMSGA+CT <BR>
 * 3: RMSGA+SRC+CT <BR>
 */
class GlauberGridThick : public GlauberThickGridCommon{
public:
  /*! Constructor
   * \param r_grid gridsize in r 
   * \param cth_grid gridsize in costheta
   * \param phi_grid gridsizein phi
   * \param pnuclthick pointer to an instance of MeanFieldNucleusThick
   * \param prec precision you want in the integrations
   * \param integrator which integrator: numint::IntegrationMethod::Romberg (numint::cube_romb) or numint::IntegrationMethod::Adaptive (numint::cube_adaptive); the integers 1 and 2 convert with nucqcd::integration_method()
   * \param dir std::string that contains dir with all input, should be the ./share subdir of the project!
   */
  GlauberGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick *pnuclthick, 
		   double prec, IntegrationMethod integrator, std::string dir);
  virtual ~GlauberGridThick();/*!< Destructor */

  virtual void printFsi_grid();/*!< Prints the FSI grid for a certain situation*/
  virtual void printFsi_ct_grid(); /*!< Prints the FSI+CT grid for a certain situation */
  virtual void printFsi_src_grid();/*!< Prints the FSI+SRC grid for a certain situation*/
  virtual void printFsi_src_ct_grid();/*!< Prints the FSI+SRC+CT grid for a certain situation*/
  /*!gives you a glauberphase at a certain point (all 4 RMSGA variants)
   * \param phases result
   * \param rhit [fm] radial coordinate of hard interaction point
   * \param costhetahit cos of polar angle of hard interaction point
   * \param phihit phi angle of hard interaction point
   */
  void getFsiphaseAll(std::vector<std::complex<double> > &phases, double rhit, double costhetahit, double phihit);
  /*!gives you a glauberphase at a certain point (only 2 CT+RMSGA variants)
   * \param phases result
   * \param rhit [fm] radial coordinate of hard interaction point
   * \param costhetahit cos of polar angle of hard interaction point
   * \param phihit phi angle of hard interaction point
   */
  void getFsiphaseCt(std::vector<std::complex<double> > &phases, double rhit, double costhetahit, double phihit);



  
protected:
  
  
private:

  
  virtual void setFilenames(std::string dir); /*!< set filenames of the grids \param dir dir where all input/output is located */ 

  
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi
   */
  virtual void calcGlauberphasesBoth(const int i, const int j, const int k); 
  /*! calculates the glauberphases for one gridpoint (only FSI+CT)
   * \param i grid index in r
   * \param j grid index in costheta
   * \param k grid index in phi
   */
  virtual void calcGlauberphasesCt(const int i, const int j, const int k);
  /*! calculates the glauberphases for one previously set point (both FSI and FSI+CT)
   * \param protonphases glauberphases (for proton scattering)
   * \param neutronphases glauberphases (for neutron scattering)
   */
  virtual void calcGlauberphasesBoth(std::vector< std::complex<double> > &protonphases, std::vector< std::complex<double> > &neutronphases); 
  /*! calculates the glauberphases for one previously set point (only FSI+CT)
   * \param protonphases glauberphases (for proton scattering)
   * \param neutronphases glauberphases (for neutron scattering)
   */
  virtual void calcGlauberphasesCt(std::vector<std::complex<double> > &protonphases, std::vector<std::complex<double> > &neutronphases);
  
  /*! struct that is used for integrators (dirty ones)*/
  struct Ftor_all {
    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_all &p = * (Ftor_all *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.proton);
    }
    GlauberGridThick *grid; /*!< pointer to the grid where the integration is performed */
    int proton; /*!< proton glauberphase or not */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param proton proton glauberphase or not
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, GlauberGridThick & grid, int proton);
  };
 /*! struct that is used for integrators (clean ones)*/
   struct Ftor_bound {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
      Ftor_bound &p = * (Ftor_bound *) param;
      p.f(ret,x[0],x[1],x[2],*p.grid,p.proton);
    }
    GlauberGridThick *grid;/*!< pointer to the grid where the integration is performed */
    int proton;/*!< proton glauberphase or not */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param grid the grid instance
    * \param proton proton glauberphase or not
    */
    void (*f)(numint::vector_d &, double x, double y, double z, GlauberGridThick &, int proton);
  };

  /*! integrandum function  (dirty ones)*/
  static void klaas_int_all(numint::vector_z & ret, double r, double costheta, double phi, GlauberGridThick & grid, int proton);
  /*! integrandum function (clean ones)*/
  static void klaas_int_bound(numint::vector_d &, double b, double z, double phi, GlauberGridThick & grid, int proton);
  /*! integrandum function  (dirty ones), only CT*/
  static void klaas_int_all_ct(numint::vector_z & ret, double r, double costheta, double phi, GlauberGridThick & grid, int proton);
  /*! integrandum function (clean ones), only CT*/
  static void klaas_int_bound_ct(numint::vector_d &, double b, double z, double phi, GlauberGridThick & grid, int proton);

  /*! integrandum in the frame of the particle (see GlauberThickGridCommon::integrateTilted), results as klaas_int_bound: fsi, fsi+src, fsi+ct, fsi+ct+src */
  static void klaas_int_tilted(numint::vector_d &, double b, double z, double phi, GlauberThickGridCommon & grid, int proton, const TiltedFrame &frame);
  /*! integrandum in the frame of the particle, only CT: fsi+ct, fsi+ct+src */
  static void klaas_int_tilted_ct(numint::vector_d &, double b, double z, double phi, GlauberThickGridCommon & grid, int proton, const TiltedFrame &frame);
};


}  // namespace nucqcd

/** @} */
#endif