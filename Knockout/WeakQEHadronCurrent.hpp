/*! \file WeakQEHadronCurrent.hpp 
 * \brief Contains declaration of class WeakQEHadronCurrent, used to compute A(nu,lN) amplitudes
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */

#ifndef WEAKQEHADRONCURRENT_HPP
#define WEAKQEHADRONCURRENT_HPP

#include <Tools/Core/FourVector.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Glauber/AbstractFsiCTGrid.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Glauber/NucleonWeakOperator.hpp>
#include <Tools/Core/Vector3.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/OneGlauberGrid.hpp>
#include <Glauber/ROMEAGrid.hpp>
#include <Integrator/numint/numint.hpp>
#include <Tools/Core/Spinor.hpp>

#include <string>
#include <cstdarg>

namespace nucqcd {




/*! \brief A class WeakQEHadronCurrent, used to compute A(nu,lN) amplitudes */
class WeakQEHadronCurrent{
public:
    /*! Constructor
   * \param pnucleus pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator (1:numint::cube_romb, 2:numint::cube_adaptive (0, the old rombergerN, was removed)
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   * \param max_Eval max # of function evaluations in the coordinate integration
   * \param charged neutral weak current [0] or charged weak current [1]
   * \param M_A_in [MeV] value of axial mass
   * \param user_sigma does the user want to change sigma?
   * \param gA_s [] strange contrib to G_A
   * \param r_s2 [fm^2] strange contrib to F1weak
   * \param mu_s [] strange contrib to F2weak
   * \param sigma_screening [%] screening change of sigma
   * \param enable_ROMEA enable ROMEA FSI for slow nucleons...
   */
  WeakQEHadronCurrent(MeanFieldNucleusThick *pnucleus, double prec, int integrator, std::string dir, int max_Eval,
	bool charged, double M_A_in, bool user_sigma, bool enable_ROMEA, double gA_s=-0.19, double r_s2=0., 
	double mu_s=0., double sigma_screening=0.);
  ~WeakQEHadronCurrent(); /*!< Destructor */
  
  
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) <BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1, 3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param CT do you want CT effects in the Glauber FSI?
   * \param pw no FSI?
   * \param SRC SRC or not in the FSI
   * \param thick thickness in the glauber
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param vector_or_axial [1] calc vector part of the current [0] calc axial vector part of the current
   */
  void getMatrixEl(Kinematics2to2 &tk, Matrix<2,4> &results, int shellindex, int m, int CT, int pw, int current,
    int SRC, int thick, bool vector_or_axial);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1, 3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   */
  void getAllMatrixElMult(Kinematics2to2 &tk, Matrix<2,4> *results, int shellindex, int m, int current, int thick);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all four vector boson polarizations (0,-1,+1,z) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2} contains amplitudes<BR>
   * (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) if thickness <BR>
   * (0=RMSGA,1=+SRC,2=plane-wave) if no thickness <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1,3 is z) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   * \param medium medium modifications (0=none, 1=QCM, 2=CQSM)
   */
  void getAllMatrixEl(Kinematics2to2 &tk, Matrix<2,4> *results, int shellindex, int m, int current, int thick, int medium);
   /*! Computes the on-shell amplitude \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu u(0,m_i) \f$ in the nucleon rest fram, q vector along the z axis
   * \param Q2 [MeV^2] virtual photon Q-squared
   * \param proton [1] proton or [0] neutron
   * \param charged [0] Z boson [1] W boson
   * \param M_A [MeV] axial mass
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param spinin spin \f$ m_i \f$ of the initial nucleon times TWO!!! (-1 or +1)
   * \param spinout spin \f$ m_f \f$ of the outgoing nucleon times TWO!!! (-1 or +1)
   * \param photonpol polarization of the photon (-1,0, +1, or z(3))
   * \return [dimensionless] on-shell amplitude
   */
  static std::complex<double> getFreeMatrixEl(double Q2, bool proton, bool charged, double M_A, int current, int spinin, int spinout, int photonpol);
  double getPrec() const{return prec;} /*!< returns precision in the integrations */
  bool getUsersigma() const{return usersigma;} /*!< returns 1 if user has changed sigma with some screening */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] returns screening change to sigma */
  AbstractFsiGrid *getGrid() const{return grid;} /*!< Glauber grid used in the current calculation (with or without thickness) */
  GlauberGridThick *getThickGrid() {return &gridthick;} /*!< Glauber grid with thickness (SRC), see gridthick */
  OneGlauberGrid *getOneGrid() {return &onegrid;} /*!< Glauber grid without thickness, see onegrid */
  ROMEAGrid *getROMEAGrid() {return &romgrid;} /*!< ROMEA grid for slow nucleons, see romgrid */
  int getShell() const{return shell;} /*!< shell index of the knocked-out nucleon in the current calculation */
  int getM() const{return mm;} /*!< \f$ 2m_j \f$ of the knocked-out nucleon in the current calculation */
  int getIntegrator() const{return integrator;} /*!< choice of integrator */
  bool getEnable_ROMEA() const{return enable_ROMEA;} /*!< ROMEA instead of Glauber FSI for slow nucleons */
  MeanFieldNucleusThick* getPnucleus() {return pnucl;} /*!< the nucleus */
  const Matrix<1,4> & getBarcontract0up() const{return barcontract0up;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, polarization vector (1,0,0,0): time component J^0 */
  const Matrix<1,4> & getBarcontractminup() const{return barcontractminup;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, transverse polarization, helicity -1 with respect to q */
  const Matrix<1,4> & getBarcontractplusup() const{return barcontractplusup;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, transverse polarization, helicity +1 with respect to q */
  const Matrix<1,4> & getBarcontractzup() const{return barcontractzup;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, spatial component J^z along q */
  const Matrix<1,4> & getBarcontract0down() const{return barcontract0down;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, polarization vector (1,0,0,0): time component J^0 */
  const Matrix<1,4> & getBarcontractmindown() const{return barcontractmindown;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, transverse polarization, helicity -1 with respect to q */
  const Matrix<1,4> & getBarcontractplusdown() const{return barcontractplusdown;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, transverse polarization, helicity +1 with respect to q */
  const Matrix<1,4> & getBarcontractzdown() const{return barcontractzdown;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, spatial component J^z along q */
  const Matrix<1,4> & getBarcontract() const{return barcontract;} /*!< contraction of the outgoing-nucleon bar spinor with the current (the last one computed) */
  const NucleonWeakOperator & getJ() const{return *J;} /*!< nucleon current operator of the current calculation */
  nucqcd::Vector3 & getPm() {return pm;} /*!< [MeV] missing momentum of the current calculation */
  
private:
  double prec; /*!< precision in the integrations */
  int integrator; /*!< choice of integrator */
  AbstractFsiGrid *grid; /*!< pointer to Glauber Grid */
  MeanFieldNucleusThick *pnucl; /*!<  pointer to nucleus */
  NucleonWeakOperator *J; /*!<  pointer to nucleon formfactor instance */
  Matrix<1,4> barcontract; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0up; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractminup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractzup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0down; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractmindown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractzdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  nucqcd::Vector3 pm; /*!< missing momentum */
  std::string homedir; /*!< contains the share dir with all input */
  bool charged; /*!< neutral weak current [0] or charged [1] */
  bool usersigma; /*!< has the user set some screening to sigma */
  double sigmascreening; /*!< [%] screening effect to sigma */
  bool enable_ROMEA; /*!< enable ROMEA FSI for slow nucleons instead of Glauber */
  GlauberGridThick gridthick; /*!< grid for calc with thickness */
  OneGlauberGrid onegrid; /*!< grid for calc without thickness */
  ROMEAGrid romgrid; /*!< grid for slow nucleons (ROMEA formalism) */
  int maxEval; /*!< max Evals in integrations */
  double gA_s; /*!< [] strange contrib to G_A */
  double mu_s; /*!< [] strange contrib to F2_weak */
  double r_s2; /*!< [fm^2] strange contrib to F1_weak */
  double M_A; /*!< [MeV] axial mass value */
  
  int shell; /*!< what shell has the ejected nucleon */
  int mm; /*!< \f$ m_j \f$ quantum number of ejected nucleon */

  
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.SRC,p.pw);
    }
    WeakQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int SRC;/*!< SRC or not */
    int pw; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model WeakQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, WeakQEHadronCurrent & model, int SRC, int pw);
  };

  /*! struct that is used for integrators (medium modif ones)*/
  struct Ftor_medium {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_medium &p = * (Ftor_medium *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.thick,p.total,p.spinup, p.spindown, p.medium, p.q, p.pi,p.pf, 
	  p.polmin, p.polplus, p.pol0, p.polz, p.current);
    }
    WeakQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int thick;/*!< thickness or not */
    int total; /*!< number of different FSI? */
    Matrix<1,4> spinup; /*!< bar spinor of the outgoing nucleon, spin up */
    Matrix<1,4> spindown; /*!< bar spinor of the outgoing nucleon, spin down */
    int medium; /*!< which medium modif */
    FourVector<double> q; /*!< [MeV] four-momentum of the virtual boson */
    FourVector<double> pi; /*!< [MeV] four-momentum of the initial nucleon */
    FourVector<double> pf; /*!< [MeV] four-momentum of the outgoing nucleon */
    FourVector<std::complex<double> > polmin; /*!< polarization four-vector of the virtual boson, transverse, helicity -1 with respect to q */
    FourVector<std::complex<double> > polplus; /*!< polarization four-vector of the virtual boson, transverse, helicity +1 with respect to q */
    FourVector<std::complex<double> > pol0; /*!< polarization four-vector (1,0,0,0) of the virtual boson: picks the time component J^0 */
    FourVector<std::complex<double> > polz; /*!< unit vector along q: picks the spatial component J^z along q */
    int current; /*!<which current operator */
    
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model WeakQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, WeakQEHadronCurrent & model, int thick, int total,
	      Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
	      const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0, const FourVector<std::complex<double> > &polz,
				   int current);
    
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int SRC, int pw);
  /*! integrandum function (clean ones)*/
  static void klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int thick, int total);
  /*! integrandum function (medium modification)*/
  static void klaas_all_amp_medium(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model,
				   int thick, int total,Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
				   const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0,const FourVector<std::complex<double> > &polz,
				   int current);
  /*! integrandum function (clean ones)*/
  static void klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, WeakQEHadronCurrent & model, int SRC, int pw);
   

};

}  // namespace nucqcd

/** @} */
#endif
