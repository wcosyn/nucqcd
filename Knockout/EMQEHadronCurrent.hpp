/*! \file EMQEHadronCurrent.hpp 
 * \brief Contains declaration of class EMQEHadronCurrent, used to compute A(e,e'N) amplitudes
 * \author Wim Cosyn
 * \date 28/08/2012
 * 
 * \addtogroup Knockout
 * @{
 */

#ifndef EMQEHADRONCURRENT_HPP
#define EMQEHADRONCURRENT_HPP

#include <Tools/Core/FourVector.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Glauber/AbstractFsiCTGrid.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Glauber/NucleonEMOperator.hpp>
#include <Tools/Core/Vector3.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/OneGlauberGrid.hpp>
#include <Integrator/numint/numint.hpp>
#include <Tools/Core/Spinor.hpp>

#include <string>
#include <cstdarg>

namespace nucqcd {




/*! \brief A class EMQEHadronCurrent, used to compute A(e,e'N) amplitudes */
class EMQEHadronCurrent{
public:
    /*! Constructor
   * \param pnucleus pointer to a MF nucleus instance
   * \param prec precision you want in the integrations
   * \param integrator which integrator: numint::IntegrationMethod::Romberg (numint::cube_romb) or numint::IntegrationMethod::Adaptive (numint::cube_adaptive); the integers 1 and 2 convert with nucqcd::integration_method()
   * \param dir string that contains dir with all input, should be the ./share subdir of the project!
   * \param max_Eval max # of function evaluations in the coordinate integration
   * \param user_sigma does the user want to change sigma?
   * \param sigma_screening [%] screening change of sigma
   */
  EMQEHadronCurrent(MeanFieldNucleusThick *pnucleus, double prec, IntegrationMethod integrator, std::string dir, int max_Eval,
	bool user_sigma, double sigma_screening=0.);
  ~EMQEHadronCurrent(); /*!< Destructor */
//   void setSRC(int setSRC) {SRC=setSRC;} /*!< sets SRC in the Glauber FSI */
//   std::complex<double> getMatrixEl(Kinematics2to2 &tk,int spinout, int photonpol, int shellindex, int m, int CT, 
// 				   int pw, int current, int SRC, int thick);
  
  
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) <BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param CT do you want CT effects in the Glauber FSI?
   * \param pw no FSI?
   * \param SRC SRC or not in the FSI
   * \param thick thickness in the glauber
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param lc_mod [] modification of coherence length in CT 
   * \param nkt_mod [] modification of initial CT sigma value
   */
  void getMatrixEl(Kinematics2to2 &tk, Matrix<2,3> &results, int shellindex, int m, int CT, int pw, int current,
    int SRC, int thick, double lc_mod, double nkt_mod);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2}] contains amplitudes (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   * \param lc_mod [] modification of coherence length in CT 
   * \param nkt_mod [] modification of initial CT sigma value
   */
  void getAllMatrixElMult(Kinematics2to2 &tk, Matrix<2,3> *results, int shellindex, int m, int current, int thick, double lc_mod, double nkt_mod);
  /*! Computes amplitudes \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu \phi^{D}_{\alpha}(\vec{p}_m,m) \f$
   * for all three photon polarizations (0,-1,+1) and both final nucleon helicities (-1,+1) and all glauber varieties<BR>
   * Computed in the frame where the z-axis lies along the ejected nucleon!!!!!!
   * \param tk contains the hadron kinematics
   * \param [out] results [ fm^{3/2} contains amplitudes<BR>
   * (0=RMSGA,1=+SRC,2=+CT,3+CT+SRC,4=plane-wave) if thickness <BR>
   * (0=RMSGA,1=+SRC,2=plane-wave) if no thickness <BR>
   * First index is final nucleon helicity (0 is down, 1 is up) <BR>
   * Second index is photon polarization (0 is 0, 1 is -1, 2 is +1) <BR>
   * \param shellindex which \f$ \alpha \f$ shell do we eject from
   * \param m \f$ m_j \f$ times TWO!!! of the initial nucleon
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param thick thickness in the glauber
   * \param medium medium modifications (0=none, 1=QCM, 2=CQSM)
   * \param lc_mod [] modification of coherence length in CT 
   * \param nkt_mod [] modification of initial CT sigma value
   */
  void getAllMatrixEl(Kinematics2to2 &tk, Matrix<2,3> *results, int shellindex, int m, 
        int current, int thick, int medium, double lc_mod, double nkt_mod);
   /*! Computes the off-shell amplitude \f$ \bar{u}(\vec{p}_f,m_f)\Gamma^{\mu}\epsilon_\mu u(\vec{p}_m,m_i) \f$
   * \param tk contains the hadron kinematics
   * \param current selects the current operator [1=CC1, 2=CC2, 3=CC3], see T. de Forest, Nucl. Phys. A 392, 232 (1983).
   * \param spinin spin \f$ m_i \f$ of the initial nucleon times TWO!!! (-1 or +1)
   * \param spinout spin \f$ m_f \f$ of the outgoing nucleon times TWO!!! (-1 or +1)
   * \param photonpol polarization of the photon (-1,0 or +1)
   * \return [dimensionless] off-shell amplitude
   */
  std::complex<double> getFreeMatrixEl(Kinematics2to2 &tk, int current, int spinin, int spinout, int photonpol);
  double getPrec() const{return prec;} /*!< returns precision in the integrations */
  bool getUsersigma() const{return usersigma;} /*!< returns 1 if user has changed sigma with some screening */
  double getSigmascreening() const{return sigmascreening;} /*!< [%] returns screening change to sigma */
  AbstractFsiCTGrid *getGrid() const{return grid;} /*!< Glauber grid used in the current calculation (with or without thickness) */
  GlauberGridThick* getThickGrid() {return &gridthick;} /*!< Glauber grid with thickness (SRC), see gridthick */
  int getShell() const{return shell;} /*!< shell index of the knocked-out nucleon in the current calculation */
  int getM() const{return mm;} /*!< \f$ 2m_j \f$ of the knocked-out nucleon in the current calculation */
  IntegrationMethod getIntegrator() const{return integrator;} /*!< choice of integrator */
  MeanFieldNucleusThick* getPnucleus() {return pnucl;} /*!< the nucleus */
  const Matrix<1,4> & getBarcontract0up() const{return barcontract0up;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, polarization vector (1,0,0,0): time component J^0 */
  const Matrix<1,4> & getBarcontractminup() const{return barcontractminup;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, transverse polarization, helicity -1 with respect to q */
  const Matrix<1,4> & getBarcontractplusup() const{return barcontractplusup;} /*!< bar spinor of the outgoing nucleon (spin up) contracted with the current, transverse polarization, helicity +1 with respect to q */
  const Matrix<1,4> & getBarcontract0down() const{return barcontract0down;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, polarization vector (1,0,0,0): time component J^0 */
  const Matrix<1,4> & getBarcontractmindown() const{return barcontractmindown;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, transverse polarization, helicity -1 with respect to q */
  const Matrix<1,4> & getBarcontractplusdown() const{return barcontractplusdown;} /*!< bar spinor of the outgoing nucleon (spin down) contracted with the current, transverse polarization, helicity +1 with respect to q */
  const Matrix<1,4> & getBarcontract() const{return barcontract;} /*!< contraction of the outgoing-nucleon bar spinor with the current (the last one computed) */
  const NucleonEMOperator & getJ() const{return *J;} /*!< nucleon current operator of the current calculation */
  nucqcd::Vector3 & getPm() {return pm;} /*!< [MeV] missing momentum of the current calculation */
  
private:
  double prec; /*!< precision in the integrations */
  IntegrationMethod integrator; /*!< choice of integrator */
  AbstractFsiCTGrid *grid; /*!< pointer to Glauber Grid */
  MeanFieldNucleusThick *pnucl; /*!<  pointer to nucleus */
  NucleonEMOperator *J; /*!<  pointer to nucleon formfactor instance */
  Matrix<1,4> barcontract; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0up; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractminup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusup; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontract0down; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractmindown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  Matrix<1,4> barcontractplusdown; /*!< intermediate contraction of bar spinor with 4*4 current  */
  nucqcd::Vector3 pm; /*!< missing momentum */
  std::string homedir; /*!< contains the share dir with all input */
  bool usersigma; /*!< has the user set some screening to sigma */
  double sigmascreening; /*!< [%] screening effect to sigma */
  GlauberGridThick gridthick; /*!< Glauber grid with thickness, used with SRC or thickness */
  OneGlauberGrid onegrid; /*!< Glauber grid without thickness (one particle along z) */
  int maxEval; /*!< max Evals in integrations */
  
  int shell; /*!< what shell has the ejected nucleon */
  int mm; /*!< \f$ m_j \f$ quantum number of ejected nucleon */
//     /*! function that gets integrated over r, integral over r of matrix element
//    * \param r [fm] radial coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJR(const double r, std::complex<double> *results, va_list ap);
//   /*! function that gets integrated over cos(theta), integral of matrix element
//    * \param costheta [] cos of theta coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJCosTheta(const double costheta, std::complex<double> *results, va_list ap);
//   /*! function that gets integrated over phi, integral of matrix element
//    * \param phi [rad] phi coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJPhi(const double phi, std::complex<double> *results, va_list ap);
//     /*! function that gets integrated over r, computes 12 different amplitudes:
//      * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
//    * \param r [fm] radial coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJR12(const double r, std::complex<double> *results, va_list ap);
//   /*! function that gets integrated over cos(theta), computes 12 different amplitudes:
//      * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
//    * \param costheta [] cos of theta coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJCosTheta12(const double costheta, std::complex<double> *results, va_list ap);
//   /*! function that gets integrated over phi,computes 12 different amplitudes:
//      * 3 photon polarizations, 2 final nucleon helicities, including and excluding CT
//    * \param phi [rad] phi coordinate
//    * \param results result: contains the glauberphases (fsi and fsi+ct) for a gridpoint
//    * \param ap variable parameter list
//    */
//   void intJPhi12(const double phi, std::complex<double> *results, va_list ap);
  
  /*! struct that is used for integrators (clean ones)*/
  struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.SRC,p.pw);
    }
    EMQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int SRC;/*!< SRC or not */
    int pw; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model EMQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, EMQEHadronCurrent & model, int SRC, int pw);
  };
  /*! struct that is used for integrators (radial ones)*/
  struct Ftor_r {

    /*! integrandum function */
    static void exec(const numint::array<double,1> &x, void *param, numint::vector_z &ret) {
      Ftor_r &p = * (Ftor_r *) param;
      p.f(ret,x[0],*p.model,p.thick,p.total);
    }
    EMQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int thick;/*!< SRC or not */
    int total; /*!< plane-wave? */
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param model EMQEHadronCurrent instance
    * \param thick thick in FSI?
    * \param total number of FSI situations?
    */
    void (*f)(numint::vector_z & res, double r, EMQEHadronCurrent & model, int thick, int total);
  };
  /*! struct that is used for integrators (angular ones)*/
  struct Ftor_angles {

    /*! integrandum function */
    static void exec(const numint::array<double,2> &x, void *param, numint::vector_z &ret) {
      Ftor_angles &p = * (Ftor_angles *) param;
      p.f(ret,x[0],x[1],*p.model,p.thick,p.total,p.r);
    }
    EMQEHadronCurrent *model;/*!< pointer to model instance that contains all */
    int thick;/*!< SRC or not */
    int total; /*!< plane-wave? */
    double r; /*!<radial coordinate*/
    /*! integrandum 
    * \param res results
    * \param r first integration variable
    * \param model EMQEHadronCurrent instance
    * \param thick thick in FSI?
    * \param total number of FSI situations?
    */
    void (*f)(numint::vector_z & res, double costheta, double phi, EMQEHadronCurrent & model, int thick, int total, double r);
  };
  /*! struct that is used for integrators (medium modif ones)*/
  struct Ftor_medium {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_medium &p = * (Ftor_medium *) param;
      p.f(ret,x[0],x[1],x[2],*p.model,p.thick,p.total,p.spinup, p.spindown, p.medium, p.q, p.pi,p.pf, 
	  p.polmin, p.polplus, p.pol0, p.current);
    }
    EMQEHadronCurrent *model;/*!< pointer to model instance that contains all */
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
    int current; /*!<which current operator */
    
    /*! integrandum 
    * \param res results
    * \param x first integration variable
    * \param y second integration variable
    * \param z third integration variable
    * \param model EMQEHadronCurrent instance
    * \param SRC SRC in FSI?
    * \param pw plane-wave?
    */
    void (*f)(numint::vector_z & res, double x, double y, double z, EMQEHadronCurrent & model, int thick, int total,
	      Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
	      const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0,
				   int current);
    
  };
  /*! integrandum function (clean ones)*/
  static void klaas_one_amp(numint::vector_z & results, double r, double costheta, double phi, EMQEHadronCurrent & model, int SRC, int pw);
  /*! integrandum function (clean ones)*/
  static void klaas_all_amp(numint::vector_z & results, double r, double costheta, double phi, EMQEHadronCurrent & model, int thick, int total);
  /*! integrandum function (medium modification)*/
  static void klaas_all_amp_medium(numint::vector_z & results, double r, double costheta, double phi, EMQEHadronCurrent & model,
				   int thick, int total,Matrix<1,4> &spinup, Matrix<1,4> &spindown, const int medium,
				   const FourVector<double> &q, const FourVector<double> &pi, const FourVector<double> &pf,
				   const FourVector<std::complex<double> > &polmin, const FourVector<std::complex<double> > &polplus, 
				   const FourVector<std::complex<double> > &pol0,
				   int current);
  /*! integrandum function (clean ones)*/
  static void klaas_mult_amp(numint::vector_z & results, double r, double costheta, double phi, EMQEHadronCurrent & model, int SRC, int pw);
  /*! integrandum function (radial ones)*/
  static void klaas_all_radial(numint::vector_z & results, double r, EMQEHadronCurrent & model, int SRC, int pw);
  /*! integrandum function (angular ones)*/
  static void klaas_all_angular(numint::vector_z & results, double costheta, double phi, EMQEHadronCurrent & model, int SRC, int pw, double r);
  

};

}  // namespace nucqcd

/** @} */
#endif