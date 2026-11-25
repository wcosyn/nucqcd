/*! \file Resonance.hpp 
 * \brief Contains declaration of class Resonance, has particle and scattering properties
 * \author Wim Cosyn
 * \date 29/08/2012
 * 
 * \addtogroup DIS
 * @{
 */
#ifndef RESONANCE_HPP
#define RESONANCE_HPP

#include <complex>

namespace nucqcd {


/*! \brief A class for a certain resonance with scattering parameters and others
 *
 * Used for the produced hadronic states X in the final-state interactions of DIS on the deuteron:
 * the rescattering of X off the spectator nucleon is described with a diffractive amplitude
 * \f$ \sigma (i + \epsilon) e^{\beta t/2} \f$, with off-shell corrections. See W. Cosyn and
 * M. Sargsian, Phys. Rev. C 84, 014601 (2011), and W. Cosyn, W. Melnitchouk and M. Sargsian,
 * Phys. Rev. C 89, 014612 (2014). The constructor takes beta, betaoff and lambda in GeV units, the
 * getters return them in MeV units.
 */
class Resonance{
public:
  /*! Constructor
   * \param m mass [MeV]
   * \param coeff coefficient with which the resonance gets weighted in the DIS amplitude
   * \param sigma0 [mb] base value of sigma0
   * \param sigmaslope [GeV] slope with which sigma varies in function of Q^2
   * \param betain [GeV^-2] slope parameter
   * \param epsilonin ratio of real to imaginary part of the scattering amplitude
   * \param betaoffin [GeV^-2] off-shell beta parameter
   * \param lambdain [GeV^2] lambda cutoff off-shell parameter
   */
  Resonance(double m, std::complex<double> coeff, double sigma0=40., double sigmaslope=0.,
	    double betain=6.,double epsilonin=-0.5,double betaoffin=8., double lambdain=1.2);
  Resonance(const Resonance &copy); /*!< Copy Constructor */
  Resonance& operator=(const Resonance&); /*!< copy assignment */
  ~Resonance(); /*!< Destructor */
  /*! Gives you sigma as a function of Q^2 (inferred from semi-inclusive DIS deuteron fits)
   * increases with W but levels off for W>2.4 GeV
   * \f$ \sigma(Q^2) = \sigma_0 + \sigma_{slope} *(M-M_p)/Q^2\f$
   * \param Q2 [MeV^2] Q^2
   * \return [MeV^-2] sigma scattering parameter
   */
  double getSigma(double Q2) const;
  double getBeta(double Q2=0.) const{return beta;} /*!< [MeV^-2] slope parameter of the rescattering amplitude \param Q2 [MeV^2] unused (no Q^2 dependence) */
  double getBetaoff() const{return betaoff;} /*!< [MeV^-2] off-shell slope parameter */
  double getSigma0() const{return sigma0;} /*!< [mb] base value of the total cross section sigma, see getSigma */
  double getSigmaslope() const{return sigmaslope;} /*!< [GeV] slope of sigma in (M-M_p)/Q^2, see getSigma */
  double getEpsilon() const{return epsilon;} /*!< ratio of real to imaginary part of the rescattering amplitude */
  double getLambda() const{return lambda;} /*!< [MeV^2] cutoff of the off-shell suppression */
  double getMass2() const{return mass2;} /*!< [MeV^2] returns mass squared */
  double getMass() const{return mass;} /*!< [MeV] mass of the resonance */
  std::complex<double> getCoeff() const{return coeff;} /*!< coefficient with which the resonance is weighted in the DIS amplitude */
  
  void setCoeff(std::complex<double> coeff_in){coeff=coeff_in;} /*!< change the coefficient for the DIS amplitude */
  
private:
  double mass; /*!< mass [MeV]*/
  double mass2; /*!< [MeV^2] mass squared*/
  double sigma0; /*!< [mb] base value of sigma0*/
  double sigmaslope; /*!< [GeV] slope with which sigma varies in function of Q^2*/
  double beta; /*!< [MeV^-2] slope parameter*/
  double betaoff; /*!< [MeV^-2] off-shell beta parameter*/
  double epsilon; /*!<  real part of scattering amplitude*/
  double lambda; /*!< [MeV^2] lambda cutoff off-shell parameter*/
  std::complex<double> coeff; /*!< coefficient with which the resonance gets weighted in the DIS amplitude */  
  
};


}  // namespace nucqcd

/*! @} */
#endif