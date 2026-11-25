/*! \file GrossWavefunction.hpp
 * \brief Relativistic deuteron wave functions of Gross and collaborators (IIb, L4, WJC1, WJC2)
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on February 12, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef GROSSWAVEFUNCTION_H
#define GROSSWAVEFUNCTION_H

#include <Tools/WaveFunctions/WavefunctionImplementation.hpp>

namespace nucqcd {


/*! \brief Relativistic (Gross) deuteron wave function: S, D and the two P waves, parametrized with
 * functions f_{L,i} of masses \f$ m_{L,i} = \alpha_L + (i-1) m_0 \f$ (see GetF, GetG and GetK for the
 * expansion and the coefficients that impose the correct small-r behaviour)
 */
class GrossWavefunction : public WavefunctionImplementation
{
 private:
  /// coordinate or momentum space
  enum ESpace {
    kR = 1,     ///< configuration space
    kP = 2      ///< momentum space
  };
  
  /// type of wave
  enum EWave {
    kU = 0,    ///< S wave
    kW = 1,    ///< D wave
    kVT = 2,   ///< triplet P wave
    kVS = 3    ///< singlet P wave
  };

 public:
  /*! constructor
   * \param order number of terms N
   * \param alpha [MeV] mass parameters alpha_L (array of 3, for L = 0, 1, 2)
   * \param m0 [MeV] second mass parameter, m_{L,i} = alpha_L + (i-1) m0
   * \param parU,parW,parVT,parVS [MeV^1/2] parameters b_{L,n} of the S, D, triplet P and singlet P waves */
  GrossWavefunction(int order, double *alpha, double m0, 
	    double *parU, double *parW, double *parVT, double *parVS);
  GrossWavefunction(const GrossWavefunction&); ///< copy constructor
  GrossWavefunction& operator=(const GrossWavefunction&); ///< assignment
  virtual ~GrossWavefunction(); ///< destructor

  virtual GrossWavefunction *Clone(const char* ="") const {return new GrossWavefunction(*this); } ///< virtual copy constructor

  virtual double GetUr(double r) const; ///< [fm^-1/2] S wave u (L=0) in coordinate space \param r [fm]
  virtual double GetWr(double r) const; ///< [fm^-1/2] D wave w (L=2) in coordinate space \param r [fm]
  virtual double GetVTr(double r) const; ///< [fm^-1/2] spin-triplet P wave (L=1) in coordinate space \param r [fm]
  virtual double GetVSr(double r) const; ///< [fm^-1/2] spin-singlet P wave (L=1) in coordinate space \param r [fm]

  virtual double GetUp(double p) const; ///< [MeV^-3/2] S wave u (L=0) in momentum space \param p [MeV]
  virtual double GetWp(double p) const; ///< [MeV^-3/2] D wave w (L=2) in momentum space \param p [MeV]
  virtual double GetVTp(double p) const; ///< [MeV^-3/2] spin-triplet P wave (L=1) in momentum space \param p [MeV]
  virtual double GetVSp(double p) const; ///< [MeV^-3/2] spin-singlet P wave (L=1) in momentum space \param p [MeV]
  virtual double GetUpoff(const nucqcd::Vector3& p) const; ///< [MeV^-3/2] off-shell S-wave part, see Deuteron::Wavefunction::DeuteronPStateOff \param p [MeV]
  virtual double GetWpoff1(const nucqcd::Vector3& p) const; ///< [MeV^-3/2] first off-shell D-wave part \param p [MeV]
  virtual double GetWpoff2(double pperp2) const; ///< [MeV^-3/2] second off-shell D-wave part \param pperp2 [MeV^2]

 protected:
  double Mass(int L, int i) const; ///< [MeV] mass m_{L,i} = alpha_L + (i-1) m0
  double GetK(int L, int term, int index) const; ///< dimensionless coefficient K^n_{L,i} that fixes the small-r behaviour
  double GetG(ESpace space, int L, int term, double q) const; ///< expansion function g_{L,n}(q) = f_{L,n} - sum_i K^n_{L,i} f_{L,N+i}(q); q in fm (r space) or MeV (p space)
  double GetF(ESpace space, int L, int index, double q) const; ///< basis function f_{L,i}(q) (Yukawa-type in r space); q in fm or MeV

 private:
  /// number of terms N
  int     fOrder; // Number of terms in parametrization
  /// [MeV] first mass parameter for each L
  double *fAlpha; //[3] first mass parameter for each L-wave [MeV]
  /// [MeV] second mass parameter
  double  fM0;    // second mass parameter [MeV]
  /// [MeV^1/2] parameters of the four waves
  double *fPar;   //[4*fOrder] parameters [MeV^1/2]


};


}  // namespace nucqcd

/** @} */
#endif
