/*! \file TYukawaPWF.h
 * \brief Deuteron wave functions as a sum of Yukawa functions (CD-Bonn, Paris, AV18 parametrizations)
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on February 5, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TYUKAWAPWF_H
#define TYUKAWAPWF_H

#include "TWavefunctionImplementation.h"

/*! \brief Non-relativistic deuteron wave function as a sum of N Yukawa terms with masses m_i
 *
 * \f$ u(r) = \sum_i C_i e^{-m_i r} \f$ and w(r) the corresponding D-wave sum with coefficients D_i
 * (the parametrization used for Paris, CD-Bonn and AV18), with masses \f$ m_i = \alpha + (i-1) m_0 \f$
 * or given explicitly. The P waves are zero.
 */
class TYukawaPWF : public TWavefunctionImplementation
{ 
 public:
  TYukawaPWF(int order, double alpha, double m0, double *c, double *d); ///< constructor: order N, masses m_i = alpha + (i-1) m0 with alpha, m0 [fm^-1], coefficients c (S wave), d (D wave) [fm^1/2]
  TYukawaPWF(int order, double *c, double *d, double *m); ///< constructor: order N, coefficients c, d [fm^1/2] and masses m [fm^-1]
  TYukawaPWF(const TYukawaPWF&); ///< copy constructor
  TYukawaPWF& operator=(const TYukawaPWF&); ///< assignment
  virtual ~TYukawaPWF(); ///< destructor

  virtual TYukawaPWF *Clone(const char* ="") const { return new TYukawaPWF(*this); } ///< virtual copy constructor

  virtual double GetUr(double r) const; ///< [fm^-1/2] S wave u (L=0) in coordinate space \param r [fm]
  virtual double GetWr(double r) const; ///< [fm^-1/2] D wave w (L=2) in coordinate space \param r [fm]
  virtual double GetVTr(double r) const; ///< [fm^-1/2] spin-triplet P wave (L=1) in coordinate space (always 0: non-relativistic) \param r [fm]
  virtual double GetVSr(double r) const; ///< [fm^-1/2] spin-singlet P wave (L=1) in coordinate space (always 0: non-relativistic) \param r [fm]

  virtual double GetUp(double p) const; ///< [MeV^-3/2] S wave u (L=0) in momentum space \param p [MeV]
  virtual double GetWp(double p) const; ///< [MeV^-3/2] D wave w (L=2) in momentum space \param p [MeV]
  virtual double GetVTp(double p) const; ///< [MeV^-3/2] spin-triplet P wave (L=1) in momentum space (always 0: non-relativistic) \param p [MeV]
  virtual double GetVSp(double p) const; ///< [MeV^-3/2] spin-singlet P wave (L=1) in momentum space (always 0: non-relativistic) \param p [MeV]
  virtual double GetUpoff(const TVector3& p) const; ///< [MeV^-3/2] off-shell S-wave part, see TDeuteron::Wavefunction::DeuteronPStateOff \param p [MeV]
  virtual double GetWpoff1(const TVector3& p) const; ///< [MeV^-3/2] first off-shell D-wave part \param p [MeV]
  virtual double GetWpoff2(double pperp2) const; ///< [MeV^-3/2] second off-shell D-wave part \param pperp2 [MeV^2]
  double getResidu() const; ///< [MeV^1/2] residue of the pole of the wave function (from the first term), for the on-shell extrapolation
  
 protected:
  double Mass(int i) const; ///< [fm^-1] mass m_i of term i (1..N)

 private:
  /// number of terms N
  int     fOrder;  // Number of terms in parametrization
  /// [fm^-1] mass parameter alpha
  double  fAlpha;  // mass parameter in parametrization [fm^-1]
  /// [fm^-1] mass parameter m0
  double  fM0;     // second mass parameter in parametrization [fm^-1]
  /// [fm^1/2] S-wave coefficients
  double *fC;      //[fOrder] parameters for s-wave [fm^1/2]
  /// [fm^1/2] D-wave coefficients
  double *fD;      //[fOrder] parameters for d-wave [fm^1/2]
  /// [fm^-1] masses (when given explicitly)
  double *fM;	//[fOrder] parameters for masses [fm^-1]	

  
};




/** @} */
#endif
