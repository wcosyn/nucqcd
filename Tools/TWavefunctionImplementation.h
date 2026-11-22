/*! \file TWavefunctionImplementation.h
 * \brief Interface of the deuteron wave function implementations
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on February 5, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TWAVEFUNCTIONIMPLEMENTATION_H
#define TWAVEFUNCTIONIMPLEMENTATION_H

#include <TVector3.h>

/*! \brief Interface of a deuteron wave function implementation (used by TDeuteron::Wavefunction)
 *
 * Radial functions with the conventions of TDeuteron::Wavefunction: reduced u, w (and the
 * relativistic P waves) in coordinate space [fm^-1/2], U, W in momentum space [MeV^-3/2]; see there
 * for the normalization and the sign of the D wave.
 */
class TWavefunctionImplementation
{
 public:
  static const double fgHbarc; ///< [MeV fm] hbar c

 public:
  TWavefunctionImplementation() = default; ///< default constructor (the derived classes used the former ROOT I/O constructor)
  virtual ~TWavefunctionImplementation(); ///< destructor

  // pure virtual copy constructor
  virtual TWavefunctionImplementation *Clone(const char* ="") const =0; ///< virtual copy constructor

  // L=0 wave in configuration space (pure virtual)
  virtual double GetUr(double r) const =0; ///< [fm^-1/2] S wave u (L=0) in coordinate space \param r [fm]

  // L=2 wave in configuration space (pure virtual)
  virtual double GetWr(double r) const =0; ///< [fm^-1/2] D wave w (L=2) in coordinate space \param r [fm]

  // L=1 triplet wave in configuration space (pure virtual)
  virtual double GetVTr(double r) const =0; ///< [fm^-1/2] spin-triplet P wave (L=1) in coordinate space \param r [fm]

  // L=1 singlet wave in configuration space (pure virtual)
  virtual double GetVSr(double r) const =0; ///< [fm^-1/2] spin-singlet P wave (L=1) in coordinate space \param r [fm]

  // L=0 wave in momentum space (pure virtual)
  virtual double GetUp(double p) const =0; ///< [MeV^-3/2] S wave u (L=0) in momentum space \param p [MeV]

  // L=2 wave in momentum space (pure virtual)
  virtual double GetWp(double p) const =0; ///< [MeV^-3/2] D wave w (L=2) in momentum space \param p [MeV]

  // L=1 triplet wave in momentum space (pure virtual)
  virtual double GetVTp(double p) const =0; ///< [MeV^-3/2] spin-triplet P wave (L=1) in momentum space \param p [MeV]

  // L=1 singlet wave in momentum space (pure virtual)
  virtual double GetVSp(double p) const =0; ///< [MeV^-3/2] spin-singlet P wave (L=1) in momentum space \param p [MeV]

  // L=0 off-shell wave in momentum space (pure virtual)
  virtual double GetUpoff(const TVector3& p) const =0; ///< [MeV^-3/2] off-shell S-wave part, see TDeuteron::Wavefunction::DeuteronPStateOff \param p [MeV]

  // L=2 off-shell wave in momentum space (pure virtual)
  virtual double GetWpoff1(const TVector3& p) const =0; ///< [MeV^-3/2] first off-shell D-wave part \param p [MeV]

  // L=2 off-shell wave in momentum space (pure virtual)
  virtual double GetWpoff2(double pperp2) const =0; ///< [MeV^-3/2] second off-shell D-wave part \param pperp2 [MeV^2]


  virtual double Radial_p(int l, double p) const; ///< [MeV^-3/2] U(p) for l=0, W(p) for l=2 \param l orbital angular momentum \param p [MeV]
  virtual double Radial_r(int l, double r) const; ///< [fm^-1/2] u(r) for l=0, w(r) for l=2 \param l orbital angular momentum \param r [fm]


};

/** @} */
#endif



