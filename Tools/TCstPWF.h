/*! \file TCstPWF.h
 * \brief Deuteron wave functions in the covariant spectator theory (CST)
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on July 23, 2010
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TCSTPWF_H
#define TCSTPWF_H

#include "TWavefunctionImplementation.h"

/*! \brief Covariant spectator theory (CST) deuteron wave function: S, D and the two P waves,
 * expanded in functions G with masses \f$ m_{z,i} = \alpha_z + (i-1) m^z_x \f$ for each wave z (see GetG)
 */
class TCstPWF : public TWavefunctionImplementation
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
  /*! constructor; the first five arguments are arrays of size 4 (u, w, vt, vs)
   * \param order number of terms N
   * \param leadingMass [MeV] leading mass alpha_z \param stepMass [MeV] step mass m_x \param tailMass [MeV] tail mass
   * \param exponent exponent n_z in the expansion functions G
   * \param parU,parW,parVT,parVS [MeV^-3/2] expansion coefficients of the S, D, triplet P and singlet P waves */
  TCstPWF(int *order, double *leadingMass, double *stepMass, double *tailMass,
	  int *exponent, double *parU, double *parW, double *parVT, double *parVS);
  TCstPWF(TRootIOCtor*); ///< ROOT I/O constructor
  TCstPWF(const TCstPWF&); ///< copy constructor
  TCstPWF& operator=(const TCstPWF&); ///< assignment
  virtual ~TCstPWF(); ///< destructor

  virtual TCstPWF *Clone(const char* ="") const {return new TCstPWF(*this); } ///< virtual copy constructor

  virtual double GetUr(double r) const; ///< [fm^-1/2] S wave u (L=0) in coordinate space \param r [fm]
  virtual double GetWr(double r) const; ///< [fm^-1/2] D wave w (L=2) in coordinate space \param r [fm]
  virtual double GetVTr(double r) const; ///< [fm^-1/2] spin-triplet P wave (L=1) in coordinate space \param r [fm]
  virtual double GetVSr(double r) const; ///< [fm^-1/2] spin-singlet P wave (L=1) in coordinate space \param r [fm]

  virtual double GetUp(double p) const; ///< [MeV^-3/2] S wave u (L=0) in momentum space \param p [MeV]
  virtual double GetWp(double p) const; ///< [MeV^-3/2] D wave w (L=2) in momentum space \param p [MeV]
  virtual double GetVTp(double p) const; ///< [MeV^-3/2] spin-triplet P wave (L=1) in momentum space \param p [MeV]
  virtual double GetVSp(double p) const; ///< [MeV^-3/2] spin-singlet P wave (L=1) in momentum space \param p [MeV]
  virtual double GetUpoff(const TVector3& p) const; ///< [MeV^-3/2] off-shell S-wave part, see TDeuteron::Wavefunction::DeuteronPStateOff \param p [MeV]
  virtual double GetWpoff1(const TVector3& p) const; ///< [MeV^-3/2] first off-shell D-wave part \param p [MeV]
  virtual double GetWpoff2(double pperp2) const; ///< [MeV^-3/2] second off-shell D-wave part \param pperp2 [MeV^2]

 protected:
  static int AngularMomentum(EWave wave); ///< orbital angular momentum of the wave (u: 0, vt and vs: 1, w: 2)
  double Mass(EWave wave, int i) const; ///< [MeV] mass m_{z,i} = alpha_z + (i-1) m_x of wave z
  double GetG(ESpace space, EWave wave, int term, double q) const; ///< expansion function G of wave z, term i; q in fm (r space) or MeV (p space)
 
 private:
  /// number of terms for each wave
  int     *fOrder; //[4] Number of terms in parametrization
  /// [MeV^-3/2] expansion coefficients for each wave
  double **fCoeff; //[4][fOrder[i]] expansion coefficients [MeV^-3/2]
  /// [MeV] leading mass for each wave
  double  *fAlpha; //[4] leading mass [MeV]
  /// [MeV] step mass for each wave
  double  *fMx;    //[4] step mass [MeV]
  /// [MeV] tail mass for each wave
  double  *fMn;    //[4] tail mass [MeV]
  /// exponent in the expansion functions for each wave
  int     *fExp;   //[4] exponent in expansion functions
  
  ClassDef(TCstPWF,1); //  Covariant spectator theory parametrized wavefunction

}; // class TCstPWF

/** @} */
#endif // TCSTPWF_H
