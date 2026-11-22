/*! \file TInterpolatingWavefunction.h
 * \brief Deuteron wave function interpolated from tables
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on February 16, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TINTERPOLATINGWAVEFUNCTION_H
#define TINTERPOLATINGWAVEFUNCTION_H

#include <map>
#include "TWavefunctionImplementation.h"
#include <complex>

/*! \brief Deuteron wave function interpolated linearly from tabulated values (e.g. Nijmegen III, or a
 * copy of another wave function on a grid, as in Poldeut and Deut_Conv_EMTLF)
 *
 * DeuteronPState and DeuteronRState have the same conventions as those of TDeuteron::Wavefunction
 * (D wave with a minus sign in momentum space, a plus sign in coordinate space, RState returns
 * r psi(r)), and the S and D waves can be switched off with setS and setD.
 */
class TInterpolatingWavefunction : public TWavefunctionImplementation
{
 public:
  TInterpolatingWavefunction(); ///< constructor (empty tables)
  TInterpolatingWavefunction(const TInterpolatingWavefunction&); ///< copy constructor
  TInterpolatingWavefunction& operator=(const TInterpolatingWavefunction&); ///< assignment
  virtual ~TInterpolatingWavefunction(); ///< destructor

  /// virtual copy constructor
  virtual TInterpolatingWavefunction *Clone(const char* ="") const
  { return new TInterpolatingWavefunction(*this); }

  void           SetUr(int n,double *r, double *list); ///< adds n tabulated values of the S wave u (L=0) in coordinate space: r [fm], values [fm^-1/2]
  void           SetWr(int n,double *r, double *list); ///< adds n tabulated values of the D wave w (L=2) in coordinate space: r [fm], values [fm^-1/2]
  void           SetVTr(int n,double *r, double *list); ///< adds n tabulated values of the spin-triplet P wave (L=1) in coordinate space: r [fm], values [fm^-1/2]
  void           SetVSr(int n,double *r, double *list); ///< adds n tabulated values of the spin-singlet P wave (L=1) in coordinate space: r [fm], values [fm^-1/2]
  
  void           SetUp(int n,double *p, double *list); ///< adds n tabulated values of the S wave u (L=0) in momentum space: p [MeV], values [MeV^-3/2]
  void           SetWp(int n,double *p, double *list); ///< adds n tabulated values of the D wave w (L=2) in momentum space: p [MeV], values [MeV^-3/2]
  void           SetVTp(int n,double *p, double *list); ///< adds n tabulated values of the spin-triplet P wave (L=1) in momentum space: p [MeV], values [MeV^-3/2]
  void           SetVSp(int n,double *p, double *list); ///< adds n tabulated values of the spin-singlet P wave (L=1) in momentum space: p [MeV], values [MeV^-3/2]

  void           AddUr(double r, double value); ///< adds one value of the S wave u (L=0) in coordinate space: r [fm], value [fm^-1/2]
  void           AddWr(double r, double value); ///< adds one value of the D wave w (L=2) in coordinate space: r [fm], value [fm^-1/2]
  void           AddVTr(double r, double value); ///< adds one value of the spin-triplet P wave (L=1) in coordinate space: r [fm], value [fm^-1/2]
  void           AddVSr(double r, double value); ///< adds one value of the spin-singlet P wave (L=1) in coordinate space: r [fm], value [fm^-1/2]
  
  void           AddUp(double p, double value); ///< adds one value of the S wave u (L=0) in momentum space: p [MeV], value [MeV^-3/2]
  void           AddWp(double p, double value); ///< adds one value of the D wave w (L=2) in momentum space: p [MeV], value [MeV^-3/2]
  void           AddVTp(double p, double value); ///< adds one value of the spin-triplet P wave (L=1) in momentum space: p [MeV], value [MeV^-3/2]
  void           AddVSp(double p, double value); ///< adds one value of the spin-singlet P wave (L=1) in momentum space: p [MeV], value [MeV^-3/2]

  void           RemoveUr(); ///< clears the table of the S wave u (L=0) in coordinate space
  void           RemoveWr(); ///< clears the table of the D wave w (L=2) in coordinate space
  void           RemoveVTr(); ///< clears the table of the spin-triplet P wave (L=1) in coordinate space
  void           RemoveVSr(); ///< clears the table of the spin-singlet P wave (L=1) in coordinate space

  void           RemoveUp(); ///< clears the table of the S wave u (L=0) in momentum space
  void           RemoveWp(); ///< clears the table of the D wave w (L=2) in momentum space
  void           RemoveVTp(); ///< clears the table of the spin-triplet P wave (L=1) in momentum space
  void           RemoveVSp(); ///< clears the table of the spin-singlet P wave (L=1) in momentum space

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
  virtual std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,int nucleon1Pol, const TVector3& p) const; ///< [MeV^-3/2] momentum-space wave function, same conventions and arguments as TDeuteron::Wavefunction::DeuteronPState
  virtual std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,int nucleon1Pol, const TVector3& r) const; ///< [fm^-1/2] r times the coordinate-space wave function, same conventions and arguments as TDeuteron::Wavefunction::DeuteronRState

  void setS(const bool S){ inclS=S;} ///< switches the S wave on (default) or off in DeuteronPState and DeuteronRState
  void setD(const bool D){ inclD=D;} ///< switches the D wave on (default) or off in DeuteronPState and DeuteronRState
  
  
 protected:
  double Interpolate(std::map<double,double>* list, double q) const; ///< linear interpolation in the ordered table list at q

 private:
  /// table of the S wave, r space
  std::map<double,double> *fUr; // list of L=0 wave in r-space
  /// table of the D wave, r space
  std::map<double,double> *fWr; // list of L=2 wave in r-space
  /// table of the triplet P wave, r space
  std::map<double,double> *fVTr; // list of L=1 triplet wave in r-space
  /// table of the singlet P wave, r space
  std::map<double,double> *fVSr; // list of L=1 singlet wave in r-space

  /// table of the S wave, p space
  std::map<double,double> *fUp; // list of L=0 wave in p-space
  /// table of the D wave, p space
  std::map<double,double> *fWp; // list of L=2 wave in p-space
  /// table of the triplet P wave, p space
  std::map<double,double> *fVTp; // list of L=1 triplet wave in p-space
  /// table of the singlet P wave, p space
  std::map<double,double> *fVSp; // list of L=1 singlet wave in p-space

  ClassDef(TInterpolatingWavefunction,1); // Interpolated wavefunction

  bool inclS; ///< [1] turn on S-wave [0] off
  bool inclD; ///< [1] turn on D-wave [0] off
  
};

/** @} */
#endif
