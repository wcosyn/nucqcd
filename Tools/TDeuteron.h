/*! \file TDeuteron.h
 * \brief Deuteron wave functions (TDeuteron::Wavefunction) and polarization states (TDeuteron::Polarization)
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on February 5, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TDEUTERON_H
#define TDEUTERON_H

#include "TWavefunctionImplementation.h"
//#include "TKinematics2to3.h"
#include <TNamed.h>
#include <TString.h>
#include "FourVector.h"
#include <complex>
#include "TSpinor.h"

// Forward declarations
//class TStrangePolarization;
class TLorentzVector;

/*! \brief Namespace-like class for the deuteron: wave functions and polarization states */
class TDeuteron
{
  ClassDef(TDeuteron,0); // Deuteron namespace

 public:
  // Forward declare the polarization class
  class Polarization;

  /*! \brief A deuteron wave function
   *
   * Create one with CreateWavefunction(name) or a named constructor; the caller owns it. Names
   * (case-insensitive): "cdbonn", "cdbonn87", "paris", "av18", "av18b", "nijm3" (non-relativistic),
   * "grossiib", "grossl4", "grosswjc1", "grosswjc2" (relativistic, Gross), "test". CD-Bonn,
   * CD-Bonn87, Paris, AV18 and AV18b are Yukawa-type parametrizations (TYukawaPWF).
   *
   * \par Radial functions and normalization
   * GetUr() = u(r) and GetWr() = w(r) are the \e reduced S- and D-wave functions [fm^-1/2] (the
   * radial wave functions are u(r)/r and w(r)/r), normalized as \f$ \int_0^\infty (u^2+w^2)\,dr = 1 \f$.
   * GetUp() = U(p) and GetWp() = W(p) are the momentum-space radial functions [MeV^-3/2], normalized
   * as \f$ \int_0^\infty p^2 (U^2+W^2)\,dp = 1 \f$, with
   * \f[ U(p) = \sqrt{2/\pi}\,(\hbar c)^{-3/2} \int_0^\infty u(r)\, j_0(pr/\hbar c)\, r\, dr, \qquad
   *     W(p) = \sqrt{2/\pi}\,(\hbar c)^{-3/2} \int_0^\infty w(r)\, j_2(pr/\hbar c)\, r\, dr, \f]
   * i.e. W(p) is the plain Bessel transform of w(r): it does \e not contain the factor
   * \f$ (-i)^2=-1 \f$ of the Fourier transform (checked numerically for Paris, AV18 and CD-Bonn).
   *
   * \par Sign of the D wave
   * w(r) has the same sign as u(r) at large r (positive asymptotic D/S ratio), so that the
   * quadrupole moment \f$ Q = \frac{1}{20}\int r^2 (\sqrt{8}\,uw - w^2)\,dr \f$ is positive; W(p) is
   * then positive at low p. Both are tested in unit_deuteron_wf.
   *
   * \par Spin states
   * DeuteronPState() is the momentum-space wave function
   * \f[ \psi^{m_D}_{m_1 m_2}(\vec p) = \sum_{L=0,2} (-i)^L R_L(p) \sum_{m_L,m_S}
   *     \langle L m_L\, 1 m_S | 1 m_D\rangle \langle \tfrac12 m_1\, \tfrac12 m_2 | 1 m_S\rangle Y_{L m_L}(\hat p), \f]
   * with \f$ R_0=U,\ R_2=W \f$: the D wave enters with a minus sign. In operator form (checked
   * numerically to machine precision)
   * \f[ \psi^{m_D}(\vec p) = \frac{1}{\sqrt{4\pi}} \left[ U(p) - \frac{W(p)}{\sqrt 8}\, S_{12}(\hat p) \right] \chi_{1 m_D},
   *     \qquad S_{12}(\hat p) = 3(\vec\sigma_1\cdot\hat p)(\vec\sigma_2\cdot\hat p) - \vec\sigma_1\cdot\vec\sigma_2, \f]
   * with W(p) > 0 at low p; in coordinate space the corresponding form is
   * \f$ [u(r) + w(r) S_{12}(\hat r)/\sqrt 8]\, \chi_{1 m_D}/(\sqrt{4\pi}\, r) \f$. It is normalized as
   * \f$ \sum_{m_1 m_2}\int d^3p\, |\psi^{m_D}_{m_1 m_2}|^2 = 1 \f$ [MeV^-3/2]. DeuteronRState() is the same
   * coupling in coordinate space with u(r) and w(r) and no phase factor (the D wave enters with a plus
   * sign); it returns \f$ r\,\psi(\vec r) \f$ [fm^-1/2], so \f$ \sum\int d\Omega\, dr\, |\cdot|^2 = 1 \f$.
   * Particle 1 is the nucleon that interacts with the photon, particle 2 the spectator; the coupling
   * order is \f$ \langle \tfrac12 m_1 \tfrac12 m_2|1 m_S\rangle \f$, but note the argument order
   * (deuteronPol, nucleon2Pol, nucleon1Pol). All spin projections are passed times two. The
   * quantization axis is z; Clebsch-Gordan coefficients and spherical harmonics use the
   * Condon-Shortley phase convention.
   */
  class Wavefunction : public TNamed
  {
  private:
    Wavefunction(const TString& name="", const TString& title=""); ///< private: use the named constructors
  public:
    Wavefunction(TRootIOCtor*); ///< ROOT I/O constructor
    Wavefunction(const Wavefunction&); ///< copy constructor (clones the implementation)
    Wavefunction& operator=(const Wavefunction&); ///< assignment (clones the implementation)
    virtual ~Wavefunction();

    // named constructors
    static Wavefunction *CreateCDBonnWavefunction(); ///< CD-Bonn (Yukawa parametrization); the caller owns it
    static Wavefunction *CreateCDBonn87Wavefunction(); ///< CD-Bonn 1987 (Yukawa parametrization); the caller owns it
    static Wavefunction *CreateParisWavefunction(); ///< Paris (Yukawa parametrization); the caller owns it
    static Wavefunction *CreateAV18Wavefunction(); ///< Argonne v18 (Yukawa parametrization); the caller owns it
    static Wavefunction *CreateAV18bWavefunction(); ///< Argonne v18, second parametrization (Yukawa); the caller owns it
    static Wavefunction *CreateGrossIIbWavefunction(); ///< Gross IIb (relativistic); the caller owns it
    static Wavefunction *CreateGrossL4Wavefunction(); ///< Gross L4 (relativistic); the caller owns it
    static Wavefunction *CreateGrossWJC1Wavefunction(); ///< Gross WJC-1 (relativistic); the caller owns it
    static Wavefunction *CreateGrossWJC2Wavefunction(); ///< Gross WJC-2 (relativistic); the caller owns it
    static Wavefunction *CreateNijmegenWavefunction(); ///< Nijmegen III, "nijm3" (interpolated tables); the caller owns it
    static Wavefunction *CreateTestWavefunction(); ///< test wave function (interpolated tables); the caller owns it
    /*! creates the wave function with the given name (see the class description); the caller owns it
     * \return the wave function \throw std::invalid_argument for an unknown name */
    static Wavefunction *CreateWavefunction(const TString& name);
  
    double GetUr(double r) const { return fImplementation->GetUr(r); } ///< [fm^-1/2] reduced S wave u(r) \param r [fm]
    double GetWr(double r) const { return fImplementation->GetWr(r); } ///< [fm^-1/2] reduced D wave w(r), same sign as u(r) at large r \param r [fm]
    double GetVTr(double r) const { return fImplementation->GetVTr(r); } ///< [fm^-1/2] spin-triplet P wave (relativistic wave functions only) \param r [fm]
    double GetVSr(double r) const { return fImplementation->GetVSr(r); } ///< [fm^-1/2] spin-singlet P wave (relativistic wave functions only) \param r [fm]
    double GetUp(double p) const { return fImplementation->GetUp(p); } ///< [MeV^-3/2] S wave U(p) in momentum space \param p [MeV]
    double GetWp(double p) const { return fImplementation->GetWp(p); } ///< [MeV^-3/2] D wave W(p) in momentum space, Bessel transform of w(r) without the factor -1 (see the class description) \param p [MeV]
    double GetUpoff(const TVector3& p) const { return fImplementation->GetUpoff(p); } ///< [MeV^-3/2] off-shell S-wave part, used in DeuteronPStateOff (Yukawa parametrizations) \param p [MeV] relative momentum
    double GetWpoff1(const TVector3& p) const { return fImplementation->GetWpoff1(p); } ///< [MeV^-3/2] first off-shell D-wave part, used in DeuteronPStateOff \param p [MeV] relative momentum
    double GetWpoff2(double pperp2) const { return fImplementation->GetWpoff2(pperp2); } ///< [MeV^-3/2] second off-shell D-wave part, used in DeuteronPStateOff \param pperp2 [MeV^2] transverse momentum squared
    double GetVTp(double p) const { return fImplementation->GetVTp(p); } ///< [MeV^-3/2] spin-triplet P wave in momentum space (relativistic only) \param p [MeV]
    double GetVSp(double p) const { return fImplementation->GetVSp(p); } ///< [MeV^-3/2] spin-singlet P wave in momentum space (relativistic only) \param p [MeV]
    double Radial_r(int l, double r) const { return fImplementation->Radial_r(l,r); } ///< [fm^-1/2] u(r) for l=0, w(r) for l=2 \param l orbital angular momentum \param r [fm]
    double Radial_p(int l, double p) const { return fImplementation->Radial_p(l,p); } ///< [MeV^-3/2] U(p) for l=0, W(p) for l=2 \param l orbital angular momentum \param p [MeV]
    double getResidu() const; ///< [MeV^1/2] residue of the pole of the wave function, for the on-shell extrapolation (Yukawa parametrizations only; throws std::logic_error for the others)
  
    bool   IsRelativistic() const { return fRelativistic; } ///< true for the relativistic (Gross) wave functions

    /*! covariant Blankenbecler-Cook d->np vertex with nucleon 2 on its mass shell, see Eq. (46) of
     * W. W. Buck and F. Gross, Phys. Rev. D 20, 2361 (1979)
     * \param p2 [MeV] four-momentum of nucleon 2 \param pD [MeV] four-momentum of the deuteron */
    FourVector<GammaStructure> Vertex(const FourVector<double>& p2,
				      const FourVector<double>& pD) const;
    /*! see Vertex(const FourVector<double>&, const FourVector<double>&) */
    FourVector<GammaStructure> Vertex(const TLorentzVector& p2,
				      const TLorentzVector& pD) const { return Vertex(ToFourVector(p2),ToFourVector(pD)); }
    /*! momentum-space wave function, see the class description for the phase and spin conventions
     * \param deuteronPol 2 m_D \param nucleon2Pol 2 m_2 (spectator) \param nucleon1Pol 2 m_1 (nucleon hit by the photon)
     * \param p [MeV] relative momentum (momentum of nucleon 1 in the deuteron rest frame)
     * \return [MeV^-3/2] amplitude */
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& p) const; 
    /*! DeuteronPState with the spatial part of p1 [MeV] */
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& p1) const; 
    /*! DeuteronPState with the spatial part of p1 [MeV] */
    std::complex<double> DeuteronPState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& p1) const; 
    /*! off-shell part of the momentum-space wave function (Yukawa parametrizations): the same spin
     * coupling with the off-shell radial parts GetUpoff, GetWpoff1, GetWpoff2, multiplied by
     * \f$ -i p_z/\hbar c \f$; zero for p_z = 0; not for relativistic wave functions
     * \param deuteronPol 2 m_D \param nucleon2Pol 2 m_2 (spectator) \param nucleon1Pol 2 m_1 (nucleon hit by the photon)
     * \param p [MeV] relative momentum */
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& p) const; 
    /*! DeuteronPStateOff with the spatial part of p1 [MeV] */
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& p1) const; 
    /*! with the spatial part of p1 [MeV]; note: this overload calls DeuteronPState, not DeuteronPStateOff */
    std::complex<double> DeuteronPStateOff(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TVector3& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const FourVector<double>& p1) const; 
//     std::complex<double> DeuteronPState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TLorentzVector& p1) const;
    /*! coordinate-space wave function times r, see the class description; not for relativistic wave functions
     * \param deuteronPol 2 m_D \param nucleon2Pol 2 m_2 (spectator) \param nucleon1Pol 2 m_1 (nucleon hit by the photon)
     * \param r [fm] relative position
     * \return [fm^-1/2] \f$ r\,\psi(\vec r) \f$ */
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TVector3& r) const; 
    /*! DeuteronRState with the spatial part of r1 [fm] */
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const FourVector<double>& r1) const; 
    /*! DeuteronRState with the spatial part of r1 [fm] */
    std::complex<double> DeuteronRState(int deuteronPol, int nucleon2Pol,
					int nucleon1Pol, 
					const TLorentzVector& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TVector3& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const FourVector<double>& r1) const; 
//     std::complex<double> DeuteronRState(const TStrangePolarization& pol,
// 				       const TSpinor::Polarization& pol1, 
// 				       const TLorentzVector& r1) const;
    static std::complex<double> SphericalHarmonic(unsigned int l, int m, double theta, double phi); ///< normalized spherical harmonic Y_lm (Condon-Shortley phase) \param l orbital angular momentum \param m magnetic quantum number \param theta [rad] \param phi [rad]
    static std::complex<double> SphericalHarmonicCos(unsigned int l, int m, double costheta, double phi); ///< as SphericalHarmonic, with cos(theta) as argument \param l orbital angular momentum \param m magnetic quantum number \param costheta cos of the polar angle \param phi [rad]
    static double ClebschGordan(int two_j1,int two_m1,int two_j2,int two_m2,int two_j,int two_m); ///< Clebsch-Gordan coefficient <j1 m1 j2 m2|j m> (Condon-Shortley), all arguments times two
    
  private:
    bool                         fRelativistic;   ///< relativistic wave function?
    TWavefunctionImplementation *fImplementation; ///< implementation of the wave function (owned)
    bool yukawa; ///< Yukawa-type parametrization (TYukawaPWF)?
    
    static const double          kMd;             ///< [MeV] deuteron mass used in Vertex

    ClassDef(Wavefunction,1); // Strangecalc wavefunction
    
  }; // class TDeuteron::Wavefunction


  /*! \brief Polarization state of a deuteron
   *
   * State m = -1, 0, +1 along the quantization axis \f$ \hat n = (-\sin\theta, 0, \cos\theta) \f$
   * (in the x-z plane); the angle phi only enters as a phase \f$ e^{-im\phi} \f$ of the
   * polarization vector, see GetPolarizationVector.
   */
  class Polarization : public TObject
  {
  public:   
    Polarization(double,double,int); ///< constructor: polar angle theta [rad] in [0,pi], phase angle phi [rad], state -1, 0 or +1
    Polarization(TRootIOCtor*); ///< ROOT I/O constructor
    Polarization(const Polarization&); ///< copy constructor
    virtual ~Polarization();
    Polarization& operator=(const Polarization&); ///< assignment

    bool operator==(const Polarization&); ///< same angles (within STRANGEUFLOW) and state
    operator int() const; ///< the state times two (2m)
    //TString HashName() const;

    void SetTheta(double theta); ///< sets the polar angle of the quantization axis \param theta [rad] in [0,pi]
    void SetPhi(double phi); ///< sets the phase angle \param phi [rad]
    void SetState(int polState); ///< sets the state \param polState -1, 0 or +1

    Polarization& operator++(); ///< raises the state by one (stays at +1)
    Polarization& operator--(); ///< lowers the state by one (stays at -1)

    FourVector<std::complex<double> > GetPolarizationVector() const; ///< polarization four-vector in the deuteron rest frame: spherical basis vector of state m around the axis n, times exp(-i m phi)

  private:
    double fTheta;             ///< [rad] polar angle of the quantization axis
    double fPhi;               ///< [rad] phase angle, see GetPolarizationVector
    int    fPolarizationState; ///< polarization state -1, 0 or +1

    ClassDef(Polarization,1); // Polarization state of a deuteron

  }; // class TDeuteron::Polarization

}; // class TDeuteron

/** @} */
#endif
