/*! \file TKinematics2to3.h
 * \brief Kinematics of the 2 -> 3 reaction gamma + D -> K + Y + N
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * \addtogroup MePhys
 * @{
 */


#ifndef TKINEMATICS2TO3_H
#define TKINEMATICS2TO3_H

#include <cmath>
#include <string>
#include "TSpinor.h"
#include "TKinematics2to2.h"
#include "FourVector.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

// Forward declarations
class TString;
class TVector3;
class TLorentzQuaternion;

// Global functions for conversion 
// from ROOT 4vectors to strangecalc 4vectors.
// template<typename T>
// FourVector<T> operator*(const TLorentzRotation&, const FourVector<T>&);
// 
// FourVector<double> ToFourVector(const TLorentzVector&);
// TLorentzVector ToLorentzVector(const FourVector<double>&);

/*! \brief Kinematics of gamma + D -> K + Y + N, for a grid of six independent variables
 *
 * The names follow strangeness photoproduction on the deuteron: photon (G), deuteron (D), kaon (K),
 * hyperon (Y) and nucleon (N). With isospin channel 7 the masses are set by the user, and the class
 * is used for other 2 -> 3 reactions (for example A(e,e'pN), with N the proton, K the other nucleon
 * and Y the residual nucleus, see TKinematics2to3WithLabAngles).
 *
 * The lab frame is the deuteron rest frame. The z axis is always along the photon momentum; one of
 * the final-state particles fixes the x and y axes, chosen with EAxisConvention (internally the y
 * axis is along p_gamma x p_N). Energies and momenta are in MeV, Q^2 and Mandelstam variables in
 * MeV^2, angles in rad unless the getter name ends in _deg. The independent variables are given by
 * the format string (see SetFormat); each can be given a range (SetVarRange) and the grid is
 * stepped through with Next() or GoTo().
 */
class TKinematics2to3
{
  std::string fName;  ///< name (formerly from the ROOT base class TNamed)
  std::string fTitle; ///< title

 public:
  const char *GetName() const { return fName.c_str(); }  ///< name
  const char *GetTitle() const { return fTitle.c_str(); } ///< title
  void SetName(const char *name) { fName = name; }        ///< sets the name
  void SetTitle(const char *title) { fTitle = title; }    ///< sets the title

  friend class TKinematics2to3WithLabAngles;

 public:
  /*! axis conventions: z along the photon momentum, x = y x z, with the y axis along
   * p_gamma x p_K (kK), p_gamma x p_Y (kY) or p_gamma x p_N (kN, the internal one); kS is the
   * standard input/output convention stored in fAC */
  enum EAxisConvention { 
    // code to specify the axis convention.
    // z ~ p_gamma and x ~ y x z
    kS = 0, // standard input/output axis convention stored in fAC
    kK = 1, // y ~ p_gamma x p_K
    kY = 2, // y ~ p_gamma x p_Y
    kN = 3  // y ~ p_gamma x p_N (internal convention)
  };

  /*! constructor
   * \param name,title ROOT name and title
   * \param isospin isospin channel: (1) K+ Lambda n, (2) K+ Sigma0 n, (3) K0 Sigma+ n, (4) K0 Lambda p,
   * (5) K0 Sigma0 p, (6) K+ Sigma- p, (7) masses given by Md, Mn, Mk, My
   * \param ac standard axis convention (kK, kY or kN; kS is not allowed)
   * \param format format string with the six independent variables, see SetFormat
   * \param var1,var2,var3,var4,var5,var6 initial values of the
   * independent variables, in the order of the format string
   * \param Md,Mn,Mk,My [MeV] masses of deuteron, nucleon, kaon and hyperon (isospin channel 7)
   */
  TKinematics2to3(const char *name, const char *title,
		  int isospin, EAxisConvention ac, 
		  const char *format,
		  double var1, double var2, double var3,
		  double var4, double var5, double var6,
		  double Md=0., double Mn=0., double Mk=0., double My=0.);
  TKinematics2to3(const TKinematics2to3&); ///< Copy constructor
  virtual ~TKinematics2to3();
  TKinematics2to3& operator=(const TKinematics2to3&); ///< assignment

  void           SetVarRange(int varNum,double low,double high,int steps); ///< Give independent variable number 'varNumber' a range in the interval [lowerLimit,upperLimit].
  void           FixVariable(int varNum); ///< Unset the range of variable number 'varNumber'. Set it to its current value.
  void           FixVariables(); ///< Calls FixVariable(int) on all variables.
  int            Next(); ///< Proceed to the next step in the grid of independent variables. Return value is zero when the end of the grid is reached.
  void           GoTo(int totalStep); ///< Proceed to the 'step'th point in the grid of independent variables.
  /*! proceeds to the point in the grid of independent variables given by the step of each variable
   * \param var1step,var2step,var3step,var4step,var5step,var6step step of variable 1-6 */
  void           GoTo(int var1step,int var2step,int var3step,
		      int var4step,int var5step,int var6step);
  int            VariableInfo() const; ///< Print information on the type of independent variables and their status.
  static void    Help(); ///< prints help on the format strings and variables
  bool           operator==(const TKinematics2to3&) const; ///< same masses, format and values of the independent variables
//   TString        HashName() const;

  void           SetIsospin(int isospin); ///< changes the isospin channel and updates the kinematics \param isospin channel, see the constructor
  void           SetFormat(const char *format); ///< sets the format string (the six independent variables), see ReadFormatString \param format format string
  void           SetVar(int varNum,double value); ///< sets independent variable varNum (1-6) to value
  void           SetVar(double,double,double,double,double,double); ///< sets the six independent variables, in the order of the format string
  
  const char     *GetFormat() const         { return fFormat; } ///< format string that lists the six independent variables, see SetFormat
  EAxisConvention GetAxisConvention() const { return fAC; } ///< the standard axis convention (fAC)
  bool           IsPhysical() const     { return fIsPhysical; } ///< does the current kinematic point lie in the physical region?
  int            GetIsospin() const     { return fIsospin; } ///< isospin channel (1-6, or 7 for masses set by the user)
  double         GetMk() const          { return fMk; } ///< [MeV] mass of the kaon (K)
  double         GetMy() const          { return fMy; } ///< [MeV] mass of the hyperon (Y)
  double         GetMn() const          { return fMn; } ///< [MeV] mass of the nucleon (N)
  double         GetMd() const          { return fMd; } ///< [MeV] mass of the deuteron
  TVector3       GetK3(EAxisConvention =kS) const; ///< [MeV] momentum three-vector of the kaon (K) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TVector3       GetY3(EAxisConvention =kS) const; ///< [MeV] momentum three-vector of the hyperon (Y) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TVector3       GetN3(EAxisConvention =kS) const; ///< [MeV] momentum three-vector of the nucleon (N) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TVector3       GetG3(EAxisConvention =kS) const; ///< [MeV] momentum three-vector of the photon in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TVector3       GetD3(EAxisConvention =kS) const; ///< [MeV] momentum three-vector of the deuteron in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector GetK4(EAxisConvention =kS) const; ///< [MeV] four-momentum of the kaon (K) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector GetY4(EAxisConvention =kS) const; ///< [MeV] four-momentum of the hyperon (Y) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector GetN4(EAxisConvention =kS) const; ///< [MeV] four-momentum of the nucleon (N) in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector GetG4(EAxisConvention =kS) const; ///< [MeV] four-momentum of the photon in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector GetD4(EAxisConvention =kS) const; ///< [MeV] four-momentum of the deuteron in the lab frame (deuteron rest frame); the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPk() const          { return fKvector->P(); } ///< [MeV] lab momentum of the kaon (K)
  double         GetPy() const          { return fYvector->P(); } ///< [MeV] lab momentum of the hyperon (Y)
  double         GetPn() const          { return fNvector->P(); } ///< [MeV] lab momentum of the nucleon (N)
  double         GetKlab() const        { return fGvector->P(); } ///< [MeV] lab momentum of the photon
  double         GetEk() const          { return fKvector->E(); } ///< [MeV] lab energy of the kaon (K)
  double         GetEy() const          { return fYvector->E(); } ///< [MeV] lab energy of the hyperon (Y)
  double         GetEn() const          { return fNvector->E(); } ///< [MeV] lab energy of the nucleon (N)
  double         GetWlab() const        { return fGvector->E(); } ///< [MeV] lab energy of the photon
  double         GetThk_deg() const; ///< polar angle of the kaon (K) with respect to the photon, in the lab [deg, 0-180]
  double         GetThy_deg() const; ///< polar angle of the hyperon (Y) with respect to the photon, in the lab [deg, 0-180]
  double         GetThn_deg() const; ///< polar angle of the nucleon (N) with respect to the photon, in the lab [deg, 0-180]
  double         GetThk() const; ///< polar angle of the kaon (K) with respect to the photon, in the lab [rad, 0-pi]
  double         GetThy() const; ///< polar angle of the hyperon (Y) with respect to the photon, in the lab [rad, 0-pi]
  double         GetThn() const; ///< polar angle of the nucleon (N) with respect to the photon, in the lab [rad, 0-pi]
  double         GetCosthk() const      { return fKvector->CosTheta(); } ///< cosine of the polar angle of the kaon (K) with respect to the photon, in the lab
  double         GetCosthy() const      { return fYvector->CosTheta(); } ///< cosine of the polar angle of the hyperon (Y) with respect to the photon, in the lab
  double         GetCosthn() const      { return fNvector->CosTheta(); } ///< cosine of the polar angle of the nucleon (N) with respect to the photon, in the lab
  double         GetPhik_deg(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the lab [deg, 0-360[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPhiy_deg(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the lab [deg, 0-360[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPhin_deg(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the lab [deg, 0-360[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPhik(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the lab [rad, 0-2pi[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPhiy(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the lab [rad, 0-2pi[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetPhin(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the lab [rad, 0-2pi[; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetThky_deg() const; ///< angle between the photon and the total momentum of the K-Y pair, in the lab [deg, 0-180]
  double         GetThkn_deg() const; ///< angle between the photon and the total momentum of the K-N pair, in the lab [deg, 0-180]
  double         GetThyn_deg() const; ///< angle between the photon and the total momentum of the Y-N pair, in the lab [deg, 0-180]
  double         GetThky() const; ///< angle between the photon and the total momentum of the K-Y pair, in the lab [rad, 0-pi]
  double         GetThkn() const; ///< angle between the photon and the total momentum of the K-N pair, in the lab [rad, 0-pi]
  double         GetThyn() const; ///< angle between the photon and the total momentum of the Y-N pair, in the lab [rad, 0-pi]
  double         GetCosthky() const     { return std::cos(GetThky()); } ///< cosine of the angle between the photon and the total momentum of the K-Y pair
  double         GetCosthkn() const     { return std::cos(GetThkn()); } ///< cosine of the angle between the photon and the total momentum of the K-N pair
  double         GetCosthyn() const     { return std::cos(GetThyn()); } ///< cosine of the angle between the photon and the total momentum of the Y-N pair
  double         GetQsquared() const    { return -1*fGvector->Mag2(); } ///< [MeV^2] Q^2 = -q^2 of the photon
  double         GetWtot() const        { return std::sqrt(GetStot()); } ///< [MeV] total invariant mass W
  double         GetWky() const         { return fKYsystem->GetW(); } ///< [MeV] invariant mass of the K-Y subsystem
  double         GetWkn() const         { return fKNsystem->GetW(); } ///< [MeV] invariant mass of the K-N subsystem
  double         GetWyn() const         { return fYNsystem->GetW(); } ///< [MeV] invariant mass of the Y-N subsystem
  double         GetStot() const; ///< [MeV^2] Mandelstam s = (q + P_D)^2 = M_D^2 + 2 M_D omega - Q^2
  double         GetSky() const         { return fKYsystem->GetS(); } ///< [MeV^2] Mandelstam s of the K-Y subsystem
  double         GetSkn() const         { return fKNsystem->GetS(); } ///< [MeV^2] Mandelstam s of the K-N subsystem
  double         GetSyn() const         { return fYNsystem->GetS(); } ///< [MeV^2] Mandelstam s of the Y-N subsystem
  double         GetTkg() const; ///< [MeV^2] momentum transfer t = (p_K - p_gamma)^2 between the photon and the kaon (K)
  double         GetTyg() const; ///< [MeV^2] momentum transfer t = (p_Y - p_gamma)^2 between the photon and the hyperon (Y)
  double         GetTng() const; ///< [MeV^2] momentum transfer t = (p_N - p_gamma)^2 between the photon and the nucleon (N)
  double         GetTkd() const; ///< [MeV^2] momentum transfer t = (p_K - p_D)^2 between the deuteron and the kaon (K)
  double         GetTyd() const; ///< [MeV^2] momentum transfer t = (p_Y - p_D)^2 between the deuteron and the hyperon (Y)
  double         GetTnd() const; ///< [MeV^2] momentum transfer t = (p_N - p_D)^2 between the deuteron and the nucleon (N)
  double         GetKYsystem_Pk() const { return fKYsystem->GetPk(); } ///< [MeV] momentum of the kaon (K) in the KY center-of-mass frame
  double         GetKYsystem_Py() const { return fKYsystem->GetPk(); } ///< [MeV] momentum of the hyperon (Y) in the KY center-of-mass frame
  double         GetKYsystem_Ek() const; ///< [MeV] energy of the kaon (K) in the KY center-of-mass frame
  double         GetKYsystem_Ey() const; ///< [MeV] energy of the hyperon (Y) in the KY center-of-mass frame
  double         GetKYsystem_S() const  { return fKYsystem->GetS(); } ///< [MeV^2] Mandelstam s of the KY subsystem
  double         GetKYsystem_W() const  { return fKYsystem->GetW(); } ///< [MeV] invariant mass of the KY subsystem
  double         GetKYsystem_costhkcm() const; ///< cosine of the polar angle of the kaon (K) in the KY center-of-mass frame
  double         GetKYsystem_costhycm() const; ///< cosine of the polar angle of the hyperon (Y) in the KY center-of-mass frame
  double         GetKYsystem_thkcm() const; ///< polar angle of the kaon (K) in the KY center-of-mass frame [rad]
  double         GetKYsystem_thycm() const; ///< polar angle of the hyperon (Y) in the KY center-of-mass frame [rad]
  double         GetKYsystem_phikcm(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the KY center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKYsystem_phiycm(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the KY center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKYsystem_phikcm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the KY center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKYsystem_phiycm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the KY center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKNsystem_Pk() const { return fKNsystem->GetPk(); } ///< [MeV] momentum of the kaon (K) in the KN center-of-mass frame
  double         GetKNsystem_Pn() const { return fKNsystem->GetPk(); } ///< [MeV] momentum of the nucleon (N) in the KN center-of-mass frame
  double         GetKNsystem_Ek() const; ///< [MeV] energy of the kaon (K) in the KN center-of-mass frame
  double         GetKNsystem_En() const; ///< [MeV] energy of the nucleon (N) in the KN center-of-mass frame
  double         GetKNsystem_S() const  { return fKNsystem->GetS(); } ///< [MeV^2] Mandelstam s of the KN subsystem
  double         GetKNsystem_W() const  { return fKNsystem->GetW(); } ///< [MeV] invariant mass of the KN subsystem
  double         GetKNsystem_costhkcm() const; ///< cosine of the polar angle of the kaon (K) in the KN center-of-mass frame
  double         GetKNsystem_costhncm() const; ///< cosine of the polar angle of the nucleon (N) in the KN center-of-mass frame
  double         GetKNsystem_thkcm() const; ///< polar angle of the kaon (K) in the KN center-of-mass frame [rad]
  double         GetKNsystem_thncm() const; ///< polar angle of the nucleon (N) in the KN center-of-mass frame [rad]
  double         GetKNsystem_phikcm(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the KN center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKNsystem_phincm(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the KN center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKNsystem_phikcm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the kaon (K) in the KN center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetKNsystem_phincm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the KN center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetYNsystem_Py() const { return fYNsystem->GetPk(); } ///< [MeV] momentum of the hyperon (Y) in the YN center-of-mass frame
  double         GetYNsystem_Pn() const { return fYNsystem->GetPk(); } ///< [MeV] momentum of the nucleon (N) in the YN center-of-mass frame
  double         GetYNsystem_Ey() const; ///< [MeV] energy of the hyperon (Y) in the YN center-of-mass frame
  double         GetYNsystem_En() const; ///< [MeV] energy of the nucleon (N) in the YN center-of-mass frame
  double         GetYNsystem_S() const  { return fYNsystem->GetS(); } ///< [MeV^2] Mandelstam s of the YN subsystem
  double         GetYNsystem_W() const  { return fYNsystem->GetW(); } ///< [MeV] invariant mass of the YN subsystem
  double         GetYNsystem_costhycm() const; ///< cosine of the polar angle of the hyperon (Y) in the YN center-of-mass frame
  double         GetYNsystem_costhncm() const;  ///< cosine of the polar angle of the nucleon (N) in the YN center-of-mass frame
  double         GetYNsystem_thycm() const; ///< polar angle of the hyperon (Y) in the YN center-of-mass frame [rad]
  double         GetYNsystem_thncm() const; ///< polar angle of the nucleon (N) in the YN center-of-mass frame [rad]
  double         GetYNsystem_phiycm(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the YN center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetYNsystem_phincm(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the YN center-of-mass frame [rad]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetYNsystem_phiycm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the hyperon (Y) in the YN center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetYNsystem_phincm_deg(EAxisConvention =kS) const; ///< azimuthal angle of the nucleon (N) in the YN center-of-mass frame [deg]; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostToKYcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the LAB frame with the axis convention ac to the CM frame of the kaon-hyperon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostToKNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the LAB frame with the axis convention ac to the CM frame of the kaon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostToYNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the LAB frame with the axis convention ac to the CM frame of the hyperon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostToGNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the LAB frame with the axis convention ac to the CM frame of the photon and the struck nucleon inside the deuteron. The z-axis lies along the photon and the y-axis is defined as being perpendicular to the gN and KY plane. This is the reference frame where the matrix elements get evaluated in strangecalc.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostToKYcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the kaon-hyperon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostToKNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the kaon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostToYNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the hyperon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostToGNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the photon and the struck nucleon inside the deuteron. The z-axis lies along the photon and the y-axis is defined as being perpendicular to the gN and KY plane. This is the reference frame where the matrix elements get evaluated in strangecalc.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostFromKYcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the CM frame of the kaon-hyperon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostFromKNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the CM frame of the kaon-nucleon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostFromYNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the CM frame of the hyperon-nucleon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector BoostFromGNcm(const TLorentzVector&, EAxisConvention =kS) const; ///< Boost the 4vector 'vector' defined in the CM frame of the photon and the struck nucleon inside the deuteron to the LAB frame with the axis convention ac.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostFromKYcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the kaon-hyperon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostFromKNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the kaon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostFromYNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the hyperon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation BoostFromGNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the photon and the struck nucleon inside the deuteron to the LAB frame with the axis convention ac.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostToKYcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the LAB frame with the axis convention ac to the CM frame of the kaon-hyperon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostToKNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the LAB frame with the axis convention ac to the CM frame of the kaon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostToYNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the LAB frame with the axis convention ac to the CM frame of the hyperon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostToGNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the LAB frame with the axis convention ac to the CM frame of the photon and the struck nucleon inside the deuteron. The z-axis lies along the photon and the y-axis is defined as being perpendicular to the gN and KY plane. This is the reference frame where the matrix elements get evaluated in strangecalc.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorToKYcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the LAB frame with the axis convention ac to the CM frame of the kaon-hyperon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorToKNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the LAB frame with the axis convention ac to the CM frame of the kaon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorToYNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the LAB frame with the axis convention ac to the CM frame of the hyperon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorToGNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the LAB frame with the axis convention ac to the CM frame of the photon and the struck nucleon inside the deuteron. The z-axis lies along the photon and the y-axis is defined as being perpendicular to the gN and KY plane. This is the reference frame where the matrix elements get evaluated in strangecalc.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostFromKYcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the CM frame of the kaon-hyperon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostFromKNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the CM frame of the kaon-nucleon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostFromYNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the CM frame of the hyperon-nucleon subsystem with the axis convention ac to the LAB frame .; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor        BoostFromGNcm(const TSpinor&, EAxisConvention =kS) const; ///< Boost the Dirac spinor 'spinor' defined in the CM frame of the photon and the struck nucleon inside the deuteron to the LAB frame with the axis convention ac.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorFromKYcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the CM frame of the kaon-hyperon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorFromKNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the CM frame of the kaon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorFromYNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the CM frame of the hyperon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation BoostSpinorFromGNcm(EAxisConvention =kS) const; ///< Object to boost a Dirac spinor defined in the CM frame of the photon and the struck nucleon inside the deuteron to the LAB frame with the axis convention ac.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionToKYcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the kaon-hyperon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionToKNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the kaon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionToYNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the hyperon-nucleon subsystem.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionToGNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the LAB frame with the axis convention ac to the CM frame of the photon and the struck nucleon inside the deuteron. The z-axis lies along the photon and the y-axis is defined as being perpendicular to the gN and KY plane. This is the reference frame where the matrix elements get evaluated in strangecalc.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionFromKYcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the kaon-hyperon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionFromKNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the kaon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionFromYNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the hyperon-nucleon subsystem with the axis convention ac to the LAB frame.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzQuaternion BoostQuaternionFromGNcm(EAxisConvention =kS) const; ///< Object to boost a 4vector defined in the CM frame of the photon and the struck nucleon inside the deuteron to the LAB frame with the axis convention ac.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation RotateToKac(EAxisConvention =kS) const; ///< Object to rotate a 4vector defined in the LAB frame from the axis convention ac to the 'kK' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation RotateToYac(EAxisConvention =kS) const; ///< Object to rotate a 4vector defined in the LAB frame from the axis convention ac to the 'kY' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzRotation RotateToNac(EAxisConvention =kS) const; ///< Object to rotate a 4vector defined in the LAB frame from the axis convention ac to the 'kN' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector RotateToKac(const TLorentzVector&, EAxisConvention =kS) const; ///< Rotate the 4vector 'vector' defined in the LAB frame from the axis convention ac to the 'kK' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector RotateToYac(const TLorentzVector&, EAxisConvention =kS) const; ///< Rotate the 4vector 'vector' defined in the LAB frame from the axis convention ac to the 'kY' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TLorentzVector RotateToNac(const TLorentzVector&, EAxisConvention =kS) const; ///< Rotate the 4vector 'vector' defined in the LAB frame from the axis convention ac to the 'kN' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation RotateSpinorToKac(EAxisConvention =kS) const; ///< Object to rotate a Dirac spinor defined in the LAB frame from the axis convention ac to the 'kK' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation RotateSpinorToYac(EAxisConvention =kS) const; ///< Object to rotate a Dirac spinor defined in the LAB frame from the axis convention ac to the 'kY' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor::LorentzRotation RotateSpinorToNac(EAxisConvention =kS) const; ///< Object to rotate a Dirac spinor defined in the LAB frame from the axis convention ac to the 'kN' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor RotateToKac(const TSpinor&, EAxisConvention =kS) const; ///< Rotate the Dirac spinor 'spinor' defined in the LAB frame from the axis convention ac to the 'kK' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor RotateToYac(const TSpinor&, EAxisConvention =kS) const; ///< Rotate the Dirac spinor 'spinor' defined in the LAB frame from the axis convention ac to the 'kY' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  TSpinor RotateToNac(const TSpinor&, EAxisConvention =kS) const; ///< Rotate the Dirac spinor 'spinor' defined in the LAB frame from the axis convention ac to the 'kN' axis convention.; the argument is the axis convention (default kS: the standard one, fAC)
  double         GetVar(int) const; ///< Returns the value of the independent variable number 'varNumber' as specified in the format string.
  double         GetVar(const TString& varName) const; ///< Returns the value of the variable given by name. The variable can effectively by any one for which their exists a getter. It need not be listed in the info for SetFormat(const char*), eg. kyphikcm_deg.
  virtual double *GetVarArray(const TString& varName) const; ///< Return an array of the variable of type 'variable'. The variable can effectively by any one for which their exists a getter. It need not be listed in the info for SetFormat(const char*), eg. kyphikcm_deg.
  virtual double *GetPhysicalVarArray(const TString& varName) const; ///< Return an array of the variable of type 'variable'. The variable can effectively by any one for which their exists a getter. It need not be listed in the info for SetFormat(const char*), eg. kyphikcm_deg.
  TString        GetVarName(int varNum) const; ///< name of independent variable varNum (1-6) \param varNum variable number
  int            GetNumberOfSteps(int varNum) const; ///< number of steps of independent variable varNum
  int            GetNumberOfSteps() const; ///< number of steps of independent variable varNum
  virtual int    GetNumberOfPhysicalSteps() const; ///< number of steps in the grid where the kinematics are physical
  int            GetNumberOfVariables() const; ///< Determine how many independent variables are not fixed
  double         GetStepSize(int varNum) const; ///< step size of independent variable varNum
  int            GetStep(int varNum) const; ///< current step of independent variable varNum
  int            GetStep() const { return fStep; } ///< current step in the grid of independent variables
  bool           IsFixed() const; ///< are all independent variables fixed?
  bool           IsFixed(int varNum) const; ///< is independent variable varNum fixed (no range)?

  double GetPhotonThreshold() const; ///< Calculate the threshold photon lab energy for the current isospin channel.
  void   GetBoundsMomentumParticle1(double th1, double& lower, double& upper) const; ///< For the current invariant mass 'w' and a given scattering angle 'th1' of the particle specified by the axis convention, the lab momentum |p_1| has an upper and lower bound.
  double GetLowerBoundMomentumParticle1(double th1) const; ///< For the current invariant mass 'w' and a given scattering angle 'th1' of the particle specified by the axis convention, the lab momentum |p_1| has a lower bound.
  double GetUpperBoundMomentumParticle1(double th1) const; ///< For the current invariant mass 'w' and a given scattering angle 'th1' of the particle specified by the axis convention, the lab momentum |p_1| has an upper bound.
  double GetMaximumScatterAngleParticle1() const; ///< For the current invariant mass 'w' and lab momentum 'p_1' of the particle specified by the axis convention, the lab scattering angle 'theta_1' has a maximum.

 protected:
  virtual void   UpdateKinematics(); ///< Solve the kinematics!
  void           SetMd(double md) { fMd = md; } ///< sets the deuteron mass \param md [MeV]
  void           SetMk(double mk) { fMk = mk; } ///< sets the kaon mass \param mk [MeV]
  void           SetMy(double my) { fMy = my; } ///< sets the hyperon mass \param my [MeV]
  void           SetMn(double mn) { fMn = mn; } ///< sets the nucleon mass \param mn [MeV]

 private:
  virtual void   ReadFormatString(const char*); ///< Process the format string.
  void           ChangeIsospinChannel(int); ///< Set the particle masses according to the isospin channel
  void           InitializeSubsystems(); ///< Create and set the KY-, KN- and YN-subsystems.
  void           UpdateSubsystems(); ///< Recalculate the KY-, KN- and YN-subsystems.
  void           UpdateVariables(); ///< After UpdateKinematics() not all data member are up-to-date. This function makes sure they all are.
  double         GetVarByGetter(const TString& variable) const; ///< Return the value of the variable of type 'variable'. The variable can effectively by any one for which their exists a getter. It need not be listed in the info for SetFormat(const char*), eg. kyphikcm_deg.
  void           SetAxisConvention(EAxisConvention); ///< Change the axis convention. Be sure you know what you're doing before you call this function. This function doesn't check whether what you are doing makes sense.; the argument is the axis convention (default kS: the standard one, fAC)
  const char    *GetSolutionType() const { return fKinSolution; } ///< how the kinematics are solved, set by ReadFormatString

  // Wrap kinematics getters that take an argument
  double         GetPhik_deg_kS() const            { return GetPhik_deg(); } ///< GetPhik_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetPhiy_deg_kS() const            { return GetPhiy_deg(); } ///< GetPhiy_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetPhin_deg_kS() const            { return GetPhin_deg(); } ///< GetPhin_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetPhik_kS() const                { return GetPhik(); } ///< GetPhik with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetPhiy_kS() const                { return GetPhiy(); } ///< GetPhiy with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetPhin_kS() const                { return GetPhin(); } ///< GetPhin with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKYsystem_phikcm_kS() const     { return GetKYsystem_phikcm(); } ///< GetKYsystem_phikcm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKYsystem_phiycm_kS() const     { return GetKYsystem_phiycm(); } ///< GetKYsystem_phiycm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKYsystem_phikcm_deg_kS() const { return GetKYsystem_phikcm_deg(); } ///< GetKYsystem_phikcm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKYsystem_phiycm_deg_kS() const { return GetKYsystem_phiycm_deg(); } ///< GetKYsystem_phiycm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKNsystem_phikcm_kS() const     { return GetKNsystem_phikcm(); } ///< GetKNsystem_phikcm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKNsystem_phincm_kS() const     { return GetKNsystem_phincm(); } ///< GetKNsystem_phincm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKNsystem_phikcm_deg_kS() const { return GetKNsystem_phikcm_deg(); } ///< GetKNsystem_phikcm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetKNsystem_phincm_deg_kS() const { return GetKNsystem_phincm_deg(); } ///< GetKNsystem_phincm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetYNsystem_phiycm_kS() const     { return GetYNsystem_phiycm(); } ///< GetYNsystem_phiycm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetYNsystem_phincm_kS() const     { return GetYNsystem_phincm(); } ///< GetYNsystem_phincm with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetYNsystem_phiycm_deg_kS() const { return GetYNsystem_phiycm_deg(); } ///< GetYNsystem_phiycm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)
  double         GetYNsystem_phincm_deg_kS() const { return GetYNsystem_phincm_deg(); } ///< GetYNsystem_phincm_deg with the standard axis convention (argument-free wrapper, used by GetVarByGetter)

  // Data Members
  // ------------
 private:
  /// underflow parameter (1e-7)
  static 
    const double   kUflow;         // underflow parameter (=1e-7)
  /// Format string
  char            *fFormat;        //[100] Format string
  /// How to solve kinematics
  char            *fKinSolution;   //[4] How to solve kinematics
  /// Isospin channel [1-6]
  int              fIsospin;       // isospin channel [1-6]
  /// Standard input/output axis convention
  EAxisConvention  fAC;            // Standard input/output axis convention
  /// Deuteron mass
  double           fMd;            // deuteron mass
  /// Kaon mass
  double           fMk;            // kaon mass
  /// Hyperon mass
  double           fMy;            // hyperon mass
  /// Nucleon mass
  double           fMn;            // nucleon mass
  
  // Kinematics variables (in the deuteron rest frame)
  /// Does kinematic point lie in physical plane
  bool             fIsPhysical;    // Does kinematic point lie in physical plane
  /// Q^2
  double           fQsquared;      // Q^2
  /// Photon energy
  double           fWlab;          // photon energy
  /// Photon momentum
  double           fKlab;          // photon momentum
  /// Kaon Energy
  double           fEk;            // Kaon Energy
  /// Hyperon Energy
  double           fEy;            // Hyperon Energy
  /// Nucleon Energy
  double           fEn;            // Nucleon Energy
  /// Kaon Momentum
  double           fPk;            // Kaon Momentum
  /// Hyperon Momentum
  double           fPy;            // Hyperon Momentum
  /// Nucleon Momentum
  double           fPn;            // Nucleon Momentum
  /// Kaon Polar angle [rad]
  double           fThk;           // Kaon Polar angle [rad]
  /// Hyperon Polar angle [rad]
  double           fThy;           // Hyperon Polar angle [rad]
  /// Nucleon Polar angle [rad]
  double           fThn;           // Nucleon Polar angle [rad]
  /// Cosine of kaon polar angle
  double           fCosthk;        // Cosine of kaon polar angle
  /// Cosine of hyperon polar angle
  double           fCosthy;        // Cosine of hyperon polar angle
  /// Cosine of nucleon polar angle
  double           fCosthn;        // Cosine of nucleon polar angle
  /// Kaon azimuthal angle [rad]
  double           fPhik;          // Kaon azimuthal angle [rad]
  /// Hyperon azimuthal angle [rad]
  double           fPhiy;          // Hyperon azimuthal angle [rad]
  /// Nucleon azimuthal angle [rad]
  double           fPhin;          // Nucleon azimuthal angle [rad]
  /// Total invariant mass
  double           fWtot;          // Total invariant mass
  /// Invariant mass KY-system
  double           fWky;           // Invariant mass KY-system
  /// Invariant mass KN-system
  double           fWkn;           // Invariant mass KN-system
  /// Invariant mass YN-system
  double           fWyn;           // Invariant mass YN-system
  /// Mandelstam s
  double           fStot;          // Mandelstam s
  /// Mandelstam s KY-system
  double           fSky;           // Mandelstam s KY-system
  /// Mandelstam s KN-system
  double           fSkn;           // Mandelstam s KN-system
  /// Mandelstam s YN-system
  double           fSyn;           // Mandelstam s YN-system
  /// Momentum transfer photon to kaon
  double           fTgk;           // Momentum transfer photon to kaon
  /// Momentum transfer photon to hyperon
  double           fTgy;           // Momentum transfer photon to hyperon
  /// Momentum transfer photon to nucleon
  double           fTgn;           // Momentum transfer photon to nucleon
  /// Momentum transfer deuteron to kaon
  double           fTdk;           // Momentum transfer deuteron to kaon
  /// Momentum transfer deuteron to hyperon
  double           fTdy;           // Momentum transfer deuteron to hyperon
  /// Momentum transfer deuteron to nucleon
  double           fTdn;           // Momentum transfer deuteron to nucleon
  /// Kaon scatt. angle in KY CM [rad]
  double           fKYthkcm;       // kaon scatt. angle in KY CM [rad]
  /// Hyperon scatt. angle in KY CM [rad]
  double           fKYthycm;       // hyperon scatt. angle in KY CM [rad]
  /// Kaon scatt. angle in KY CM
  double           fKYcosthkcm;    // kaon scatt. angle in KY CM
  /// Hyperon scatt. angle in KY CM
  double           fKYcosthycm;    // hyperon scatt. angle in KY CM
  /// Kaon azim. angle in KY CM [rad]
  double           fKYphikcm;      // kaon azim. angle in KY CM [rad]
  /// Hyperon azim. angle in KY CM [rad]
  double           fKYphiycm;      // hyperon azim. angle in KY CM [rad]
  /// Kaon scatt. angle in KN CM [rad]
  double           fKNthkcm;       // kaon scatt. angle in KN CM [rad]
  /// Nucleon scatt. angle in KN CM [rad]
  double           fKNthncm;       // nucleon scatt. angle in KN CM [rad]
  /// Kaon scatt. angle in KN CM
  double           fKNcosthkcm;    // kaon scatt. angle in KN CM
  /// Nucleon scatt. angle in KN CM
  double           fKNcosthncm;    // nucleon scatt. angle in KN CM
  /// Kaon azim. angle in KN CM [rad]
  double           fKNphikcm;      // kaon azim. angle in KN CM [rad]
  /// Nucleon azim. angle in KN CM [rad]
  double           fKNphincm;      // nucleon azim. angle in KN CM [rad]
  /// Hyperon scatt. angle in YN CM [rad]
  double           fYNthycm;       // hyperon scatt. angle in YN CM [rad]
  /// Nucleon scatt. angle in YN CM [rad]
  double           fYNthncm;       // nucleon scatt. angle in YN CM [rad]
  /// Hyperon scatt. angle in YN CM
  double           fYNcosthycm;    // hyperon scatt. angle in YN CM
  /// Nucleon scatt. angle in YN CM
  double           fYNcosthncm;    // nucleon scatt. angle in YN CM
  /// Hyperon azim. angle in YN CM [rad]
  double           fYNphiycm;      // hyperon azim. angle in YN CM [rad]
  /// Nucleon azim. angle in YN CM [rad]
  double           fYNphincm;      // nucleon azim. angle in YN CM [rad]
     
  // Subsystems
  /// Kaon-hyperon subkinematics
  TKinematics2to2 *fKYsystem;      //-> kaon-hyperon subkinematics
  /// Kaon-nucleon subkinematics
  TKinematics2to2 *fKNsystem;      //-> kaon-nucleon subkinematics
  /// Hyperon-nucleon subkinematics
  TKinematics2to2 *fYNsystem;      //-> hyperon-nucleon subkinematics

  // 4vectors for all particles (stored in internal reference frame)
  /// Kaon 4vector
  TLorentzVector  *fKvector;       //-> kaon 4vector
  /// Hyperon 4vector
  TLorentzVector  *fYvector;       //-> hyperon 4vector
  /// Nucleon 4vector
  TLorentzVector  *fNvector;       //-> nucleon 4vector
  /// Photon 4vector
  TLorentzVector  *fGvector;       //-> photon 4vector
  /// Deuteron 4vector
  TLorentzVector  *fDvector;       //-> deuteron 4vector
  
  // 6 independent variables
  /// Pointer to independent kinematic variable
  double         **fVar;           //! pointer to independent kinematic variable
  /// Points to independent variable
  double          *fE1Var;         //! points to independent variable
  /// Points to independent variable
  double          *fTh1Var;        //! points to independent variable
  /// Points to independent variable
  double          *fE2Var;         //! points to independent variable
  /// Points to independent variable
  double          *fTh2Var;        //! points to independent variable
  /// Points to independent variable
  double          *fPhi2Var;       //! points to independent variable
  /// Points to independent variable
  double          *fEVar;          //! points to independent variable
  /// Is non-constant kinematic variable?
  bool            *fIsVar;         //[6] is non-constant kinematic variable?
  /// Lower limit non-constant kinematic variable
  double          *fLowerLimit;    //[6] Lower limit non-constant kinematic variable
  /// Upper limit non-constant kinematic variable
  double          *fUpperLimit;    //[6] Upper limit non-constant kinematic variable
  /// Step size non-constant kinematic variable
  double          *fStepSize;      //[6] Step size non-constant kinematic variable
  /// Number of steps non-constant kinematic variable
  int             *fNumberOfSteps; //[6] Number of steps non-constant kinematic variable
  /// step counter (0 or higher)
  int              fStep;          // Step counter (0 or higher)
  /// number of physical steps (-1 if not yet determined)
  mutable int      fNrOfPhysicals; // Number of physical steps (-1 if not determined)
  
  
    
};

// //_____________________________________________________________________
// template<typename T>
// FourVector<T> operator*(const TLorentzRotation& rot, const FourVector<T>& v)
// {
//   // Perform a Lorentz transformation on a strangecalc FourVector object.
//   FourVector<T> rotated;
// 
//   // In principle we need a simple matrix multiplication:
//   // rotated(i) = sum(j) rot(i,j) * v(j)
//   //
//   // Unfortunately, FourVector and ROOT label vectors differently
//   //    | FourVector  |  ROOT  |  (ROOT+1)%4
//   //    | ------------|--------|-------------
//   //  T |     0       |    3   |      0
//   //  X |     1       |    0   |      1
//   //  Y |     2       |    1   |      2
//   //  Z |     3       |    2   |      3
// 
//   for(int rooti=0; rooti<4; ++rooti) {
//     for(int rootj=0; rootj<4; ++rootj) {
//       rotated[(rooti+1)%4] += rot(rooti,rootj) * v[(rootj+1)%4];
//     }
//   }
//   
//   return rotated;
// }

/** @} */
#endif
