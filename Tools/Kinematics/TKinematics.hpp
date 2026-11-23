/*! \file TKinematics.hpp
 * \brief Kinematics of photon- or electron-induced 2 -> 2 reactions on a nucleon, for a grid of points
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@UGent.be)
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TKINEMATICS_H
#define TKINEMATICS_H

#include <string>
#include <string>
// #include <Structures.h>
#include <Tools/Core/constants.hpp>
#include <vector>
#include <cmath>


/*! \brief Kinematics of gamma(*) + N -> K + Y (2 -> 2), for a grid of independent variables
 *
 * The names follow meson production on a nucleon: photon (real or virtual, Q^2), nucleon N,
 * meson ("kaon") K and baryon ("hyperon") Y. Isospin channels (SetIsospin): (1) g p -> K+ Lambda,
 * (2) g p -> K+ Sigma0, (3) g p -> K0 Sigma+, (4) g n -> K0 Lambda, (5) g n -> K0 Sigma0,
 * (6) g n -> K+ Sigma-, (7) g p -> pi0 p, (8) g p -> pi+ n, (9) g n -> pi0 n, (10) g n -> pi- p,
 * (11) K- p -> g Lambda, (12) K- p -> g Sigma0 (radiative capture). TKinematics2to2 sets the masses
 * directly; in this code base it is used for example for A(e,e'p) with N the target nucleus, K the
 * residual nucleus and Y the knocked-out nucleon (see Knockout/Model).
 *
 * "lab" is the rest frame of the target N, "cm" the center-of-mass frame of photon and N; the polar
 * angles are with respect to the photon momentum. Energies and momenta in MeV, Q^2 and the Mandelstam
 * variables in MeV^2, angles in rad. The independent variables are given by the format string (see
 * SetFormat); each can be given a range (SetVarRange) and the grid is stepped through with Next() or
 * GoTo().
 */
class TKinematics
{
  std::string fName;  ///< name (formerly from the ROOT base class TNamed)
  std::string fTitle; ///< title

 public:
  const char *GetName() const { return fName.c_str(); }  ///< name
  const char *GetTitle() const { return fTitle.c_str(); } ///< title
  void SetName(const char *name) { fName = name; }        ///< sets the name
  void SetTitle(const char *title) { fTitle = title; }    ///< sets the title

 public:
  /*! constructor with 2 independent variables
   * \param name,title ROOT name and title \param isospin isospin channel, see SetIsospin
   * \param formatString format string with the 2 independent variables, see SetFormat
   * \param var1,var2 initial values of the independent variables, in the order of the format string */
  TKinematics(const char *name,const char *title,
	      int isospin,const char *formatString,
	      double var1,double var2);
  /*! constructor with 3 independent variables
   * \param name,title ROOT name and title \param isospin isospin channel, see SetIsospin
   * \param formatString format string with the 3 independent variables, see SetFormat
   * \param var1,var2,var3 initial values of the independent variables, in the order of the format string */
  TKinematics(const char *name,const char *title,
	      int isospin,const char *formatString,
	      double var1,double var2,double var3);
  /*! constructor with 4 independent variables
   * \param name,title ROOT name and title \param isospin isospin channel, see SetIsospin
   * \param formatString format string with the 4 independent variables, see SetFormat
   * \param var1,var2,var3,var4 initial values of the independent variables, in the order of the format string */
  TKinematics(const char *name,const char *title,
	      int isospin,const char *formatString,
	      double var1,double var2,double var3,double var4);
  TKinematics(const TKinematics&); ///< copy constructor
  virtual TKinematics* Clone(const char *newname ="") const; ///< virtual copy constructor \param newname name of the copy (optional)
  virtual ~TKinematics(); ///< destructor
  TKinematics& operator=(const TKinematics&); ///< assignment
  bool          operator==(const TKinematics&) const; ///< same isospin channel, format and independent variables

  void          SetVarRange(int varNum,double lowerLimit,double upperLimit,int steps); ///< gives independent variable varNum a range [lowerLimit,upperLimit] with steps points
  void          FixVariable(int varNumber); ///< removes the range of variable varNumber and fixes it at its current value
  void          FixVariables(); ///< calls FixVariable on all variables
  int           Next(); ///< proceeds to the next point of the grid; returns zero at the end of the grid
  void          GoTo(int totalStep); ///< proceeds to point totalStep of the grid
  void          GoTo(int* steps); ///< proceeds to the grid point with the given step for each variable (array of size GetNumberOfIndependentVariables())
  void          GoTo(int stepVar1, int stepVar2); ///< proceeds to the grid point with the given step for each variable
  void          GoTo(int stepVar1, int stepVar2, int stepVar3); ///< proceeds to the grid point with the given step for each variable
  void          GoTo(int stepVar1, int stepVar2, int stepVar3, int stepVar4); ///< proceeds to the grid point with the given step for each variable
  virtual void  Print(const char* option="") const; ///< prints the object and the current kinematics
  int           VariableInfo() const; ///< prints the ranges given to the variables
  static void   Help(); ///< prints help on the format strings and variables

  void          SetIsospin(int isospin); ///< changes the isospin channel (see the class description) and updates the kinematics
  void          SetFormat(const char *formatString); ///< sets the format string that lists the independent variables (Q^2 and phi are zero when not listed)
  void          SetVar(int varNum,double value); ///< sets independent variable varNum (not allowed when it has a range)
  void          SetPhiLimits(double phiMin, double phiMax); ///< sets phi = phiMax and phiMin, and fixes phi (for observables integrated over phi)

  int           GetIsospin() const     { return fIsospin; } ///< isospin channel, see the class description
  const char   *GetFormat() const      { return fFormat; } ///< format string, see SetFormat
  double        GetNucleonMass() const { return fMp; } ///< [MeV] mass of the target N
  double        GetMesonMass() const   { return fMk; } ///< [MeV] mass of the meson K
  double        GetHyperonMass() const { return fMy; } ///< [MeV] mass of the baryon Y
  double        GetWlab() const        { return fWlab; } ///< [MeV] photon energy in the lab
  double        GetWcm() const         { return fWcm; } ///< [MeV] photon energy in the cm
  double        GetKlab() const        { return fKlab; } ///< [MeV] photon momentum in the lab
  double        GetKcm() const         { return fKcm; } ///< [MeV] photon momentum in the cm
  double        GetPk() const          { return fPk; } ///< [MeV] meson momentum in the cm
  double        GetPkcm() const        { return fPk; } ///< [MeV] meson momentum in the cm (same as GetPk)
  double        GetPklab() const       { return fPklab; } ///< [MeV] meson momentum in the lab
  double        GetCosthkcm() const    { return fCosthkcm; } ///< cosine of the meson angle with the photon in the cm
  double        GetCosthklab() const    { return fCosthklab; } ///< cosine of the meson angle with the photon in the lab
  double        GetPYlab() const       { return fPYlab; } ///< [MeV] baryon momentum in the lab
  double        GetCosthYlab() const    { return fCosthYlab; } ///< cosine of the baryon angle with the photon in the lab
  double        GetEYlab() const {return sqrt(fMy*fMy+fPYlab*fPYlab);} ///< [MeV] baryon energy in the lab
  double        GetEklab() const {return sqrt(fMk*fMk+fPklab*fPklab);} ///< [MeV] meson energy in the lab
  double        GetQsquared() const    { return fQsquared; } ///< [MeV^2] Q^2 of the photon
  double        GetW() const           { return fW; } ///< [MeV] invariant mass
  double        GetS() const           { return fS; } ///< [MeV^2] Mandelstam s
  double        GetT() const           { return fT; } ///< [MeV^2] Mandelstam t = (q - p_K)^2
  double        GetTMin() const; ///< [MeV^2] minimum value of -t
  double        GetU() const           { return fU; } ///< [MeV^2] Mandelstam u = (q - p_Y)^2
  double        GetXb() const          { return fXb; } ///< Bjorken x = Q^2/(2 nucqcd::mass_n omega_lab)
  double        GetPhi() const         { return fPhi; } ///< [rad] angle between the lepton and the hadron plane
  double        GetPhiMin() const      { return fPhiMin; } ///< [rad] lower limit of phi for observables integrated over phi, see SetPhiLimits
  bool          IsPhysical() const     { return fIsPhysical; } ///< does the current point lie in the physical region?
  int           GetNumberOfIndependentVariables() const { return fNVars; } ///< number of independent variables in the format string (2 to 4)
  double        GetVar(int) const; ///< current value of independent variable varNum (1 to GetNumberOfIndependentVariables())
  double       *GetVarArray(const std::string& varName) const; ///< values of the named variable (any variable with a getter) at all grid points; the caller owns the array (size GetNumberOfSteps())
  double       *GetPhysicalVarArray(const std::string& varName) const; ///< values of the named variable at the physical grid points; the caller owns the array (size GetNumberOfPhysicalSteps())
  std::string   GetVarName(int varNumber) const; ///< name of independent variable varNumber in the format string
  int           GetNumberOfSteps(int varNumber) const; ///< number of grid points of variable varNumber
  int           GetNumberOfSteps() const; ///< number of points in the grid
  int           GetNumberOfPhysicalSteps() const; ///< number of physical points in the grid
  int           GetNumberOfVariables() const; ///< number of independent variables that have a range
  double        GetStepSize(int varNumber) const; ///< step size of independent variable varNumber
  int           GetStep(int varNumber) const; ///< current grid position of variable varNumber
  int           GetStep() const; ///< current position in the grid
  bool          IsFixed() const; ///< true when none of the variables has a range
  bool          IsFixed(int varNumber) const; ///< true when variable varNumber has no range

 protected:
  void          SetMn(double mn) { fMp = mn; } ///< sets the nucleon (target) mass \param mn [MeV]
  void          SetMk(double mk) { fMk = mk; } ///< sets the meson mass \param mk [MeV]
  void          SetMy(double my) { fMy = my; } ///< sets the baryon mass \param my [MeV]
  void          UpdateKinematics(); ///< solves the kinematics for the current values of the independent variables
  // Underflow parameter for TKinematics
  virtual 
    double      Underflow() const { return nucqcd::strange_underflow; } ///< underflow parameter (nucqcd::strange_underflow)

 private:
  void          InitializeArrays(); ///< resizes the member vectors to the number of independent variables, all constant
  void          ReadFormatString(const std::string&); ///< processes the format string and points the independent variables at the right members
  void          ChangeIsospinChannel(int); ///< sets the masses for the isospin channel
  double        EnergyConservation(double); ///< deviation from energy conservation in the cm for meson momentum pk [MeV]
  
  // Data Members
  // ------------
  /// format string
  char    *fFormat;        //[100] Format string
  /// isospin channel
  int      fIsospin;       // isospin channel [1-6]
  /// number of independent variables
  int      fNVars;         // number of independent variables
  /// [MeV] nucleon (target) mass
  double   fMp;            // nucleon mass
  /// [MeV] meson mass
  double   fMk;            // meson mass
  /// [MeV] baryon mass
  double   fMy;            // hyperon mass

  // Kinematic variables
  /// [MeV] photon energy (lab)
  double   fWlab;          // photon energy (LAB)
  /// [MeV] photon energy (cm)
  double   fWcm;           // photon energy (CM)
  /// [MeV] photon momentum (lab)
  double   fKlab;          // photon momentum (LAB)
  /// [MeV] photon momentum (cm)
  double   fKcm;           // photon momentum (CM)
  /// [MeV] meson momentum (cm)
  double   fPk;            // kaon momentum (CM)
  /// [MeV] meson momentum (lab)
  double   fPklab;         // kaon momentum (LAB)
  /// cosine of the meson angle (cm)
  double   fCosthkcm;      // cos of kaon scattering angle in CM
  /// cosine of the meson angle (lab)
  double fCosthklab;
  /// [MeV] baryon momentum (lab)
  double   fPYlab;	//hyperon mom in lab
  /// cosine of the baryon angle (lab)
  double fCosthYlab;
  /// [MeV^2] Q^2
  double   fQsquared;      // Q^2
  /// Bjorken x
  double   fXb;            // Bjorken x
  /// [MeV] invariant mass
  double   fW;             // invariant mass
  /// [MeV^2] Mandelstam s
  double   fS;             // s (Mandelstam variable)
  /// [MeV^2] Mandelstam t
  double   fT;             // t (Mandelstam variable)
  /// [MeV^2] Mandelstam u
  double   fU;             // u (Mandelstam variable)
  /// [rad] angle between the lepton plane and the hadron plane
  double   fPhi;
  /// [rad] lower limit of phi, used only when the virtual-photon cross section is integrated over phi ("diff_int" observable)
  double   fPhiMin;
  /// is the current point physical?
  bool     fIsPhysical;    // Does kinematic point lie in physical plane
  /// number of physical points (-1 when not yet determined)
  mutable 
    int    fNrOfPhysicals; // Number of physical points (-1 when not determined)

  // 3 independent variables
  /// pointers to the independent variables
  std::vector<double*> fVar;           //! pointer to independent kinematic variables
  /// the independent energy variable
  double              *fEVar;          //! points to independent energy variable
  /// the independent angular variable
  double              *fAVar;          //! points to independent angular variable
  /// does the variable have a range?
  std::vector<bool>    fIsVar;         // is non-constant kinematic variable?
  /// lower limits of the ranges
  std::vector<double>  fLowerLimit;    // Lower limit non-constant kinematic variable
  /// upper limits of the ranges
  std::vector<double>  fUpperLimit;    // Upper limit non-constant kinematic variable
  /// step sizes of the ranges
  std::vector<double>  fStepSize;      // Step size non-constant kinematic variable
  /// numbers of steps of the ranges
  std::vector<int>     fNumberOfSteps; // Number of steps non-constant kinematic variable
  /// current step of each variable
  std::vector<int>     fStepVar;       // Step counter for each variable

};


/** @} */
#endif
