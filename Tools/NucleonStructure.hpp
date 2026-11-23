/*! \file NucleonStructure.hpp 
 * \brief Contains declaration of class NucleonStructure, contains a variety of nucleon structure function parametrizations
 * \author Wim Cosyn
 * \date 12/09/2016
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef NUCLEONSTRUCTURE_HPP
#define NUCLEONSTRUCTURE_HPP

#include <string>
#include <cmath>
#include <memory>
namespace nucqcd { class LHAPDFStructure; struct LHAPDFOptions; }

/*! \brief A class for a variety of nucleon structure parametrizations: Christy & Bosted, SLAC (Bodek),
 * Alekhin (leading twist), parton-model F2 from CTEQ, HMRS, MSTW or any LHAPDF6 set, and combinations
 * of two of them (one below and one above a W^2 boundary).
 */
class NucleonStructure{
  
public:
  NucleonStructure(); /*!< default constructor, SLAC parametrization taken */
  /*! Constructor
   * \param name Name of the parametriation. Possibilities: <BR>
   * "CB": Christy & Bosted parametrization (see F1F209.f file); valid up to W^2 = 25 GeV^2, throws
   * std::out_of_range above (no switch to SLAC here, unlike NuclStructure, which uses SLAC above
   * W^2 = 9 GeV^2)<BR>
   * "SLAC": SLAC paramtetrization from Bodek <BR>
   * "Alekhin": leading twist parametrization by Alekhin [see PRD 68,014002], also see alekhin.f file <BR>
   * "CTEQ": F2 based on the pdf's from CTEQ (code from Misak, see cteq.f file) <BR>
   * "HMRS": F2 from pdfs from unknown source (Shunzo Kumano is the source, not enough info in his code file, see hmrs-b.f) <BR>
   * "MSTW": F2 from the MSTW 2008 LO PDFs, an alias of "LHAPDF:MSTW2008lo68cl" (central member; needs
   * nucqcd built with LHAPDF and the set installed). Until 2026-11 nucqcd had its own copy of the MSTW
   * grid code and grid (share/mstw_grids); LHAPDF interpolates the same grid differently: the same at
   * the grid nodes, up to ~1e-3 relative between them, ~0.1% below Q^2 = 1 GeV^2 (extrapolation).<BR>
   * "LHAPDF:\<set\>" or "LHAPDF:\<set\>:\<polarized set\>": leading-order parton model from LHAPDF6 sets,
   * see nucqcd::LHAPDFStructure (default options: F1 with R1998 and the target-mass factor, as for MSTW);
   * a set name may end in "/\<member\>" (default member 0), the unpolarized set may be empty
   * ("LHAPDF::NNPDFpol11_100" for g1 only). Needs nucqcd built with LHAPDF. <BR>
   * "\<low\>+\<high\>" or "\<low\>+\<high\>@\<W2\>": the parametrization \<low\> up to W^2 = \<W2\> GeV^2 (default 9)
   * and \<high\> above it, e.g. "CB+SLAC" (what NuclStructure uses for "CB") or "CB+LHAPDF:CT18LO@4";
   * see the constructor with three arguments. <BR>
   * \throw std::invalid_argument for a malformed "LHAPDF:" or combination name (other unknown names
   * only throw when a structure function is asked for)
   */
  NucleonStructure(const std::string name);
  /*! Constructor for LHAPDF sets with options (number of flavours, R, ...); getName() is
   * "LHAPDF:\<unpolarized set\>:\<polarized set\>"
   * \param options the sets and options, see nucqcd::LHAPDFOptions */
  explicit NucleonStructure(const nucqcd::LHAPDFOptions &options);
  /*! Constructor for a combination of two parametrizations: low for W^2 <= Wsq_switch, high above.
   * getName() is "\<low\>+\<high\>@\<Wsq_switch in GeV^2\>".
   * \param low name of the parametrization at low W^2 (e.g. "CB")
   * \param high name of the parametrization at high W^2 (e.g. "SLAC", "LHAPDF:CT18LO")
   * \param Wsq_switch [MeV^2] the boundary in W^2 */
  NucleonStructure(const std::string &low, const std::string &high, double Wsq_switch);
  NucleonStructure(const NucleonStructure&); /*!< Copy Constructor */
  NucleonStructure& operator=(const NucleonStructure&); /*!< assignment operator */
  ~NucleonStructure(); /*!< Destructor */
  
  
  /*! Returns the F1 and F2 structure functions 
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_xQ(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
  /*! Returns the F1 and F2 structure functions for general (off-shell) kinematics, with x and W^2
   * independent: Christy-Bosted and SLAC (Bodek) use W^2, the other parametrizations x. With the
   * on-shell W^2 = Q^2(1/x-1) + M^2 this is getF_xQ. Used by NuclStructure.
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \param Wsq [MeV^2] invariant mass squared of the produced X
   */
  void getF_xQW(double &F1, double &F2, const bool proton, const double x, const double Q2, const double Wsq) const;
  /*! Returns the F1 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF1_xQ(const bool proton, const double x, const double Q2) const; 
  /*! Returns the F2 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF2_xQ(const bool proton, const double x, const double Q2) const; 
  /*! Returns the F1 and F2 structure functions 
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_WQ(double &F1, double &F2, const bool proton, const double Wsq, const double Q2) const;
  /*! Returns the F1 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF1_WQ(const bool proton, const double Wsq, const double Q2) const; 
  /*! Returns the F2 structure functions starting from bjorken x and Q^2
   * \param proton 0 neutron, 1 proton
   * \param Wsq [MeV^2] invariant mass produced X
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 structure function
   */
  double getF2_WQ(const bool proton, const double Wsq, const double Q2) const; 
  const std::string getName() const{return name;} /*!< returns the name of the chosen parametrization */
  /*! \param Wsq [MeV^2] invariant mass squared W^2
   * \return the name of the parametrization used at W^2 (for a combination its low or high part, else getName()) */
  std::string getName(double Wsq) const;
  /*! Whether the parametrization used at W^2 is a fit of F2 only, so that F1 is not fitted: SLAC (Bodek),
   * where F1 follows from F2 and R1998. DeuteronStructure then builds F1 from F2 itself (with R = 0.18
   * and the off-shell kinematics). False for the others, also for the parton-model ones (CTEQ, MSTW,
   * LHAPDF), whose F1 comes from F2 and R1998 too.
   * \param Wsq [MeV^2] invariant mass squared W^2 */
  bool isF2Fit(double Wsq) const;

  /*! g1 for an LHAPDF parametrization with a polarized set (see nucqcd::LHAPDFStructure)
   * \param proton 1=proton 0=neutron \param x Bjorken x \param Q2 [MeV^2] Q^2
   * \throw std::logic_error for other parametrizations (GRSV: the static getG1_grsv2000) */
  double getG1(const bool proton, const double x, const double Q2) const;
  /*! g1+g2 (Wandzura-Wilczek) for an LHAPDF parametrization with a polarized set
   * \param proton 1=proton 0=neutron \param x Bjorken x \param Q2 [MeV^2] Q^2
   * \throw std::logic_error for other parametrizations (GRSV: the static getG1plusg2_grsv2000) */
  double getG1plusg2(const bool proton, const double x, const double Q2) const;
  
  /*! returns the g1 structure function from the grsv2000 parametrization (leading order,  hep-ph/0011215, Tools/grsv2000pdf_g1.f )
   * \param proton 1=proton 0=neutron
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return g1n structure functions
   */
  static double getG1_grsv2000(const bool proton, const double x, const double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1990 parametrization 0.1 < x < 0.9, 0.6 < Q^2 < 20 GeV^2
   * L. W. Withlow et al., Phys. Lett. B250 194(1990)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1990(const double x, const double Q2);
  /*! returns the sigma_L/sigma_T ratio, from SLAC fits (code from Shunzo Kumano, see R199x.f)
   * This is the 1998 parametrization 0.005 < x < 0.86, 0.5 < Q^2 < 130 GeV^2
   * (E143) K. Abe et al., Phys. Lett B452 194 (1999)
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return R=sigma_L/sigma_T (indicates scaling violations from Callan-Gross LT relation)
   */
  static double getr1998(const double x, const double Q2);


  /*! returns the g1+g2 structure function using the Wandura-Wilczek relation [g1+g2=int dy/y g1(y)] from the grsv2000 parametrization (leading order,  hep-ph/0011215, Tools/grsv2000pdf_g1.f )
   * This sum of structure functions is related to the transversity pdfs (and higher twist).
   * \param proton 1=proton 0=neutron
   * \param x Bjorken x nucleon
   * \param Q2 [MeV^2] virtual photon four-momentum sq
   * \return g1+g2 structure functions
   */
  static double getG1plusg2_grsv2000(const bool proton, double x, double Q2);
  
private:
  /*! F1 and F2 in the SLAC parametrization with an independent W^2 [MeV^2], see getF_xQW */
  void getF_SLAC_W(double &F1, double &F2, const bool proton, const double x, const double Q2, const double Wsq) const;
  /*! F2 in the SLAC parametrization (Bodek) with an independent W^2 [MeV^2] */
  double getF2_SLAC_W(const bool proton, const double x, const double Q2, const double Wsq) const;
  /*! F1 and F2 in the Christy & Bosted parametrization at Q^2 [MeV^2] and W^2 [MeV^2] (throws above W^2 = 25 GeV^2) */
  void getF_CB_W(double &F1, double &F2, const bool proton, const double Q2, const double Wsq) const;
  std::string name; /*!< name of the paramtetrization */
  std::string dir; /*!< share dir to read input files */
  std::shared_ptr<const nucqcd::LHAPDFStructure> lhapdf; ///< LHAPDF sets (shared, see LHAPDFStructure::get), or null
  std::unique_ptr<NucleonStructure> low;  ///< combination: parametrization for W^2 <= Wsq_switch, or null
  std::unique_ptr<NucleonStructure> high; ///< combination: parametrization for W^2 > Wsq_switch, or null
  double Wsq_switch = 0.; ///< [MeV^2] combination: W^2 boundary between low and high
  
  /*! proton F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of proton
   */
  double f2p_b(const double massi, const double xp, const double q2, const double fm) const;
  /*! neutron F2 SLAC parametrization
   * \param massi [GeV] initial mass, off-shell
   * \param xp bjorken x
   * \param q2 [GeV^2] Q^2 of virtual photon
   * \param fm [GeV] invariant mass W
   * \return F2 of neutron
   */
  double f2n_b(const double massi, const double xp, const double q2, const double fm) const;
  
  /*! SLAC parametrization helper function
   * \param wm [GeV] invariant mass
   * \param qsq [GeV^2] Q^2 of virtual photon
   * \return something :p
   */
  double bodek(const double wm, const double qsq) const;
  
  
  /*! returns F1 and F2 in the alekhin parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_Alekhin(double &F1, double &F2, const bool proton, const double x, const double Q2)const;
  /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 in Alekhin parametriation
   */
  double getF1_Alekhin(const bool proton, const double x, const double Q2) const; 
  /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F2 in Alekhin parametriation
   */
  double getF2_Alekhin(const bool proton, const double x, const double Q2) const; 
  /*! returns F1 and F2 in the SLAC parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_SLAC(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F2 in SLAC parametriation
   */
  double getF2_SLAC(const bool proton, const double x, const double Q2) const; 
  /*! returns F1 and F2 in the Christy & Bosted parametrization
   * \param[out] F1 F1 structure function, determined through Callan-Gross relation incl scaling violation
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_CB(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F1 in Christy & Bosted parametriation
   */
  double getF1_CB(const bool proton, const double x, const double Q2) const;
    /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F2 in Christy & Bosted parametriation
   */
double getF2_CB(const bool proton, const double x, const double Q2) const;
  /*! returns F1 and F2 in the CTEQ parametrization
   * \param[out] F1 F1 structure function
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_CTEQ(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
    /*! \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   * \return F2 in CTEQ parametriation
   */
double getF2_CTEQ(const bool proton, const double x, const double Q2) const;/*!< return F2 structure function in the CTEQ parametrization */
  
  /*! returns F1 and F2 in the HRMS-B parametrization
   * \param[out] F1 F1 structure function, determined through Callan-Gross relation incl scaling violation
   * \param[out] F2 F2 structure function
   * \param proton 0 neutron, 1 proton
   * \param x Bjorken x
   * \param Q2 [MeV^2] fourmomentum sq
   */
  void getF_HMRS(double &F1, double &F2, const bool proton, const double x, const double Q2) const;
  

  
};
/*! @} */
#endif