/*! \file LHAPDFStructure.hpp
 * \brief Leading-order parton-model nucleon structure functions from LHAPDF6 sets (class nucqcd::LHAPDFStructure)
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * \addtogroup MePhys
 * @{
 */
#ifndef LHAPDFSTRUCTURE_HPP
#define LHAPDFSTRUCTURE_HPP

#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace LHAPDF { class PDF; }

namespace nucqcd {

/*! \brief Options of an LHAPDFStructure: the LHAPDF sets and how F1 is built */
struct LHAPDFOptions {
  /*! unpolarized set for F1 and F2, as LHAPDF names it: "CT18LO", or "CT18LO/3" for member 3 (default
   * member 0); empty for none */
  std::string unpolarized;
  std::string polarized; ///< polarized set for g1 and g1+g2 (e.g. "NNPDFpol11_100"), same format; empty for none
  int nf_unpol = 5;      ///< number of quark flavours in the F2 sum
  int nf_pol = 3;        ///< number of quark flavours in the g1 sum
  /*! R = sigma_L/sigma_T used for F1; negative: the SLAC R1998 fit (NucleonStructure::getr1998), as for
   * the MSTW and CTEQ parametrizations of NucleonStructure */
  double R = -1.;
  /*! include the target-mass factor (1 + 4 M^2 x^2/Q^2) in F1 = F2 (1 + 4 M^2 x^2/Q^2) / (2x (1+R)), as
   * the MSTW and CTEQ parametrizations of NucleonStructure do (poldeut's LHAPDFStructure does not) */
  bool target_mass = true;
  int n_ww_x = 300;   ///< number of x nodes of the Wandzura-Wilczek g1+g2 table
  int n_ww_Q2 = 30;   ///< number of Q^2 nodes of the Wandzura-Wilczek g1+g2 table

  /*! the options as a string (also the key of the cache in LHAPDFStructure::get) */
  std::string key() const;
};

/*! \brief Leading-order parton-model structure functions from LHAPDF6 sets
 *
 * \f[ F_2 = \sum_q e_q^2\, x\,(q + \bar q), \qquad g_1 = \tfrac12 \sum_q e_q^2\, (\Delta q + \Delta\bar q), \f]
 * summed over nf_unpol (default 5) and nf_pol (default 3) flavours, with u <-> d interchanged for the
 * neutron (isospin symmetry). F1 = F2 (1 + 4 M^2 x^2/Q^2) / (2x (1+R)) (see LHAPDFOptions for R and
 * the target-mass factor). g1+g2 from the Wandzura-Wilczek relation \f$ g_1+g_2 = \int_x^1 dy/y\, g_1(y) \f$,
 * tabulated on (x, Q^2) nodes the first time it is needed (log-spaced: x from max(xMin, 1e-5) to 0.95
 * plus x = 1, Q^2 in the set's range within [1, 1e5] GeV^2; Gauss-Legendre with 40 points between
 * neighbouring x nodes) and interpolated bilinearly in (ln x, ln Q^2), clamped to the table edges in Q^2.
 * All structure functions are 0 outside 0 < x < 1.
 *
 * Outside a set's (x, Q^2) range LHAPDF's own extrapolation is used (by default "continuation" below
 * Q^2_min and at small x).
 *
 * The C++ counterpart of LHAPDFStructure in taggedpd/nucleon.py of the poldeut repository (there F1
 * uses a constant R = 0.18 and no target-mass factor by default). Use it through NucleonStructure with
 * a name "LHAPDF:\<unpolarized set\>[:\<polarized set\>]", or get one with get(). Loading a set takes a
 * while, so get() keeps every LHAPDFStructure it made (one per set of options) for the rest of the
 * program, and NucleonStructure objects share them.
 *
 * Without LHAPDF (nucqcd built without it, see CLAUDE.md) available() is false and the constructor
 * throws std::runtime_error.
 */
class LHAPDFStructure {
public:
  /*! Loads the sets
   * \param options the sets and the options
   * \throw std::runtime_error if nucqcd was built without LHAPDF, std::invalid_argument if both set
   * names are empty or a set is not installed (see installedSet) */
  explicit LHAPDFStructure(const LHAPDFOptions &options);
  ~LHAPDFStructure();
  LHAPDFStructure(const LHAPDFStructure &) = delete;
  LHAPDFStructure &operator=(const LHAPDFStructure &) = delete;

  /*! The LHAPDFStructure with these options: made at the first call, later calls (from any
   * NucleonStructure) return the same object */
  static std::shared_ptr<const LHAPDFStructure> get(const LHAPDFOptions &options);
  /*! \return whether nucqcd was built with LHAPDF */
  static bool available();
  /*! \param set set name, optionally with "/member" \return whether the set is installed (false without LHAPDF) */
  static bool installedSet(const std::string &set);

  const LHAPDFOptions &options() const { return opts; } ///< the options
  bool hasUnpolarized() const { return unpol != nullptr; } ///< an unpolarized set is loaded (F1, F2)
  bool hasPolarized() const { return pol != nullptr; }     ///< a polarized set is loaded (g1, g1+g2)

  /*! F1 and F2
   * \param[out] F1 F1 structure function \param[out] F2 F2 structure function
   * \param proton proton (true) or neutron (false) \param x Bjorken x \param Q2 [MeV^2] Q^2
   * \throw std::logic_error without an unpolarized set */
  void getF(double &F1, double &F2, bool proton, double x, double Q2) const;
  /*! g1
   * \param proton proton (true) or neutron (false) \param x Bjorken x \param Q2 [MeV^2] Q^2
   * \throw std::logic_error without a polarized set */
  double getG1(bool proton, double x, double Q2) const;
  /*! g1+g2 from the Wandzura-Wilczek relation (interpolated in the table, see the class description)
   * \param proton proton (true) or neutron (false) \param x Bjorken x \param Q2 [MeV^2] Q^2
   * \throw std::logic_error without a polarized set */
  double getG1plusg2(bool proton, double x, double Q2) const;
  /*! \f$ \int_x^1 dy/y\, g_1(y) \f$ by direct integration (no table; slow), to check the table
   * \param proton proton (true) or neutron (false) \param x Bjorken x \param Q2 [MeV^2] Q^2 */
  double getG1plusg2_direct(bool proton, double x, double Q2) const;

  /*! \f$ \sum_q e_q^2 x (q + \bar q) \f$ over nf flavours, straight from LHAPDF
   * \param polarized the polarized set (then \f$ x \sum_q e_q^2 (\Delta q + \Delta\bar q) \f$) or the unpolarized one
   * \param proton proton (true) or neutron (false: u <-> d) \param x Bjorken x \param Q2 [GeV^2] Q^2
   * \param nf number of flavours \return the sum (0 outside 0 < x < 1) */
  double partonSum(bool polarized, bool proton, double x, double Q2, int nf) const;

private:
  void buildWWTable() const; ///< fills the Wandzura-Wilczek table (once, thread safe)
  LHAPDFOptions opts; ///< the options
  std::unique_ptr<LHAPDF::PDF> unpol; ///< unpolarized PDF (may be null)
  std::unique_ptr<LHAPDF::PDF> pol;   ///< polarized PDF (may be null)
  double pol_xmin = 0.;  ///< smallest x of the polarized set
  double pol_q2min = 0.; ///< [GeV^2] smallest Q^2 of the polarized set
  double pol_q2max = 0.; ///< [GeV^2] largest Q^2 of the polarized set
  mutable std::once_flag ww_once;          ///< the WW table is built once
  mutable std::vector<double> ww_lx;       ///< ln x nodes of the WW table
  mutable std::vector<double> ww_lq;       ///< ln Q^2 [GeV^2] nodes of the WW table
  mutable std::vector<double> ww_p;        ///< WW table of the proton g1+g2, [ix*nq + iq]
  mutable std::vector<double> ww_n;        ///< WW table of the neutron g1+g2, [ix*nq + iq]
};

}  // namespace nucqcd
/*! @} */
#endif
