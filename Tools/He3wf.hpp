/*! \file He3wf.hpp 
 * \brief Contains declaration of class He3wf, has a variety of CD-Bonn and AV18 he3 wave functions
 * uses fortran code obtained from Misak Sargsian [getwavemod_tm.f,clebsch.f,spline.f].  This class is basically a wrapper for that code.
 * \author Wim Cosyn
 * \date 19/9/2014
 * 
 * \addtogroup MePhys
 * @{
 */
#ifndef HE3WF_HPP
#define HE3WF_HPP

#include <string>
#include <TVector3.h>
#include <complex>

/*! \brief A class that returns the momentum space he3 wave function for a variety of parametrizations
 *
 * The wave function is the sum of the three Faddeev components (pair 12 with spectator 3 and its two
 * cyclic permutations) computed by the Fortran code. It is normalized to 1:
 * \f$ \sum_{m_{s1} m_{s2} m_{s3}} \sum_{m_{t1} m_{t2}} \int d^3p_1\, d^3p_2\, |\Psi_{m_A}|^2 = 1 \f$
 * (0.9998-0.9999 for the five parametrizations, checked in Jacobi coordinates, unit_he3wf).
 *
 * \note The spins of the three nucleons do not have to add up to m_A: the D and P waves (L >= 1) also
 * have configurations with ms1+ms2+ms3 = mA -+ 2 (and +-4). Until 2026-11 getWF took no ms3 and set
 * ms3 = mA - ms1 - ms2, returning 0 when that was not +-1, so those components were missing: the norm
 * was 0.967-0.975, and the spin-summed density was not rotation invariant (by up to 30% for some
 * momenta). He3Cross summed over ms1 and ms2 only. The Fortran code was right.
 */
class He3wf{
public:
  /*!Constructor
   * \param wf_name parametrization of the wave function.  Possibilities (Case sensitive!!):
   * "AV18", "AV18-UIX", "AV18-TML", "CDBonn", "CDBonn-TML"
   * \param dir dir where the data files are located
   * \throw std::invalid_argument for an unknown name
   */
  He3wf(std::string wf_name, std::string dir); 
  He3wf(); /*!< Default Constructor*/
  He3wf(const He3wf&); /*!< Copy Constructor */
  He3wf& operator=(const He3wf&); /*!< assignment operator */
  /*! Calculate the he3 wf for a certain momentum, spin and isospin configuration (in the He3 rest
   * frame: p3 = -p1-p2)
   * \param p1 [MeV] momentum of particle 1
   * \param p2 [MeV] momentum of particle 2
   * \param ms1 [-1,1] twice the spin projection of particle 1
   * \param ms2 [-1,1] twice the spin projection of particle 2
   * \param ms3 [-1,1] twice the spin projection of particle 3 (independent of the others and of mA)
   * \param mt1 [-1,1] twice the isospin projection of particle 1 (1 proton, -1 neutron)
   * \param mt2 [-1,1] twice the isospin projection of particle 2 (1 proton, -1 neutron)
   * \param mA [-1,1] twice the spin projection of the nucleus
   * \return [MeV^-3] returns the complex valued He3 wave function (0 if mt3 = 1-mt1-mt2 is not +-1)
   * \throw std::invalid_argument for spin or isospin projections other than +-1
   */
  std::complex<double> getWF(const TVector3 &p1, const TVector3 &p2, int ms1, int ms2, int ms3, int mt1, int mt2, int mA);
  /*! Calculate the he3 wf for a certain momentum, spin and isospin configuration, see the other getWF
   * \param p1 [MeV] momentum of particle 1
   * \param p2 [MeV] momentum of particle 2
   * \param ms1 [-1,1] twice the spin projection of particle 1
   * \param ms2 [-1,1] twice the spin projection of particle 2
   * \param ms3 [-1,1] twice the spin projection of particle 3
   * \param mt1 [-1,1] twice the isospin projection of particle 1 (1 proton, -1 neutron)
   * \param mt2 [-1,1] twice the isospin projection of particle 2 (1 proton, -1 neutron)
   * \param mA [-1,1] twice the spin projection of the nucleus
   * \return [MeV^-3] returns the complex valued He3 wave function
   */
  std::complex<double> getWF(const double p1[3], const double p2[3], int ms1, int ms2, int ms3, int mt1, int mt2, int mA);
};
/** @} */
#endif