#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <complex>
#include <vector>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {


/*! \brief Older implementation of the non-relativistic deuteron wave function (Paris, AV18, CD-Bonn, AV18b)
 *
 * deuteronwf gives the same amplitudes as TDeuteron::Wavefunction::DeuteronPState (with nucleon 1 the
 * proton; checked numerically for AV18). The radial functions absorb factors: U = U_TDeuteron/sqrt(4 pi)
 * and W = -W_TDeuteron/sqrt(4 pi), so the wave function is written as U + W S12/sqrt(8), and
 * \f$ 4\pi\int p^2 (U^2+W^2)\,dp = 1 \f$.
 */
class OldDeuteron{

public:

/**
 * @brief Construct a new Old Deuteron object
 * 
 * @param wf_name parmatrization.  Possibilities "Paris", "AV18", "CDBonn", "AV18b"
 */
OldDeuteron(const std::string &wf_name);

/**
 * @brief computes non-relativistic deuteron wf at certain kinematics and spin states
 * 
 * @param dspin  two times deuteron mj
 * @param proton active nucleon proton or not?
 * @param spinp  two times spin ms of proton
 * @param spinn two times spin ms of neutron
 * @param p  [MeV] norm of momentum 
 * @param theta theta angle of momentum
 * @param phi  azimuthal angle of momentum
 * @return std::complex< double>  [MeV^{-3/2}] wave function value
 */
std::complex< double> deuteronwf(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi);
/**
 * @brief computes radial S-wave
 * 
 * @param p [MeV]
 * @return double [MeV^-3/2] S-wave, U_TDeuteron/sqrt(4 pi)
 */
double U(const double p);

/**
 * @brief computes radial D-wave
 * 
 * @param p [MeV]
 * @return double [MeV^-3/2] D-wave, -W_TDeuteron/sqrt(4 pi)
 */
double W(const double p);


private:

std::vector<double> mass; /*!< mass parameters of the Yukawa terms */
std::vector<double> S_pre; /*!< [fm^-1/2] S-wave coefficients */
std::vector<double> D_pre; /*!< [fm^-1/2] D-wave coefficients */
std::string wf_name; ///< name of the parametrization
double norm; ///< normalization factor

/**
 * @brief CG factor with S-wave
 * 
 * @param M  two times deuteron Mj
 * @param spinp two times active nucleon ms
 * @param spinn two times spectator nucleon ms
 * @return double CG factor
 */
double Ufront(int M, int spinp, int spinn);

/**
 * @brief computes D-wave spherical harmonic + CG coefficient
 * 
 * @param dspin two times deuteron mj
 * @param spinp two times proton ms
 * @param spinn two times neuttron ms
 * @param theta angle theta
 * @param phi azimuthal angle
 * @return std::complex<double> D-wave coefficient
 */
std::complex<double> get_stensor(int dspin, int spinp, int spinn, double theta, double phi);


};

}  // namespace nucqcd

