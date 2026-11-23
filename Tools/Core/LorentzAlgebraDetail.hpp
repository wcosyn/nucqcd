/*! \file LorentzAlgebraDetail.hpp
 * \brief Internal helpers of nucqcd::Vector2/Vector3/LorentzVector/Rotation/LorentzRotation: the
 * constants and the atan2 of ROOT's TMath, so that the results are the same as ROOT's.
 */
#ifndef NUCQCD_LORENTZALGEBRADETAIL_HPP
#define NUCQCD_LORENTZALGEBRADETAIL_HPP

#include <cmath>

namespace nucqcd {
namespace detail {

constexpr double pi = 3.14159265358979323846;  ///< TMath::Pi()
constexpr double two_pi = 2.0 * pi;            ///< ROOT's kTWOPI (TVector2)

/*! TMath::ATan2: atan2(y, x), except for x = 0 (also -0), where it returns 0 for y = 0 and +-pi/2
 * otherwise (std::atan2 gives +-0 or +-pi depending on the signs of the zeros). */
inline double root_atan2(double y, double x) {
  if (x != 0) return std::atan2(y, x);
  if (y == 0) return 0;
  if (y > 0) return pi / 2;
  return -pi / 2;
}

}  // namespace detail
}  // namespace nucqcd

#endif
