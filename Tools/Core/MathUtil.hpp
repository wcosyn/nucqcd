/*! \file MathUtil.hpp
 * \brief Small numerical helpers that replace ROOT's MathMore and matrix classes (roadmap step 9, phase 2)
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * The functions reproduce ROOT's results bit for bit (checked against ROOT 6.40 on millions of
 * random inputs, see tests/unit/unit_mathutil.cpp), so that replacing ROOT does not change any values:
 * - quadratic_real_roots: ROOT::Math::Polynomial(a,b,c).FindRealRoots() for a quadratic
 * - solve_3x3: TDecompLU (Crout with implicit partial pivoting) followed by Solve (in MathUtil.cpp)
 *
 * \addtogroup MePhys
 * @{
 */
#ifndef MATHUTIL_HPP
#define MATHUTIL_HPP

#include <Tools/Core/ErrorMessage.hpp>

#include <gsl/gsl_poly.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace nucqcd {

/*! Real roots of a x^2 + b x + c = 0, as ROOT::Math::Polynomial(a,b,c).FindRealRoots()
 *
 * For a != 0 the roots come from gsl_poly_complex_solve_quadratic and only the roots with an imaginary
 * part that is exactly zero are kept, in the order GSL returns them (ascending real part; a double root
 * appears twice). For a == 0 the single root -c/b of the linear equation is returned.
 * \param a coefficient of x^2
 * \param b coefficient of x
 * \param c constant term
 * \return the real roots (0, 1 or 2 of them)
 * \throw std::invalid_argument if a == b == 0 (ROOT crashes in that case)
 */
inline std::vector<double> quadratic_real_roots(double a, double b, double c) {
  std::vector<double> roots;
  if (a == 0.) {
    if (b == 0.) throw std::invalid_argument(nucqcd::ErrorMessage() << "quadratic_real_roots: a = b = 0");
    roots.push_back(-c / b);
    return roots;
  }
  gsl_complex z0, z1;
  if (gsl_poly_complex_solve_quadratic(a, b, c, &z0, &z1) != 2) return roots;
  if (GSL_IMAG(z0) == 0.) roots.push_back(GSL_REAL(z0));
  if (GSL_IMAG(z1) == 0.) roots.push_back(GSL_REAL(z1));
  return roots;
}

/*! Solves the 3x3 linear system A x = b, with the same algorithm (and so the same rounding) as ROOT's
 * TDecompLU: Crout LU decomposition with implicit (row-scaled) partial pivoting, then forward and
 * backward substitution.
 * \param a the matrix A, row major: a[3*i+j] = A_ij
 * \param b the right-hand side
 * \return the solution x
 * \throw std::runtime_error if A is singular (a zero pivot)
 *
 * Defined in MathUtil.cpp, which is compiled with -ffp-contract=off: with fused multiply-adds
 * (-march=native) the rounding, and so the Yukawa D-wave coefficients of YukawaWavefunction, would differ from
 * ROOT's (its libMatrix is built without them).
 */
std::array<double, 3> solve_3x3(std::array<double, 9> a, std::array<double, 3> b);

}  // namespace nucqcd
/*! @} */
#endif
