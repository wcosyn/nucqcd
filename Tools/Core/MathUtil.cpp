// Compiled with -ffp-contract=off (see Tools/Core/CMakeLists.txt and MathUtil.hpp), so that the rounding
// is the same as that of ROOT's TDecompLU.
#include "MathUtil.hpp"

#include <cmath>
#include <stdexcept>

namespace nucqcd {

std::array<double, 3> solve_3x3(std::array<double, 9> a, std::array<double, 3> b) {
  constexpr int n = 3;
  double scale[n];
  int index[n];
  for (int i = 0; i < n; i++) {
    double max = 0.;
    for (int j = 0; j < n; j++) {
      const double t = std::abs(a[i * n + j]);
      if (t > max) max = t;
    }
    scale[i] = (max == 0. ? 0. : 1. / max);
  }
  for (int j = 0; j < n; j++) {
    // elements of U in column j above the diagonal
    for (int i = 0; i < j; i++) {
      double r = a[i * n + j];
      for (int k = 0; k < i; k++) r -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = r;
    }
    // residuals on and below the diagonal, pivot = largest scaled one
    double max = 0.;
    int imax = 0;
    for (int i = j; i < n; i++) {
      double r = a[i * n + j];
      for (int k = 0; k < j; k++) r -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = r;
      const double t = scale[i] * std::abs(r);
      if (t >= max) {
        max = t;
        imax = i;
      }
    }
    if (j != imax) {
      for (int k = 0; k < n; k++) {
        const double t = a[imax * n + k];
        a[imax * n + k] = a[j * n + k];
        a[j * n + k] = t;
      }
      scale[imax] = scale[j];
    }
    index[j] = imax;
    if (a[j * n + j] == 0.) throw std::runtime_error(nucqcd::ErrorMessage() << "solve_3x3: singular matrix");
    if (j != n - 1) {
      const double t = 1. / a[j * n + j];
      for (int i = j + 1; i < n; i++) a[i * n + j] *= t;
    }
  }
  // forward substitution (L has a unit diagonal), skipping leading zeros of b as ROOT does
  int nonzero = -1;
  for (int i = 0; i < n; i++) {
    const int lp = index[i];
    double r = b[lp];
    b[lp] = b[i];
    if (nonzero >= 0)
      for (int j = nonzero; j < i; j++) r -= a[i * n + j] * b[j];
    else if (r != 0.)
      nonzero = i;
    b[i] = r;
  }
  // backward substitution
  for (int i = n - 1; i >= 0; i--) {
    double r = b[i];
    for (int j = i + 1; j < n; j++) r -= a[i * n + j] * b[j];
    b[i] = r / a[i * n + i];
  }
  return b;
}

}  // namespace nucqcd
