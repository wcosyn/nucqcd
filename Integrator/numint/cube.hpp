/**
 * @file cube.hpp
 * @brief numint::cube: one call for numint::cube_romb and numint::cube_adaptive, chosen by IntegrationOptions
 */

#ifndef NUMINT_CUBE_HPP
#define NUMINT_CUBE_HPP 1

#include <stdexcept>
#include <type_traits>

#include <Integrator/numint/integration_method.hpp>
#include <Integrator/numint/numint.hpp>

OPEN_NUMINT_NAMESPACE

/*! Integrates \p f over the box [a,b] with the method of \p options: numint::cube_romb for
 * IntegrationMethod::Romberg (epsabs, epsrel), numint::cube_adaptive for IntegrationMethod::Adaptive (epsabs,
 * epsrel, min_eval, max_eval). The routines are called with exactly these arguments, so the result is the same as
 * from calling them directly. This form is for functions with a vector result (vector_d or vector_z), for which the
 * routines have the arguments \p dirty and \p compare; the form for double and dcomplex results is below.
 * \param options method and precisions
 * \param f the functor to integrate
 * \param a lower limits
 * \param b upper limits
 * \param result the integral value
 * \param neval value incremented with the total number of function evaluations
 * \param dirty check convergence on the first dirty values
 * \param compare use another function to determine convergence
 * \return the return value of the routine that was called
 * \throw std::invalid_argument for a method that is not in IntegrationMethod
 */
template <typename T, unsigned N,
          typename = std::enable_if_t<std::is_same_v<T, vector_d> || std::is_same_v<T, vector_z>>>
int cube(const IntegrationOptions &options, const mdfunction<T,N> &f, const numint::array<double,N> &a,
         const numint::array<double,N> &b, T &result, unsigned &neval, int dirty,
         typename Convergence<T>::type compare = NULL) {
  switch (options.method) {
    case IntegrationMethod::Romberg:
      return cube_romb(f, a, b, options.epsabs, options.epsrel, result, neval, dirty, compare);
    case IntegrationMethod::Adaptive:
      return cube_adaptive(f, a, b, options.epsabs, options.epsrel, options.min_eval, options.max_eval, result, neval,
                           dirty, compare);
  }
  throw std::invalid_argument("numint::cube: unknown integration method");
}

/*! As the form above, for functions with a double or dcomplex result: numint::cube_romb has no arguments \p dirty and
 * \p compare for these, and numint::cube_adaptive is called with dirty = 0 and the default \p compare.
 */
template <typename T, unsigned N,
          typename = std::enable_if_t<std::is_same_v<T, double> || std::is_same_v<T, dcomplex>>>
int cube(const IntegrationOptions &options, const mdfunction<T,N> &f, const numint::array<double,N> &a,
         const numint::array<double,N> &b, T &result, unsigned &neval) {
  switch (options.method) {
    case IntegrationMethod::Romberg:
      return cube_romb(f, a, b, options.epsabs, options.epsrel, result, neval);
    case IntegrationMethod::Adaptive:
      return cube_adaptive(f, a, b, options.epsabs, options.epsrel, options.min_eval, options.max_eval, result, neval, 0);
  }
  throw std::invalid_argument("numint::cube: unknown integration method");
}

CLOSE_NUMINT_NAMESPACE

#endif
