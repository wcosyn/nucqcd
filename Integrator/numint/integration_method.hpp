/**
 * @file integration_method.hpp
 * @brief The choice between the multidimensional integrators numint::cube_romb and numint::cube_adaptive
 *
 * Header-only and without dependencies, so that classes that only store the choice do not have to include
 * numint.hpp. The dispatch itself is numint::cube in Integrator/numint/cube.hpp.
 */

#ifndef NUMINT_INTEGRATION_METHOD_HPP
#define NUMINT_INTEGRATION_METHOD_HPP 1

#include <Integrator/numint/macros.hpp>

OPEN_NUMINT_NAMESPACE

/*! Which multidimensional integrator to use. The values are the integers that the classes took as an `integrator`
 * option before 2026-11; they are part of the names of the cached Glauber grids (`intgr2`, `integr2`), so they must
 * not change. 0, the old nested Romberg of `rombergerN`, was removed in 2026-11.
 */
enum class IntegrationMethod {
  Romberg = 1,   /*!< numint::cube_romb, nested Romberg */
  Adaptive = 2   /*!< numint::cube_adaptive, adaptive cubature (what the tests and the published results use) */
};

/*! The integer that stands for a method, as used in the cache file names of the grids
 * \param method the method
 * \return 1 for IntegrationMethod::Romberg, 2 for IntegrationMethod::Adaptive
 */
constexpr int method_flag(IntegrationMethod method) { return static_cast<int>(method); }

/*! The arguments of a call of numint::cube: the method and the precisions of the integral */
struct IntegrationOptions {
  IntegrationMethod method;  /*!< which integrator */
  double epsabs;             /*!< absolute precision */
  double epsrel;             /*!< relative precision */
  int min_eval;              /*!< minimum number of function evaluations, only used by IntegrationMethod::Adaptive */
  int max_eval;              /*!< maximum number of function evaluations (0 = no limit), only used by IntegrationMethod::Adaptive */
};

CLOSE_NUMINT_NAMESPACE

#endif
