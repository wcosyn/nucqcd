/*! \file IntegratorFlag.hpp
 * \brief Validation of the `integrator` option of the Glauber and Knockout classes
 */
#ifndef INTEGRATORFLAG_HPP
#define INTEGRATORFLAG_HPP

#include <Tools/Core/ErrorMessage.hpp>

#include <stdexcept>

namespace nucqcd {

/*! The classes with an `integrator` option accept 1 (numint::cube_romb) and 2 (numint::cube_adaptive,
 * what the tests and the published results use). 0, the old `rombergerN` (the "Wim's romberg fubini"),
 * was removed in 2026-11: in the Glauber grids its calls were already commented out, so it silently gave
 * grids without FSI.
 * \param integrator the option given to a constructor
 * \return \p integrator
 * \throw std::invalid_argument for any other value
 */
inline int checked_integrator_flag(int integrator) {
  if (integrator != 1 && integrator != 2)
    throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator must be 1 (numint::cube_romb) or 2 (numint::cube_adaptive), got "
                                                       << integrator << " (0, the old rombergerN, was removed)");
  return integrator;
}

}  // namespace nucqcd

#endif
