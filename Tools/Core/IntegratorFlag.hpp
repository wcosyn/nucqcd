/*! \file IntegratorFlag.hpp
 * \brief Validation of the `integrator` option of the Glauber and Knockout classes
 */
#ifndef INTEGRATORFLAG_HPP
#define INTEGRATORFLAG_HPP

#include <Integrator/numint/integration_method.hpp>
#include <Tools/Core/ErrorMessage.hpp>

#include <cctype>
#include <stdexcept>
#include <string>

namespace nucqcd {

using numint::IntegrationMethod;

/*! The integrator of the Glauber and Knockout classes: 1 (numint::cube_romb) or 2 (numint::cube_adaptive, what the
 * tests and the published results use). 0, the old `rombergerN` (the "Wim's romberg fubini"), was removed in 2026-11:
 * in the Glauber grids its calls were already commented out, so it silently gave grids without FSI.
 * \param flag the integer that the command lines of the programs give
 * \return the method
 * \throw std::invalid_argument for any other value
 */
inline IntegrationMethod integration_method(int flag) {
  if (flag != 1 && flag != 2)
    throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator must be 1 (numint::cube_romb) or 2 (numint::cube_adaptive), got "
                                                       << flag << " (0, the old rombergerN, was removed)");
  return static_cast<IntegrationMethod>(flag);
}

/*! As integration_method(int), for a name or a number given as a string
 * \param name "romberg", "adaptive" (any case), or "1" or "2"
 * \throw std::invalid_argument for any other string
 */
inline IntegrationMethod integration_method(const std::string &name) {
  std::string lower;
  for (char c : name) lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (lower == "romberg" || lower == "1") return IntegrationMethod::Romberg;
  if (lower == "adaptive" || lower == "2") return IntegrationMethod::Adaptive;
  throw std::invalid_argument(nucqcd::ErrorMessage() << "integration method must be \"romberg\" (1) or \"adaptive\" (2), got \"" << name << "\"");
}

/*! Validation of an `integrator` option that is still an integer (the Knockout classes, until they take an
 * IntegrationMethod). \param integrator the option \return \p integrator \throw std::invalid_argument as integration_method(int) */
inline int checked_integrator_flag(int integrator) {
  return numint::method_flag(integration_method(integrator));
}

}  // namespace nucqcd

#endif
