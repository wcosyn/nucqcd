/**
 * @file adaptive_stats.hpp
 * @brief Convergence statistics of numint::cube_adaptive (header-only, standard library only, so
 * that nucqcd::RunLog can report them without linking the integrators).
 */

#ifndef NUMINT_ADAPTIVE_STATS
#define NUMINT_ADAPTIVE_STATS 1

#include <atomic>
#include <cmath>
#include <cstdio>
#include <ostream>

namespace numint {

/**
 * @brief Convergence statistics of all calls to numint::cube_adaptive in this process.
 *
 * The cubature code returns 0 (success) also when it stops because it reached @c maxEval without
 * meeting the requested precision, and it only gives an error *estimate*. cube_adaptive therefore
 * checks the estimate against the request (the same test as the cubature code: for every component
 * either the absolute error <= @c epsabs or the relative error <= @c epsrel) and counts the calls
 * for which that failed ("capped" calls). The return value of cube_adaptive is unchanged.
 * The first capped call prints a warning to stderr; nucqcd::RunLog writes the totals as comment
 * lines at the end of a run. Counts are atomic, so the statistics are safe with threads.
 */
struct AdaptiveStats {
  std::atomic<unsigned long> calls{0};         ///< number of cube_adaptive calls
  std::atomic<unsigned long> capped{0};        ///< calls that stopped without reaching the requested precision
  std::atomic<double> worst_rel_error{0.};     ///< largest estimated relative error (err/|val|, worst component) among capped calls
  std::atomic<unsigned> max_eval_seen{0};      ///< largest maxEval of a capped call
  std::atomic<bool> warned{false};             ///< the warning has been printed
  /// Set to false to silence the warning on stderr (the counts are kept).
  std::atomic<bool> warn{true};
  /// Forget all counts (the warning can be printed again).
  void reset() { calls = 0; capped = 0; worst_rel_error = 0.; max_eval_seen = 0; warned = false; }
};

/// The process-wide statistics of cube_adaptive.
inline AdaptiveStats &adaptive_stats() {
  static AdaptiveStats stats;
  return stats;
}

/**
 * @brief Write a summary of adaptive_stats() to @p out, one line per item, each starting with @p prefix.
 * @return true if something was written (at least one call), false if cube_adaptive was never used.
 */
inline bool adaptive_report(std::ostream &out, const char *prefix = "# ") {
  const AdaptiveStats &s = adaptive_stats();
  const unsigned long n = s.calls, c = s.capped;
  if (n == 0) return false;
  out << prefix << "cube_adaptive: " << n << " integrals, " << c
      << " stopped at maxEval without reaching the requested precision";
  if (c > 0) out << " (worst estimated relative error " << s.worst_rel_error.load()
                 << ", largest maxEval " << s.max_eval_seen.load() << ")";
  out << "\n";
  return true;
}

/// Record the outcome of one adaptive integration (used by cube_adaptive).
inline void adaptive_record(unsigned fdim, const double *val, const double *err, double epsabs, double epsrel,
                            unsigned maxEval) {
  AdaptiveStats &s = adaptive_stats();
  ++s.calls;
  bool converged = true;
  double worst = 0.;
  for (unsigned j = 0; j < fdim; ++j) {
    const double rel = (val[j] == 0. ? HUGE_VAL : std::fabs(err[j] / val[j]));
    if (!(err[j] <= epsabs || rel <= epsrel)) {  // negated so that NaN counts as not converged
      converged = false;
      if (!(rel <= worst)) worst = rel;
    }
  }
  if (converged) return;
  ++s.capped;
  double old = s.worst_rel_error.load();
  while (worst > old && !s.worst_rel_error.compare_exchange_weak(old, worst)) {}
  unsigned oldm = s.max_eval_seen.load();
  while (maxEval > oldm && !s.max_eval_seen.compare_exchange_weak(oldm, maxEval)) {}
  if (s.warn && !s.warned.exchange(true)) {
    std::fprintf(stderr,
                 "warning: numint::cube_adaptive stopped at maxEval=%u without reaching the requested precision "
                 "(epsrel=%g, estimated relative error %g). The result is only as accurate as that estimate. "
                 "This warning is printed once; totals are in numint::adaptive_stats() and in the RunLog footer.\n",
                 maxEval, epsrel, worst);
  }
}

}  // namespace numint

#endif
