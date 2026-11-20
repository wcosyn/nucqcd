/*! \file test_support.hpp
 * \brief Helpers for the nucqcd test suite: sandbox share directory and golden-value store.
 *
 * Regression tests compare against reference values stored in tests/reference/<name>.ref
 * (plain text, one "key value" pair per line). These values were captured from the code on
 * master when the test suite was introduced. A mismatch means the physics output changed: either
 * a new bug, or a deliberate fix of an old one (the reference values may themselves contain
 * long-standing bugs). It is a lead to investigate, not a verdict.
 *
 * Re-baselining: run with NUCQCD_UPDATE_REFERENCE=1, then review `git diff tests/reference`.
 */
#ifndef NUCQCD_TEST_SUPPORT_HPP
#define NUCQCD_TEST_SUPPORT_HPP

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <map>
#include <string>

namespace nucqcd_test {

/*! Sandbox share directory (copy/symlinks of <repo>/share with an empty grids/), to be passed as
 * the dir/homedir argument of the physics classes. Overridable with NUCQCD_TEST_SHARE. */
std::string share_dir();

/*! true if NUCQCD_UPDATE_REFERENCE=1: reference values are (re)written instead of checked */
bool updating_references();

/*! Result of a comparison with a stored reference value */
struct RefResult {
  bool ok;
  std::string message;
};

/*! Golden-value store for one test executable: tests/reference/<name>.ref */
class RefStore {
public:
  explicit RefStore(const std::string &name);
  /*! Compare value with the stored reference for key: |value-ref| <= max(abs_tol, rel_tol*|ref|).
   * In update mode the value is stored instead (and the check passes). */
  RefResult compare(const std::string &key, double value, double rel_tol = 1e-8, double abs_tol = 0.);

private:
  std::string path;
  std::map<std::string, std::string> values;  // key -> value as written in the file
  void store(const std::string &key, double value);
};

}  // namespace nucqcd_test

/*! CHECK a value against the golden value stored under key (optional args: rel_tol, abs_tol) */
#define CHECK_REF(store, key, ...)                                   \
  do {                                                               \
    const auto nucqcd_ref_result_ = (store).compare(key, __VA_ARGS__); \
    INFO(nucqcd_ref_result_.message);                              \
    CHECK(nucqcd_ref_result_.ok);                                  \
  } while (0)

#endif
