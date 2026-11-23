# Test suite

Catch2 v3 tests registered with CTest. They protect the main physics functionality while code
downstream is changed or refactored. Catch2 is downloaded by CMake (FetchContent, pinned tag) on the
first configure, so that configure needs network access.

## Running

```bash
cd build
cmake ..                          # BUILD_TESTING is ON by default; -DBUILD_TESTING=OFF skips the tests
make -j8                          # builds libraries and test executables (tests/bin/ in the build dir)
ctest -LE slow -j8                # default run: unit + regression tests, ~30 s
ctest -L slow                     # Glauber-grid tests (computes FSI grids from scratch, minutes)
ctest -R reg_poldeut --output-on-failure   # one executable; -V for full output
```

A single executable can also be run directly from `build/tests/work` (its working directory), with
Catch2 options, e.g. `../bin/reg_poldeut "[poldeut]" -s`.

## Labels

| label        | content |
|--------------|---------|
| `unit`       | exact or analytic results: integrals with known values, Dirac algebra, kinematics, wave-function normalizations, D-state probabilities, published constants |
| `regression` | golden values captured from master (`reference/*.ref`) |
| `slow`       | regression tests that compute Glauber grids |

Every `TEST_CASE` runs as its own process, which isolates the static and Fortran COMMON-block state
of the libraries.

## Regression tests and golden values

`reference/<executable>.ref` stores `key value` pairs. These were captured from master when the
suite was introduced (2026-11-20). A regression failure means the output changed by more than the
tolerance (default 1e-8 relative). That can be a new bug, but it can also be a deliberate fix of an
old one, because the reference values may contain long-standing bugs. **Treat a failure as a lead to
investigate.** The failure message prints the computed value, the reference value and the difference.

When a change of a reference value is understood and intended, re-baseline and review the diff:

```bash
NUCQCD_UPDATE_REFERENCE=1 ctest -L "regression|slow"
git diff tests/reference/
```

New keys are added the same way, so a new `CHECK_REF` has no value until it is first run in update mode.

## Isolation from your data

Tests do not use `share/` directly. The fixture `share_sandbox` (re)creates `build/tests/share/`:

- read-only data directories are symlinked
- `input/` is copied, because the code writes into it
- `grids/` is emptied on every ctest run, so Glauber grids are always recomputed by the code under
  test, and your `share/grids/` cache is never read or written (the test `main()` unsets
  `NUCQCD_GRID_DIR`, so a cache directory set in your environment is not used either)

Pass `nucqcd_test::share_dir()` wherever a class takes a `dir`/`homedir`. Classes that use the
compiled-in `nucqcd::share_dir` (GRSV) only read from `share/`.

## Known bugs

Tests tagged `[known-bug]` document bugs found while writing the suite:

- `[!mayfail]` tests run, are reported, and do not fail the suite.
- Hidden `[.]` tests crash the process, so they are not registered with CTest. Run them by hand, e.g.
  `../bin/reg_deuteron_tagged "[known-bug]"`.

When a bug is fixed, remove the tags (and capture reference values where needed).

## Adding tests

- Add `unit/<name>.cpp` or `regression/<name>.cpp`.
- Register it in `CMakeLists.txt`:
  `nucqcd_add_test(<name> DIR unit|regression LABEL unit|regression|slow LINK <highest library> cuba)`.
- Regression values: `nucqcd_test::RefStore refs("<name>");` then
  `CHECK_REF(refs, "key", value[, rel_tol[, abs_tol]]);`, then capture with `NUCQCD_UPDATE_REFERENCE=1`.
- Keep each test case at a few seconds or less, or label the executable `slow`.
