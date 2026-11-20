// Common main() of all test executables.
#include <catch2/catch_session.hpp>

#include <cstdlib>

int main(int argc, char *argv[]) {
  // Cuba forks worker processes by default; run it serially so results are reproducible
  // (an explicit CUBACORES in the environment wins).
  setenv("CUBACORES", "0", 0);
  return Catch::Session().run(argc, argv);
}
