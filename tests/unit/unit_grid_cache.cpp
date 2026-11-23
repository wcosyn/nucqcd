// Tools/Core/GridCache.hpp: location of the grid cache (NUCQCD_GRID_DIR or <share>/grids/) and atomic
// writes of cache files (temporary file + rename, so no truncated cache file is ever left behind).
#include <test_support.hpp>

#include <GridCache.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace fs = std::filesystem;

namespace {
// fresh scratch directory under the sandbox share dir
fs::path scratch(const std::string &name) {
  fs::path p = fs::path(nucqcd_test::share_dir()) / "grid_cache_test" / (name + "." + std::to_string(::getpid()));
  fs::remove_all(p);
  return p;
}

std::string contents(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

int files_in(const fs::path &dir) {
  int n = 0;
  for ([[maybe_unused]] const auto &e : fs::directory_iterator(dir)) n++;
  return n;
}
}  // namespace

TEST_CASE("grid_cache_dir: <share>/grids/ unless NUCQCD_GRID_DIR is set", "[gridcache]") {
  unsetenv("NUCQCD_GRID_DIR");
  CHECK(nucqcd::grid_cache_dir("/some/share") == "/some/share/grids/");
  setenv("NUCQCD_GRID_DIR", "", 1);  // empty counts as unset
  CHECK(nucqcd::grid_cache_dir("/some/share") == "/some/share/grids/");
  setenv("NUCQCD_GRID_DIR", "/scratch/cache", 1);  // a trailing slash is added
  CHECK(nucqcd::grid_cache_dir("/some/share") == "/scratch/cache/");
  setenv("NUCQCD_GRID_DIR", "/scratch/cache/", 1);
  CHECK(nucqcd::grid_cache_dir("/some/share") == "/scratch/cache/");
  unsetenv("NUCQCD_GRID_DIR");
}

TEST_CASE("write_file_atomically writes the file and creates missing directories", "[gridcache]") {
  const fs::path dir = scratch("write");
  const fs::path file = dir / "sub" / "grid.dat";
  REQUIRE(nucqcd::write_file_atomically(file.string(), [](std::ofstream &out) { out << "first"; }));
  CHECK(contents(file) == "first");
  // an existing file is replaced
  REQUIRE(nucqcd::write_file_atomically(file.string(), [](std::ofstream &out) { out << "second"; }));
  CHECK(contents(file) == "second");
  CHECK(files_in(dir / "sub") == 1);  // no temporary files left
  fs::remove_all(dir);
}

TEST_CASE("write_file_atomically leaves nothing behind when the writer throws", "[gridcache]") {
  const fs::path dir = scratch("throw");
  fs::create_directories(dir);
  const fs::path file = dir / "grid.dat";
  CHECK_THROWS_AS(nucqcd::write_file_atomically(file.string(),
                                                [](std::ofstream &out) {
                                                  out << "half a grid";
                                                  throw std::runtime_error("interrupted");
                                                }),
                  std::runtime_error);
  CHECK_FALSE(fs::exists(file));
  CHECK(files_in(dir) == 0);
  // an older complete file survives a failed rewrite
  REQUIRE(nucqcd::write_file_atomically(file.string(), [](std::ofstream &out) { out << "complete"; }));
  CHECK_THROWS(nucqcd::write_file_atomically(file.string(), [](std::ofstream &out) {
    out << "partial";
    throw std::runtime_error("interrupted");
  }));
  CHECK(contents(file) == "complete");
  CHECK(files_in(dir) == 1);
  fs::remove_all(dir);
}

TEST_CASE("write_file_atomically returns false when the directory cannot be created", "[gridcache]") {
  const fs::path dir = scratch("fail");
  fs::create_directories(dir);
  std::ofstream(dir / "plainfile") << "x";
  // a regular file is in the way of the directory
  CHECK_FALSE(nucqcd::write_file_atomically((dir / "plainfile" / "grid.dat").string(),
                                            [](std::ofstream &out) { out << "data"; }));
  CHECK(files_in(dir) == 1);
  fs::remove_all(dir);
}
