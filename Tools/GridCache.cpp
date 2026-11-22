#include "GridCache.hpp"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <system_error>
#include <unistd.h>

namespace fs = std::filesystem;

namespace nucqcd {

std::string grid_cache_dir(const std::string &share_dir) {
  const char *env = std::getenv("NUCQCD_GRID_DIR");
  std::string dir = (env != nullptr && *env != '\0') ? std::string(env) : share_dir + "/grids/";
  if (dir.back() != '/') dir += '/';
  return dir;
}

bool write_file_atomically(const std::string &path, const std::function<void(std::ofstream &)> &writer) {
  const fs::path target(path);
  std::error_code ec;
  if (target.has_parent_path()) {
    fs::create_directories(target.parent_path(), ec);
    if (ec) {
      std::cerr << "could not create directory " << target.parent_path() << " for " << path << ": "
                << ec.message() << std::endl;
      return false;
    }
  }
  // unique per process and per call, in the same directory so that the rename is atomic
  static std::atomic<unsigned long> counter{0};
  const std::string tmp =
      path + ".tmp." + std::to_string(::getpid()) + "." + std::to_string(counter++);

  bool ok = false;
  try {
    std::ofstream out(tmp, std::ios::out | std::ios::binary | std::ios::trunc);
    if (out.is_open()) {
      writer(out);
      out.close();
      ok = !out.fail();
    }
  } catch (...) {
    fs::remove(tmp, ec);
    throw;
  }
  if (!ok) {
    std::cerr << "could not write cache file " << path << std::endl;
    fs::remove(tmp, ec);
    return false;
  }
  fs::rename(tmp, target, ec);
  if (ec) {
    std::cerr << "could not rename " << tmp << " to " << path << ": " << ec.message() << std::endl;
    fs::remove(tmp, ec);
    return false;
  }
  return true;
}

}  // namespace nucqcd
