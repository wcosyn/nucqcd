#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace nucqcd_test {

namespace {

std::string env_or(const char *name, const char *fallback) {
  const char *v = std::getenv(name);
  return (v && *v) ? std::string(v) : std::string(fallback);
}

std::map<std::string, std::string> read_ref_file(const std::string &path) {
  std::map<std::string, std::string> result;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    std::string key, value;
    if (ls >> key >> value) result[key] = value;
  }
  return result;
}

std::string format_double(double v) {
  std::ostringstream os;
  os << std::setprecision(17) << v;
  return os.str();
}

// RAII exclusive lock on <path>.lock: CTest may run several test cases of the same executable
// in parallel, and in update mode they all rewrite the same .ref file.
class FileLock {
public:
  explicit FileLock(const std::string &path) : fd(::open((path + ".lock").c_str(), O_RDWR | O_CREAT, 0644)) {
    if (fd >= 0) ::flock(fd, LOCK_EX);
  }
  ~FileLock() {
    if (fd >= 0) {
      ::flock(fd, LOCK_UN);
      ::close(fd);
    }
  }
  FileLock(const FileLock &) = delete;
  FileLock &operator=(const FileLock &) = delete;

private:
  int fd;
};

}  // namespace

std::string share_dir() { return env_or("NUCQCD_TEST_SHARE", NUCQCD_TEST_SHARE_DEFAULT); }

bool updating_references() { return env_or("NUCQCD_UPDATE_REFERENCE", "0") == "1"; }

RefStore::RefStore(const std::string &name)
    : path(env_or("NUCQCD_TEST_REFDIR", NUCQCD_TEST_REFDIR_DEFAULT) + "/" + name + ".ref"),
      values(read_ref_file(path)) {}

void RefStore::store(const std::string &key, double value) {
  FileLock lock(path);
  auto current = read_ref_file(path);  // re-read: other processes may have added keys
  current[key] = format_double(value);
  const std::string tmp = path + ".tmp";
  {
    std::ofstream out(tmp);
    out << "# Golden values for the nucqcd regression tests (key value).\n"
        << "# Captured from the code on master; regenerate with NUCQCD_UPDATE_REFERENCE=1.\n";
    for (const auto &kv : current) out << kv.first << ' ' << kv.second << '\n';
  }
  std::rename(tmp.c_str(), path.c_str());
  values = current;
}

RefResult RefStore::compare(const std::string &key, double value, double rel_tol, double abs_tol) {
  std::ostringstream msg;
  msg << std::setprecision(17) << "reference '" << key << "' in " << path << ": computed " << value;
  if (updating_references()) {
    store(key, value);
    msg << " (stored, update mode)";
    return {std::isfinite(value), msg.str()};
  }
  auto it = values.find(key);
  if (it == values.end()) {
    msg << " -- no reference value stored; run with NUCQCD_UPDATE_REFERENCE=1 to create it";
    return {false, msg.str()};
  }
  const double ref = std::stod(it->second);
  const double diff = std::fabs(value - ref);
  const double tol = std::max(abs_tol, rel_tol * std::fabs(ref));
  const bool ok = std::isfinite(value) && diff <= tol;
  msg << ", reference " << ref << ", |diff| " << diff << ", tolerance " << tol;
  if (!ok)
    msg << "\n  -> The output differs from the value captured from master. This can be a new bug,"
           "\n     or a deliberate change/fix (the old value may itself have been wrong): investigate,"
           "\n     and if the new value is right re-baseline with NUCQCD_UPDATE_REFERENCE=1.";
  return {ok, msg.str()};
}

}  // namespace nucqcd_test
