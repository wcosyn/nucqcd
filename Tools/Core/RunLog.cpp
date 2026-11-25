#include <Tools/Core/RunLog.hpp>

#include <Integrator/numint/adaptive_stats.hpp>

#include "nucqcd_build_config.inc"  // NUCQCD_BUILD_TYPE, ... (generated at configure time)
#include "nucqcd_git_info.inc"      // NUCQCD_GIT_DESCRIBE, NUCQCD_GIT_DATE (generated at build time)

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>

namespace nucqcd {

namespace build_info {
std::string git_describe() { return NUCQCD_GIT_DESCRIBE; }
std::string git_date() { return NUCQCD_GIT_DATE; }
std::string build_type() { return NUCQCD_BUILD_TYPE; }
std::string compiler() { return NUCQCD_COMPILER; }
std::string flags() { return NUCQCD_FLAGS; }
std::string libraries() { return NUCQCD_LIBRARIES; }
}  // namespace build_info

namespace {
// RunLogs that still have to write their end at exit() (exit does not run the destructors of main's locals).
// The registry and its mutex are never destroyed: finish_at_exit runs during exit(), interleaved with the
// destructors of static objects in reverse order of registration, and the atexit handler is registered
// before a function-local static set would be constructed, so such a set would already be gone (seen as
// a crash at exit() in builds without ROOT).
struct Registry {
  std::mutex mutex;
  std::set<RunLog *> logs;
};
Registry &registry() {
  static Registry *r = new Registry;  // intentionally leaked
  return *r;
}
void finish_at_exit() {
  Registry &r = registry();
  std::lock_guard<std::mutex> lock(r.mutex);
  for (RunLog *log : r.logs) log->finish();
}

std::string now_string() {
  const std::time_t t = std::time(nullptr);
  std::tm tm{};
  localtime_r(&t, &tm);
  char buf[64];
  std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S %z (%Z)", &tm);
  return buf;
}

std::string quote(const std::string &arg) {
  if (!arg.empty() && arg.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+=.,/:@%^") ==
                          std::string::npos)
    return arg;
  std::string q = "'";
  for (char c : arg) q += (c == '\'') ? std::string("'\\''") : std::string(1, c);
  return q + "'";
}
}  // namespace

std::string RunLog::quotedCommandLine(int argc, char *argv[]) {
  std::string s;
  for (int i = 0; i < argc; i++) s += (i ? " " : "") + quote(argv[i] ? argv[i] : "");
  return s;
}

int RunLog::mpiRank() {
  for (const char *var : {"OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK", "MV2_COMM_WORLD_RANK", "SLURM_PROCID"}) {
    if (const char *v = std::getenv(var)) return std::atoi(v);
  }
  return 0;
}

RunLog::RunLog(int argc, char *argv[], std::vector<std::string> arg_names, std::ostream *stream)
    : out(mpiRank() == 0 ? (stream ? stream : &std::cout) : nullptr), start(std::chrono::steady_clock::now()) {
  {
    Registry &r = registry();
    std::lock_guard<std::mutex> lock(r.mutex);
    static const bool registered = (std::atexit(finish_at_exit) == 0);
    (void)registered;
    r.logs.insert(this);
  }
  if (!out) return;
  std::ostream &o = *out;
  o << "# nucqcd run log\n";
  o << "# command: " << quotedCommandLine(argc, argv) << "\n";
  const int nargs = argc > 0 ? argc - 1 : 0;
  const int n = std::max<int>(nargs, arg_names.size());
  for (int i = 0; i < n; i++) {
    o << "#   " << (i < int(arg_names.size()) ? arg_names[i] : "argument " + std::to_string(i + 1)) << " = "
      << (i < nargs ? std::string(argv[i + 1]) : std::string("(missing)")) << "\n";
  }
  o << "# started: " << now_string() << "\n";
  char host[256] = "unknown";
  gethostname(host, sizeof host - 1);
  char cwd[4096] = "unknown";
  if (!getcwd(cwd, sizeof cwd)) std::snprintf(cwd, sizeof cwd, "unknown");
  o << "# host: " << host << ", working directory: " << cwd << "\n";
  for (const char *var : {"NUCQCD_GRID_DIR"})
    if (const char *v = std::getenv(var)) o << "# environment: " << var << "=" << v << "\n";
  o << "# nucqcd: git " << build_info::git_describe() << " (" << build_info::git_date() << "), "
    << build_info::build_type() << " build, " << build_info::compiler() << ", flags: " << build_info::flags() << "\n";
  o << "# libraries: " << build_info::libraries() << "\n";
  o << std::flush;
}

void RunLog::finish() {
  if (finished) return;
  finished = true;
  if (!out) return;
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::ostringstream wall;
  wall << std::fixed << std::setprecision(1) << seconds;
  numint::adaptive_report(*out, "# ");  // only if cube_adaptive was used
  *out << "# finished: " << now_string() << ", wall time " << wall.str() << " s" << std::endl;
}

RunLog::~RunLog() {
  finish();
  Registry &r = registry();
  std::lock_guard<std::mutex> lock(r.mutex);
  r.logs.erase(this);
}

}  // namespace nucqcd
