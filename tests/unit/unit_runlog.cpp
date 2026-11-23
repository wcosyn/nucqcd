// MePhys: run bookkeeping (Tools/Core/RunLog, roadmap step 11): comment marker, quoting, argument names, end
// of the log (also at exit()), MPI rank, build information.
#include <test_support.hpp>

#include <Tools/Core/RunLog.hpp>

#include <sys/wait.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::string> lines(const std::string &s) {
  std::vector<std::string> v;
  std::istringstream in(s);
  for (std::string l; std::getline(in, l);) v.push_back(l);
  return v;
}
bool contains(const std::string &s, const std::string &part) { return s.find(part) != std::string::npos; }
}  // namespace

TEST_CASE("RunLog: comment lines, quoted command line, named arguments, end with wall time", "[runlog]") {
  char a0[] = "./prog", a1[] = "AV18", a2[] = "two words", a3[] = "it's";
  char *argv[] = {a0, a1, a2, a3, nullptr};
  std::ostringstream out;
  {
    nucqcd::RunLog log(4, argv, {"wave function", "label"}, &out);
    CHECK_FALSE(contains(out.str(), "# finished"));
  }
  const std::string s = out.str();
  for (const std::string &l : lines(s)) CHECK(l.rfind("# ", 0) == 0);
  CHECK(contains(s, "# command: ./prog AV18 'two words' 'it'\\''s'\n"));
  CHECK(contains(s, "#   wave function = AV18\n"));
  CHECK(contains(s, "#   label = two words\n"));
  CHECK(contains(s, "#   argument 3 = it's\n"));  // more arguments than names
  CHECK(contains(s, "# started: "));
  CHECK(contains(s, "# host: "));
  CHECK(contains(s, "# nucqcd: git " + nucqcd::build_info::git_describe()));
  CHECK(contains(s, "# libraries: ROOT "));
  CHECK(contains(s, ", wall time "));
  CHECK(s.find("# finished") == s.rfind("# finished"));  // written once
}

TEST_CASE("RunLog: missing arguments are marked, finish() writes the end once", "[runlog]") {
  char a0[] = "prog";
  char *argv[] = {a0, nullptr};
  std::ostringstream out;
  {
    nucqcd::RunLog log(1, argv, {"Q2", "x"}, &out);
    log.finish();
  }
  const std::string s = out.str();
  CHECK(contains(s, "#   Q2 = (missing)\n"));
  CHECK(contains(s, "#   x = (missing)\n"));
  CHECK(s.find("# finished") == s.rfind("# finished"));
}

TEST_CASE("RunLog: the end is written when the program calls exit()", "[runlog]") {
  const std::string file = "runlog_exit_test.txt";
  const pid_t pid = fork();
  REQUIRE(pid >= 0);
  if (pid == 0) {  // child: exit() does not run the destructor of the RunLog
    // silence the child's copy of the Catch2 session, which reports at exit() too
    if (std::FILE *null = std::fopen("/dev/null", "w")) {
      dup2(fileno(null), STDOUT_FILENO);
      dup2(fileno(null), STDERR_FILENO);
    }
    // a stream that is never destroyed: it must still exist when the exit handler writes to it
    std::ofstream &out = *new std::ofstream(file);
    char a0[] = "prog";
    char *argv[] = {a0, nullptr};
    nucqcd::RunLog *log = new nucqcd::RunLog(1, argv, {}, &out);  // never destroyed either
    (void)log;
    std::exit(0);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  REQUIRE(WIFEXITED(status));
  std::ifstream in(file);
  const std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  CHECK(contains(s, "# command: prog"));
  CHECK(contains(s, "# finished: "));
  std::remove(file.c_str());
}

TEST_CASE("RunLog: only MPI rank 0 writes", "[runlog]") {
  char a0[] = "prog";
  char *argv[] = {a0, nullptr};
  setenv("OMPI_COMM_WORLD_RANK", "1", 1);
  CHECK(nucqcd::RunLog::mpiRank() == 1);
  std::ostringstream out;
  { nucqcd::RunLog log(1, argv, {}, &out); }
  CHECK(out.str().empty());
  unsetenv("OMPI_COMM_WORLD_RANK");
  CHECK(nucqcd::RunLog::mpiRank() == 0);
}

TEST_CASE("RunLog: build information", "[runlog]") {
  CHECK_FALSE(nucqcd::build_info::git_describe().empty());
#ifdef NUCQCD_TEST_HAVE_GIT_CHECKOUT
  CHECK(nucqcd::build_info::git_describe().rfind("unknown", 0) != 0);
#endif
  CHECK_FALSE(nucqcd::build_info::build_type().empty());
  CHECK_FALSE(nucqcd::build_info::compiler().empty());
  CHECK(contains(nucqcd::build_info::libraries(), "GSL "));
}
