/*! \file RunLog.hpp
 * \brief Run bookkeeping for the programs: how a result was produced (class nucqcd::RunLog, namespace nucqcd::build_info)
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * \addtogroup MePhys
 * @{
 */
#ifndef RUNLOG_HPP
#define RUNLOG_HPP

#include <chrono>
#include <iosfwd>
#include <string>
#include <vector>

namespace nucqcd {

/*! \brief Information about the build of nucqcd, fixed when the library was built (not at run time) */
namespace build_info {
std::string git_describe();  ///< `git describe --always --dirty` of the source tree at build time ("-dirty": uncommitted changes)
std::string git_date();      ///< commit date (ISO 8601) of the built commit
std::string build_type();    ///< CMake build type (Release, Debug, ...)
std::string compiler();      ///< C++ compiler and version
std::string flags();         ///< compiler flags of the build type
std::string libraries();     ///< versions of ROOT, GSL and (if used) LHAPDF
}  // namespace build_info

/*! \brief Writes how a program was run, as comment lines starting with "# "
 *
 * Construct it as the first statement of main:
 * \code
 * int main(int argc, char *argv[]) {
 *   nucqcd::RunLog runlog(argc, argv, {"wave function", "Q2 [GeV^2]"});
 *   ...
 * \endcode
 * It writes the full command line (quoted so that it can be pasted back into a shell), the arguments
 * with the given names, date and time with time zone, host, working directory, the nucqcd environment
 * variables that are set (NUCQCD_GRID_DIR), and the build (git describe with "-dirty" for uncommitted
 * changes, build type, compiler, flags, library versions; see build_info). When the program ends (end
 * of main, or exit()) it writes the end time and the wall time.
 *
 * The lines go to stdout by default, where most programs write their data, so the provenance stays in
 * the data file; every line starts with "# ", which gnuplot, numpy.loadtxt and pandas (comment='#') skip.
 * With MPI only rank 0 writes (the rank is read from the environment of the MPI launcher).
 *
 * Argument names are optional: fewer names than arguments prints the others as "argument <i>", more names
 * than arguments marks the missing ones. (It replaces the former Bookkeep function of Utilfunctions,
 * which read past its array of names when a program was called with extra arguments, ran git at run
 * time, and wrote without a comment marker.)
 */
class RunLog {
public:
  /*! Writes the start of the log
   * \param argc argc of main \param argv argv of main
   * \param arg_names names of the arguments argv[1], argv[2], ... (optional)
   * \param out stream to write to (default std::cout); if the program may end with exit(), it must
   * still exist then (std::cout does) */
  RunLog(int argc, char *argv[], std::vector<std::string> arg_names = {}, std::ostream *out = nullptr);
  /*! Writes the end of the log (end time, wall time) */
  ~RunLog();
  RunLog(const RunLog &) = delete;
  RunLog &operator=(const RunLog &) = delete;

  /*! \return the command line with each argument quoted for a POSIX shell where needed */
  static std::string quotedCommandLine(int argc, char *argv[]);
  /*! \return the MPI rank from the launcher's environment (Open MPI, MPICH/Intel MPI, Slurm), 0 without MPI */
  static int mpiRank();

  void finish(); ///< writes the end of the log now (once; called by the destructor and at exit())

private:
  std::ostream *out;  ///< stream, or null when this process does not write (MPI rank > 0)
  std::chrono::steady_clock::time_point start; ///< start time, for the wall time
  bool finished = false; ///< the end has been written
};

}  // namespace nucqcd
/*! @} */
#endif
