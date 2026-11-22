/*! \file GridCache.hpp
 * \brief Location of the grid cache and atomic writes of cache files
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * \details Glauber FSI grids, momentum-distribution grids and classical SCX grids are expensive to
 * compute, so they are written to a cache directory and read back on later runs. By default the
 * cache is the grids/ subdirectory of the share directory passed to the classes (usually HOMEDIR).
 * The environment variable NUCQCD_GRID_DIR overrides it, so the cache can live elsewhere (for
 * example on scratch space, or per job) while the input data are still read from share/.
 *
 * Cache files are written atomically: the data go to a temporary file in the same directory,
 * which is renamed to the final name only when it is complete. A crash or a full disk then never
 * leaves a truncated cache file, and concurrent jobs that need the same grid never read one that
 * is still being written (the last writer wins, with identical contents).
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef GRIDCACHE_HPP
#define GRIDCACHE_HPP

#include <fstream>
#include <functional>
#include <string>

namespace nucqcd {

/*! Directory of the grid cache, with a trailing '/'.
 * \param share_dir the share directory passed to the physics classes (usually HOMEDIR)
 * \return $NUCQCD_GRID_DIR if that is set and not empty, otherwise share_dir + "/grids/"
 */
std::string grid_cache_dir(const std::string &share_dir);

/*! Writes a file atomically: writer fills a temporary file in the same directory, which is then
 * renamed to path. Missing parent directories are created.
 * \param path final file name
 * \param writer writes the contents to the stream it is given (opened in binary mode)
 * \return true on success. On failure (directory or file cannot be created, write error) a
 * message goes to std::cerr, no file is left behind, and false is returned. Exceptions thrown by
 * writer are passed on, also without leaving a file behind.
 */
bool write_file_atomically(const std::string &path, const std::function<void(std::ofstream &)> &writer);

}  // namespace nucqcd

/** @} */
#endif
