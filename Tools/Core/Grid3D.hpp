/*! \file Grid3D.hpp
 * \brief A three-dimensional array that owns its (contiguous) storage, with the index syntax of a C array
 * \author Wim Cosyn
 * \date 26/11/2026
 *
 * \addtogroup Tools
 * @{
 */
#ifndef GRID3D_HPP
#define GRID3D_HPP

#include <cstddef>
#include <vector>

namespace nucqcd {

/*! The elements `[j][k]` of one `i` of a Grid3D: `plane[j]` is a pointer to the row of fixed `j`. Valid as long as the
 * Grid3D it came from is not resized or destroyed. \tparam T `T` for a plane that can be written, `const T` for read only */
template <class T>
class Grid3DPlane {
public:
  /*! \param first the element `[0][0]` of the plane \param n2 size of the last dimension */
  Grid3DPlane(T *first, std::size_t n2) : first_(first), n2_(n2) {}
  T *operator[](std::size_t j) const { return first_ + j * n2_; }  /*!< the row of fixed i and j */
private:
  T *first_;  /*!< the element [0][0] of the plane */
  std::size_t n2_;  /*!< size of the last dimension */
};

/*! \brief A three-dimensional array with `grid[i][j][k]` indexing, stored in one contiguous std::vector
 *
 * Replaces the `T ***` arrays of the Glauber grids that were allocated with a `new` for every row and freed with
 * matching loops: the elements are value-initialized (a complex number is 0), copying and moving work, nothing leaks if
 * a constructor throws, and the elements are contiguous. The storage order is that of a C array: the last index runs
 * fastest, so the iterators visit the elements in the order of the loops `for i, for j, for k` (the order of the cache
 * files of the grids). `grid[i][j]` is a pointer to the `k` row. There is no bounds checking, as with the pointers.
 */
template <class T>
class Grid3D {
public:
  Grid3D() : n0_(0), n1_(0), n2_(0) {}  /*!< an empty array */
  /*! \param n0,n1,n2 sizes of the three dimensions; the elements are value-initialized */
  Grid3D(std::size_t n0, std::size_t n1, std::size_t n2) : n0_(n0), n1_(n1), n2_(n2), data_(n0 * n1 * n2) {}

  /*! the elements of fixed i */
  Grid3DPlane<T> operator[](std::size_t i) { return Grid3DPlane<T>(data_.data() + i * n1_ * n2_, n2_); }
  /*! the elements of fixed i, read only */
  Grid3DPlane<const T> operator[](std::size_t i) const { return Grid3DPlane<const T>(data_.data() + i * n1_ * n2_, n2_); }

  std::size_t size0() const { return n0_; }  /*!< size of the first dimension */
  std::size_t size1() const { return n1_; }  /*!< size of the second dimension */
  std::size_t size2() const { return n2_; }  /*!< size of the third dimension */

  typename std::vector<T>::iterator begin() { return data_.begin(); }  /*!< first element, in the order of the loops over i, j, k */
  typename std::vector<T>::iterator end() { return data_.end(); }  /*!< after the last element */
  typename std::vector<T>::const_iterator begin() const { return data_.begin(); }  /*!< first element */
  typename std::vector<T>::const_iterator end() const { return data_.end(); }  /*!< after the last element */

private:
  std::size_t n0_;  /*!< size of the first dimension */
  std::size_t n1_;  /*!< size of the second dimension */
  std::size_t n2_;  /*!< size of the third dimension */
  std::vector<T> data_;  /*!< the elements, last index fastest */
};

}  // namespace nucqcd

/*! @} */
#endif
