/*! \file GlauberThickGridCommon.hpp
 * \brief The part that GlauberGridThick and GlauberDecayGridThick share: the grid arrays and everything that is done to them
 * \author Wim Cosyn
 * \date 25/11/2026
 *
 * \addtogroup Glauber
 * @{
 */
#ifndef GLAUBERTHICKGRIDCOMMON_H
#define GLAUBERTHICKGRIDCOMMON_H

#include <Glauber/AbstractFsiCTGrid.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>

#include <complex>
#include <fstream>
#include <string>

namespace nucqcd {

/*! \brief The grids of a Glauber FSI grid in the thickness approximation, with and without color transparency
 *
 * GlauberGridThick (variants: RMSGA, RMSGA+SRC) and GlauberDecayGridThick (RMSGA, +SRC, +decay, +SRC+decay)
 * had the same code for everything except the integrals of a grid point: they were copies of each other
 * (738 of 1005 lines identical). This class holds that code, for a number of "variants" given to the
 * constructor, in two arrays fsi_grid and fsi_ct_grid (without and with color transparency) with the index
 * `[variant][nucleon: 0 neutron, 1 proton][r][cos theta][phi]`:
 * - allocation, construction of all grids or of the CT grids only (the loops over the grid points, the
 *   shortcuts for r=0 and cos theta = +-1, and the threshold beyond which a point is one),
 * - the interpolation of one variant, the interpolation by grid number and the printing,
 * - reading and writing the cache files.
 * The derived classes provide the integral of one grid point (calcGlauberphasesBoth, calcGlauberphasesCt),
 * the interpolation functions that name the variants, and the file names.
 *
 * Grid numbers (getFsiGridN_interp, print_grid): 0 .. variants-1 are the variants without CT, then the same
 * variants with CT (fsi_grid[g], fsi_ct_grid[g - variants]).
 *
 * The class derives virtually from AbstractFsiCTGrid, so the derived class has one AbstractFsiGrid and one
 * AbstractFsiCTGrid, whichever abstract class it also derives from (AbstractFsiCTGridThick or
 * AbstractFsiCTDecayGridThick), and must initialize both of them in its constructor.
 */
class GlauberThickGridCommon : public virtual AbstractFsiCTGrid{
public:
  virtual ~GlauberThickGridCommon(); /*!< Destructor, frees the arrays */
  /*! the interpolated value of grid number grid (see the class description) at the point set with setinterp
   * \throw std::invalid_argument for a grid number out of range */
  virtual std::complex<double> getFsiGridN_interp(int grid);
  /*! prints grid number grid (see the class description) on all grid points */
  virtual void print_grid(int grid);

protected:
  /*! Constructor
   * \param nvariants number of variants in the first index of the grid arrays (2 or 4)
   * \param r_grid,cth_grid,phi_grid,pnuclthick,prec,integrator,dir as for GlauberGridThick */
  GlauberThickGridCommon(int nvariants, int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnuclthick,
                         double prec, int integrator, std::string dir);

  /*! interpolated value of one variant of the grid at the point set with setinterp: the proton and neutron
   * factors to the powers of the numbers of nucleons that are not knocked out
   * \param array fsi_grid or fsi_ct_grid \param variant index in its first dimension
   * \throw std::logic_error if the grids are not filled */
  std::complex<double> interpolateVariant(std::complex<double> *****array, int variant);

  virtual void constructAllGrids(); /*!< construct both the fsi and the fsi+ct grids */
  virtual void constructCtGrid(); /*!< construct only the fsi+ct grids */
  virtual void readinFsiGrid(std::ifstream &infile); /*!< read in both the fsi and the fsi+ct grids */
  virtual void readinFsiCtGrid(std::ifstream &infile); /*!< read in only the fsi+ct grids */
  virtual void writeoutFsiGrid(std::ofstream &outfile); /*!< write out both the fsi and the fsi+ct grids */
  virtual void writeoutFsiCtGrid(std::ofstream &outfile); /*!< write out only the fsi+ct grids */
  /*! calculates the glauberphases for one gridpoint (both FSI and FSI+CT), for all variants
   * \param i grid index in r \param j grid index in costheta \param k grid index in phi */
  virtual void calcGlauberphasesBoth(const int i, const int j, const int k)=0;
  /*! calculates the glauberphases for one gridpoint (only FSI+CT), for all variants
   * \param i grid index in r \param j grid index in costheta \param k grid index in phi */
  virtual void calcGlauberphasesCt(const int i, const int j, const int k)=0;

  const int variants; /*!< number of variants: size of the first index of the arrays */
  std::complex<double> *****fsi_grid; /*!< grid that contains the fsi factor, `[variant][nucleon: 0 neutron, 1 proton][r][cos theta][phi]` */
  std::complex<double> *****fsi_ct_grid; /*!< grid that contains the fsi+ct factor, same indices */
  int ** treshold; /*!< array that checks if calculated glauberphases are close to one, then doesn't compute them for larger r, to save computing time*/
  double error; /*!< integration error estimate of the cubature integrators */

private:
  /*! the loop over the grid points of constructAllGrids and constructCtGrid
   * \param ctOnly only the fsi+ct grids */
  void constructGrids(bool ctOnly);
  /*! true when the fsi+ct (and, if not ctOnly, the fsi) factors at a grid point are all within 1e-6 of one
   * (to the power of the numbers of neutrons and protons) */
  bool closeToOne(bool ctOnly, int i, int j, int k) const;
  GlauberThickGridCommon(const GlauberThickGridCommon &) = delete;
  GlauberThickGridCommon &operator=(const GlauberThickGridCommon &) = delete;
};

}  // namespace nucqcd

/** @} */
#endif
