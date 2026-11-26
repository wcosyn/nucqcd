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
#include <Glauber/FsiCorrelator.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Tools/Core/Grid3D.hpp>

#include <complex>
#include <fstream>
#include <string>
#include <vector>

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
 * It also holds what the thickness approximation adds to the nucleus: the MeanFieldNucleusThick and the FsiCorrelator
 * (the SRC grid). The grids of a class are numbered as in getFsiGridN_interp, see getGridNumber for the variants.
 *
 * Until 2026-11 there was a layer of abstract classes between AbstractFsiCTGrid and the two grid classes (AbstractFsiGridThick,
 * AbstractFsiDecayGridThick, AbstractFsiCTGridThick, AbstractFsiCTDecayGrid, AbstractFsiCTDecayGridThick: one named
 * interface per variant, repeated, with virtual inheritance); the single chain is now
 * AbstractFsiGrid <- AbstractFsiCTGrid <- GlauberThickGridCommon <- GlauberGridThick or GlauberDecayGridThick.
 */
class GlauberThickGridCommon : public AbstractFsiCTGrid{
public:
  virtual ~GlauberThickGridCommon(); /*!< Destructor, frees the arrays */
  /*! the interpolated value of grid number grid (see the class description) at the point set with setinterp
   * \throw std::invalid_argument for a grid number out of range */
  virtual std::complex<double> getFsiGridN_interp(int grid);
  /*! prints grid number grid (see the class description) on all grid points */
  virtual void print_grid(int grid);
  /*! the grid number of a variant: the variants without CT are numbered 0 .. variants-1 (bit 0 of the variant: SRC, bit 1:
   * decay), the same with CT follow; -1 for a variant that needs more variants than the class has (decay for
   * GlauberGridThick) */
  virtual int getGridNumber(GridVariant variant) const;

  const FsiCorrelator& getFsiCorrelator() const {return fsicorrelator;} /*!< the FsiCorrelator (SRC grid) */
  FsiCorrelator& getFsiCorrelator() {return fsicorrelator;} /*!< the FsiCorrelator (SRC grid) */
  MeanFieldNucleusThick *getPnucleusthick() const {return pnucleusthick;} /*!< the nucleus */

protected:
  /*! Constructor
   * \param nvariants number of variants in the first index of the grid arrays (2 or 4)
   * \param r_grid,cth_grid,phi_grid,pnuclthick,prec,integrator,dir as for GlauberGridThick */
  GlauberThickGridCommon(int nvariants, int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnuclthick,
                         double prec, IntegrationMethod integrator, std::string dir);

  /*! interpolated value of one variant of the grid at the point set with setinterp: the proton and neutron
   * factors to the powers of the numbers of nucleons that are not knocked out
   * \param array fsi_grid or fsi_ct_grid \param variant index in its first dimension
   * \throw std::logic_error if the grids are not filled */
  std::complex<double> interpolateVariant(const std::vector<std::vector<Grid3D<std::complex<double> > > > &array, int variant);

  /*! set the file names of the grids: the grid cache directory with the prefix "Thick." added \param dir cache directory plus the prefix of the derived class */
  virtual void setFilenames(std::string dir);

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
  MeanFieldNucleusThick *pnucleusthick; /*!< the nucleus */
  FsiCorrelator fsicorrelator; /*!< the SRC correlation grid */
  std::vector<std::vector<Grid3D<std::complex<double> > > > fsi_grid; /*!< grid that contains the fsi factor, `[variant][nucleon: 0 neutron, 1 proton][r][cos theta][phi]` */
  std::vector<std::vector<Grid3D<std::complex<double> > > > fsi_ct_grid; /*!< grid that contains the fsi+ct factor, same indices */
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
