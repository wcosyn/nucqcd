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
#include <Integrator/numint/numint.hpp>
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

  /*! use the coordinates around the momentum of the particle (TiltedFrame) for a single particle off the z axis (default),
   * or the spherical coordinates of the lab for it like for several particles (the code until 2026-11: up to 100 times slower
   * off axis, the cap on the number of evaluations is reached, and errors of a few 1e-3 at precision 1e-3; for comparisons) */
  void setUseTiltedFrame(bool use) {useTiltedFrame = use;}
  bool getUseTiltedFrame() const {return useTiltedFrame;} /*!< see setUseTiltedFrame */
  unsigned getEvalCount() const {return evalCount;} /*!< number of integrand evaluations of the last integrateTilted call (diagnostics) */
  const FsiCorrelator& getFsiCorrelator() const {return fsicorrelator;} /*!< the FsiCorrelator (SRC grid) */
  FsiCorrelator& getFsiCorrelator() {return fsicorrelator;} /*!< the FsiCorrelator (SRC grid) */
  MeanFieldNucleusThick *getPnucleusthick() const {return pnucleusthick;} /*!< the nucleus */

public:
  /*! The coordinates for a single particle in any direction: cylindrical (b, z, phi) around the momentum of the
   * particle, with phi measured from the direction of the b vector of the hit point. In these the step of the
   * Glauber phase (only the nucleons the particle still has to pass count) is a limit of the z integral and the
   * profile function a Gaussian in b, where in the spherical coordinates of the lab it is a cut along a plane and a
   * narrow tilted tube. (Only the A(e,e'p)-like case of one particle; for several particles there is no common frame.) */
  struct TiltedFrame {
    double u[3]; /*!< unit vector along the b vector of the hit point (any perpendicular one if that is zero) */
    double w[3]; /*!< unit vector perpendicular to the momentum and to u */
  };
  /*! what the integrands in the frame of the particle have in common */
  struct TiltedPoint {
    double dens; /*!< density times the Jacobian (b / r^2), without the profile function */
    double fsi; /*!< dens times the Gaussian profile: the integrand of the plain FSI factor, without the scattering amplitude */
    double src; /*!< the SRC factor at the point (correlation function times the correlation grid) */
  };
  /*! computes TiltedPoint at b, z, phi (frame coordinates) for the single particle of the grid
   * \return false if the point is outside the nucleus (r larger than its range): then the integrand is zero */
  static bool tiltedPoint(GlauberThickGridCommon &grid, double b, double z, double phi, int proton, const TiltedFrame &frame, TiltedPoint &point);
  /*! the FSI factor from the results of integrateTilted
   * \param proton proton or neutron \param ret results of integrateTilted \param index which result (odd: an SRC variant, which uses N and the SRC factor at the hit point) \param src the SRC factor at the hit point
   * (getCorrGrid_interp there) */
  std::complex<double> tiltedPhase(int proton, const std::vector<double> &ret, int index, double src);
  /*! integrand of the SRC normalization N (see integrateTilted) */
  static void srcNormIntegrand(numint::vector_d &results, double b, double z, double phi, GlauberThickGridCommon &grid, int proton, const TiltedFrame &frame);
  /*! integrand function in the frame of the particle: results (the number is that of the integrand) at the point b, z, phi */
  typedef void (*TiltedIntegrand)(numint::vector_d &results, double b, double z, double phi, GlauberThickGridCommon &grid, int proton, const TiltedFrame &frame);
  /*! integrates a TiltedIntegrand over the nucleus for the single particle of the grid and the point that is set (the hit point)
   * \param proton proton or neutron \param f integrand \param nresults number of results of f
   * \param maxeval maximum number of evaluations
   * \param ret [out] nresults+1 values: the integrals of f, then the SRC normalization N = integral of dens * src over the
   * whole nucleus (z from -range to range). The SRC variants of the FSI factor are (N - amplitude * integral of f) * src(hit point), the
   * others 1 - amplitude * integral of f. N is not exactly 1/src(hit point), and the case of a particle along the z axis
   * (klaas_int_bound) takes the term to be 1.
   * The z limits of f run from the hit point to the edge of the nucleus (to the other side for an incoming particle). */
  void integrateTilted(int proton, TiltedIntegrand f, int nresults, int maxeval, std::vector<double> &ret);

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
  bool useTiltedFrame = true; /*!< see setUseTiltedFrame */
  unsigned evalCount = 0; /*!< number of integrand evaluations of the last integrateTilted call */

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
