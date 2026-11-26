/*! \file Utilfunctions.hpp
 * \brief Some useful utilities functions, mostly math related
 * \author Wim Cosyn
 * \date 16/08/2011
 */

#ifndef UTILFUNCTIONS_H
#define UTILFUNCTIONS_H


#define SIGN(x) ((x>0)-(x<0)) /*!< sign of x: -1, 0 or 1 */

#include <cstdarg>
#include <Tools/Core/ErrorMessage.hpp>
#include <sstream>
#include <complex>
#include <Tools/Core/Grid3D.hpp>

namespace nucqcd {


// Bookkeep(argc, argv, arg_names) was replaced by nucqcd::RunLog (RunLog.hpp) in 2026-11.

/*! interpolation of an array
 * \param array array to interpolate
 * \param r interpolation value
 * \param rstep gridsize of the array
 * \param lines size of the array
 * \param offset units of rstep where the grid starts
 * \return interpolation value of the grid at r
 */
template <class T> T interpolate(const T *array, double r, double rstep,int lines, int offset){
  //interpolation points is double of order
  //rstep is set distance between known points
  //offset indicates times rstep from where the interpolation array starts
  //lines is number of elements in array

  if(array==NULL) {
    throw std::logic_error(nucqcd::ErrorMessage() << "Array has not been initialized!!!!" );
  }
  int index=int(floor(r/rstep))-offset;

  if(index==lines-1) index--; // if index is at the end of the array.
  if(index<0) index=0;
  double a = (r-rstep*(index+1+offset))/-rstep;
  return a*(array[index]-array[index+1]) + array[index+1];
}


double power(double x,int n);/*!< returns x^n */
int power(int x,int n);/*!< returns x^n */
/*! returns the distance between points r1 and r2
 * \param r1 r value of point1
 * \param costheta1 cos(theta) value of point1
 * \param sintheta1 sin(theta) value of point1
 * \param cosphi1 cos(phi) value of point1
 * \param sinphi1 sin(phi) value of point1
 * \param r2 r value of point2
 * \param costheta2 cos(theta) value of point2
 * \param sintheta2 sin(theta) value of point2
 * \param cosphi2 cos(phi) value of point2
 * \param sinphi2 sin(phi) value of point2
 */
double normr(double r1, double costheta1, double sintheta1, double cosphi1, double sinphi1,
	     double r2, double costheta2, double sintheta2, double cosphi2, double sinphi2);
/*! solves matrix equation a*x=b with LU solve algorithm and store solution in b, from numerical recipes p46
 * \param a square matrix a
 * \param b column matrix b
 * \param order number of equations
 */
void LUsolve(double **a, double *b, int order);



// void sincos(const double x, double * sin, double * cos);

/*! performs a fit to data with a set of fitfunctions, also computes the chi squared of the fit
 * \param xarray data x array
 * \param yarray data y array
 * \param sigma data y errors
 * \param lower lower index of the data array for the fit
 * \param upper upper index of the data array for the fit
 * \param order order of the fit (in fitfunctions)
 * \param chi2 chi squard of the fit
 */
double* fitdata(double *xarray, double *yarray, double *sigma, int lower, int upper, int order, double & chi2);
/*!set of base functions used for the fit
 * \param x x value
 * \param i order index of the base function set
 * \return value of f_i(x)
 */
double fitfunction(double x, int i);
double **matrix(int n, int m); /*!< create a n*m matrix on the heap*/
void freematrix(double **matrix, int n); /*!< free a n*? matrix on the heap*/


/*! 3d interpolation of an array
 */
// template <class T> T Interp3d(T ***grid, double s, double t, double u, 
// 			  double comps, double compt, double compu, int sindex, int tindex, int uindex){
//   return compu*(compt*(comps*grid[sindex][tindex][uindex]
// 	+ s*grid[sindex+1][tindex][uindex]) 
// 	+ t*(comps*grid[sindex][tindex+1][uindex] 
// 	+ s*grid[sindex+1][tindex+1][uindex]))
// 	+  u*(compt*(comps*grid[sindex][tindex][uindex+1]
// 	+ s*grid[sindex+1][tindex][uindex+1]) 
// 	+ t*(comps*grid[sindex][tindex+1][uindex+1]
// 	+ s*grid[sindex+1][tindex+1][uindex+1]));
// }
/// trilinear interpolation in a 3D grid: weights s, t, u (and comps = 1-s, ...) between the points with indices (sindex, tindex, uindex) and the next ones
std::complex<double> Interp3d(const Grid3D<std::complex<double> > &grid, double s, double t, double u, 
			  double comps, double compt, double compu, int sindex, int tindex, int uindex);
/// trilinear interpolation in a real 3D grid, see the complex version
double Interp3d(double ***grid, double s, double t, double u, 
			  double comps, double compt, double compu, int sindex, int tindex, int uindex);
/*!convert a certain type to a string
 * \tparam T type to convert to string
 * \return returns a string
 */
template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

}  // namespace nucqcd

#endif

