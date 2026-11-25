#ifndef CLASS_GRIDTHICK_SCX_HPP
#define CLASS_GRIDTHICK_SCX_HPP

#include <vector>
#include <cmath>
#include <Glauber/FastParticle.hpp>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_integration.h>
#include <fstream>
#include <string>

namespace nucqcd {


/** A class to calculate charge exchange
 *  using the mean field approximation.
 *  Using a semiclassical formalism (bye bye glauber)
 *  although it resembles the glauber formalism.
 *  The grid contains the probability for interaction.
 *  (The stored grid does NOT, only storing the energy
 *   independent integral).
 *  
 *  NOTE: For the same density function but different normalisations
 *        (e.g. normed to A,A-1,A-2,...) you should not calculate different
 *        grids constructed with different density function (_densf)'s but 
 *        simply set the _dens_fctr accordingly via set_dens_fctr().
 *        
 *        For example MeanFieldNucleusThick densities are normed to 1. in Wims Code.
 *        To get for example a density normed to total protons do
 *           $> get_dens_fctr( MeanFieldNucleusThick::getZ() );
 * 
 *  Class deliberately kept fairly modular. Only uses FastParticle.
 *  Please keep it that way. No hard coded densities etc. plz, thx.
 * 
 *  \author Camille Colle (Ghent University)
 *  \date 21/09/2015
 *
 *  */
class ClassGridThick_SCX {
	public:
		/** make sure both bpoints and zpoints are always larger than 2 plz **/
		ClassGridThick_SCX(double (*densf)(double,void*),void* p,const char* densName,double rmax, FastParticle* fp,int bpoints, int zpoints,const char* SHAREDIR);
		void constructGrid();
		void printGrid(); /**< prints the grid and its error estimates (b index, z index, value, error), for debugging */
		void printDensityGrid(); /**< prints the density at the grid points (b index, z index, density), for debugging */
		void calcP(double b,double z,double& ret,double& err,double (*)(double,void*) );
		double getInterp(double[],double& res,double& err); // remember x is r,cos\theta and \phi
		
		double get_dens_fctr(){ return _dens_fctr; } /**< correction factor for the normalization of the density */
        /** sets the correction factor for the normalization of the density; call it before constructGrid \param d factor */
        void   set_dens_fctr(double d){ assert(_grid.empty()); _dens_fctr = d; };
		
		
	private:
		double (*_densf)(double r,void* p); /**< density function rho(r, parameters), spherical symmetry assumed */
		void* _dens_param;    /**< density parameters (void*) pointer */
		double _rmax;         /**< radius of the system */
		FastParticle* _fp;    /**< fastparticle: particle on which to calculate SCX probabilities */
		int _bpoints; /**< number of grid points in b */
		int _zpoints; /**< number of grid points in z */
		double _bstep; /**< [fm] step size in b of the interpolation grid */
		double _zstep; /**< [fm] step size in z of the interpolation grid */
		double _dens_fctr;    /**< correction for density normalization. Density should be normed to A for A particles...*/
		std::vector<std::vector<double>> _grid;       /**< SCX probability data goes here */
		std::vector<std::vector<double>> _errorGrid;  /**< integrator error estimations go here */
		std::string _filename;  /**< input/output file for the grid (in the grid cache directory) */
};

}  // namespace nucqcd



#endif // CLASS_GRIDTHICK_SCX_HPP
