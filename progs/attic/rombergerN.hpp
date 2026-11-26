/*! \file rombergerN.hpp
 * \brief The nested 1D Romberg integration routine `rombergerN`, kept for the archived programs
 *
 * Until 2026-11 this lived in Tools/Core/Utilfunctions.hpp and was the integrator behind the `integrator
 * == 0` option of the Glauber and Knockout classes (removed in #69, the option had silently stopped working
 * in the Glauber grids) and behind a few hard-coded calls. The last users are the archived programs
 * (progs/attic/rhodeuteron: RhoDeuteron and rhodint; progs/attic/threeres: InclusiveCrossRes, not built). The
 * integrators to use are numint::cube_romb and numint::cube_adaptive (see docs/integrator_benchmark.md).
 */
#ifndef ROMBERGERN_ATTIC_HPP
#define ROMBERGERN_ATTIC_HPP

#include <cmath>
#include <cstdarg>
#include <complex>
#include <iostream>

#include <Tools/Core/Utilfunctions.hpp>  // power()

namespace nucqcd {

/*! template function for romberg integration of N functions (slightly different preferably) T f(x)
 * \tparam T template return type of the function you're integrating
 * \param function function you want to integrate, must return N values (through array)
 * \param a lower limit
 * \param b upper limit
 * \param N number of functions you want to integrate
 * \param results array with integration results
 * \param acc relative accuracy that has to be reached
 * \param min min number of steps in the Romberg algorithm
 * \param max max number of steps in the Romberg algorithm
 * \param estimate pointer to a variable containing an upper limit on nested integration, to speed up things
 */
template <class T> void rombergerN(void (*function)(double, T*,  va_list), double a, double b, 
		       int N, T* results, double acc, int min, int max, double *estimate, ...)
{

  if(std::abs(a-b)<1e-06){for(int i=0;i<N;i++) results[i]=0.;return;}
  T **DN1= new T*[N];
  T **DN2= new T*[N];
  double x,h=b-a;
  T sum[N], value[N];
  double dev[N];
  
  va_list ap;
  va_start(ap,estimate);
  
  for(int i=0;i<N;i++) {
    DN1[i] = new T[1];
  }

  function(a,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++) DN1[i][0] = value[i]; 
  va_start(ap,estimate);
  function(b,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++){
    DN1[i][0] += value[i]; 
  DN1[i][0] *= 0.5*(b-a);
  }
  

  
  for (int n = 1; n <= max ; n++) {
    for(int i=0;i<N;i++) {
      sum[i] = 0.;
      //store new results
      DN2[i] = new T[n+1];
    }
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,estimate);
      function(x,value,ap);
      va_end(ap);
      for(int i=0;i<N;i++) sum[i] += value[i];
    }
    for(int i=0;i<N;i++) DN2[i][0] = 0.5*DN1[i][0] + power(0.5,n)*h*sum[i];
    int p=4;
    for (int m = 1; m <= n; m++) {
      for(int i=0;i<N;i++) DN2[i][m] = (double(p)*DN2[i][m-1] - DN1[i][m-1])/(p - 1.);
      p*=4;
    }
    if (n >= min) {
      double deviation=0.;
      for(int i=0;i<N;i++) {
	dev[i] = (DN2[i][n] == std::complex<double>(0.,0.)) ? 0. : std::abs((DN2[i][n]-DN1[i][n-1]))/std::abs(DN2[i][n]);
	if(dev[i]>deviation) deviation=dev[i];
      }
     if (((deviation < acc ) ) || ((std::abs(DN2[0][n]-DN1[0][n-1]) <acc*1e-07 ))) {
       for(int i=0;i<N;i++) {
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];
	 if(std::abs(results[i])>(*estimate)) (*estimate)=std::abs(results[i]);
       }
       delete [] DN1; delete [] DN2;
       return;
     }
     double dummy = 0.;
     for(int i=0;i<N;i++) dummy +=std::abs(DN2[i][n]);
     if(dummy <(*estimate)*1e-05){
       for(int i=0;i<N;i++){
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];	
       }
       delete [] DN1; delete [] DN2;
       return;
     }
    }
    for(int i=0;i<N;i++){
      delete [] DN1[i];
      DN1[i] = DN2[i];
    }
  }
    
  for(int i=0;i<N;i++){
    results[i] = DN2[i][max];
    delete [] DN2[i];
  }
  delete [] DN1; delete [] DN2;
  va_end(ap);
  return;  
  
}


/*! template function for romberg integration of N member functions (slightly different preferably) T F->f(x)
 * \tparam T template return type of the function you're integrating
 * \tparam F template class of which the integration functions are member function
 * \param object pointer to class object that contains the memberfunction
 * \param function member function you want to integrate, must return N values (through array)
 * \param a lower limit
 * \param b upper limit
 * \param N number of functions you want to integrate
 * \param results array with integration results
 * \param acc relative accuracy that has to be reached
 * \param min min number of steps in the Romberg algorithm
 * \param max max number of steps in the Romberg algorithm
 * \param estimate pointer to a variable containing an upper limit on nested integration, to speed up things
 */
template <class T, class F> void rombergerN(F* object, void (F::*function)(double, T*,  va_list), double a, double b, 
		       int N, T* results, double acc, int min, int max, double *estimate, ...)
{

  if(std::abs(a-b)<1e-06){for(int i=0;i<N;i++) results[i]=0.;return;}
  T **DN1= new T*[N];
  T **DN2= new T*[N];
  double x,h=b-a;
  T sum[N], value[N];
  double dev[N];
  
  va_list ap;
  va_start(ap,estimate);
  
  for(int i=0;i<N;i++) {
    DN1[i] = new T[1];
  }

  (object->*function)(a,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++) DN1[i][0] = value[i]; 
  va_start(ap,estimate);
  (object->*function)(b,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++){
    DN1[i][0] += value[i]; 
  DN1[i][0] *= 0.5*(b-a);
  }
  

  
  for (int n = 1; n <= max ; n++) {
    for(int i=0;i<N;i++) {
      sum[i] = 0.;
      //store new results
      DN2[i] = new T[n+1];
    }
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,estimate);
      (object->*function)(x,value,ap);
      va_end(ap);
      for(int i=0;i<N;i++) sum[i] += value[i];
    }
    for(int i=0;i<N;i++) DN2[i][0] = 0.5*DN1[i][0] + power(0.5,n)*h*sum[i];
    int p=4;
    for (int m = 1; m <= n; m++) {
      for(int i=0;i<N;i++) DN2[i][m] = (double(p)*DN2[i][m-1] - DN1[i][m-1])/(p - 1.);
      p*=4;
    }
    if (n >= min) {
      double deviation=0.;
      for(int i=0;i<N;i++) {
	dev[i] = (DN2[i][n] == std::complex<double>(0.,0.)) ? 0. : std::abs((DN2[i][n]-DN1[i][n-1]))/std::abs(DN2[i][n]);
	if(dev[i]>deviation) deviation=dev[i];
      }
     if (((deviation < acc ) ) || ((std::abs(DN2[0][n]-DN1[0][n-1]) <acc*1e-07 ))) {
       for(int i=0;i<N;i++) {
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];
	 if(std::abs(results[i])>(*estimate)) (*estimate)=std::abs(results[i]);
       }
       delete [] DN1; delete [] DN2;
       return;
     }
     double dummy = 0.;
     for(int i=0;i<N;i++) if(dummy<std::abs(DN2[i][n])) dummy=std::abs(DN2[i][n]);
     if(dummy <(*estimate)*1e-05){
       for(int i=0;i<N;i++){
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];	
       }
       delete [] DN1; delete [] DN2;
       return;
     }
    }
    for(int i=0;i<N;i++){
      delete [] DN1[i];
      DN1[i] = DN2[i];
    }
  }
    
  for(int i=0;i<N;i++){
    results[i] = DN2[i][max];
    delete [] DN2[i];
  }
  delete [] DN1; delete [] DN2;
  va_end(ap);
  return;  
  
}

}  // namespace nucqcd

#endif
