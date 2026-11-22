/*! \file GammaStructure.h
 * \brief General complex 4x4 matrix as a combination of the 16 Dirac basis matrices (Dirac representation)
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef GAMMASTRUCTURE_H
#define GAMMASTRUCTURE_H

#include <iostream>
#include <complex>
#include "Matrix.h"

// ********************************************************** //
// First we define the gamma matrices (Dirac representation) //
// ********************************************************** //

/// unit 4x4 matrix (basis code 0)
const Matrix<4,4> U(1.0, 0.0, 0.0, 0.0,
		    0.0, 1.0, 0.0, 0.0,
		    0.0, 0.0, 1.0, 0.0,
		    0.0, 0.0, 0.0, 1.0);

/// gamma^0 (basis code 1), Dirac representation
const Matrix<4,4> G0(1.0, 0.0,  0.0,  0.0,
		     0.0, 1.0,  0.0,  0.0,
		     0.0, 0.0, -1.0,  0.0,
		     0.0, 0.0,  0.0, -1.0);

/// gamma^1 (code 2)
const Matrix<4,4> G1( 0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0,
		      0.0,-1.0, 0.0, 0.0,
		     -1.0, 0.0, 0.0, 0.0);

/// gamma^2 (code 3)
const Matrix<4,4> G2(0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0,
		     0.0,std::complex<double>(0,1.0),0.0,0.0,
		     std::complex<double>(0,-1.0),0.0,0.0,0.0);

/// gamma^3 (code 4)
const Matrix<4,4> G3( 0.0, 0.0, 1.0, 0.0,
		      0.0, 0.0, 0.0,-1.0,
		     -1.0, 0.0, 0.0, 0.0,
		      0.0, 1.0, 0.0, 0.0);

/// gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3 (code 5)
const Matrix<4,4> G5(0.0, 0.0, 1.0, 0.0,
		     0.0, 0.0, 0.0, 1.0,
		     1.0, 0.0, 0.0, 0.0,
		     0.0, 1.0, 0.0, 0.0); // i * G0*G1*G2*G3

/// gamma^5 gamma^0 (code 6)
const Matrix<4,4> G5G0(0.0, 0.0,-1.0, 0.0,
		       0.0, 0.0, 0.0,-1.0,
		       1.0, 0.0, 0.0, 0.0,
		       0.0, 1.0, 0.0, 0.0); // G5*G0

/// gamma^5 gamma^1 (code 7)
const Matrix<4,4> G5G1( 0.0,-1.0, 0.0, 0.0,
		       -1.0, 0.0, 0.0, 0.0,
		        0.0, 0.0, 0.0, 1.0,
		        0.0, 0.0, 1.0, 0.0); // G5*G1

/// gamma^5 gamma^2 (code 8)
const Matrix<4,4> G5G2(0.0,std::complex<double>(0,1.0),0.0,0.0,
		       std::complex<double>(0,-1.0),0.0,0.0,0.0,
		       0.0,0.0,0.0,std::complex<double>(0,-1.0),
		       0.0,0.0,std::complex<double>(0,1.0),0.0); // G5*G2

/// gamma^5 gamma^3 (code 9)
const Matrix<4,4> G5G3(-1.0, 0.0, 0.0, 0.0,
		        0.0, 1.0, 0.0, 0.0,
		        0.0, 0.0, 1.0, 0.0,
		        0.0, 0.0, 0.0,-1.0); // G5*G3

/// gamma^0 gamma^1 (code 10)
const Matrix<4,4> A1( 0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0,
		      0.0, 1.0, 0.0, 0.0,
		      1.0, 0.0, 0.0, 0.0); // G0*G1

/// gamma^0 gamma^2 (code 11)
const Matrix<4,4> A2(0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0,
		     0.0,std::complex<double>(0,-1.0),0.0,0.0,
		     std::complex<double>(0,1.0),0.0,0.0,0.0); // G0*G2

/// gamma^0 gamma^3 (code 12)
const Matrix<4,4> A3( 0.0, 0.0, 1.0, 0.0,
		      0.0, 0.0, 0.0,-1.0,
		      1.0, 0.0, 0.0, 0.0,
		      0.0,-1.0, 0.0, 0.0); // G0*G3

/// gamma^5 gamma^0 gamma^1 (code 13)
const Matrix<4,4> S1( 0.0, 1.0, 0.0, 0.0,
		      1.0, 0.0, 0.0, 0.0,
		      0.0, 0.0, 0.0, 1.0,
		      0.0, 0.0, 1.0, 0.0); // G5*G0*G1

/// gamma^5 gamma^0 gamma^2 (code 14)
const Matrix<4,4> S2(0.0,std::complex<double>(0,-1.0),0.0,0.0,
		     std::complex<double>(0,1.0),0.0,0.0,0.0,
		     0.0,0.0,0.0,std::complex<double>(0,-1.0),
		     0.0,0.0,std::complex<double>(0,1.0),0.0); // G5*G0*G2

/// gamma^5 gamma^0 gamma^3 (code 15)
const Matrix<4,4> S3(1.0, 0.0, 0.0, 0.0,
		     0.0,-1.0, 0.0, 0.0,
		     0.0, 0.0, 1.0, 0.0,
		     0.0, 0.0, 0.0,-1.0); // G5*G0*G3


// *********************************************** //
// Our basis of lineary independent gamma matrices //
// *********************************************** //
/// the 16 basis matrices, indexed by their code
const Matrix<4,4> gammaBasis[16] =
  { U, G0, G1, G2, G3, G5, G5G0, G5G1, G5G2, G5G3,
    A1, A2, A3, S1, S2, S3 };


// ******************************** //
// This is the GammaStructure class //
// ******************************** //

/*! \brief Complex 4x4 matrix as a combination \f$ \sum_c a_c \Gamma_c \f$ of the 16 basis matrices
 *
 * Basis codes: 0 unity, 1-4 gamma^0..gamma^3, 5 gamma^5, 6-9 gamma^5 gamma^mu, 10-12 gamma^0 gamma^i,
 * 13-15 gamma^5 gamma^0 gamma^i, all in the Dirac representation (the Matrix<4,4> constants above).
 * Only components with a nonzero coefficient are stored.
 *
 * \note The constructors take the coefficients in the order (unity, gamma^5, gamma^0, gamma^1,
 * gamma^2, gamma^3, gamma^5 gamma^0, gamma^5 gamma^1, gamma^5 gamma^2, gamma^5 gamma^3, A_1, A_2,
 * A_3, S_1, S_2, S_3): gamma^5 comes second, which differs from the basis codes.
 */
class GammaStructure
{  
 public:
  // Constructors
  // ------------
  // Syntax of GammaStructure constructor is
  // GammaStructure(unity,gamma^5,gamma^0,gamma^1,
  //                gamma^2,gamma^3,gamma^5*gamma^0,gamma^5*gamma^1,
  //                gamma^5*gamma^2,gamma^5*gamma^3,
  //                A_1,A_2,A_3,S_1,S_2,S_3)
  GammaStructure(const std::complex<double>& =0); ///< constructor: a times unity
  GammaStructure(const std::complex<double>&, const std::complex<double>&); ///< constructor: coefficients of unity and gamma^5
  /// constructor: coefficients of unity, gamma^5, gamma^0, gamma^1, gamma^2, gamma^3
  GammaStructure(const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0);
  /// constructor with all 16 coefficients, in the order of the class description
  GammaStructure(const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>&,
		 const std::complex<double>&, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0,
		 const std::complex<double>& =0, const std::complex<double>& =0);
  GammaStructure(const GammaStructure&); ///< copy constructor

  // Destructor
  // ----------
  ~GammaStructure(); ///< destructor

  // Assignment
  // ----------
  GammaStructure& operator=(const GammaStructure&); ///< assignment
  
  // Evaluate the structure and return the 4x4 matrix
  // ------------------------------------------------
  Matrix<4,4> value() const; ///< the 4x4 matrix
  operator Matrix<4,4>() const { return value(); } ///< the 4x4 matrix

  // Addition
  // --------
  GammaStructure operator+(const GammaStructure&) const; ///< sum
  GammaStructure operator-(const GammaStructure&) const; ///< difference
  
  GammaStructure& operator+=(const GammaStructure&); ///< addition
  GammaStructure& operator-=(const GammaStructure&); ///< subtraction

  /* // Addition with scalar multiplication */
  /* // ----------------------------------- */
  /* GammaStructure& multiAdd(double,const GammaStructure&, */
  /* 			   double,const GammaStructure&, */
  /* 			   double,const GammaStructure&); */
  /* GammaStructure multiAdd(double,const GammaStructure&, */
  /* 			  double,const GammaStructure&, */
  /* 			  double,const GammaStructure&) const; */
  /* static */
  /*   GammaStructure multiAdd(double,const GammaStructure&, */
  /* 			    double,const GammaStructure&, */
  /* 			    double,const GammaStructure&, */
  /* 			    double,const GammaStructure&); */

  // Multiplication
  // --------------
  GammaStructure operator*(double) const; ///< multiplication by a number
  GammaStructure operator*(const std::complex<double>&) const; ///< multiplication by a number
  GammaStructure operator*(const GammaStructure&) const; ///< product, with the multiplication table of the basis
  /// product with a 4xK matrix
  template<int K>
  Matrix<4,K>    operator*(const Matrix<4,K>& m) const
    { return (Matrix<4,4>)*this * m; }

  /// number times gamma structure
  friend
    GammaStructure operator*(double, const GammaStructure&);
  /// number times gamma structure
  friend
    GammaStructure operator*(const std::complex<double>&, const GammaStructure&);
  /// product of a Kx4 matrix with a gamma structure
  template<int K> friend
    Matrix<K,4>    operator*(const Matrix<K,4>& m, const GammaStructure& g)
    { return m * (Matrix<4,4>)g; }

  GammaStructure& operator*=(double); ///< multiplication by a number
  GammaStructure& operator*=(const std::complex<double>&); ///< multiplication by a number
  GammaStructure& operator*=(const GammaStructure&); ///< right multiplication

  // Matrix operations
  // -----------------
  GammaStructure Hermitian() const; ///< Hermitian conjugate
  GammaStructure Conjugate() const; ///< Dirac conjugate \f$ \gamma^0 \Gamma^\dagger \gamma^0 \f$
  
  // Output
  // ------
  void print(void) const; ///< prints the nonzero components
  /// prints the nonzero components
  friend 
    std::ostream& operator<<(std::ostream&, const GammaStructure&); // Overwrite <<

 private:
  // Data members
  // ------------
  int nrComp; ///< number of nonzero components
  int* comp; ///< basis codes of the nonzero components
  std::complex<double>* compCoeff; ///< their coefficients

  // multiplication tables
  static const int multTable[16][16]; ///< basis code of the product of two basis matrices
  static const std::complex<double> multTableCoeff[16][16]; ///< coefficient of that product

  // Private member funtions
  // -----------------------
  GammaStructure(int,int*,std::complex<double>*); ///< constructor from the members (takes ownership of the arrays), used by the products
  void crop(); ///< removes the components with a zero coefficient
};


/* The 4vector gamma_{mu} (GMU) containing 4 GammaStructures is defined in
 * ../FourVector/FourVector.h for nesting reasons */

/** @} */
#endif
