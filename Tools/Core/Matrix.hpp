/*! \file Matrix.hpp
 * \brief Complex K x L matrix with a fixed size
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <Tools/Core/ErrorMessage.hpp>
#include <complex>
#include <cstdlib>

namespace nucqcd {


//-------------------------------------------------------------------------
//---- Certain classes and functions need to be pre-declared
//---- See C++ FAQ lite 35.16 as to why

/// forward declaration
template<int K,int L> class Matrix;  // pre-declare Matrix<K,L>
// pre-declare friend template functions of Matrix<K,L>
template<int K,int L> Matrix<K,L> operator*(double,const Matrix<K,L>&);
template<int K,int L> Matrix<K,L> operator*(const std::complex<double>&,
					    const Matrix<K,L>&);
template<int K,int L> std::ostream& operator<<(std::ostream&,const Matrix<K,L>&);

//-------------------------------------------------------------------------
//---- class Matrix<K,L> declaration

/*! \brief Complex K x L matrix; indices run from 0, the components are stored row by row */
template< int K, int L >
class Matrix
{
 private:
  // Data member
  // -----------
  std::complex<double>* m; ///< the K*L components, row by row

 public:
  // Constructors
  // ------------
  Matrix(); ///< default constructor: all components zero
  //Matrix(std::complex<double>);               // Ambiguous
  /// constructor for a 2x2, 4x1 or 1x4 matrix, components row by row
  Matrix(std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>);
  /// constructor for a 4x4, 8x2, 2x8, 1x16 or 16x1 matrix, components row by row
  Matrix(std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>,
	 std::complex<double>,std::complex<double>);
  Matrix(const Matrix&); ///< copy constructor

  // Destructor
  // ----------
  ~Matrix(); ///< destructor

  // Assignment
  // ----------
  Matrix& operator=(const Matrix&); ///< assignment

  // Access individual components
  // ----------------------------
  std::complex<double>& operator()(int,int); ///< component (row, column), from 0; throws std::out_of_range outside the matrix
  const std::complex<double> operator()(int,int) const; ///< component (row, column), from 0; throws std::out_of_range outside the matrix

  // Addition
  // --------
  Matrix& operator+=(const Matrix&); ///< addition
  Matrix& operator-=(const Matrix&); ///< subtraction

  Matrix operator+(const Matrix&) const; ///< sum
  Matrix operator-(const Matrix&) const; ///< difference
  
  // Scalar multiplication
  // ---------------------
  Matrix& operator*=(double); ///< multiplication by a number
  Matrix& operator*=(const std::complex<double>&); ///< multiplication by a number

  friend Matrix operator*<>(double,const Matrix&); ///< number times matrix
  friend Matrix operator*<>(const std::complex<double>&,const Matrix&); ///< number times matrix
  
  Matrix operator*(double) const; ///< matrix times number
  Matrix operator*(const std::complex<double>&) const; ///< matrix times number

  // Addition with Scalar multiplication
  // -----------------------------------
  Matrix& multiAdd(double,const Matrix&); ///< adds a number times a matrix: *this += f * m
  Matrix& multiAdd(const std::complex<double>&,const Matrix&); ///< adds a number times a matrix: *this += f * m

  // Matrix multiplication
  // ---------------------
  /// matrix product; throws std::invalid_argument when the sizes do not match
  template< int K1, int L1, int K2, int L2 >
    friend Matrix<K1,L2> operator*(const Matrix<K1,L1>& left,
				   const Matrix<K2,L2>& right);
  /// product of a row and a column: a number
  template< int L1 , int K2 >
    friend std::complex<double> operator*(const Matrix<1,L1>& left,
				     const Matrix<K2,1>& right);
  /// right multiplication by a square matrix
  template < int K1, int L1 >
    friend Matrix<K1,L1>& operator*=(Matrix<K1,L1>& left,
				     const Matrix<L1,L1>& right);

  // Matrix manipulations
  // --------------------
  Matrix<L,K> T() const; ///< transpose
  Matrix<L,K> H() const; ///< Hermitian conjugate
  friend Matrix<2,2> sqrt(const Matrix<2,2>&); ///< square root of a Hermitian 2x2 matrix (by diagonalization)
  /// trace of a square matrix
  template< int K1 >
    friend std::complex<double> Trace(const Matrix<K1,K1>&); // Trace

  // Output
  // ------
  void print(int, int) const; ///< prints component (row, column)
  void print() const; ///< prints all components
  void printReal() const; ///< prints the real parts of all components
  friend std::ostream& operator<< <>(std::ostream&, const Matrix&); ///< prints all components
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//---- class Matrix<K,L> definitions
//---- We do not have a separate source file
//---- See C++ FAQ lite 35.12 for more info

// Constructors
// ------------

// Default constructor
template< int K, int L >
Matrix<K,L>::Matrix()
  : m(new std::complex<double>[K*L])
{
  static std::complex<double> zero(0,0);

  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      m[i*L+j] = zero;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 1x1 matrix
/*
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> number)
{
  if(K==1 && L==1) { // check dimensions of Matrix
      
    m = new std::complex<double>[1*1];
    m[0*1+0] = number; // Assign value
  }
  else {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Error in Matrix::Matrix(std::complex,std::complex,std::complex,std::complex): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 1 argument!" );
  }  
}*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 2x2 matrix
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> v11,std::complex<double> v12,
		    std::complex<double> v21,std::complex<double> v22)
{
  if( K*L == 4 ) { // check dimensions of Matrix
      
    m = new std::complex<double>[2*2];
    // Assign values
    m[0*2+0] = v11;
    m[0*2+1] = v12;
    m[1*2+0] = v21;
    m[1*2+1] = v22;
  }
  else {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Error in Matrix::Matrix(complex,complex,complex,complex): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 4 arguments!" );
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Constructor for 4x4 matrix
template< int K, int L >
Matrix<K,L>::Matrix(std::complex<double> v11,std::complex<double> v12,
		    std::complex<double> v13,std::complex<double> v14,
		    std::complex<double> v21,std::complex<double> v22,
		    std::complex<double> v23,std::complex<double> v24,
		    std::complex<double> v31,std::complex<double> v32,
		    std::complex<double> v33,std::complex<double> v34,
		    std::complex<double> v41,std::complex<double> v42,
		    std::complex<double> v43,std::complex<double> v44)
{
  if( K*L == 16 ) { // check dimensions of Matrix
    
    m = new std::complex<double>[4*4];
    // Assign values
    m[0*4+0] = v11;
    m[0*4+1] = v12;
    m[0*4+2] = v13;
    m[0*4+3] = v14;
    m[1*4+0] = v21;
    m[1*4+1] = v22;
    m[1*4+2] = v23;
    m[1*4+3] = v24;
    m[2*4+0] = v31;
    m[2*4+1] = v32;
    m[2*4+2] = v33;
    m[2*4+3] = v34;
    m[3*4+0] = v41;
    m[3*4+1] = v42;
    m[3*4+2] = v43;
    m[3*4+3] = v44;
  }
  else {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Error in Matrix::Matrix(complex(16times)): a ("
	 << K << "x" << L << ") matrix can not be initialized "
	 << "with 16 arguments!" );
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Copy constructor
template< int K, int L >
Matrix<K,L>::Matrix(const Matrix<K,L>& toCopy)
  : m(new std::complex<double>[K*L])
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      m[i*L+j] = toCopy.m[i*L+j];
}

//-------------------------------------------------------------------------

// Destructor
// ----------

template< int K, int L >
Matrix<K,L>::~Matrix()
{
  delete[] m;
}

//-------------------------------------------------------------------------

// Assignment
// ----------

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator=(const Matrix<K,L>& toCopy)
{
  if(this != &toCopy) // avoid self-assignment
    {
      for(int i=0; i<K; i++)
	for(int j=0; j<L; j++)
	  m[i*L+j] = toCopy.m[i*L+j];
    }

  return *this;
}

//-------------------------------------------------------------------------

// Access individual components
// ----------------------------

template< int K, int L >
std::complex<double>& Matrix<K,L>::operator()(int k, int l)
{
  if( !(k>=0 && k<K && l>=0 && l<L) ) {
    throw std::out_of_range(nucqcd::ErrorMessage() << "Error in Matrix::operator()(int,int): "
	 << "Out of bounds!" );
  }
  return m[k*L+l];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
const std::complex<double> Matrix<K,L>::operator()(int k, int l) const
{
  if( !(k>=0 && k<K && l>=0 && l<L) ) {
    throw std::out_of_range(nucqcd::ErrorMessage() << "Error in Matrix::operator()(int,int): "
	 << "Out of bounds!" );
  }
  return m[k*L+l];
}

//-------------------------------------------------------------------------

// Addition
// --------

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator+=(const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += that.m[i*L+j];
     
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator-=(const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] -= that.m[i*L+j];
     
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator+(const Matrix<K,L>& that) const
{
  Matrix<K,L> sum = *this;
  sum += that;
  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator-(const Matrix<K,L>& that) const
{
  Matrix<K,L> sum = *this;
  sum -= that;
  return sum;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::multiAdd(double factor, const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += factor * that.m[i*L+j];
  
  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template< int K, int L >
Matrix<K,L>& Matrix<K,L>::multiAdd(const std::complex<double>& factor, 
				   const Matrix<K,L>& that)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] += factor * that.m[i*L+j];
  
  return *this;
}

//-------------------------------------------------------------------------

// Scalar multiplication 
// ---------------------

// (matrix*double)
template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator*=(double factor)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*complex)
template< int K, int L >
Matrix<K,L>& Matrix<K,L>::operator*=(const std::complex<double>& factor)
{
  for(int i=0; i<K; i++)
    for(int j=0; j<L; j++)
      this->m[i*L+j] *= factor;

  return *this;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*double)
template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator*(double factor) const
{
  Matrix<K,L> product = *this;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (matrix*complex)
template< int K, int L >
Matrix<K,L> Matrix<K,L>::operator*(const std::complex<double>& factor) const
{
  Matrix<K,L> product = *this;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (double*matrix)
/// number times matrix
template< int K, int L >
Matrix<K,L> operator*(double factor, const Matrix<K,L>& matrix)
{
  Matrix<K,L> product = matrix;
  product *= factor;
  
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// (complex*matrix)
/// number times matrix
template< int K, int L >
Matrix<K,L> operator*(const std::complex<double>& factor,
		      const Matrix<K,L>& matrix)
{
  Matrix<K,L> product = matrix;
  product *= factor;
  
  return product;
}

//-------------------------------------------------------------------------

// Matrix multiplication
// ---------------------

/// matrix product; throws std::invalid_argument when the sizes do not match
template< int K1, int L1, int K2, int L2 >
Matrix<K1,L2> operator*(const Matrix<K1,L1>& left, const Matrix<K2,L2>& right)
{
  if(L1 != K2) { // Check whether multiplication is meaningful
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Error in operator*(Matrix,Matrix): "
	 << "Matrix sizes are not compatible for multiplication");
  }
  
  Matrix<K1,L2> product;
//   double realpart[K1][L2];
//   double imagpart[K1][L2];
  for(int i=0; i<K1; i++)
    for(int j=0; j<L2; j++){
//       realpart[i][j]=0.; imagpart[i][j]=0.;
      for(int x=0; x<K2; x++){
 	product.m[i*L2+j] += (left.m[i*L1+x]*right.m[x*L2+j]);
// 	 realpart[i][j] += real(left.m[i*L1+x])*real(right.m[x*L2+j])-imag(left.m[i*L1+x])*imag(right.m[x*L2+j]);
// 	 imagpart[i][j]+=    real(left.m[i*L1+x])*imag(right.m[x*L2+j])+imag(left.m[i*L1+x])*real(right.m[x*L2+j]);
      }
//       product.m[i*L2+j]=std::complex<double>(realpart[i][j],imagpart[i][j]);
    }
  return product;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Matrix multiplication resulting in a scalar
/// product of a row and a column: a number
template< int L1 , int K2 >
std::complex<double> operator*(const Matrix<1,L1>& left, const Matrix<K2,1>& right)
{
  if(L1 != K2) { // Check whether multiplication is meaningful
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Error in operator*(Matrix,Matrix): "
	 << "Matrix sizes are not compatible for multiplication");
  }

  std::complex<double> scalar(0,0);
  for(int i=0; i<L1; i++)
    scalar += left.m[0*L1+i]*right.m[i*1+0];
  return scalar;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// right multiplication by a square matrix
template < int K1, int L1 >
Matrix<K1,L1>& operator*=(Matrix<K1,L1>& left,const Matrix<L1,L1>& right)
{
  left = left * right;
  
  return left;
}

//-------------------------------------------------------------------------

// Transpose
// ---------
template< int K, int L >
Matrix<L,K> Matrix<K,L>::T() const
{
  Matrix<L,K> transposed;
  for(int i=0; i<L; i++)
    for(int j=0; j<K; j++)
      transposed(i,j) = m[j*L+i];

  return transposed;
}

//-------------------------------------------------------------------------

// Hermitian conjugate
// -------------------
template< int K, int L >
Matrix<L,K> Matrix<K,L>::H() const
{
  Matrix<L,K> conjugate;
  for(int i=0; i<L; i++)
    for(int j=0; j<K; j++)
      conjugate(i,j) = conj(m[j*L+i]);

  return conjugate;
}

//-------------------------------------------------------------------------

// Trace of square matrices
// ------------------------
/// trace of a square matrix
template< int K>
std::complex<double> Trace(const Matrix<K,K>& m)
{
  std::complex<double> trace = 0.;
  
  for(int i=0; i<K; ++i) trace += m.m[i*K+i];

  return trace;
}

//-------------------------------------------------------------------------

// Print the (x,y) component
template< int K, int L >
void Matrix<K,L>::print(int x, int y) const
{
  std::cout << (*this)(x,y) << std::endl;
}

//-------------------------------------------------------------------------

// Print all components
template< int K, int L >
void Matrix<K,L>::print() const
{
  std::cout << "{";
  for(int i=0; i<K; i++)
    {
      std::cout << "{";
      for(int j=0; j<L; j++)
	{
	  std::cout << m[i*L+j]
	       << ( j<L-1 ? "," : "" );
	}
      std::cout << "}";
    }
  std::cout << "}" << std::endl;
}

//-------------------------------------------------------------------------

// Print the real part of all components
template< int K, int L >
void Matrix<K,L>::printReal() const
{
  std::cout << "{";
  for(int i=0; i<K; i++)
    {
      std::cout << "{";
      for(int j=0; j<L; j++)
	{
	  std::cout << real(m[i*L+j])
	       << ( j<L-1 ? "," : "" );
	}
      std::cout << "}";
    }
  std::cout << "}" << std::endl;
}

//-------------------------------------------------------------------------

// Enables output using << operator
/// prints all components
template< int K, int L >
std::ostream& operator<<(std::ostream& output, const Matrix<K,L>& matrix)
{
  output << "{";
  for(int i=0; i<K; i++)
    {
      output << "{";
      for(int j=0; j<L; j++)
	{
	  output << matrix.m[i*L+j]
	       << ( j<L-1 ? "," : "" );
	}
      output << "}";
    }
  output << "}";

  return output;
}




}  // namespace nucqcd

/** @} */
#endif
