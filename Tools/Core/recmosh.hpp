#ifndef RECMOSH_H
#define RECMOSH_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <cmath>
using std::fabs;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <map>
using std::map;
#include <array>

#include <string>

namespace nucqcd {

using namespace std;

/**
 * Moshinsky brackets <(n l, N Lambda) L | (n1 l1, n2 l2) L> for equal-mass particles, computed
 * with the recursion relations in n1 and n2 down to the closed form for n1 = n2 = 0.
 * Intermediate brackets are memoized, so the full set for a pair (n1 l1, n2 l2) costs well under a
 * second even for the 5th oscillator shell. (Earlier versions read a cache from share/recmosh/.)
 */
class RecMosh
{
private: 	
	int n1; ///< n1 of the first particle
	int l1; ///< l1 of the first particle
	int n2; ///< n2 of the second particle
	int l2; ///< l2 of the second particle
	map< std::array<int,9>, double > memo; ///< memoized calculate() results
	string path; ///< unused (the former file cache)
	double Wc( int a, int b, int c, int d, int e, int f); ///< Racah W coefficient
	double A( int l1, int l, int l2, int Lambda, int x ); ///< coefficient A of the closed form for n1 = n2 = 0
	double calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L); ///< bracket <(n l, N Lambda) L | (n1 l1, n2 l2) L>, memoized
	double calculateUncached( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L); ///< bracket by the recursion relations (used by calculate)
	double getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f); ///< matrix element of the recursion relation
public:
	/** path is unused; it is kept for source compatibility with the former file cache */
	RecMosh(int n1, int l1, int n2, int l2, string path);
	/// destructor
	~RecMosh();
	double getCoefficient( int n, int l, int N, int Lambda, int L ); ///< Moshinsky bracket <(n l, N Lambda) L | (n1 l1, n2 l2) L> for the pair of this object
// 	RecMosh(int n1, int l1, int n2, int l2, int n, int l, int N, int Lambda, int L, bool reversed);
};

}  // namespace nucqcd


#endif // RECMOSH_H
