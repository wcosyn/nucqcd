#include <Tools/Core/recmosh.hpp>

namespace nucqcd {

RecMosh::RecMosh(int n1, int l1, int n2, int l2, string path)
	: n1(n1),
	  l1(l1),
	  n2(n2),
	  l2(l2),
	  path(path)
{
}

RecMosh::~RecMosh()
{
}

double RecMosh::getCoefficient( int n, int l, int N, int Lambda, int L )
{
	if( 2*n+l+2*N+Lambda != 2*n1+l1+2*n2+l2)
	{
		return 0;
	}
	// memoized in calculate(); the old int key 100000*n + 10000*l + ... collided once l >= 10
	return calculate(n,l,N,Lambda,n1,l1,n2,l2,L);
}

double RecMosh::calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L)
{
	const std::array<int,9> key{ n, l, N, Lambda, n1, l1, n2, l2, L };
	auto it = memo.find(key);
	if( it != memo.end() ) return it->second;
	double result = calculateUncached(n, l, N, Lambda, n1, l1, n2, l2, L);
	memo[key] = result;
	return result;
}

double RecMosh::calculateUncached( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L)
{
	if( L < fabs(l-Lambda) || L > l+Lambda)
		return 0;
	if( L < fabs(l1-l2) || L > l1+l2 )
		return 0;
	if( 2*n+l + 2*N+Lambda != 2*n1+l1+2*n2+l2) return 0;
	if( n1 > 0 )
	{
// 		cout << "n1>0" << endl;
		double factor = 1./sqrt((n1)*(n1+l1+0.5));
		double sum = 0.;
		for( int na = n-1; na <= n ; na++ )
		{
			if( na<0) continue;
			for( int la = l-1; la <= l+1; la++)
			{
				if( la<0) continue;
				for( int Na = N-1; Na <=N; Na++ )
				{
					if( Na<0) continue;
					for( int Lambdaa = Lambda-1; Lambdaa <= Lambda+1; Lambdaa++)
					{
						if( Lambdaa<0) continue;
						double me = getMatrixElement(n,l,N,Lambda,na,la,Na,Lambdaa, L, 1);
						double moshbr = calculate(na,la,Na,Lambdaa,n1-1,l1,n2,l2,L);
						sum += me*moshbr;
					}
				}
			}
		}
		return factor*sum;
	}
	else if( n2 > 0 )
	{
// 		cout << "n2>0" << endl;
		double factor = 1./sqrt((n2)*(n2+l2+0.5));
		double sum = 0.;
		for( int na = n-1; na <= n; na++ )
		{
			if( na<0) continue;
			for( int la = l-1; la <= l+1; la++)
			{
				if( la<0) continue;
				for( int Na = N-1; Na <=N; Na++ )
				{
					if( Na<0) continue;
					for( int Lambdaa = Lambda-1; Lambdaa <= Lambda+1; Lambdaa++)
					{
						if( Lambdaa<0) continue;
						double me = getMatrixElement(n,l,N,Lambda,na,la,Na,Lambdaa, L, 2);
						double moshbr = calculate(na,la,Na,Lambdaa,n1,l1,n2-1,l2,L);
						sum += me*moshbr;
					}
				}
			}
		}
		return factor*sum;
	}
	else if( n1 == 0 && n2 == 0)
	{
		double factor = gsl_sf_fact(l1) * gsl_sf_fact(l2)/ gsl_sf_fact(2*l1) / gsl_sf_fact(2*l2);
		factor *= (2*l+1)*(2*Lambda+1)/pow(2.,l+Lambda);
		factor *= gsl_sf_fact(n+l)/gsl_sf_fact(n) / gsl_sf_fact(2*n+2*l+1);
		factor *= gsl_sf_fact(N+Lambda)/gsl_sf_fact(N)/ gsl_sf_fact(2*N+2*Lambda+1);
		double sign = pow( -1., n+l+Lambda-L);
		double sum = 0;
// 	  	for( int x = 0; x <= 4; x++)
		for( int x = fabs(l-l1); x <= l+l1; x++)
		{
// 			cout << "x allowed ?" << x << endl;
			if( x < fabs(Lambda - l2) || x > Lambda+l2 )
				continue;
			double W = Wc( l, Lambda, l1, l2, L, x);
// 			cout << x << " " << W << endl;
			if( W == 0 ) 
			{
	// 			cout << "NOOOOOOOOOOO" << endl;
				continue;
			}
// 			cout << "x allowed: " << x << endl;
// 			cout << "-------" << endl;
			double term = 2*x+1;
			term *= A( l1, l, l2, Lambda, x );
	// 		cout << " A" << term << endl;
			term *= W;
			sum += term;
	// 		cout << " ------E " << endl;
		}
// 		cout << "<(" << n << l << N << Lambda << ")" << L;
// 		cout << "|(0" << l1 << 0 << l2 << ")" << L << "> = " << sqrt(factor)*sign*sum << endl;
		return sqrt(factor)*sign*sum;
	}
	
	cerr << " IMPOSSIBLE " << endl;
	return 0.;
	
}

double RecMosh::getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f)
{
	if( na == n-1 )
	{
		if( la == l)
		{
			if( Na == N)
			{
				if( Lambdaa == Lambda) return 0.5*sqrt(n*(n+l+0.5));
				else return 0;
			}
			else return 0;
		}
		else if( la == l+1)
		{
			if( Na == N-1 )
			{
				if ( Lambdaa == Lambda+1)
					return pow(-1.,f+1)*sqrt( n*N*(l+1)*(Lambda+1))*pow(-1.,L+Lambda+l)*Wc(l,l+1,Lambda,Lambda+1,1,L);
				else return 0;
			}
			else if (Na == N)
			{
				if( Lambdaa == Lambda-1)
					return pow(-1.,f+1)*sqrt(n*(N+Lambda+0.5)*(l+1)*Lambda)*pow(-1.,L+Lambda+l)*Wc(l,l+1,Lambda,Lambda-1,1,L);
				else return 0;
			}
			else return 0;
		}
		else return 0;
	}
	else if ( na == n)
	{
		if (la == l)
		{
			if( Na == N-1)
			{
				if( Lambdaa == Lambda) return 0.5*sqrt(N*(N+Lambda+0.5));
				else return 0;
			}
			else return 0;
		}
		else if ( la == l-1)
		{
			if( Na == N-1)
			{
				if( Lambdaa == Lambda+1)
					return pow(-1.,f+1)*sqrt((n+l+0.5)*N*l*(Lambda+1))*pow(-1.,L+Lambda+l)*Wc(l,l-1,Lambda,Lambda+1,1,L);
				else return 0;
			}
			else if ( Na == N)
			{
				if( Lambdaa == Lambda-1)
					return pow(-1.,f+1)*sqrt((n+l+0.5)*(N+Lambda+0.5)*l*Lambda)*pow(-1.,L+Lambda+l)*Wc(l,l-1,Lambda,Lambda-1,1,L);
				else return 0;
			}
			else return 0;
		}
		else return 0;
	}
	else return 0;
		
}
				

double RecMosh::Wc( int a, int b, int c, int d, int e, int f)
{
	double sign = pow(-1.,a+b+c+d);
	double sixj = gsl_sf_coupling_6j( 2*a, 2*b, 2*e, 2*d, 2*c, 2*f);
	return sign*sixj;
}

double RecMosh::A( int l1, int l, int l2, int Lambda, int x )
{
	double factor = gsl_sf_fact(l1+l+x+1) ;
	factor*= gsl_sf_fact(l1+l-x);
	factor*= gsl_sf_fact(l1+x-l);
	factor *= gsl_sf_fact(l2+Lambda+x+1)*gsl_sf_fact(l2+Lambda-x)*gsl_sf_fact(l2+x-Lambda);
	factor *= 1./gsl_sf_fact(l+x-l1)/gsl_sf_fact(Lambda+x-l2);
	double sum = 0.;
	for( int q = x; q <= l1+l && q<= l2+Lambda; q++ )
	{
// 		cout << "q " << q << endl;
		if( ((l+q-l1)%2) || ((Lambda+q-l2)%2) )
			continue;
// 		cout << ((l+q-l1)%2) << " " <<  ((Lambda+q-l2)%2)  <<  q << endl;
		double sign = pow(-1., 0.5*(l+q-l1) );
		double term = gsl_sf_fact(l+q-l1)/gsl_sf_fact((l+q-l1)/2)/gsl_sf_fact((l+l1-q)/2);
		term *= 1./gsl_sf_fact(q-x)/gsl_sf_fact(q+x+1);
		term*= gsl_sf_fact(Lambda+q-l2)/gsl_sf_fact((Lambda+q-l2)/2)/gsl_sf_fact((Lambda+l2-q)/2);
		term *= sign;
// 		cout << " t " << term << endl;
		sum += term;
	}
	return sqrt(factor)*sum;
	
	
}

}  // namespace nucqcd


