#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__
/*
 *  *    Function calculates Legendre Polynomials Pn(x)
 *   */
namespace Legendre
{
	// n = 0
	template <class T> inline double P0(const T& x)
	{
		return static_cast<T>(1);
	}

	// n = 1
	template <class T> inline double P1(const T& x)
	{
		return x;
	}

	// n = 2
	template <class T> inline double P2(const T& x)
	{
		return ((static_cast<T>(3) * x*x) - static_cast<T>(1)) /
			static_cast<T>(2);
	}

	/*
	 *      Pn(x)
	 */
	template <class T> inline double Pn(unsigned int n, const T& x)
	{
		switch(n)
		{
			case 0:
				return P0<T>(x);

			case 1:
				return P1<T>(x);

			case 2:
				return P2<T>(x);

			default:
				break;
		}

		/*  We could simply do this:
			 return (static_cast<T>(((2 * n) - 1)) * x * Pn(n - 1, x) -
			 (static_cast<T>(n - 1)) * Pn(n - 2, x)) / static_cast<T>(n);
			 but it could be slow for large n */

		double pnm1(P2<T>(x));
		double pnm2(P1<T>(x));
		T pn;

		for (double m = 3u ; m <= n ; ++m)
		{ 
			pn = ((static_cast<T>((2 * m) - 1)) * x * pnm1
					- (static_cast<T>(m - 1) * pnm2)) / static_cast<T>(m);
			pnm2 = pnm1;
			pnm1 = pn;
		}

		return pn;
	}
}
#endif
