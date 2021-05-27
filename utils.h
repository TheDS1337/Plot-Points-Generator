#pragma once

#define M_PI           3.1415926535897932384

constexpr long double factorial(int n)
{
	if( n > 65 )
	{
		return exp(n * log(n) - n + (log(n) / 2) - log(2 * M_PI));
	}

	return n <= 0 ? 1 : n * factorial(n - 1);
}

constexpr long double permutation(int n, int p, bool inversed = false)
{
	return inversed ? factorial(p) / factorial(n) : factorial(n) / factorial(p);
}

constexpr long double combination(int n, int p, bool inversed = false)
{
	return inversed ? factorial(p) / permutation(n, p) : permutation(n, p) / factorial(p);
}

constexpr long double deformed_factorial(int n, long double (*pF) (int, long double), long double q)
{
	long double f = pF(n, q), p = 1.0;

	if( f != 0.0 )
	{
		p = pow(pF(n, q), 2);
	}

	return n <= 0 ? 1 : n * p * deformed_factorial(n - 1, pF, q);
}

constexpr long double f_factorial(int n, long double (*pF) (int, long double), long double q)
{
	long double f = pF(n, q), p = 1.0;

	if( f != 0.0 )
	{
		p = pow(pF(n, q), 2);
	}

	return n <= 0 ? 1 : p * f_factorial(n - 1, pF, q);
}

bool IsEvenNumber(int n)
{
	return n & 1;
}