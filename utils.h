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

constexpr long double permutation(int n, int p)
{
	return factorial(n) / factorial(p);
}

constexpr long double combination(int n, int p)
{
	return permutation(n, p) / factorial(p);
}

constexpr long double deformed_factorial(int n, long double (*pF) (float))
{
	return n <= 0.0 ? 1.0 : n * pF((float) n) * deformed_factorial(n - 1, pF);
}

bool IsEvenNumber(int n)
{
	return n & 1 ? true : false;
}