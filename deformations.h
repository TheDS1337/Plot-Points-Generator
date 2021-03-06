#pragma once

#include "utils.h"

long double f_jackson(int x, long double q)
{
	long double epsl = 5 * LDBL_EPSILON;

	// Harmonious
	if( q <= epsl )
	{
		return 1 / sqrt(x);
	}
	// Canonical
	else if( 1 - q <= epsl )
	{
		return 1;
	}

	return sqrt((1 - pow(q, x)) / ((1 - q) * x));
}

long double f_macfarlane(int x, long double q)
{
	return sqrt((pow(q, -x) - pow(q, x)) / ((1 / q - q) * x));
}

long double f_tsallis(int x, long double q)
{
	long double epsl = 5 * LDBL_EPSILON;

	// Harmonious
	if( 2 - q <= epsl )
	{
		return 1 / sqrt(x);
	}
	// Canonical
	else if( q - 1 <= epsl )
	{		
		return 1;
	}

	return 1 / sqrt(1 + (q - 1) * (x - 1.0));
}

long double deformed_exp(long double x, long double (*pF) (int, long double), long double q, long double exponent = 1.0)
{
	if( pF == &f_jackson )
	{
		long double value = 0.0;

		for( auto n = 1; n <= SUM_INFTY; n++ )
		{
			value += pow((1 - q) * x, n) / ((1 - pow(q, n)) * n);
		}

		return exp(exponent * value);
	}
	else if( pF == &f_macfarlane )
	{
		long double value = 0.0;

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			value += pow(x, n) / factorial(n, pF, q);
		}

		return pow(value, exponent);
	}
	else if( pF == &f_tsallis )
	{
		long double y = 1 + (1 - q) * x;

		if( y <= 0.0 )
		{
			return 0.0;
		}

		return pow(y, exponent / (1 - q));
	}

	return exp(x * exponent);
}

long double normalization_factor(long double (*pF) (int, long double), long double q, long double alpha_squared)
{
	return deformed_exp(alpha_squared, pF, q, -1.0);
}

long double normalization_factor_superposition(long double (*pF) (int, long double), long double q, long double alpha_squared, long double theta)
{
	long double sum = 0.0;

	for( auto l = 0; l <= SUM_INFTY; l++ )
	{
		int sign = IsEvenNumber(l) ? 1 : -1;
		sum += pow(alpha_squared + LDBL_EPSILON, l) * (1 + sign * cos(theta)) / factorial(l, pF, q);
	}

	return 1 / sum;
}

long double alpha_domain(long double (*pF) (int, long double), long double q)
{
	return pF == &f_macfarlane ? 1 / sqrt(1 - pow(q, 2)) : 1 / sqrt(abs(1 - q));
}

bool exponential_convergence_check(long double (*pF) (int, long double), long double q, long double alpha)
{
	return alpha < alpha_domain(pF, q) ? true : false;
}