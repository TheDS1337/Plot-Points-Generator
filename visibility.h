#pragma once

#include "utils.h"

void visibility_superposition_n(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double t, long double dn = 0.0)
{
	if( !exponential_convergence_check(pF, q, alpha) )
	{
		std::cout << "alpha is out of range." << std::endl;
		return;
	}

	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << '\n';

	long double alpha_squared = pow(alpha, 2);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, pF, q) * factorial(k);

			p += pow(d, k) / e;
			c += pow(-d, k) / e;
		}

		long double v = abs(c) / p;

		dataFile << n + dn << ' ' << v << '\n';
	}

	dataFile.close();

	std::cout << "Quantum visibility (n as a variable)... OKAY!" << std::endl;
}

void visibility_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, int n)
{
	if( !exponential_convergence_check(pF, q, alpha) )
	{
		std::cout << "alpha is out of range." << std::endl;
		return;
	}

	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_T / 500.0;

	for( auto t = LDBL_EPSILON; t <= 10.0; t += interval )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, pF, q) * factorial(k);

			p += pow(d, k) / e;
			c += pow(-d, k) / e;
		}

		long double v = abs(c) / p;

		dataFile << t << ' ' << v << '\n';
	}

	dataFile.close();

	std::cout << "Quantum visibility (t as a variable)... OKAY!" << std::endl;
}