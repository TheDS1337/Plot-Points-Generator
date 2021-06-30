#pragma once

#include "utils.h"

void photon_statistics_n(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double dn = 0.0)
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

	dataFile << "## q = " << q << ", alpha = " << alpha << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double normalization = normalization_factor(pF, q, alpha_squared);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double sum = normalization * pow(alpha_squared, n) / factorial(n, pF, q);

		dataFile << n + dn << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Photon statistics (n as a variable)... OKAY!" << std::endl;
}

void photon_statistics_superposition_n(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double t, long double theta = 0.0, long double dn = 0.0)
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

	dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << ", theta = " << theta << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		int sign = IsEvenNumber(n) ? 1 : -1;
		long double sum2 = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign2 = IsEvenNumber(k) ? 1 : -1;
			sum2 += pow(alpha_squared, n + k) * (1 + (sign * sign2 * cos(theta))) * pow(1 - exp(-t) + LDBL_EPSILON, k) / (f_factorial(n + k, pF, q) * factorial(k));
		}

		long double sum = normalization * sum2 * exp(-n * t) / factorial(n);

		dataFile << n + dn << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Photon statistics (n as a variable)... OKAY!" << std::endl;
}

void photon_statistics_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, int n, long double theta = 0.0)
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


	dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << ", theta = " << theta << '\n';

	int sign = IsEvenNumber(n) ? 1 : -1;
	long double alpha_squared = pow(alpha, 2);
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		long double sum2 = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign2 = IsEvenNumber(k) ? 1 : -1;
			sum2 += pow(alpha_squared, n + k) * (1 + (1 + (sign * sign2 * cos(theta)))) * pow(1 - exp(-t) + LDBL_EPSILON, k) / (f_factorial(n + k, pF, q) * factorial(k));
		}

		long double sum = normalization * sum2 * exp(-n * t) / factorial(n);

		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Photon statistics (t as a variable)... OKAY!" << std::endl;
}

void photon_statistics_superposition(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double t)
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
	long double interval = INTERVAL_THETA / MAX_POINTS;

	for( auto theta = LDBL_EPSILON; theta <= INTERVAL_THETA; theta += interval )
	{
		long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta);

		for( auto n = 0; n <= 30; n++ )
		{
			int sign = IsEvenNumber(n) ? 1 : -1;
			long double sum2 = 0.0;

			for( auto k = 0; k <= SUM_INFTY; k++ )
			{
				int sign2 = IsEvenNumber(k) ? 1 : -1;
				sum2 += pow(alpha_squared, n + k) * (1 + (sign * sign2 * cos(theta))) * pow(1 - exp(-t) + LDBL_EPSILON, k) / (f_factorial(n + k, pF, q) * factorial(k));
			}

			long double sum = normalization * sum2 * exp(-n * t) / factorial(n);

			dataFile << n << ' ' << theta << ' ' << sum << '\n';
		}
	}

	dataFile.close();

	std::cout << "Photon statistics... OKAY!" << std::endl;
}