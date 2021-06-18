#pragma once

#include "utils.h"

void photon_statistics_superposition_n(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double t, long double theta = 0.0)
{
	long double domain = 1 / sqrt(abs(1 - q));

	if( alpha >= domain )
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
	long double norm_sum = 0.0;

	for( auto l = 0; l <= SUM_INFTY; l++ )
	{
		int sign_l = IsEvenNumber(l) ? 1 : -1;
		norm_sum += pow(alpha_squared, l) * (1 + sign_l * cos(theta)) / factorial(l, pF, q);
	}

	long double normalization = 1 / norm_sum;

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

		dataFile << n << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Photon statistics (n as a variable)... OKAY!" << std::endl;
}

void photon_statistics_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, int n, long double theta = 0.0)
{
	long double domain = 1 / sqrt(abs(1 - q));

	if( alpha >= domain )
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
	long double norm_sum = 0.0;

	for( auto l = 0; l <= SUM_INFTY; l++ )
	{
		int sign_l = IsEvenNumber(l) ? 1 : -1;
		norm_sum += pow(alpha_squared, l) * (1 + sign_l * cos(theta)) / factorial(l, pF, q);
	}

	long double normalization = 1 / norm_sum;
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