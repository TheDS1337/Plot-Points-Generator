#pragma once

#include "utils.h"

void second_order_correlation_parameter_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
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

	dataFile << "## q = " << q << ", alpha = " << alpha << ", theta = " << theta << '\n';

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
		long double n_sum = 0.0, n_squared_sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign = IsEvenNumber(k) ? 1 : -1;
			long double d = pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k);

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign2 = IsEvenNumber(n) ? 1 : -1;
				long double e = d * pow(alpha_squared, n + k) * (1 + sign * sign2 * cos(theta)) * exp(-n * t) / (factorial(n) * f_factorial(n + k, pF, q));

				n_sum += e * n;
				n_squared_sum += e * pow(n + LDBL_EPSILON, 2);
			}
		}

		long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Second order correlation parameter (t as a variable)... OKAY!" << std::endl;
}

void second_order_correlation_parameter_superposition_alpha(std::string filename, long double (*pF) (int, long double), long double q, long double t, long double theta = 0.0)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", t = " << t << ", theta = " << theta << '\n';

	long double domain = 1 / sqrt(abs(1 - q));

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double norm_sum = 0.0;

		for( auto l = 0; l <= SUM_INFTY; l++ )
		{
			int sign_l = IsEvenNumber(l) ? 1 : -1;
			norm_sum += pow(alpha_squared, l) * (1 + sign_l * cos(theta)) / factorial(l, pF, q);
		}

		long double normalization = 1 / norm_sum;
		long double n_sum = 0.0, n_squared_sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign = IsEvenNumber(k) ? 1 : -1;
			long double d = pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k);

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign2 = IsEvenNumber(n) ? 1 : -1;
				long double e = d * pow(alpha_squared, n + k) * (1 + sign * sign2 * cos(theta)) * exp(-n * t) / (factorial(n) * f_factorial(n + k, pF, q));

				n_sum += e * n;
				n_squared_sum += e * pow(n + LDBL_EPSILON, 2);
			}
		}

		long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}
	

	dataFile.close();

	std::cout << "Second order correlation parameter (alpha as a variable)... OKAY!" << std::endl;
}