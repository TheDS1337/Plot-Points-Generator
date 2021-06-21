#pragma once

#include "utils.h"

void second_order_correlation_parameter_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
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

	dataFile << "## q = " << q << ", alpha = " << alpha << ", theta = " << theta << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		long double n_sum = 0.0, n_squared_sum = 0.0;

		for( auto n = 1; n <= SUM_INFTY; n++ )
		{
			int sign = IsEvenNumber(n) ? 1 : -1;
			long double d = pow(alpha_squared, n) * (1 + sign * cos(theta)) / f_factorial(n, pF, q);

			n_sum += d * exp(-t) / factorial(n - 1);

			if( n >= 2 )
			{
				n_squared_sum += d * exp(-2 * t) / factorial(n - 2);
			}
		}

		n_squared_sum += n_sum;
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

	long double domain = alpha_domain(pF, q);

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta);
		long double n_sum = 0.0, n_squared_sum = 0.0;

		for( auto n = 1; n <= SUM_INFTY; n++ )
		{
			int sign = IsEvenNumber(n) ? 1 : -1;
			long double d = pow(alpha_squared, n) * (1 + sign * cos(theta)) / f_factorial(n, pF, q);

			n_sum += d * exp(-t) / factorial(n - 1);

			if( n >= 2 )
			{
				n_squared_sum += d * exp(-2 * t) / factorial(n - 2);
			}
		}

		n_squared_sum += n_sum;
		long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}
	

	dataFile.close();

	std::cout << "Second order correlation parameter (alpha as a variable)... OKAY!" << std::endl;
}