#pragma once

#include "utils.h"

void second_order_correlation_parameter_alpha(std::string filename, long double (*pF) (int, long double), long double q)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << '\n';

	long double domain = alpha_domain(pF, q);

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double normalization = normalization_factor(pF, q, alpha_squared);
		long double n_sum = 0.0, n_squared_sum = 0.0;

		for( auto n = 1; n <= SUM_INFTY; n++ )
		{
			long double d = pow(alpha_squared, n) / f_factorial(n, pF, q);

			n_sum += d / factorial(n - 1);

			if( n >= 2 )
			{
				n_squared_sum += d / factorial(n - 2);
			}
		}

		n_squared_sum += n_sum;
		long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}


	dataFile.close();

	std::cout << "Second order correlation parameter (alpha as a variable)... OKAY!" << std::endl;
}

void second_order_correlation_parameter_superposition_alpha(std::string filename, long double (*pF) (int, long double), long double q, long double theta = 0.0)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", theta = " << theta << '\n';

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

			n_sum += d / factorial(n - 1);

			if( n >= 2 )
			{
				n_squared_sum += d / factorial(n - 2);
			}
		}

		n_squared_sum += n_sum;
		long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}
	

	dataFile.close();

	std::cout << "Second order correlation parameter (alpha as a variable)... OKAY!" << std::endl;
}

void second_order_correlation_parameter_superposition(std::string filename, long double (*pF) (int, long double), long double q)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << '\n';

	long double domainAlpha = alpha_domain(pF, q);

	if( domainAlpha > INTERVAL_ALPHA )
	{
		domainAlpha = INTERVAL_ALPHA;
	}

	long double intervalTheta = INTERVAL_THETA / MAX_POINTS;

	for( auto theta = 0.0; theta < INTERVAL_THETA; theta += intervalTheta )
	{
		long double intervalAlpha = domainAlpha / MAX_POINTS;

		for( auto alpha = 10 * intervalAlpha; alpha < domainAlpha; alpha += intervalAlpha )
		{
			long double alpha_squared = pow(alpha, 2);
			long double normalization = normalization_factor_superposition(pF, q, alpha_squared, intervalAlpha);
			long double n_sum = 0.0, n_squared_sum = 0.0;

			for( auto n = 1; n <= SUM_INFTY; n++ )
			{
				int sign = IsEvenNumber(n) ? 1 : -1;
				long double d = pow(alpha_squared, n) * (1 + sign * cos(theta)) / f_factorial(n, pF, q);

				n_sum += d / factorial(n - 1);

				if( n >= 2 )
				{
					n_squared_sum += d / factorial(n - 2);
				}
			}

			n_squared_sum += n_sum;
			long double sum = (n_squared_sum - n_sum) * pow(n_sum, -2) / normalization;

			dataFile << theta << ' ' << alpha << ' ' << sum << '\n';
		}
	}


	dataFile.close();

	std::cout << "Second order correlation parameter... OKAY!" << std::endl;
}