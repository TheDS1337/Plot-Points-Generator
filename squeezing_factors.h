#pragma once

#include "utils.h"

void squeezing_factors_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
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
		long double x_sum = 0.0, x_squared_sum = 0.0;
		long double p_sum = 0.0, p_squared_sum = 0.0;

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			int sign = IsEvenNumber(n) ? 1 : -1;

			p_sum += pow(alpha_squared, n) * sign * sin(theta) * exp(-0.5 * t) / (factorial(n) * pF(n + 1, q) * f_factorial(n, pF, q));

			if( n >= 1 )
			{
				long double d = pow(alpha_squared, n) * exp(-t) / (factorial(n - 1) * f_factorial(n, pF, q));
				long double e = 1 + sign * cos(theta);
				long double f = (1 - sign * cos(theta)) * pF(n, q) / pF(n + 1, q);

				x_squared_sum += d * (e + f);
				p_squared_sum += d * (e - f);
			}
		}

		x_sum = 0.0;
		x_squared_sum = 0.5 + normalization * x_squared_sum;

		p_sum *= -sqrt(2) * alpha * normalization;
		p_squared_sum = 0.5 + normalization * p_squared_sum;

		long double i1 = (x_squared_sum - pow(x_sum + LDBL_EPSILON, 2)) - 0.5;
		long double i2 = (p_squared_sum - pow(p_sum + LDBL_EPSILON, 2)) - 0.5;

		dataFile << t << ' ' << i1 << ' ' << i2 << '\n';
	}

	dataFile.close();

	std::cout << "Squeezing factors (t as a variable)... OKAY!" << std::endl;
}

void squeezing_factors_superposition_alpha(std::string filename, long double (*pF) (int, long double), long double q, long double t, long double theta = 0.0)
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

		long double x_sum = 0.0, x_squared_sum = 0.0;
		long double p_sum = 0.0, p_squared_sum = 0.0;

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			int sign = IsEvenNumber(n) ? 1 : -1;

			p_sum += pow(alpha_squared, n) * sign * sin(theta) * exp(-0.5 * t) / (factorial(n) * pF(n + 1, q) * f_factorial(n, pF, q));

			if( n >= 1 )
			{
				long double d = pow(alpha_squared, n) * exp(-t) / (factorial(n - 1) * f_factorial(n, pF, q));
				long double e = 1 + sign * cos(theta);
				long double f = (1 - sign * cos(theta)) * pF(n, q) / pF(n + 1, q);

				x_squared_sum += d * (e + f);
				p_squared_sum += d * (e - f);
			}
		}

		x_sum = 0.0;
		x_squared_sum = 0.5 + normalization * x_squared_sum;

		p_sum *= -sqrt(2) * alpha * normalization;
		p_squared_sum = 0.5 + normalization * p_squared_sum;

		long double i1 = (x_squared_sum - pow(x_sum + LDBL_EPSILON, 2)) - 0.5;
		long double i2 = (p_squared_sum - pow(p_sum + LDBL_EPSILON, 2)) - 0.5;

		dataFile << alpha << ' ' << i1 << ' ' << i2 << '\n';
	}
	
	dataFile.close();

	std::cout << "Squeezing factors (alpha as a variable)... OKAY!" << std::endl;
}