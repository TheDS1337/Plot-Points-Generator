#pragma once

#pragma once

#include "utils.h"

void heisenberg_alpha(std::string filename, long double (*pF) (int, long double), long double q)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << '\n';

	long double domain = alpha_domain(pF, q);

	if( domain > 3.0 )
	{
		domain = 3.0;
	}

	long double interval = domain / 100;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double normalization = normalization_factor(pF, q, alpha_squared);

		long double sum = 0.0;
		long double x_sum = 0.0, x_squared_sum = 0.0;
		long double p_sum = 0.0, p_squared_sum = 0.0;

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			long double d = pow(alpha_squared, n) / factorial(n, pF, q);
			long double e = pow(pF(n + 1, q), 2) + alpha_squared * (1.0 - pow(pF(n + 2, q) / pF(n + 1, q), 2));

			sum += d * e;
		}

		x_sum = sqrt(2) * alpha;
		x_squared_sum = alpha_squared + 0.5 * (normalization * sum + 2 * alpha_squared);

		p_sum = 0.0;
		p_squared_sum = alpha_squared + 0.5 * (normalization * sum - 2 * alpha_squared);

		long double dx = sqrt(x_squared_sum - pow(x_sum + LDBL_EPSILON, 2));
		long double dp = sqrt(p_squared_sum - pow(p_sum + LDBL_EPSILON, 2));

		long double heisenberg = dx * dp;

		dataFile << alpha << ' ' << heisenberg << '\n';
	}

	dataFile.close();

	std::cout << "Heisenberg (alpha as a variable)... OKAY!" << std::endl;
}