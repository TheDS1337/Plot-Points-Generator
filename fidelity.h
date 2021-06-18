#pragma once

#include "utils.h"

void fidelity_superposition_t(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
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
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta) / 4;
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		long double sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign = IsEvenNumber(k) ? 1 : -1;
			long double sum2_real = 0.0, sum2_img = 0.0;
			long double d = pow(alpha_squared, k) * pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k);

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign2 = IsEvenNumber(n) ? 1 : -1;
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (f_factorial(n + k, pF, q, false) * f_factorial(n, pF, q, false) * factorial(n));

				sum2_real += e * (1 + sign + sign2 * (sign + 1) * cos(theta));
				sum2_img += e * sign2 * (sign - 1) * sin(theta);
			}

			sum += d * (pow(sum2_real + LDBL_EPSILON, 2) + pow(sum2_img + LDBL_EPSILON, 2));
		}

		sum *= normalization;

		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Fidelity (t as a variable)... OKAY!" << std::endl;
}

void fidelity_superposition_alpha(std::string filename, long double (*pF) (int, long double), long double q, long double t, long double theta = 0.0)
{
	std::fstream dataFile;
	dataFile.open(filename, std::ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", theta = " << theta << ", t = " << t << '\n';

	long double domain = 1 / sqrt(abs(1 - q));

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = LDBL_EPSILON; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta) / 4;
		long double sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign = IsEvenNumber(k) ? 1 : -1;
			long double sum2_real = 0.0, sum2_img = 0.0;
			long double d = pow(alpha_squared, k) * pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k);

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign2 = IsEvenNumber(n) ? 1 : -1;
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (f_factorial(n + k, pF, q, false) * f_factorial(n, pF, q, false) * factorial(n));

				sum2_real += e * (1 + sign + sign2 * (sign + 1) * cos(theta));
				sum2_img += e * sign2 * (sign - 1) * sin(theta);
			}

			sum += d * (pow(sum2_real + LDBL_EPSILON, 2) + pow(sum2_img + LDBL_EPSILON, 2));
		}

		sum *= normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}

	dataFile.close();

	std::cout << "Fidelity (alpha as a variable)... OKAY!" << std::endl;
}
