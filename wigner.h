#pragma once

#include "utils.h"

void wigner(std::string filename, long double (*pF) (int, long double), long double q, long double alpha)
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
	long double normalization = 2 * normalization_factor(pF, q, alpha_squared) / M_PI;
	long double domainXY = 3.0;
	long double interval = domainXY / 10;
	long double intervalX = interval;
	long double intervalY = interval;

	for( long double x = -domainXY; x < domainXY; x += intervalX )
	{
		intervalX = interval / (1.0 + 9.0 * exp(-pow(x + LDBL_EPSILON, 2)));
		
		for( long double y = -domainXY; y < domainXY; y += intervalY )
		{
			intervalY = interval / (1.0 + 9.0 * exp(-pow(y + LDBL_EPSILON, 2)));

			long double beta_squared = pow(x + LDBL_EPSILON, 2) + pow(y + LDBL_EPSILON, 2);
			long double beta = sqrt(beta_squared), phi = atan2(y, x);

			long double sum = 0.0;

			for( auto k = 0; k <= 20; k++ )
			{
				int sign = IsEvenNumber(k) ? 1 : -1;
				long double sum2_cos_1 = 0.0, sum2_sin_1 = 0.0;
				long double sum2_cos_2 = 0.0, sum2_sin_2 = 0.0;
				long double sum2_cos_3 = 0.0, sum2_cos_4 = 0.0, sum2_sin_3 = 0.0, sum2_sin_4 = 0.0;

				for( auto n = 0; n <= 20; n++ )
				{
					int sign2 = IsEvenNumber(n) ? 1 : -1;
					long double r = pow(alpha, n) * exp(-beta_squared / 2) / f_factorial(n, pF, q, false);

					if( n <= k )
					{
						long double v = sign2 * r * pow(beta, k - n) * laguerre(n, k - n, beta_squared) / sqrt(factorial(k));
						long double u = sign2 * r * pow(beta, -n) * laguerre(n, k - n, beta_squared);

						sum2_cos_1 += v * cos(n * phi);
						sum2_sin_1 += v * sin(n * phi);

						sum2_cos_3 += u * cos(n * phi);
						sum2_sin_3 += u * sin(n * phi);
					}
					else
					{
						long double v = r * pow(beta, n - k) * laguerre(k, n - k, beta_squared) * sqrt(factorial(k)) / factorial(n);
						long double u = sign * r * pow(beta, n) * laguerre(k, n - k, beta_squared) / factorial(n);

						sum2_cos_2 += v * cos(n * phi);
						sum2_sin_2 += v * sin(n * phi);

						sum2_cos_4 += u * cos(n * phi);
						sum2_sin_4 += u * sin(n * phi);
					}
				}

				sum += sign * (pow(sum2_cos_1 + LDBL_EPSILON, 2) + pow(sum2_cos_2 + LDBL_EPSILON, 2) + pow(sum2_sin_1 + LDBL_EPSILON, 2) + pow(sum2_sin_2 + LDBL_EPSILON, 2) + 2 * (sum2_cos_3 * sum2_cos_4 + sum2_sin_3 * sum2_sin_4));
			}

			sum *= normalization;

			dataFile << x << ' ' << y << ' ' << sum << '\n';
		}
	}

	dataFile.close();

	std::cout << "Wigner... OKAY!" << std::endl;
}

void wigner_superposition(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
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
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta) / M_PI;
	long double domainXY = 4.0;
	long double interval = domainXY / 10;
	long double intervalX = interval;
	long double intervalY = interval;

	for( long double x = -domainXY; x < domainXY; x += intervalX )
	{
		intervalX = interval / (1.0 + 9.0 * exp(-pow(x + LDBL_EPSILON, 2)));

		for( long double y = -domainXY; y < domainXY; y += intervalY )
		{
			intervalY = interval / (1.0 + 9.0 * exp(-pow(y + LDBL_EPSILON, 2)));

			long double beta_squared = pow(x + LDBL_EPSILON, 2) + pow(y + LDBL_EPSILON, 2);
			long double beta = sqrt(beta_squared), phi = atan2(y, x);

			long double sum = 0.0;

			for( auto k = 0; k <= 20; k++ )
			{
				int sign = IsEvenNumber(k) ? 1 : -1;
				long double sum2_cos_1_1 = 0.0, sum2_cos_1_2 = 0.0, sum2_sin_1_1 = 0.0, sum2_sin_1_2 = 0.0;
				long double sum2_cos_2_1 = 0.0, sum2_cos_2_2 = 0.0, sum2_sin_2_1 = 0.0, sum2_sin_2_2 = 0.0;
				long double sum2_cos_3_1 = 0.0, sum2_cos_3_2 = 0.0, sum2_cos_4_1 = 0.0, sum2_cos_4_2 = 0.0, sum2_sin_3_1 = 0.0, sum2_sin_3_2 = 0.0, sum2_sin_4_1 = 0.0, sum2_sin_4_2 = 0.0;

				for( auto n = 0; n <= 20; n++ )
				{
					int sign2 = IsEvenNumber(n) ? 1 : -1;
					long double r = pow(alpha, n) * exp(-beta_squared / 2) / f_factorial(n, pF, q, false);

					if( n <= k )
					{
						long double v = sign2 * r * pow(beta, k - n) * laguerre(n, k - n, beta_squared) / sqrt(factorial(k));
						long double u = sign2 * r * pow(beta, -n) * laguerre(n, k - n, beta_squared);

						sum2_cos_1_1 += v * cos(n * phi);
						sum2_cos_1_2 += sign2 * v * cos(n * phi);
						sum2_sin_1_1 += v * sin(n * phi);
						sum2_sin_1_2 += sign2 * v * sin(n * phi);

						sum2_cos_3_1 += u * cos(n * phi);
						sum2_cos_3_2 += sign2 * u * cos(n * phi);
						sum2_sin_3_1 += u * sin(n * phi);
						sum2_sin_3_2 += sign2 * u * sin(n * phi);
					}
					else
					{
						long double v = r * pow(beta, n - k) * laguerre(k, n - k, beta_squared) * sqrt(factorial(k)) / factorial(n);
						long double u = sign * r * pow(beta, n) * laguerre(k, n - k, beta_squared) / factorial(n);

						sum2_cos_2_1 += v * cos(n * phi);
						sum2_cos_2_2 += sign2 * v * cos(n * phi);
						sum2_sin_2_1 += v * sin(n * phi);
						sum2_sin_2_2 += sign2 * v * sin(n * phi);

						sum2_cos_4_1 += u * cos(n * phi);
						sum2_cos_4_2 += sign2 * u * cos(n * phi);
						sum2_sin_4_1 += u * sin(n * phi);
						sum2_sin_4_2 += sign2 * u * sin(n * phi);
					}
				}

				long double sum_1 = pow(sum2_cos_1_1 + LDBL_EPSILON, 2) + pow(sum2_cos_1_2 + LDBL_EPSILON, 2) + pow(sum2_sin_1_1 + LDBL_EPSILON, 2) + pow(sum2_sin_1_2 + LDBL_EPSILON, 2) + 2 * cos(theta) * (sum2_cos_1_1 * sum2_cos_1_2 + sum2_sin_1_1 * sum2_sin_1_2);
				long double sum_2 = pow(sum2_cos_2_1 + LDBL_EPSILON, 2) + pow(sum2_cos_2_2 + LDBL_EPSILON, 2) + pow(sum2_sin_2_1 + LDBL_EPSILON, 2) + pow(sum2_sin_2_2 + LDBL_EPSILON, 2) + 2 * cos(theta) * (sum2_cos_2_1 * sum2_cos_2_2 + sum2_sin_2_1 * sum2_sin_2_2);
				long double sum_3 = sum2_cos_3_1 * sum2_cos_4_1 + sum2_cos_3_2 * sum2_cos_4_2 + sum2_sin_3_1 * sum2_sin_4_1 + sum2_sin_3_2 * sum2_sin_4_2 + cos(theta) * (sum2_cos_3_1 * sum2_cos_4_2 + sum2_cos_3_2 * sum2_cos_4_1 + sum2_sin_3_1 * sum2_sin_4_2 + sum2_sin_3_2 * sum2_sin_4_1);

				sum += sign * (sum_1 + sum_2 + 2 * sum_3);
			}

			sum *= normalization;

			dataFile << x << ' ' << y << ' ' << sum << '\n';
		}

	}

	dataFile.close();

	std::cout << "Wigner... OKAY!" << std::endl;
}

void wigner_superposition_batouli(std::string filename, long double (*pF) (int, long double), long double q, long double alpha, long double theta = 0.0)
{
	using namespace std::complex_literals;
	std::complex<long double> imaginary = 1i;

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
	long double normalization = normalization_factor_superposition(pF, q, alpha_squared, theta) / M_PI;
	long double domainXY = 4.0;
	long double interval = domainXY / 10;
	long double intervalX = interval;
	long double intervalY = interval;

	for( long double x = -domainXY; x < domainXY; x += intervalX )
	{
		intervalX = interval / (1.0 + 9.0 * exp(-pow(x + LDBL_EPSILON, 2)));

		for( long double y = -domainXY; y < domainXY; y += intervalY )
		{
			intervalY = interval / (1.0 + 9.0 * exp(-pow(y + LDBL_EPSILON, 2)));

			std::complex<long double> beta = x + imaginary * y;
			long double beta_squared = norm(beta);

			std::complex<long double> sum = 0.0;

			for( auto n = 0; n <= 20; n++ )
			{
				int sign = IsEvenNumber(n) ? 1 : -1;
				std::complex<long double> c_n = pow(alpha, n) * (std::complex<long double>(1.0) + exp(imaginary * std::complex<long double>(theta + n * M_PI))) / sqrt(factorial(n, pF, q));

				for( auto m = 0; m <= 20; m++ )
				{
					std::complex<long double> c_m = pow(alpha, m) * (std::complex<long double>(1.0) + exp(-imaginary * std::complex<long double>(theta + m * M_PI))) / sqrt(factorial(m, pF, q));

					if( m <= n )
					{
						sum += std::complex<long double>(sign) * c_n * c_m * sqrt(factorial(m) / factorial(n)) * exp(-2 * beta_squared) * pow(-std::complex<long double>(2) * conj(beta), n - m) * laguerre(m, n - m, 4 * beta_squared);
					}
					else
					{
						sum += std::complex<long double>(sign) * c_n * c_m * sqrt(factorial(n) / factorial(m)) * exp(-2 * beta_squared) * pow(std::complex<long double>(2) * beta, m - n) * laguerre(n, m - n, 4 * beta_squared);
					}
				}
			}

			sum *= normalization;

			dataFile << x << ' ' << y << ' ' << real(sum) << '\n';
		}

	}

	dataFile.close();

	std::cout << "Wigner (BATOULI)... OKAY!" << std::endl;
}