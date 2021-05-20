#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "utils.h"

using namespace std;

#define SUM_INFTY		2000				// This is where the fractional series are truncated (65 = magic number for long long unsigned int, 171 for long doubles)
#define MAX_POINTS		100					// More points = more smooth = more processing (overleaf can't handle much)

#define Q_DEFORMER		1.0					// q = 1 is the undeformed limit, q = 0 is the harmonious limit and anything 0 < q < 1 is a case for study.

inline long double f(int x, float q)
{
	return sqrt((1 - pow(q, x)) / ((1 - q) * x));
}

inline long double f(int x)
{						
	return f(x, Q_DEFORMER);
}

inline long double deformed_exp(float x, float q)
{
	long double value = 0.0;

	for( auto n = 1; n <= SUM_INFTY; n++ )
	{
		value += pow((1 - q) * x, n) / ((1 - pow(q, x)) * x);
	}

	return exp(value);
}

inline long double deformed_exp(float x)
{
	return deformed_exp(x, Q_DEFORMER);
}

void plot_variances(float alpha)
{
	fstream dataFile;
	dataFile.open("variances.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	float alpha_squared = pow(alpha, 2);

	float interval = 1.0 / MAX_POINTS;
	float end_point = 1.0 - interval;

	for( auto q = interval; q <= end_point; q += interval )
	{
		long double x = 0.0, y = 0.0, x_squared = 0.0, y_squared = 0.0;
		long double inv_def_exp_sum = 0.0;

		for( auto k = 1; k <= SUM_INFTY; k++ )
		{
			inv_def_exp_sum += pow((1 - q) * alpha_squared, k) / (k * (1 - pow(q, k)));
		}

		long double inv_def_exp = 1 / sqrt(exp(inv_def_exp_sum));

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			long double def_factorial = deformed_factorial(n, &f, q);
			long double fn1 = f(n + 1, q), fn2 = f(n + 2, q);
			long double sum1 = alpha_squared / fn1, sum2 = alpha_squared / fn2;
			long double prod = pow(alpha_squared, n) / (fn1 * def_factorial);

			x += sqrt(2) * alpha * prod;
			x_squared += prod * (sum1 + sum2);
			y_squared += prod * (sum1 - sum2);
		}

		x *= inv_def_exp;
		x_squared *= inv_def_exp;		x_squared += 0.5;
		y_squared *= inv_def_exp;		y_squared += 0.5;

		long double x_variance = sqrt(x_squared - pow(x, 2));
		long double y_variance = sqrt(y_squared - pow(y, 2));
		long double xy_prod = x_variance * y_variance;

		dataFile << q << ' ' << xy_prod << '\n';
	}

	dataFile.close();

	cout << "Quadratic variances... OKAY!" << endl;
}
/*
void plot_fidelity()
{
	float alpha = 0.5;

	fstream dataFile;
	dataFile.open("fidelity.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	const double interval = 1.0 / MAX_POINTS;
	for( auto t = 0.0; t <= 1.0; t += interval )
	{
		long double sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k += 2 )
		{
			long double sum1 = 0.0;

			for( auto n = k; n <= SUM_INFTY; n += 2 )
			{
				sum1 += 4 * sqrt(combination(n, k)) * exp(((k - n) / 2) * t) * pow(1 - exp(-t), k / 2) * (pow(alpha, n) / sqrt(deformed_factorial(n, &f))) * (pow(alpha, n - k) / sqrt(deformed_factorial(n - k, &f)));
			}

			sum += pow(sum1, 2);
		}

		long double exp1 = deformed_exp(pow(alpha, 2));
		long double exp2 = deformed_exp(-pow(alpha, 2));
		long double normalization_constant = 1 / pow(2 * (exp1 + exp2), 2);

//				if( t == 0.0 )
//				{
//					normalization_constant = sum;
//					sum = 1.0;
//				}
//				else
//				{
//					sum /= normalization_constant;
//				}
		
		sum *= normalization_constant;
		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Fidelity... OKAY!" << endl;
}
*/

int main()
{
	plot_variances(0.5);
//	plot_fidelity();
	system("pause");
	return 0;
}