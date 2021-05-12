#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "utils.h"

using namespace std;

#define SUM_INFTY		171					// This is where the series are truncated (65 = magic number for long long unsigned int, 171 for long doubles)
#define MAX_POINTS		1000				// More points = more smooth = more processing (overleaf can't handle much)

inline long double f(float x)
{
	return 1 / sqrt(x);
}

inline long double deformed_exp(float x)
{
	long double value = 0.0;

	for( auto n = 0; n <= SUM_INFTY; n++ )
	{
		value += pow(x, n) / deformed_factorial(n, &f);
	}

	return value;
}

void plot_fidelity()
{
	float alpha = 1.0;

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

/*		if( t == 0.0 )
		{
			normalization_constant = sum;
			sum = 1.0;
		}
		else
		{
			sum /= normalization_constant;
		}
*/
		sum *= normalization_constant;
		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Fidelity... OKAY!" << endl;
}

int main()
{
	plot_fidelity();
	system("pause");
	return 0;
}