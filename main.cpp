#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#define SUM_INFTY		10		// This is where the series are truncated
#define MAX_POINTS		5000		// More points = more smooth = more processing (overleaf can't handle much)

constexpr double factorial(int n)
{
	return n <= 0.0 ? 1.0 : n * factorial(n - 1);
}

constexpr double permutation(int n, int p)
{
	return factorial(n) / factorial(p);
}

constexpr double combination(int n, int p)
{	
	return permutation(n, p) / factorial(p);
}

constexpr double deformed_factorial(float x, double (*pF) (float))
{
	return x <= 0.0 ? 1.0 : x * pF(x) * deformed_factorial(x - 1, pF);
}

inline double f(float x)
{
	return 1 / sqrt(x);
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

	float normalization_constant = 0.0;

	const float interval = 2.0 / MAX_POINTS;
	for( auto t = 0.0; t <= 2.0; t += interval )
	{
		float sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			float sum1 = 0.0;

			for( auto n = k; n <= SUM_INFTY; n++ )
			{
				sum1 += sqrt(combination(n, k)) * exp(((k - n) / 2) * t) * pow(1 - exp(-t), k / 2) * (((1 + pow(-1, n)) * pow(alpha, n)) / sqrt(deformed_factorial(n, &f))) * (((1 + pow(-1, n - k)) * pow(alpha, n - k)) / sqrt(deformed_factorial(n - k, &f)));
			}

			sum += pow(sum1, 2);
		}

		if( t == 0.0 )
		{
			normalization_constant = sum;
			sum = 1.0;
		}
		else
		{
			sum /= normalization_constant;
		}

		dataFile << t << '\t' << sum << '\n';
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