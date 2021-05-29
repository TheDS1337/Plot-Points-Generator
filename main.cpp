#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "utils.h"

using namespace std;

#define SUM_INFTY		50				// This is where the fractional series are truncated (65 = magic number for long long unsigned int, 171 for long doubles)
#define MAX_POINTS		50					// More points = more smooth = more processing (overleaf can't handle much)

#define INTERVAL_Q		1.0
#define INTERVAL_ALPHA	5.0
#define INTERVAL_T		1.0
#define INTERVAL_N		50

inline long double f(int x, long double q)
{
	return sqrt((1 - pow(q, x)) / ((1 - q) * x));
}

inline long double deformed_exp(long double x, long double q, long double exponent = 1.0)
{
	long double value = 0.0;

	for( auto n = 1; n <= SUM_INFTY; n++ )
	{
		value += pow((1 - q) * x, n) / ((1 - pow(q, n)) * n);
	}

	return exp(exponent * value);
}

void plot_variances_q(long double alpha)
{
	fstream dataFile;
	dataFile.open("variances_q.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## alpha = " << alpha << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_Q / MAX_POINTS;

	for( auto q = interval; q < INTERVAL_Q; q += interval )
	{
		long double x = 0.0, y = 0.0, x_squared = 0.0, y_squared = 0.0;

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			long double def_factorial = factorial(n, &f, q);
			long double fn1 = f(n + 1, q), fn2 = f(n + 2, q);
			long double sum1 = alpha_squared / fn1, sum2 = alpha_squared / fn2;
			long double prod = pow(alpha_squared, n) / (fn1 * def_factorial);

			x += sqrt(2) * alpha * prod;
			x_squared += prod * (sum1 + sum2);
			y_squared += prod * (sum1 - sum2);
		}

		long double inv_def_exp = deformed_exp(alpha_squared, q, -1.0);

		x *= inv_def_exp;
		x_squared *= inv_def_exp;		x_squared += 0.5;
		y_squared *= inv_def_exp;		y_squared += 0.5;

		long double x_variance = sqrt(x_squared - pow(x, 2));
		long double y_variance = sqrt(y_squared - pow(y, 2));
		long double xy_prod = x_variance * y_variance;

		dataFile << q << ' ' << xy_prod << '\n';
	}

	dataFile.close();

	cout << "Quadratic variances (q as a variable)... OKAY!" << endl;
}

void plot_variances_alpha(long double q)
{
	fstream dataFile;
	dataFile.open("variances_alpha.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << '\n';

	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double x = 0.0, y = 0.0, x_squared = 0.0, y_squared = 0.0;	

		for( auto n = 0; n <= SUM_INFTY; n++ )
		{
			long double def_factorial = factorial(n, &f, q);
			long double fn1 = f(n + 1, q), fn2 = f(n + 2, q);
			long double sum1 = alpha_squared / fn1, sum2 = alpha_squared / fn2;
			long double prod = pow(alpha_squared, n) / (fn1 * def_factorial);

			x += sqrt(2) * alpha * prod;
			x_squared += prod * (sum1 + sum2);
			y_squared += prod * (sum1 - sum2);
		}

		long double inv_def_exp = deformed_exp(alpha_squared, q, -1.0);

		x *= inv_def_exp;
		x_squared *= inv_def_exp;		x_squared += 0.5;
		y_squared *= inv_def_exp;		y_squared += 0.5;

		long double x_variance = sqrt(x_squared - pow(x, 2));
		long double y_variance = sqrt(y_squared - pow(y, 2));
		long double xy_prod = x_variance * y_variance;

		dataFile << alpha << ' ' << xy_prod << '\n';
	}

	dataFile.close();

	cout << "Quadratic variances (alpha as a variable)... OKAY!" << endl;
}

void plot_photon_statistics_n(long double q, long double alpha, long double t, bool even = true)
{
	fstream dataFile;
	dataFile.open("photon_statistics_n.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	if( even )
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << ", even\n";
	}
	else
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << ", odd\n";
	}

	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	int sign = even ? 1 : -1;
	long double alpha_squared = pow(alpha, 2);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double p = 0.0, c = 0.0;
		long double normalization = (1 / (2 * (deformed_exp(alpha_squared, q, 1.0) + sign * deformed_exp(-alpha_squared, q, 1.0)))) * (exp(-n * t) / factorial(n));

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (pow(1 - exp(-t) + LDBL_EPSILON, k) * pow(alpha_squared, n + k)) / (f_factorial(n + k, &f, q) * factorial(k));

			p += d;
			c += IsEvenNumber(n + k) ? d : -d;
		}

		long double sum = normalization * (p + sign * c);

		dataFile << n << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Photon statistics (n as a variable)... OKAY!" << endl;
}

void plot_photon_statistics_t(long double q, long double alpha, int n, bool even = true)
{
	fstream dataFile;
	dataFile.open("photon_statistics_t.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	if( even )
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << ", even\n";
	}
	else
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << ", odd\n";
	}

	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	int sign = even ? 1 : -1;
	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = interval; t < INTERVAL_T; t += interval )
	{
		long double p = 0.0, c = 0.0;
		long double normalization = (1 / (2 * (deformed_exp(alpha_squared, q, 1.0) + sign * deformed_exp(-alpha_squared, q, 1.0)))) * (exp(-n * t) / factorial(n));

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (pow(1 - exp(-t) + LDBL_EPSILON, k) * pow(alpha_squared, n + k)) / (f_factorial(n + k, &f, q) * factorial(k));

			p += d;
			c += IsEvenNumber(n + k) ? d : -d;
		}

		long double sum = normalization * (p + sign * c);

		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Photon statistics (t as a variable)... OKAY!" << endl;
}

void plot_quantum_visibility_n(long double q, long double alpha, long double t)
{
	fstream dataFile;
	dataFile.open("quantum_visibility_n.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << '\n';

	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	long double alpha_squared = pow(alpha, 2);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, &f, q) * factorial(k);

			p += pow(d, k) / e;
			c += pow(-d, k) / e;
		}

		long double v = abs(c) / p;

		dataFile << n << ' ' << v << '\n';
	}

	dataFile.close();

	cout << "Quantum visibility (n as a variable)... OKAY!" << endl;
}

void plot_quantum_visibility_t(long double q, long double alpha, int n)
{
	fstream dataFile;
	dataFile.open("quantum_visibility_t.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << '\n';

	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = interval; t < INTERVAL_T; t += interval )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, &f, q) * factorial(k);

			p += pow(d, k) / e;
			c += pow(-d, k) / e;			
		}

		long double v = abs(c) / p;

		dataFile << t << ' ' << v << '\n';
	}

	dataFile.close();

	cout << "Quantum visibility (t as a variable)... OKAY!" << endl;
}

void plot_fidelity_t(long double q, long double alpha, bool even = true)
{
	fstream dataFile;
	dataFile.open("fidelity_t.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	if( even )
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", even\n";
	}
	else
	{
		dataFile << "## q = " << q << ", alpha = " << alpha << ", odd\n";
	}

	int sign = even ? 1 : -1;
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	long double alpha_squared = pow(alpha, 2);
	long double normalization = pow(2 * (deformed_exp(alpha_squared, q, 1.0) + sign * deformed_exp(-alpha_squared, q, 1.0)), -2);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = interval; t < INTERVAL_T; t += interval )
	{
		long double sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign2 = IsEvenNumber(k) ? 1 : -1;
			long double sum2 = 0.0;
			long double d = pow(alpha, k) * sqrt(pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k));

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign3 = IsEvenNumber(n) ? 1 : -1;
				int g = 1 + sign * sign3 * (1 + sign2) + sign2;
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (sqrt(f_factorial(n + k, &f, q) * f_factorial(n, &f, q)) * factorial(n));

				sum2 += d * e * g;
			}

			sum += pow(sum2, 2);
		}
		
		sum *= normalization;

		dataFile << t << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Fidelity (t as a variable)... OKAY!" << endl;
}

void plot_fidelity_alpha(long double q, long double t, bool even = true)
{
	fstream dataFile;
	dataFile.open("fidelity_alpha.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	if( even )
	{
		dataFile << "## q = " << q << ", t = " << t << ", even\n";
	}
	else
	{
		dataFile << "## q = " << q << ", t = " << t << ", odd\n";
	}

	int sign = even ? 1 : -1;	
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( domain > INTERVAL_ALPHA )
	{
		domain = INTERVAL_ALPHA;
	}

	long double interval = domain / MAX_POINTS;

	for( auto alpha = interval; alpha < domain; alpha += interval )
	{
		long double alpha_squared = pow(alpha, 2);
		long double normalization = pow(2 * (deformed_exp(alpha_squared, q, 1.0) + sign * deformed_exp(-alpha_squared, q, 1.0)), -2);

		long double sum = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			int sign2 = IsEvenNumber(k) ? 1 : -1;
			long double sum2 = 0.0;
			long double d = pow(alpha, k) * sqrt(pow(1 - exp(-t) + LDBL_EPSILON, k) / factorial(k));

			for( auto n = 0; n <= SUM_INFTY; n++ )
			{
				int sign3 = IsEvenNumber(n) ? 1 : -1;
				int g = 1 + sign * sign3 * (1 + sign2) + sign2;
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (sqrt(f_factorial(n + k, &f, q) * f_factorial(n, &f, q)) * factorial(n));

				sum2 += d * e * g;
			}

			sum += pow(sum2, 2);
		}

		sum *= normalization;

		dataFile << alpha << ' ' << sum << '\n';
	}

	dataFile.close();

	cout << "Fidelity (alpha as a variable)... OKAY!" << endl;
}

int main()
{
//	plot_variances(0.5);
//	plot_variances2(0.01);
//	plot_photon_statistics(0.95, 1.0, 1.0);
//	plot_photon_statistics2(0.5, 1.0, 0);
//	plot_quantum_visibility(0.5, 1.0, 1.0); 
//	plot_quantum_visibility2(0.5, 1.0, 0);
//	plot_fidelity_t(0.00001, 0.9, false);
//	plot_fidelity_alpha(0.9999, 1.0, false);

	system("pause");
	return 0;
}