#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <chrono>
#include <ctime> 

#include "utils.h"

using namespace std;
using namespace std::complex_literals;

#define SUM_INFTY				50					// This is where the fractional series are truncated (65 = magic number for long long unsigned int, 171 for long doubles)
#define MAX_POINTS				50					// More points = more smooth = more processing (overleaf can't handle much)

#define INTERVAL_Q				1.0
#define INTERVAL_ALPHA			5.0
#define INTERVAL_T				1.0
#define INTERVAL_N				50

// Compiler can't handle conversion by itsellf, these needs to be here
const complex<long double> zero = 0.0;
const complex<long double> one = 1.0;
const complex<long double> two = 2.0;
const complex<long double> four = 4.0;
const complex<long double> imaginary = 1.0i;

inline long double f_jackson(int x, long double q)
{
	long double epsl = 5 * LDBL_EPSILON;

	// Harmonious
	if( q <= epsl )
	{
		return 1 / sqrt(x);
	}
	// Canonical
	else if( 1 - q <= epsl )
	{
		return 1;
	}

	return sqrt((1 - pow(q, x)) / ((1 - q) * x));
}

inline long double f_tsallis(int x, long double q)
{
	long double epsl = 5 * LDBL_EPSILON;

	// Harmonious
	if( 2 - q <= epsl )
	{
		return 1 / sqrt(x);
	}
	// Canonical
	else if( q - 1 <= epsl )
	{		
		return 1;
	}

	return 1 / sqrt(1 + (q - 1) * (x - 1.0));
}

inline long double deformed_exp_jackson(long double x, long double q, long double exponent = 1.0)
{
	long double value = 0.0;

	for( auto n = 1; n <= SUM_INFTY; n++ )
	{
		value += pow((1 - q) * x, n) / ((1 - pow(q, n)) * n);
	}

	return exp(exponent * value);
}

inline complex<long double> deformed_exp_tsallis(complex<long double> z, long double q, long double exponent = 1.0)
{
	if( z == zero )
	{
		return 0.0;
	}

	return pow(one + (one - q) * z, exponent / (1 - q));
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
			long double def_factorial = factorial(n, &f_jackson, q);
			long double fn1 = f_jackson(n + 1, q), fn2 = f_jackson(n + 2, q);
			long double sum1 = alpha_squared / fn1, sum2 = alpha_squared / fn2;
			long double prod = pow(alpha_squared, n) / (fn1 * def_factorial);

			x += sqrt(2) * alpha * prod;
			x_squared += prod * (sum1 + sum2);
			y_squared += prod * (sum1 - sum2);
		}

		long double inv_def_exp = deformed_exp_jackson(alpha_squared, q, -1.0);

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
			long double def_factorial = factorial(n, &f_jackson, q);
			long double fn1 = f_jackson(n + 1, q), fn2 = f_jackson(n + 2, q);
			long double sum1 = alpha_squared / fn1, sum2 = alpha_squared / fn2;
			long double prod = pow(alpha_squared, n) / (fn1 * def_factorial);

			x += sqrt(2) * alpha * prod;
			x_squared += prod * (sum1 + sum2);
			y_squared += prod * (sum1 - sum2);
		}

		long double inv_def_exp = deformed_exp_jackson(alpha_squared, q, -1.0);

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
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

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

	int sign = even ? 1 : -1;
	long double alpha_squared = pow(alpha, 2);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double p = 0.0, c = 0.0;
		long double normalization = (1 / (2 * (deformed_exp_jackson(alpha_squared, q, 1.0) + sign * deformed_exp_jackson(-alpha_squared, q, 1.0)))) * (exp(-n * t) / factorial(n));

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (pow(1 - exp(-t) + LDBL_EPSILON, k) * pow(alpha_squared, n + k)) / (f_factorial(n + k, &f_jackson, q) * factorial(k));

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
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

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

	int sign = even ? 1 : -1;
	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = interval; t < INTERVAL_T; t += interval )
	{
		long double p = 0.0, c = 0.0;
		long double normalization = (1 / (2 * (deformed_exp_jackson(alpha_squared, q, 1.0) + sign * deformed_exp_jackson(-alpha_squared, q, 1.0)))) * (exp(-n * t) / factorial(n));

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (pow(1 - exp(-t) + LDBL_EPSILON, k) * pow(alpha_squared, n + k)) / (f_factorial(n + k, &f_jackson, q) * factorial(k));

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
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	fstream dataFile;
	dataFile.open("quantum_visibility_n.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", t = " << t << '\n';

	long double alpha_squared = pow(alpha, 2);

	for( auto n = 0; n <= INTERVAL_N; n++ )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, &f_jackson, q) * factorial(k);

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
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	fstream dataFile;
	dataFile.open("quantum_visibility_t.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", n = " << n << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double interval = INTERVAL_T / MAX_POINTS;

	for( auto t = interval; t < INTERVAL_T; t += interval )
	{
		long double p = 0.0, c = 0.0;

		for( auto k = 0; k <= SUM_INFTY; k++ )
		{
			long double d = (1 - exp(-t) + LDBL_EPSILON) * alpha_squared;
			long double e = f_factorial(n + k, &f_jackson, q) * factorial(k);

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
	long double domain = 1 / (1 - q + LDBL_EPSILON);

	if( alpha >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

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

	long double alpha_squared = pow(alpha, 2);
	long double normalization = pow(2 * (deformed_exp_jackson(alpha_squared, q, 1.0) + sign * deformed_exp_jackson(-alpha_squared, q, 1.0)), -2);
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
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (sqrt(f_factorial(n + k, &f_jackson, q) * f_factorial(n, &f_jackson, q)) * factorial(n));

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
		long double normalization = pow(2 * (deformed_exp_jackson(alpha_squared, q, 1.0) + sign * deformed_exp_jackson(-alpha_squared, q, 1.0)), -2);

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
				long double e = pow(alpha_squared, n) * exp(-n * t / 2) / (sqrt(f_factorial(n + k, &f_jackson, q) * f_factorial(n, &f_jackson, q)) * factorial(n));

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

void plot_i1(long double q, long double psi = 0.0, long double theta = 0.0)
{
	fstream dataFile;
	dataFile.open("i1.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", psi = " << psi << ", theta = " << theta << '\n';

	long double interval = 1.0 / MAX_POINTS;

	for( complex<long double> alpha = interval; abs(alpha) < 1.0; alpha += interval )
	{
		long double alpha_length_squared = pow(abs(alpha), 2);
		complex<long double> norm_real = two * (deformed_exp_tsallis(alpha_length_squared, q) + exp(imaginary * psi) * deformed_exp_tsallis(alpha_length_squared * exp(imaginary * theta), q));
		long double normalization = pow(real(norm_real), -1);

		complex<long double> sum_a = 0.0, sum_a_dagger = 0.0;
		complex<long double> sum_a_squared = 0.0, sum_a_dagger_squared = 0.0;

		complex<long double> sum_a_dagger_a = 0.0;

		for( auto k = 0; k <= 2000; k++ )
		{
			complex<long double> c = alpha * (one + exp(imaginary * theta) + exp(-imaginary * (psi + k * theta)) + exp(imaginary * (psi + (k + 1.0) * theta)));
			sum_a += pow(alpha_length_squared, k) * c / (f_tsallis(k + 1, q) * factorial(k, &f_tsallis, q));

			complex<long double> d = pow(alpha, 2) * (one + exp(two * imaginary * theta) + exp(-imaginary * (psi + k * theta)) + exp(imaginary * (psi + (k + 2.0) * theta)));
			sum_a_squared += pow(alpha_length_squared, k) * d / (f_tsallis(k + 1, q) * f_tsallis(k + 2, q) * factorial(k, &f_tsallis, q));


			sum_a_dagger_a += two * pow(alpha_length_squared, k) * complex<long double> (1.0 + k) * (one + cos(psi + k * theta)) / factorial(k, &f_tsallis, q);
		}

		sum_a *= normalization;
		sum_a_squared *= normalization;

		sum_a_dagger = conj(sum_a); 
		sum_a_dagger_squared = conj(sum_a_squared);

		sum_a_dagger_a *= normalization;

		complex<long double> x = (sum_a + sum_a_dagger) / two;
		complex<long double> y = (sum_a - sum_a_dagger) / (two * imaginary);

		complex<long double> x_squared = (one + two * sum_a_dagger_a + sum_a_squared + sum_a_dagger_squared) / four;
		complex<long double> y_squared = (one + two * sum_a_dagger_a - sum_a_squared - sum_a_dagger_squared) / four;

		complex<long double> variance_x_squared = x_squared - pow(x, 2);
		complex<long double> variance_y_squared = y_squared - pow(y, 2);
		complex<long double> variance_prod = variance_x_squared * variance_y_squared;

		complex<long double> i1 = sum_a_squared + sum_a_dagger_squared - pow(sum_a, 2) - pow(sum_a_dagger, 2) - two * (sum_a * sum_a_dagger - sum_a_dagger_a);
		complex<long double> i2 = -sum_a_squared - sum_a_dagger_squared + pow(sum_a, 2) + pow(sum_a_dagger, 2) - two * (sum_a * sum_a_dagger - sum_a_dagger_a);

		dataFile << real(alpha) << ' ' << i1 << ' ' << i2 << ' ' << variance_x_squared << ' ' << variance_y_squared << ' ' << variance_prod << '\n';
	}

	dataFile.close();

	cout << "I1... OKAY!" << endl;
}

void plot_wigner(long double q, long double alpha, long double theta = 0.0)
{
	long double domain = 1 / sqrt(q - 1 + LDBL_EPSILON);

	if( abs(alpha) >= domain )
	{
		cout << "alpha is out of range." << endl;
		return;
	}

	fstream dataFile;
	dataFile.open("wigner.dat", ios::out);

	if( !dataFile )
	{
		return;
	}

	dataFile << "## q = " << q << ", alpha = " << alpha << ", theta = " << theta << '\n';

	long double alpha_squared = pow(alpha, 2);
	long double sum = 0.0;

	for( auto l = 0; l <= SUM_INFTY; l++ )
	{
		sum += pow(alpha_squared, l) * (1 + cos(theta + l * M_PI)) / factorial(l, &f_tsallis, q);
	}

	sum *= M_PI;

	long double normalization = 1 / sum;
	long double domainXY = 5.0;
	long double interval = domainXY / 10;
	long double intervalX = interval;
	long double intervalY = interval;

	for( long double x = intervalX - domainXY; x < domainXY; x += intervalX )
	{
		intervalX = interval / (1.0 + 20.0 * exp(-pow(x, 2)));

		for( long double y = intervalY - domainXY; y < domainXY; y += intervalY )
		{
			intervalY = interval / (1.0 + 20.0 * exp(-pow(y, 2)));

			long double phi = atan2(y, x);
			long double beta_squared = pow(x, 2) + pow(y, 2);

			long double sum2 = 0.0, sum2_img = 0.0;

			for( auto k = 0; k <= 10; k++ )
			{
				int sign = IsEvenNumber(k) ? 1 : -1;
				complex<long double> sum3 = 0.0;

				for( auto n = 0; n < 2 * k; n++ )
				{
					int sign2 = IsEvenNumber(n) ? 1 : -1;

					for( auto m = 0; m < 2 * k; m++ )
					{
						int sign3 = IsEvenNumber(m) ? 1 : -1;

						if( n > k && m > k )
						{
							sum3 += complex<long double>(sign2) * complex<long double>(sign3) * pow(alpha, n + m) * (one + complex<long double>(sign2) * exp(imaginary * theta)) * (one + complex<long double>(sign3) * exp(-imaginary * theta))  * factorial(k) * exp(-beta_squared) * pow(beta_squared, (n + m - 2 * k) / 2) * exp(imaginary * complex<long double>((n - m) * phi)) * laguerre(k, n - k, beta_squared) * laguerre(k, m - k, beta_squared) / (factorial(n) * factorial(m) * f_factorial(n, &f_tsallis, q) * f_factorial(m, &f_tsallis, q));
						}
						else if( n > k && m <= k )
						{
							sum3 += complex<long double>(sign2) * pow(alpha, n + m) * (one + complex<long double>(sign2) * exp(imaginary * theta)) * (one + complex<long double>(sign3) * exp(-imaginary * theta)) * exp(-beta_squared) * pow(beta_squared, (n - m) / 2) * exp(imaginary * complex<long double>((n - m) * phi)) * laguerre(k, n - k, beta_squared) * laguerre(m, k - m, beta_squared) / (factorial(n) * f_factorial(n, &f_tsallis, q) * f_factorial(m, &f_tsallis, q));
						}
						else if( n <= k && m > k )
						{
							sum3 += complex<long double>(sign3) * pow(alpha, n + m) * (one + complex<long double>(sign2) * exp(imaginary * theta)) * (one + complex<long double>(sign3) * exp(-imaginary * theta)) * exp(-beta_squared) * pow(beta_squared, (m - n) / 2) * exp(imaginary * complex<long double>((n - m) * phi)) * laguerre(n, k - n, beta_squared) * laguerre(k, m - k, beta_squared) / (factorial(m) * f_factorial(n, &f_tsallis, q) * f_factorial(m, &f_tsallis, q));
						}
						else if( n <= k && m <= k )
						{
							sum3 += pow(alpha, n + m) * (one + complex<long double>(sign2) * exp(imaginary * theta)) * (one + complex<long double>(sign3) * exp(-imaginary * theta)) * exp(-beta_squared) * pow(beta_squared, (2 * k - m - n) / 2) * exp(imaginary * complex<long double>((n - m) * phi)) * laguerre(n, k - n, beta_squared) * laguerre(m, k - m, beta_squared) / (factorial(k) * f_factorial(n, &f_tsallis, q) * f_factorial(m, &f_tsallis, q));
						}
					}

				}

				sum2 += sign * real(sum3);
				sum2_img += sign * imag(sum3);				
			}

			sum2 *= normalization;
			sum2_img *= normalization;

//			cout << sum2 << endl;
			long double ratio = sum2_img * 100 / sum2;
			if( ratio >= 0.5 )
			{
				cout << "Imaginary contribution: " << sum2_img << " % compared to real: " << ratio << endl;
			}

			dataFile << x << ' ' << y << ' ' << sum2 << '\n';
		}

	}

	dataFile.close();

	cout << "Wigner ... OKAY!" << endl;
}

int main()
{
//	plot_variances(0.5);
//	plot_variances2(0.01);
//	plot_photon_statistics(0.95, 1.0, 1.0);
//	plot_photon_statistics2(0.5, 1.0, 0);
//	plot_quantum_visibility(0.5, 1.0, 1.0); 
//	plot_quantum_visibility2(0.5, 1.0, 0);
//	plot_fidelity_t(0.5, 1.8, false);
//	plot_fidelity_alpha(0.9999, 1.0, false);

//	plot_i1(1.0 + LDBL_EPSILON, M_PI / 2, M_PI / 2 );
//	plot_i1(1.5 + LDBL_EPSILON);

	auto start = std::time(NULL);
	plot_wigner(1.5, 0.5, M_PI / 3);
	auto end = std::time(NULL);

	cout << "Started at " << time_stamp(localtime_xp(start), "%T") << endl;
	cout << "Finished at " << time_stamp(localtime_xp(end), "%T") << endl;

	auto duration = difftime(end, start) / 60.0;
	cout << "Duration: " << duration << " minutes" << endl;

	system("pause");
	return 0;
}