#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <ctime>

#define M_PI					3.1415926535897932384

#define SUM_INFTY				50					// This is where the fractional series are truncated (65 = magic number for long long unsigned int, 171 for long doubles)
#define MAX_POINTS				50					// More points = more smooth = more processing (overleaf can't handle much)

#define INTERVAL_Q				1.0
#define INTERVAL_ALPHA			1.5
#define INTERVAL_T				1.0
#define INTERVAL_N				50

constexpr long double factorial(int n, long double (*pF) (int, long double) = nullptr, long double q = 1.0)
{
	long double p = 1.0;

	if( pF && n > 0 )
	{
		long double f = pF(n, q);

		if( f != 0.0 && f != 1.0 )
		{
			p = pow(f, 2);
		}
	}

	return n <= 0 ? 1 : n * p * factorial(n - 1, pF, q);
}

constexpr long double permutation(int n, int p, long double (*pF) (int, long double) = nullptr, long double q = 1.0, bool inversed = false)
{
	return inversed ? factorial(n - p, pF, q) / factorial(n, pF, q) : factorial(n, pF, q) / factorial(n - p, pF, q);
}

constexpr long double combination(int n, int p, long double (*pF) (int, long double) = nullptr, long double q = 1.0, bool inversed = false)
{
	return inversed ? factorial(p, pF, q) / permutation(n, p, pF, q) : permutation(n, p, pF, q) / factorial(p, pF, q);
}

constexpr long double f_factorial(int n, long double (*pF) (int, long double), long double q, bool squared = true)
{
	long double p = 1.0;

	if( n > 0 )
	{
		long double f = pF(n, q);

		if( f != 0.0 && f != 1.0 )
		{
			if( squared )
			{
				p = pow(f, 2);
			}
			else
			{
				p = f;
			}
		}
	}

	return n <= 0 ? 1 : p * f_factorial(n - 1, pF, q, squared);
}

bool IsEvenNumber(int n)
{
	return n % 2 == 0;
}

long double laguerre(int n, int m, long double x)
{
	if( n == 0 )
	{
		return 1.0;
	}

	long double value = 0.0;

	for( auto k = 0; k <= n; k++ )
	{
		int sign = IsEvenNumber(k) ? 1 : -1;
		value += sign * combination(n + m, n - k) * pow(x, k) / factorial(k);
	}

	return value;
}

enum LevinTransform
{
	LevinTransform_TypeU,
	LevinTransform_TypeT,
	LevinTransform_TypeBatouli
};

long double levin_transformation(int k, long double (*s) (int, long double), long double x, LevinTransform type = LevinTransform_TypeU, int n = 0, int b = 1)
{
	long double A = 0.0, B = 0.0;

	switch( type )
	{
	case LevinTransform_TypeU:
		for( auto j = 0; j <= k; j++ )
		{
			int sign_j = IsEvenNumber(j) ? 1 : -1;
			int p = j + n;

			long double partial_sum = 0.0;

			for( auto m = 0; m <= p; m++ )
			{
				partial_sum += s(m, x);
			}

			long double d = sign_j * combination(k, j) * pow(j + n + b, k - 1) / ((n + b) * s(j + n, x));
			A += d * partial_sum;
			B += d;
		}
		break;

	case LevinTransform_TypeT:
		for( auto j = 0; j <= k; j++ )
		{
			int sign_j = IsEvenNumber(j) ? 1 : -1;
			int p = j + n;

			long double partial_sum = 0.0;

			for( auto m = 0; m <= p; m++ )
			{
				partial_sum += s(m, x);
			}

			long double d = sign_j * combination(k, j) * pow(j + n + b, k - 1) / s(j + n, x);
			A += d * partial_sum;
			B += d;
		}
		break;

	case LevinTransform_TypeBatouli:
		for( auto j = 0; j <= k; j++ )
		{
			int sign_j = IsEvenNumber(j) ? 1 : -1;
			int p = j + n;

			long double partial_sum = 0.0;

			for( auto m = 0; m <= p; m++ )
			{
				partial_sum += s(m, x);
			}

			long double d = sign_j * combination(k, j) * pow((1 + j) / (1 + k), k - 2) / s(j + 1, x);
			A += d * partial_sum;
			B += d;
		}
		break;
	}

	return A / B;
}

inline std::tm localtime_xp(std::time_t timer)
{
	std::tm bt{};
#if defined(__unix__)
	localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
	localtime_s(&bt, &timer);
#else
	static std::mutex mtx;
	std::lock_guard<std::mutex> lock(mtx);
	bt = *std::localtime(&timer);
#endif
	return bt;
}

// default = "YYYY-MM-DD HH:MM:SS"
inline std::string time_stamp(std::tm bt, const std::string& fmt = "%F %T")
{
	char buf[64];
	return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
}

#include "deformations.h"