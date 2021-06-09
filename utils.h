#pragma once

#define M_PI           3.1415926535897932384

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

constexpr long double f_factorial(int n, long double (*pF) (int, long double), long double q)
{
	long double p = 1.0;

	if( n > 0 )
	{
		long double f = pF(n, q);

		if( f != 0.0 && f != 1.0 )
		{
			p = pow(f, 2);
		}
	}

	return n <= 0 ? 1 : p * f_factorial(n - 1, pF, q);
}

bool IsEvenNumber(int n)
{
	return n % 2 == 0;
}

long double laguerre(int n, int m, long double x)
{
	long double value = 0.0;

	for( auto k = 0; k <= n; k++ )
	{
		int sign = IsEvenNumber(k) ? 1 : -1;
		value += sign * combination(n + m, n - k) * pow(x, k) / factorial(k);
	}

	return value;
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