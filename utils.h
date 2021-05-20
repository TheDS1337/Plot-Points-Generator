#pragma once

#define M_PI           3.1415926535897932384

constexpr long double factorial(int n)
{
	if( n > 65 )
	{
		return exp(n * log(n) - n + (log(n) / 2) - log(2 * M_PI));
	}

	return n <= 0 ? 1 : n * factorial(n - 1);
}

constexpr long double permutation(int n, int p, bool inversed = false)
{
	return inversed ? factorial(p) / factorial(n) : factorial(n) / factorial(p);
}

constexpr long double combination(int n, int p, bool inversed = false)
{
	return inversed ? factorial(p) / permutation(n, p) : permutation(n, p) / factorial(p);
}

constexpr long double deformed_factorial(int n, long double (*pF) (int, float), float q)
{
	return n <= 0 ? 1 : n * pow(pF(n, q), 2) * deformed_factorial(n - 1, pF, q);
}

bool IsEvenNumber(int n)
{
	return n & 1;
}

// From Quake 3
double Q_rsqrt(double number)
{
	long long i;
	double x2, y;
	const float threehalfs = 1.5;

	x2 = number * 0.5;
	y = number;
	i = *(long long*) &y;                       // evil floating point bit level hacking
	i = 0x5fe6eb50c7b537a9 - (i >> 1);               // what the fuck? 
	y = *(double*) &i;
	y = y * (threehalfs - (x2 * y * y));   // 1st iteration
//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}