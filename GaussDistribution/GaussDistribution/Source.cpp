#define _USE_MATH_DEFINES
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "spec_func.hpp"

using namespace std;

double CalcRandValue(double P, double v, double a)
{
	double r1, r2, r3, r4, r5;
	double x1, x2, z;

	double one_seven = 1.0 / 7.0;
	double two_seven = 2.0 / 7.0;

	double seven_one_seven = pow(7, -one_seven);
	double seven_two_seven = pow(7, two_seven);

	double five_fourteen = 5.0 / 14.0;

	r1 = (double)rand() / RAND_MAX;
	r2 = (double)rand() / RAND_MAX;
	r3 = (double)rand() / RAND_MAX;

	if (r1 <= P)
	{
		while (true)
		{
			r4 = (double)rand() / RAND_MAX;
			r5 = (double)rand() / RAND_MAX;
			z = sqrt(-2 * log(r4)) * cos(2 * M_PI * r5);
			x1 = seven_one_seven * z;
			if (abs(x1) <= v)
				if (log(r2) <= -pow(abs(x1), 7) + seven_two_seven * x1 * x1 * 0.5 - five_fourteen)
					return x1;
		}
	}

	x2 = v - (log(r3) / a);

	if (r1 < (1.0 + P) * 0.5)
		return x2;
	return -x2;
}

int main()
{
	double v = 0.9;
	double v2 = v * v;
	double v3 = v2 * v;
	double v4 = v3 * v;
	double v5 = v4 * v;
	double v6 = v5 * v;
	double v7 = v6 * v;

	double a = 7.0 * v6;
	double a2 = a * a;
	double a3 = a2 * a;
	double a4 = a3 * a;
	double a5 = a4 * a;

	double one_seven = 1.0 / 7.0;
	double two_seven = 2.0 / 7.0;
	double three_seven = 3.0 / 7.0;
	double five_seven = 5.0 / 7.0;
	double five_fourteen = 5.0 / 14.0;

	double seven_one_seven = pow(7, -one_seven);
	double seven_two_seven = pow(7, two_seven);

	double ev7 = exp(-v7);

	double K = two_seven * (igamma(one_seven, v7) + ev7 / v6);

	double variance = (2.0 / K) * ((igamma(three_seven, v7) * one_seven) + ev7 * ((2 / a3) +
		                                                                          (2 * v / a2) +
		                                                                          (v2 / a)));

	double igamma2 = (2.0 / K) * ((igamma(five_seven, v7) * one_seven) + ev7 * ((24 / a5) +
																		        (24 * v / a4) +
																		        (12 * v2 / a3) +
																		        (4 * v3 / a2) +
																		        (v4 / a)));
	igamma2 /= variance * variance;
	igamma2 -= 3.0;

	double P = 2.0 * igamma(one_seven, v7) / (7.0 * K);

	printf("v: %e\n", v);
	printf("variance: %e\n", variance);
	printf("igamma2: %e\n", igamma2);
	printf("P: %e\n", P);
	printf("K: %e\n", K);

	double result = 0;

	ofstream file;
	file.open("data.csv");

	srand(time(nullptr));

	for (int i = 0; i < 100000; i++)
	{
		
		result = CalcRandValue(P, v, a);
		file << result << ';' << endl;
		printf("%i\n", i);
	}

	file.close();

	return 0;
}