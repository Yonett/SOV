#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include "spec_func.hpp"

using namespace std;

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

	// r-i
	srand(time(nullptr));
	double r1, r2, r3;
	r1 = (double)rand() / RAND_MAX;
	r2 = (double)rand() / RAND_MAX;
	r3 = (double)rand() / RAND_MAX;

	printf("v: %e\n", v);
	printf("variance: %e\n", variance);
	printf("igamma2: %e\n", igamma2);
	printf("P: %e\n", P);
	printf("K: %e\n", K);

	printf("r1: %e\nr2: %e\nr3: %e\n", r1, r2, r3);



	return 0;
}