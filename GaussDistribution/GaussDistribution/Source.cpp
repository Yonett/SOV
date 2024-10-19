#define _USE_MATH_DEFINES

#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "spec_func.hpp"

using namespace std;

const int N = 1e7;
const double EPS = 0.15;
const double TETTA = 2.0;
const double LAMBDA = 0.3;
const double ALPHA = 0.3;
const double DELTA = 5;

double CalcPureValue(double P, double v, double a)
{
	//printf("Calculating value..\n");
	double r1, r2, r3, r4, r5;
	double x1, x2, z;

	double one_seven = 1.0 / 7.0;
	double two_seven = 2.0 / 7.0;

	double seven_one_seven = pow(7, -one_seven);
	double seven_two_seven = pow(7, two_seven);

	double five_fourteen = 5.0 / 14.0;

	do {
		r1 = (double)rand() / RAND_MAX;
		r3 = (double)rand() / RAND_MAX;
	} while (r1 == 0 || r3 == 0 || r1 == 1 || r3 == 1);

	if (r1 <= P)
	{
		while (true)
		{
			do {
				r4 = (double)rand() / RAND_MAX;
				r5 = (double)rand() / RAND_MAX;
				r2 = (double)rand() / RAND_MAX;
			} while (r4 == 0 || r5 == 0 || r2 == 0 || r4 == 1 || r5 == 1 || r2 == 1);

			z = sqrt(-2 * log(r4)) * cos(2 * M_PI * r5);
			x1 = seven_one_seven * z;
			if (abs(x1) <= v) {
				//printf("Not Done\n");
				if (log(r2) <= -pow(abs(x1), 7) + seven_two_seven * x1 * x1 * 0.5 - five_fourteen)
				{
					//printf("Done\n");
					return x1;
				}
			}
		}
	}

	x2 = v - (log(r3) / a);

	//printf("Done\n");

	if (r1 < (1.0 + P) * 0.5)
		return x2;
	return -x2;
	
}

double CalcDistributionFunc(double x, double v, double v7, double K, double a)
{
	double result = 0;
	if (abs(x) > v)
		result = exp(-a * (abs(x) - v) - v7);
	else
		result = exp(pow(-abs(x), 7));
	return result / K;
}

double CalcMLEFunc(double* values, double v, double v7, double K, double a, double tetta)
{
	double result = 0;

	for (int i = 0; i < N; i++)
		result += -log(CalcDistributionFunc((values[i] - tetta) / LAMBDA, v, v7, K, a));

	return result / N;
}

double CalcRadicalFunc(double* values, double v, double v7, double K, double a, double tetta)
{
	double result = 0;

	for (int i = 0; i < N; i++)
		result += -pow(CalcDistributionFunc((values[i] - tetta) / LAMBDA, v, v7, K, a), DELTA)
		/ //----------------------------------------------------------------------------------
				   pow(CalcDistributionFunc((values[i] - 0)     / LAMBDA, v, v7, K, a), DELTA);

	return result / N;
}

double CalcMLE(double* values, double v, double v7, double K, double a, double tetta, double eps)
{
	double gld = 1.6180339887498948482;
	double n = -4;
	double m = 20;

	double x1;
	double x2;

	while (abs(m - n) > eps)
	{
		x1 = m - (m - n) / gld;
		x2 = n + (m - n) / gld;

		if (CalcMLEFunc(values, v, v7, K, a, x1) < CalcMLEFunc(values, v, v7, K, a, x2))
			m = x2;
		else
			n = x1;
	}

	return (n + m) / 2;
}

double CalcRadical(double* values, double v, double v7, double K, double a, double tetta, double eps)
{
	double gld = 1.6180339887498948482;
	double n = -4;
	double m = 20;

	double x1;
	double x2;

	while (abs(m - n) > eps)
	{
		x1 = m - (m - n) / gld;
		x2 = n + (m - n) / gld;

		if (CalcRadicalFunc(values, v, v7, K, a, x1) < CalcRadicalFunc(values, v, v7, K, a, x2))
			m = x2;
		else
			n = x1;
	}

	return (n + m) / 2;
}

double CalcArithmeticMean(double *values)
{
	double result = 0;
	for (int i = 0; i < N; i++)
		result += values[i];
	return result / N;
}

double CalcDispersion(double* values, double y_)
{
	double result = 0;
	for (int i = 0; i < N; i++)
		result += (values[i] - y_) * (values[i] - y_);
	return result / N;
}

double CalcSampleMedian(double* values)
{
	double result = 0;
	if (N % 2 == 0)
		result = (values[N / 2] + values[(N / 2) - 1]) / 2.0;
	else
		result = values[N / 2];
	return result;
}

double CalcAsymmetryCoefficient(double* values, double y_, double D)
{
	double result = 0;
	for (int i = 0; i < N; i++)
		result += (values[i] - y_) * (values[i] - y_) * (values[i] - y_);
	return result / (N * pow(D, 1.5));
}

double CalcKurtosisCoefficient(double* values, double y_, double D)
{
	double result = 0.0;
	for (int i = 0; i < N; i++)
		result += (values[i] - y_) * (values[i] - y_) * (values[i] - y_) * (values[i] - y_);
	result /= (N * D * D);
	return result - 3.0;
}

double CalcTrimmedMean(double* values)
{
	double result = 0;
	int k = N * ALPHA;
	for (int i = k; i < N - k; i++)
		result += values[i];
	return result / (N - 2 * k);
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

	srand(time(nullptr));

	double* pure_values = new double[N];
	double* dirt_values = new double[N];
	double r, r1, r2;

	//ofstream file;
	//
	//file.open("data_pure.csv");


	for (int i = 0; i < N; i++)
	{
		pure_values[i] = CalcPureValue(P, v, a);
	//	file << pure_values[i] << ';' << endl;
	}

	//file.close();

	//file.open("data_dirt.csv");

	/*
	for (int i = 0; i < N; i++)
	{
		dirt_values[i] = CalcPureValue(P, v, a);

		r = (double)rand() / RAND_MAX;
		if (r < EPS)
		{
			dirt_values[i] += TETTA;
			dirt_values[i] *= LAMBDA;
		}

		file << dirt_values[i] << ';' << endl;
	}

	file.close();
	*/

	//double values[N];

	//values[0] = 3;
	//values[1] = 1;
	//values[2] = 5;
	//values[3] = 2;
	//values[4] = 4;
	//values[5] = 6;
	//values[6] = 8;
	//values[7] = 7;

	qsort(pure_values, N, sizeof(double),
		[](const void* a, const void* b)
		{
			const double* x = (double*)a;
			const double* y = (double*)b;

			if (*x > *y)
				return 1;
			else if (*x < *y)
				return -1;

			return 0;
		}
	);

	double y_, D, ac, kc, sm, tm, MLE;

	y_ = CalcArithmeticMean(pure_values);

	printf("Arithmetic Mean: %e\n", y_);

	D = CalcDispersion(pure_values, y_);

	printf("Dispersion: %e\n", D);

	ac = CalcAsymmetryCoefficient(pure_values, y_, D);

	printf("Assymetry Coefficient: %e\n", ac);

	kc = CalcKurtosisCoefficient(pure_values, y_, D);

	printf("Kurtosis Coefficient: %e\n", kc);

	sm = CalcSampleMedian(pure_values);

	printf("Sample Median: %e\n", sm);

	//tm = CalcTrimmedMean(values);

	//printf("Trimmed Mean: %e\n", tm);

	//MLE = CalcMLE(values, v, v7, K, a, 2, 1e-5);

	//printf("MLE: %e\n", MLE);

	return 0;
}