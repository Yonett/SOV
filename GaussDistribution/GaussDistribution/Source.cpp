#define _USE_MATH_DEFINES
//#define _WRITE

#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "spec_func.hpp"

using namespace std;

const int N = 3e5;
const double EPS = 0.30;                 // clogging probability

const double PURE_TETTA = 0.0;           // shift for pure distribution
const double PURE_LAMBDA = 1;          // scale for pure distribution

const double SYM_TETTA = PURE_TETTA;     // shift for symmetrical clogging distribution
const double SYM_LAMBDA = 3.0;           // scale for symmetrical clogging distribution

const double ASYM_TETTA = 1.5;           // shift for asymmetrical clogging distribution
const double ASYM_LAMBDA = PURE_LAMBDA + 0.25;  // scale for asymmetrical clogging distribution

double ALPHA = 0.0;                      // trimmed mean parameter
double DELTA = 0.0;                      // radical function power coef

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
	{
		double tmp = -a * (abs(x) - v) - v7;
		result = exp(tmp);
	}
	else
		result = exp(-pow(abs(x), 7));
	return result / K;
}

double CalcMLEFunc(double* values, double v, double v7, double K, double a, double lambda, double tetta)
{
	double result = 0;

	for (int i = 0; i < N; i++)
		result += -log(CalcDistributionFunc((values[i] - tetta) / lambda, v, v7, K, a));

	return result / N;
}

double CalcRadicalFunc(double* values, double v, double v7, double K, double a, double lambda, double tetta)
{
	double result = 0;

	double divider = pow(CalcDistributionFunc(0, v, v7, K, a), DELTA);

	for (int i = 0; i < N; i++)
		result += -pow(CalcDistributionFunc((values[i] - tetta) / lambda, v, v7, K, a), DELTA) / divider;

	return result / N;
}

double CalcMLE(double* values, double v, double v7, double K, double a, double lambda, double eps)
{
	double gld = 1.6180339887498948482;
	double n = -5;
	double m =  5;

	double x1;
	double x2;

	while (abs(m - n) > eps)
	{
		x1 = m - (m - n) / gld;
		x2 = n + (m - n) / gld;

		if (CalcMLEFunc(values, v, v7, K, a, lambda, x1) < CalcMLEFunc(values, v, v7, K, a, lambda, x2))
			m = x2;
		else
			n = x1;
	}

	return (n + m) / 2;
}

double CalcRadical(double* values, double v, double v7, double K, double a, double lambda, double eps)
{
	double gld = 1.6180339887498948482;
	double n = -5;
	double m =  5;

	double x1;
	double x2;

	while (abs(m - n) > eps)
	{
		x1 = m - (m - n) / gld;
		x2 = n + (m - n) / gld;

		if (CalcRadicalFunc(values, v, v7, K, a, lambda, x1) < CalcRadicalFunc(values, v, v7, K, a, lambda, x2))
			m = x2;
		else
			n = x1;
	}

	return (n + m) / 2;
}

double CalcMLEIF()
{
	double result = 0.0;

	return result;
}

double CalcRadicalIF()
{
	double result = 0.0;

	return result;
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
#pragma region Constants

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

	printf("Constants:\n");
	printf("v: %e\n", v);
	printf("variance: %e\n", variance);
	printf("igamma2: %e\n", igamma2);
	printf("P: %e\n", P);
	printf("K: %e\n\n", K);

#pragma endregion

	srand(time(nullptr));

	double* pure_values = new double[N];

	double* sym_clog_values = new double[N];
	double* sym_dirt_values = new double[N];

	double* asym_clog_values = new double[N];
	double* asym_dirt_values = new double[N];

	double r = 0.0;

	for (int i = 0; i < N; i++)
	{
		pure_values[i] = (CalcPureValue(P, v, a) + PURE_TETTA) / PURE_LAMBDA;
	}
	
	for (int i = 0; i < N; i++)
	{
		sym_clog_values[i] = (CalcPureValue(P, v, a) + SYM_TETTA) / SYM_LAMBDA;

		r = (double)rand() / RAND_MAX;
		if (r < EPS)
			sym_dirt_values[i] = sym_clog_values[i];
		else
			sym_dirt_values[i] = pure_values[i];
	}

	for (int i = 0; i < N; i++)
	{
		asym_clog_values[i] = (CalcPureValue(P, v, a) + ASYM_TETTA) / ASYM_LAMBDA;

		r = (double)rand() / RAND_MAX;
		if (r < EPS)
			asym_dirt_values[i] = asym_clog_values[i];
		else
			asym_dirt_values[i] = pure_values[i];
	}	

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

	qsort(sym_dirt_values, N, sizeof(double),
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

	qsort(asym_dirt_values, N, sizeof(double),
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

	double y_, D, ac, kc, sm, tm, MLE, RAD;

#pragma region Pure

	printf("Pure distribution (shift: %e, scale %e):\n\n", PURE_TETTA, PURE_LAMBDA);

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

	ALPHA = 0.05;
	tm = CalcTrimmedMean(pure_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.10;
	tm = CalcTrimmedMean(pure_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.15;
	tm = CalcTrimmedMean(pure_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	MLE = CalcMLE(pure_values, v, v7, K, a, PURE_LAMBDA, 1e-5);

	printf("MLE: %e\n", MLE);

	DELTA = 0.1;
	RAD = CalcRadical(pure_values, v, v7, K, a, PURE_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 0.5;
	RAD = CalcRadical(pure_values, v, v7, K, a, PURE_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 1.0;
	RAD = CalcRadical(pure_values, v, v7, K, a, PURE_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n\n\n", DELTA, RAD);

#pragma endregion


#pragma region Clean Sym Dirt

	printf("Symmetrical dirt distribution (shift: %e, scale %e):\n\n", SYM_TETTA, SYM_LAMBDA);

	y_ = CalcArithmeticMean(sym_clog_values);

	printf("Arithmetic Mean: %e\n", y_);

	D = CalcDispersion(sym_clog_values, y_);

	printf("Dispersion: %e\n", D);

	ac = CalcAsymmetryCoefficient(sym_clog_values, y_, D);

	printf("Assymetry Coefficient: %e\n", ac);

	kc = CalcKurtosisCoefficient(sym_clog_values, y_, D);

	printf("Kurtosis Coefficient: %e\n", kc);

	sm = CalcSampleMedian(sym_clog_values);

	printf("Sample Median: %e\n", sm);

	ALPHA = 0.05;
	tm = CalcTrimmedMean(sym_clog_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.10;
	tm = CalcTrimmedMean(sym_clog_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.15;
	tm = CalcTrimmedMean(sym_clog_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	MLE = CalcMLE(sym_clog_values, v, v7, K, a, SYM_LAMBDA, 1e-5);

	printf("MLE: %e\n", MLE);

	DELTA = 0.1;
	RAD = CalcRadical(sym_clog_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 0.5;
	RAD = CalcRadical(sym_clog_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 1.0;
	RAD = CalcRadical(sym_clog_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n\n\n", DELTA, RAD);

#pragma endregion



#pragma region Sym Dirt

	printf("Symmetrical dirt distribution (shift: %e, scale %e):\n\n", SYM_TETTA, SYM_LAMBDA);

	y_ = CalcArithmeticMean(sym_dirt_values);

	printf("Arithmetic Mean: %e\n", y_);

	D = CalcDispersion(sym_dirt_values, y_);

	printf("Dispersion: %e\n", D);

	ac = CalcAsymmetryCoefficient(sym_dirt_values, y_, D);

	printf("Assymetry Coefficient: %e\n", ac);

	kc = CalcKurtosisCoefficient(sym_dirt_values, y_, D);

	printf("Kurtosis Coefficient: %e\n", kc);

	sm = CalcSampleMedian(sym_dirt_values);

	printf("Sample Median: %e\n", sm);

	ALPHA = 0.05;
	tm = CalcTrimmedMean(sym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.10;
	tm = CalcTrimmedMean(sym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.15;
	tm = CalcTrimmedMean(sym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	MLE = CalcMLE(sym_dirt_values, v, v7, K, a, SYM_LAMBDA, 1e-5);

	printf("MLE: %e\n", MLE);

	DELTA = 0.1;
	RAD = CalcRadical(sym_dirt_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 0.5;
	RAD = CalcRadical(sym_dirt_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 1.0;
	RAD = CalcRadical(sym_dirt_values, v, v7, K, a, SYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n\n\n", DELTA, RAD);

#pragma endregion

#pragma region Asym Dirt

	printf("Asymmetrical dirt distribution (shift: %e, scale %e):\n\n", ASYM_TETTA, ASYM_LAMBDA);

	y_ = CalcArithmeticMean(asym_dirt_values);

	printf("Arithmetic Mean: %e\n", y_);

	D = CalcDispersion(asym_dirt_values, y_);

	printf("Dispersion: %e\n", D);

	ac = CalcAsymmetryCoefficient(asym_dirt_values, y_, D);

	printf("Assymetry Coefficient: %e\n", ac);

	kc = CalcKurtosisCoefficient(asym_dirt_values, y_, D);

	printf("Kurtosis Coefficient: %e\n", kc);

	sm = CalcSampleMedian(asym_dirt_values);

	printf("Sample Median: %e\n", sm);

	ALPHA = 0.05;
	tm = CalcTrimmedMean(asym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.10;
	tm = CalcTrimmedMean(asym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	ALPHA = 0.15;
	tm = CalcTrimmedMean(asym_dirt_values);
	printf("Trimmed Mean (ALPHA = %e): %e\n", ALPHA, tm);

	MLE = CalcMLE(asym_dirt_values, v, v7, K, a, ASYM_LAMBDA, 1e-5);

	printf("MLE: %e\n", MLE);

	DELTA = 0.1;
	RAD = CalcRadical(asym_dirt_values, v, v7, K, a, ASYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 0.5;
	RAD = CalcRadical(asym_dirt_values, v, v7, K, a, ASYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n", DELTA, RAD);

	DELTA = 1.0;
	RAD = CalcRadical(asym_dirt_values, v, v7, K, a, ASYM_LAMBDA, 1e-5);
	printf("Radical (DELTA = %e): %e\n\n\n", DELTA, RAD);

#pragma endregion

#ifdef _WRITE
#pragma region WriteToFile

	ofstream file;

	file.open("data_pure.csv");

	for (int i = 0; i < N; i++)
		file << pure_values[i] << ' ' << endl;

	file.close();

	file.open("data_sym_clog.csv");

	for (int i = 0; i < N; i++)
		file << sym_clog_values[i] << ' ' << endl;

	file.close();

	file.open("data_sym_dirt.csv");

	for (int i = 0; i < N; i++)
		file << sym_dirt_values[i] << ' ' << endl;

	file.close();

	file.open("data_asym_clog.csv");

	for (int i = 0; i < N; i++)
		file << asym_clog_values[i] << ' ' << endl;

	file.close();

	file.open("data_asym_dirt.csv");

	for (int i = 0; i < N; i++)
		file << asym_dirt_values[i] << ' ' << endl;

	file.close();

#pragma endregion
#endif

	return 0;
}