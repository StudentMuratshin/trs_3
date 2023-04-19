// trs_3.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

constexpr auto pi = 3.14159265358979323846;

double True_U(double x, double y)
{
	return sin(pi * x / 2.) * sin(pi * y) * sin(pi * y);
}

double F_right(double x, double y)
{
	return pi * pi / 4. * sin(x * pi / 2.) * (17 * cos(pi * y) * cos(pi * y) - 9.);
}

struct EllipticPDEMatrixGen {
	double alpha;
	double beta;
	double gamma;

	int N;
	int M;
	EllipticPDEMatrixGen(double alpha, double beta, double gamma, int N, int M) : alpha(alpha), beta(beta), gamma(gamma), N(N), M(M) {}

	const double& operator()(int i, int j) const {
		if (i == j) {
			return alpha;
		}
		else if (i - j == -1 && (i + 1) % M != 0) {
			return beta;
		}
		else if (i - j == 1 && (i % M != 0)) {
			return beta;
		}
		else if (abs(i - j) == M) {
			return gamma;
		}
		else {
			return 0.;
		}
	}
};

void lu_decompostion(const unsigned int n, EllipticPDEMatrixGen A, vector<double>& L, vector<double>& U) {
	for (int i = 0; i < n; ++i) {
		L[i * n] = A(i, 0);
		U[i] = A(i, 0) / L[0];
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = A(i, j);

			for (int k = 0; k < i; k++)
			{
				U[i * n + j] -= L[i * n + k] * U[k * n + j];
			}

			L[j * n + i] = A(i, j);

			for (int k = 0; k < i; k++)
			{
				L[j * n + i] -= L[j * n + k] * U[k * n + i];
			}

			L[j * n + i] /= U[i * n + i];
		}
	}
}

void backward_up(const unsigned int n, vector<double> U, vector<double> b, vector<double>& x) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= U[j + i * n] * x[j];
		}
		x[i] /= U[i + i * n];
	}
}

void backward_low(const unsigned int n, vector<double>  L, vector<double>  b, vector<double>& y) {
	for (int i = 0; i < n; ++i) {
		y[i] = b[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= L[i * n + j] * y[j];
		}
		y[i] /= L[i + i * n];
	}
}


void solve_lu(const unsigned int n, vector<double>  L, vector<double>  U, vector<double>  b, vector<double>& x) {
	vector<double> y(n);
	backward_low(n, L, b, y);
	backward_up(n, U, y, x);
}
