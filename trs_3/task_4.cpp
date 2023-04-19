#include <iostream>
#include <iomanip>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

constexpr auto pi = 3.14159265358979323846;

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

void lu_decompostion(const unsigned int n, EllipticPDEMatrixGen A, SparseMatrix<double>& L, SparseMatrix<double>& U) {
	for (int i = 0; i < n; ++i) {
		L.insert(i, 0) = A(i, 0);
		U.insert(0, i) = A(0, i) / L.coeff(0, 0);
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U.insert(i, j) = A(i, j);

			for (int k = 0; k < i; k++)
			{
				U.coeffRef(i, j) = U.coeff(i,j) - L.coeff(i, k) * U.coeff(k, j);
			}

			L.insert(j, i) = A(i, j);

			for (int k = 0; k < i; k++)
			{
				L.coeffRef(j, i) = L.coeff(j,i) - L.coeff(j, k) * U.coeff(k, i);
			}

			L.insert(j, i) /= U.coeff(i, i);
		}
	}
}

void backward_up(const unsigned int n, SparseMatrix<double>& U, vector<double> b, vector<double>& x) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= U.coeff(j, i) * x[j];
		}
		x[i] /= U.coeff(i, i);
	}
}

void backward_low(const unsigned int n, SparseMatrix<double>& L, vector<double> b, vector<double>& y) {
	for (int i = 0; i < n; ++i) {
		y[i] = b[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= L.coeff(i, j) * y[j];
		}
		y[i] /= L.coeff(i, i);
	}
}


void solve_lu(const unsigned int n, SparseMatrix<double>& L, SparseMatrix<double>& U, vector<double> b, vector<double>& x) {
	vector<double> y(n);
	backward_low(n, L, b, y);
	backward_up(n, U, y, x);
}

int main() {

	const double lx = 2., ly = 1.;

	auto f = [](double x, double y) {
		return  -pi * pi / 4. * sin(x * pi / 2.) * (17 * cos(pi * y) * cos(pi * y) - 9.);
	};

	auto phi1 = [](double x = 0, double y = 0) {
		return 0;
	};

	auto phi2 = [](double x = 0, double y = 0) {
		return 0;
	};

	auto phi3 = [](double x = 0, double y = 0) {
		return 0;
	};

	auto phi4 = [](double x = 0, double y = 0) {
		return 0;
	};


	auto exact_sol = [lx, ly](double x, double y) {
		return sin(pi * x / 2.) * sin(pi * y) * sin(pi * y);
	};


	auto getF = [phi1, phi2, phi3, phi4, f](int N, int M, double hx, double hy) {
		vector<double> F(N * M);
		{
			int i, j;
			// F äëÿ ëåâîé ãðàíèöû
			i = 0;
			j = 0;
			F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi3() / hy / hy;
			for (j = 1; j < N - 1; ++j) {
				F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx;
			}
			F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi4() / hy / hy;

			// F äëÿ âåðõíåé ãðàíèöû
			for (i = 1; i < M - 1; ++i) {
				F[j * M + i] = -f(i * hx, j * hy) - phi4() / hy / hy;
			}

			F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi4() / hy / hy;

			// F äëÿ ïðàâîé ãðàíèöû
			for (j = N - 2; j > 0; --j) {
				F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx;
			}
			F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi3() / hy / hy;

			// F  äëÿ íèæíåé ãðàíèöû
			for (i = M - 2; i > 0; --i) {
				F[j * M + i] = -f(i * hx, j * hy) - phi3() / hy / hy;
			}

			// âíóòðåííÿÿ îáëàñòü
			for (j = 1; j < N - 1; ++j) {
				for (i = 1; i < M - 1; ++i) {
					F[j * M + i] = -f(i * hx, j * hy);
				}
			}
		}
		return F;

	};

	vector<int> N = { 4,20,30,40,50,60,70,80,90,100 };
	for (auto n : N)
	{
		const int N = n;
		const int M = n;
		double hx = static_cast<double>(lx) / N;
		double hy = static_cast<double>(ly) / M;

		vector<double> U(N * M, 0.); // Íà÷àëüíîå ïðèáëèæåíèå íóëåâîé âåêòîð

		double alpha = -2. * (1. / hx / hx + 1. / hy / hy);
		double beta = 1. / hx / hx;
		double gamma = 1. / hy / hy;

		vector<double> F = getF(N, M, hx, hy);

		//SparseMatrix<double> A((N - 1) * (M - 1), (N - 1)* (M - 1));

		//for (int i = 0; i < (N - 1) * (M - 2); i++)
		//{
		//	A.insert(i,i) = alpha;
		//	A.insert(i, i + 1) = beta;
		//	A.insert(i + 1, i) = beta;

		//}
		EllipticPDEMatrixGen A(alpha, beta, gamma, N, M);
		SparseMatrix<double> Lower(N * M, N * M);
		SparseMatrix<double> Upper(N * M, N * M);

		lu_decompostion(N * M, A, Lower, Upper);

		solve_lu(N * M, Lower, Upper, F, U);

		double error = 0.;

		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < M; ++i) {
				double x = i * hx;
				double y = j * hy;
				error = std::max(error, abs(U[j * M + i] - sin(pi * x / 2.) * sin(pi * y) * sin(pi * y)));
				//cout << x << " " << y << U[j * M + i] <<endl;
			}
		}

		cout << "m, n: " << N << "   Diff: " << setprecision(8) << error << endl;

	}
    return 0;
}
