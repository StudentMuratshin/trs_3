#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <span>
#include <fstream>
#include <format>
#include <Eigen/Sparse>
#include <Eigen/Dense>

constexpr auto pi = 3.14159265358979323846;

using namespace Eigen;
using namespace std;


int main() {
	
	vector<int> N_arr = { 25 };
	const double lx = 1.;
	const double ly = 2.;

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

	auto b = [](double x, double y) {
		return sin(pi*x);
	};

	auto by = [](double x, double y) {
		return 0;
	};

	auto a = [](double x, double y) {
		return 2+sin(pi*y);
	};

	auto ax = [](double x, double y) {
		return 0;
	};

	auto c = [](double x, double y) {
		return x;
	};

	auto getF = [phi1, phi2, phi3, phi4, f](int N, int M, double hx, double hy) {
		VectorXd F(N * M);
		{
			int i, j;
			// F для левой границы
			i = 0;
			j = 0;
			F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi3() / hy / hy;
			for (j = 1; j < N - 1; ++j) {
				F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx;
			}
			F[j * M + i] = -f(i * hx, j * hy) - phi1() / hx / hx - phi4() / hy / hy;

			// F для верхней границы
			for (i = 1; i < M - 1; ++i) {
				F[j * M + i] = -f(i * hx, j * hy) - phi4() / hy / hy;
			}

			F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi4() / hy / hy;

			// F для правой границы
			for (j = N - 2; j > 0; --j) {
				F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx;
			}
			F[j * M + i] = -f(i * hx, j * hy) - phi2() / hx / hx - phi3() / hy / hy;

			// F  для нижней границы
			for (i = M - 2; i > 0; --i) {
				F[j * M + i] = -f(i * hx, j * hy) - phi3() / hy / hy;
			}

			// внутренняя область
			for (j = 1; j < N - 1; ++j) {
				for (i = 1; i < M - 1; ++i) {
					F[j * M + i] = -f(i * hx, j * hy);
				}
			}
		}
		return F;
	};

	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\out\\u4_out.csv");

	for (auto N : N_arr) {
		int M = N;
		double hx = static_cast<double>(lx) / N;
		double hy = static_cast<double>(ly) / M;
		cout << "Process: " << N << " " << M << endl;
		std::vector<Eigen::Triplet<double>> triplets;
		// заполним левые гаммы
		for (int j = 1; j < M; ++j) {
			for (int i = 0; i < N; ++i) {
				double value = -by(i * hx, j * hy) / 2. / hy + b(i * hx, j * hy) / hy / hy;
				triplets.push_back(Eigen::Triplet<double>(j * N + i, (j - 1) * N + i, value));
			}
		}

		// правые гаммы
		for (int j = 0; j < M - 1; ++j) {
			for (int i = 0; i < N; ++i) {
				double value = by(i * hx, j * hy) / 2. / hy + b(i * hx, j * hy) / hy / hy;
				triplets.push_back(Eigen::Triplet<double>(j * N + i, (j + 1) * N + i, value));
			}
		}

		// альфы
		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < N; ++i) {
				double value = -2. * (a(i * hx, j * hy) / hx / hx + b(i * hx, j * hy) / hy / hy) + c(i * hx, j * hy);
				triplets.push_back(Eigen::Triplet<double>(j * N + i, j * N + i, value));
			}
		}

		// левые беты
		for (int j = 0; j < M; ++j) {
			for (int i = 1; i < N; ++i) {
				double value = -ax(i * hx, j * hy) / 2 / hx + a(i * hx, j * hy) / hx / hx;
				triplets.push_back(Eigen::Triplet<double>(j * N + i, j * N + i - 1, value));
			}
		}

		// правые беты
		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < N - 1; ++i) {
				double value = ax(i * hx, j * hy) / 2 / hx + a(i * hx, j * hy) / hx / hx;
				triplets.push_back(Eigen::Triplet<double>(j * N + i, j * N + i + 1, value));
			}
		}

		Eigen::SparseMatrix<double> A(N * M, N * M);
		A.setFromTriplets(triplets.begin(), triplets.end());

		ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
		cg.compute(A);
		auto b = getF(N, M, hx, hy);
		auto U = cg.solve(b);


		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < N; ++i) {
				u_1 << i * hx << " " << j * hy << " " << U[j * N + i] << "\n";
			}
		}
		u_1.close();

	}

}