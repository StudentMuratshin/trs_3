#include "trs_3.h"



int main()
{
	//ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\out\\u3_out.csv");
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

	const int N = 30;
	const int M = 30;
	double hx = static_cast<double>(lx) / N;
	double hy = static_cast<double>(ly) / M;
	
	vector<double> U(N * M, 0.); // Начальное приближение нулевой вектор
	
	double alpha = -2. * (1. / hx / hx + 1. / hy / hy);
	double beta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	
	vector<double> F = getF(N, M, hx, hy);
	
	EllipticPDEMatrixGen pm(alpha, beta, gamma, N, M);
	
	vector<double> Lower(M * N * M * N);
	vector<double> Upper(M * N * M * N);
	
	lu_decompostion(N * M, pm, Lower, Upper);
	
	solve_lu(N * M, Lower, Upper, F, U);
	
	double error = 0.;
	
	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < M; ++i) {
			double x = i * hx;
			double y = j * hy;
			error = std::max(error, abs(U[j * M + i] - sin(pi * x / 2.) * sin(pi * y) * sin(pi * y)));
		}
	}
	
	cout << "LU N = " << N << " M  = " << M << " " << " error: " << error << endl;

	//u_1.close();
}