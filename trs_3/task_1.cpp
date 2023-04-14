#include "trs_3.h"

using namespace std;

constexpr auto pi = 3.14159265358979323846;

double U(double x, double y)
{
	return sin(pi * x / 2.) * sin(pi * y) * sin(pi * y);
}

double F(double x, double y)
{
	return pi * pi / 4. * sin(x * pi / 2.) * (17 * cos(pi * y) * cos(pi * y) - 9.);
}

int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\out\\u1_out.csv");
	double max_diff = 0;
	const int n = 100;
	const int m = 100;
	const int lx = 2, ly = 1;
	double hx = static_cast<double>(lx) / n;
	double hy = static_cast<double>(ly) / m;
	vector < vector <double> > u_prev(n, vector <double>(m));
	vector < vector <double> > u(n, vector <double>(m));
	vector < vector <double> > f(n, vector <double>(m));
	double betta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	double alpha = -2. / hx / hx - 2. / hy / hy;

	f[1][1] = F(hx, hy) - u_prev[0][1] / hx / hx - u_prev[1][0] / hy / hy;
	f[1][m - 2] = F(hx, (m - 1) * hy) - u_prev[0][m - 2] / hx / hx - u_prev[1][m - 1] / hy / hy;
	f[n - 2][1] = F(hx * (n - 1), hy) - u_prev[n - 1][1] / hx / hx - u_prev[n - 2][0] / hy / hy;
	f[n - 2][m - 2] = F(hx * (n - 1), hy * (m - 1)) - u_prev[n - 1][m - 2] / hx / hx - u_prev[n - 2][m - 1] / hy / hy;

	for (int j = 2; j < m - 2; j++)
	{
		f[1][j] = F(hx, j * hy) - u_prev[0][j] / hx / hx;
		f[n - 2][j] = F((n - 1) * hx, j * hy) - u_prev[n - 1][j] / hx / hx;
	}

	for (int i = 2; i < n - 2; i++)
	{
		f[i][1] = F(i * hx, hy) - u_prev[i][1] / hy / hy;
		f[i][m - 2] = F(i * hx, (m - 1) * hy) - u_prev[i][m - 2] / hy / hy;
	}

	int kol = 0;
	do
	{
		kol++;
		max_diff = -1;

		
		if (kol == 1) u_prev = f;
		else u_prev = u;

		//когда первый
		u[1][1] = f[1][1] / alpha + (-u_prev[2][1] * betta / alpha - u_prev[1][2] * gamma / alpha);
		for (int i = 2; i < n - 1; i++)
		{
			u[i][1] = f[i][1] / alpha + (-u_prev[i - 1][1] * betta / alpha - u_prev[i + 1][1] * betta / alpha - u_prev[i][2] * gamma / alpha);
		}
		u[n - 1][1] = f[n - 1][1] / alpha + (-u_prev[n - 2][1] * betta / alpha - u_prev[n - 1][2] * gamma / alpha);

		//когда между
		for (int j = 2; j < m - 1; j++)
		{
			u[1][j] = f[1][j] / alpha + (-u_prev[1][j - 1] * gamma / alpha - u_prev[2][j] * betta / alpha - u_prev[1][j + 1] * gamma / alpha);
			for (int i = 2; i < n - 1; i++)
			{
				u[i][j] = F(i * hx, j * hy) / alpha + (-u_prev[i][j - 1] * gamma / alpha - u_prev[i - 1][j] * betta / alpha - u_prev[i + 1][j] * betta / alpha - u_prev[i][j + 1] * gamma / alpha);
			}
			u[n - 1][j] = f[n - 1][j] / alpha + (-u_prev[n - 1][j - 1] * gamma / alpha - u_prev[n - 2][j] * betta / alpha - u_prev[n - 1][j + 1] * gamma / alpha);
		}

		//когда последний
		u[1][m - 1] = f[1][m - 1] / alpha + (-u_prev[2][m - 1] * betta / alpha - u_prev[1][m - 2] * gamma / alpha);
		for (int i = 2; i < n - 1; i++)
		{
			u[i][m - 1] = f[i][m - 1] / alpha + (-u_prev[i - 1][m - 1] * betta / alpha - u_prev[i + 1][m - 1] * betta / alpha - u_prev[i][m - 2] * gamma / alpha);
		}
		u[n - 1][m - 1] = f[n - 1][m - 1] / alpha + (-u_prev[n - 2][m - 1] * betta / alpha - u_prev[n - 1][m - 2] * gamma / alpha);

		for (int j = 0; j < m; j++)
		{
			for (int i = 0; i < n; i++)
			{
				double diff = abs(u_prev[i][j] - u[i][j]);
				//double diff = abs(U(i * hx, j * hy) - u_prev[i][j]);
					if (diff > max_diff)
					{
						max_diff = diff;
					}
			}
		}
		//cout << i * hx << " " << j * hy << "          " << u[i][j] << "           " << U(i * hx, j * hy) << endl;
		//cout << max_diff << endl;
	} while (max_diff > 1e-15);

	max_diff = 0;
	
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			double diff = abs(U(i * hx, j * hy) - u[i][j]);
			if (diff > max_diff)
			{
				max_diff = diff;
			}

			u_1 << i * hx << " " << j * hy << " " << u[i][j] << endl;
		}
	}
	cout << "Diff: " << max_diff << "  " << kol << endl;

}
