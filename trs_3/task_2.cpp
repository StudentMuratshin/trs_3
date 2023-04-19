#include "trs_3.h"

using namespace std;


int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\out\\u2_out.csv");
	vector<int> n = { 10,20,30,40,50,60,70,80,90,100 };
	for (auto N: n)
	{
		double max_diff = 0;
		const int n = N;
		const int m = N;
		const double w = 1.5;
		const int lx = 2, ly = 1;
		double hx = static_cast<double>(lx) / n;
		double hy = static_cast<double>(ly) / m;
		vector < vector <double> > u_prev(n, vector <double>(m));
		vector < vector <double> > u(n, vector <double>(m));
		vector < vector <double> > f(n, vector <double>(m));
		double betta = 1. / hx / hx;
		double gamma = 1. / hy / hy;
		double alpha = -2. / hx / hx - 2. / hy / hy;

		f[1][1] = F_right(hx, hy) - u_prev[0][1] / hx / hx - u_prev[1][0] / hy / hy;
		f[1][m - 2] = F_right(hx, (m - 1) * hy) - u_prev[0][m - 2] / hx / hx - u_prev[1][m - 1] / hy / hy;
		f[n - 2][1] = F_right(hx * (n - 1), hy) - u_prev[n - 1][1] / hx / hx - u_prev[n - 2][0] / hy / hy;
		f[n - 2][m - 2] = F_right(hx * (n - 1), hy * (m - 1)) - u_prev[n - 1][m - 2] / hx / hx - u_prev[n - 2][m - 1] / hy / hy;

		for (int j = 2; j < m - 2; j++)
		{
			f[1][j] = F_right(hx, j * hy) - u_prev[0][j] / hx / hx;
			f[n - 2][j] = F_right((n - 1) * hx, j * hy) - u_prev[n - 1][j] / hx / hx;
		}

		for (int i = 2; i < n - 2; i++)
		{
			f[i][1] = F_right(i * hx, hy) - u_prev[i][1] / hy / hy;
			f[i][m - 2] = F_right(i * hx, (m - 1) * hy) - u_prev[i][m - 2] / hy / hy;
		}

		int kol = 0;
		do
		{
			kol++;
			max_diff = -1;


			if (kol != 1) u_prev = u;

			//когда первый
			u[1][1] = (1 - w) * u_prev[1][1] + w / alpha * (f[1][1] - u_prev[2][1] * betta - u_prev[1][2] * gamma);
			for (int i = 2; i < n - 1; i++)
			{
				u[i][1] = (1 - w) * u_prev[i][1] + w / alpha * (f[i][1] - u[i-1][1]*betta - u_prev[i+1][1] * betta - u_prev[i][2] * gamma);
			}
			u[n - 1][1] = (1 - w) * u_prev[n - 1][1] + w / alpha * (f[n - 1][1] - u[n - 2][1] * betta - u_prev[n - 1][2] * gamma);

			//когда между
			for (int j = 2; j < m - 1; j++)
			{
				u[1][j] = (1 - w) * u_prev[1][j] + w / alpha * (f[1][j] - u[1][j - 1] * gamma - u_prev[2][j] * betta - u_prev[1][j + 1] * gamma);
				for (int i = 2; i < n - 1; i++)
				{
					u[i][j] = (1 - w) * u_prev[i][j] + w / alpha * (F_right(i * hx, j * hy) - u[i][j - 1] * gamma - u[i - 1][j] * betta - u_prev[i + 1][j] * betta - u_prev[i][j + 1] * gamma);
				}
				u[n - 1][j] = (1 - w) * u_prev[n - 1][j] + w / alpha * (f[n - 1][j] - u[n - 1][j - 1] * gamma - u[n - 2][j] * betta - u[n - 1][j + 1] * gamma);
			}

			//когда последний
			u[1][m - 1] = (1 - w) * u_prev[1][m - 1] + w / alpha * (f[1][m - 1] - u[1][m - 2] * gamma - u_prev[2][m - 1] * betta);
			for (int i = 2; i < n - 1; i++)
			{
				u[i][m - 1] = (1 - w) * u_prev[i][m - 1] + w / alpha * (f[i][m - 1] - u[i][m - 2] * gamma - u[i][m - 1] * betta - u_prev[i + 1][m - 1] * betta);
			}
			u[n - 1][m - 1] = (1 - w) * u_prev[n - 1][m - 1] + w / alpha * (f[n - 1][m - 1] - u[n - 1][m - 2] * gamma  - u[n - 2][m - 1] * betta);

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
				double diff = abs(True_U(i * hx, j * hy) - u[i][j]);
				if (diff > max_diff)
				{
					max_diff = diff;
				}

				//u_1 << i * hx << " " << j * hy << " " << u[i][j] << endl;
			}
		}
		cout << "m, n: " << N << "   Diff: " << setprecision(8) << max_diff << "   Count: " << kol << endl;

	}
	u_1.close();
}
