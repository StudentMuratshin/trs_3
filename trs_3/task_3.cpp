#include "trs_3.h"



int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\out\\u3_out.csv");
	double max_diff = 0;
	const int n = 6;
	const int m = 5;
	const double w = 1.5;
	const int lx = 2, ly = 1;
	double hx = static_cast<double>(lx) / n;
	double hy = static_cast<double>(ly) / m;
	vector < vector <double> > u_prev(n, vector <double>(m));
	vector < vector <double> > u(n, vector <double>(m));
	vector < double > f((m - 1) * (n - 1));
	vector < vector <double> > LL((m - 1) * (n - 1), vector <double>((m - 1) * (n - 1)));
	vector < vector <double> > UU((m - 1) * (n - 1), vector <double>((m - 1) * (n - 1)));
	vector < vector <double> > A((m - 1) * (n - 1), vector <double>((m - 1) * (n - 1)));
	double betta = 1. / hx / hx;
	double gamma = 1. / hy / hy;
	double alpha = -2. / hx / hx - 2. / hy / hy;

	f[0] = F(hx, hy) - u_prev[0][1] / hx / hx - u_prev[1][0] / hy / hy;
	f[1 * (m - 2)] = F(hx, (m - 1) * hy) - u_prev[0][m - 2] / hx / hx - u_prev[1][m - 1] / hy / hy;
	f[(n - 2) * 1] = F(hx * (n - 1), hy) - u_prev[n - 1][1] / hx / hx - u_prev[n - 2][0] / hy / hy;
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

	for (int j = 2; j < m - 1; j++)
	{
		for (int i = 2; i < n - 1; i++)
		{
			f[i][j] = F(i * hx, j * hy);
		}
	}
	int kol = 0;
	do
	{
		kol++;
		max_diff = -1;


		if (kol != 1) u_prev = u;

		for (int i = 0; i < (m - 2) * (n - 1); i++)
		{
			A[i][i] = alpha;
			A[i][i + 1] = betta;
			A[i + 1][i] = betta;
			A[i][n - 1 + i] = gamma;
			A[n - 1 + i][i] = gamma;
		}

		for (int i = (m - 2) * (n - 1); i < (m - 1) * (n - 1) -1; i++)
		{
			A[i][i] = alpha;
			A[i][i + 1] = betta;
			A[i + 1][i] = betta;
		}
		A[(m - 1) * (n - 1) - 1][(m - 1) * (n - 1) - 1] = alpha;
		
		for (int j = 0; j < (m - 1) * (n - 1); j++)
		{
			for (int i = 0; i < (m - 1) * (n - 1); i++)
			{
				//double diff = abs(u_prev[i][j] - u[i][j]);
				//double diff = abs(U(i * hx, j * hy) - u_prev[i][j]);
				/*if (diff > max_diff)
				{
					max_diff = diff;
				}*/
				cout << A[i][j] << " ";
			}
			cout << endl;
		}
		//cout << i * hx << " " << j * hy << "          " << u[i][j] << "           " << U(i * hx, j * hy) << endl;
		//cout << max_diff << endl;
	} while (max_diff > 1);

	max_diff = 0;

	/*for (int j = 0; j < m; j++)
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
	cout << "Diff: " << max_diff << "  Count: " << kol << endl;*/
	u_1.close();
}