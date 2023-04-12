#include "trs_3.h"

using namespace std;

int main()
{
	const int n = 200;
	const int m = 400;
	const int lx = 2, ly = 1;
	double hx = static_cast<double>(lx) / n;
	double hy = static_cast<double>(ly) / m;
	vector < vector <double> > u(n, vector <double>(m));
	vector < vector <double> > f(n, vector <double>(m));

	for (int i = 1; i < m; i++)
	{
		u[0][i] = 0;
	}
}
