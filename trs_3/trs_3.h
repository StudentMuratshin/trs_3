// trs_3.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include<iomanip>

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
// TODO: Reference additional headers your program requires here.
