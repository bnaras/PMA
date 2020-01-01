#include "GeneralFunctions.h"

double RelDif(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : Abs(a - b) / d;
}

double RelDifNoAbs(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : (a - b) / d;
}
