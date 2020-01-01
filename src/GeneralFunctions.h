#ifndef _GENERALFUNCTIONS_
#define _GENERALFUNCTIONS_

#include <limits>

#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))

const double tolerance =1.0e-8;
const double infinite = std::numeric_limits<double>::max();
const int infiniteInt = std::numeric_limits<int>::max();

double RelDif(double a, double b);
double RelDifNoAbs(double a, double b);

inline int signum(double x) {return((x>0)-(x<0));}


#endif
