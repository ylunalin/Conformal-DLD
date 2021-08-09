#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
typedef std::complex<double> Comp;

const double pi=3.1415926535897932384626433832795;
const double b_rad=1;

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
double argument(double x,double y);

inline double arg(double x,double y) {
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):
			 (x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}
#endif
