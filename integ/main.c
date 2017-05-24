#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
//#define M_PI 3.14159265358979323846
int calls = 0;

double gsl_int (double func(double), double a, double b, double epsabs, double epsrel);

double adapt_int (double func(double),double a, double b, double delta, double eps);

double f1 (double x) {
	calls++;
	return sqrt(x);
}
double f2 (double x) {
	calls++;
	return 1/sqrt(x);
}
double f3 (double x) {
	calls++;
	return log(x)/sqrt(x);
}
double f4 (double x) {
	calls++;
	return 4 * sqrt(1 - pow((1 - x),2));
}
double f5 (double x) {
	calls++;
return exp(-fabs(x));
}
double f6 (double x) {
	calls++;
return exp(-x*x);
}

int main() {
	double Q, a = 0, b = 1, delta = 1e-9, eps = 1e-9;

	printf("Comparison between adaptive integration and gsl rutines\n\n");
	calls = 0;
	Q = adapt_int (f1, a, b, delta, eps);
	printf("int (sqrt(x), %g, %g) ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.3f\n", Q);
	calls = 0;
	Q = gsl_int (f1, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.3f\n\n", Q);

	calls = 0;
	Q = adapt_int (f2, a, b, delta, eps);
	printf("int (1 / sqrt(x), %g, %g) ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.3f\n", Q);
	calls = 0;
	Q = gsl_int (f2, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.3f\n\n", Q);

	calls = 0; 
	Q = adapt_int (f3, a, b, delta, eps);
	printf("int (ln(x) / sqrt(x), %g, %g) ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.3f\n", Q);
	calls = 0; 
	Q = gsl_int (f3, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.3f\n\n", Q);

	calls = 0;
	a = -INFINITY; b = INFINITY;
	Q = adapt_int (f5, a, b, delta, eps);
	printf("int (exp(-abs(x)), %g, %g), ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.3f\n", Q);
	calls = 0;
	Q = gsl_int (f5, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.3f\n\n", Q);

	b = -2; a = INFINITY;
	calls = 0;
	Q = adapt_int (f6, a, b, delta, eps);
	printf("int (exp(-x^2), %g, %g), ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.3f\n", Q);

	calls = 0;
	Q = gsl_int (f6, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.3f\n\n", Q);

	// f4 = sqrt(1 - (1 - x)^2)
	calls = 0; delta = 1e-10; eps = 1e-10;
	a = 0; b = 1;
	printf("int (sqrt(1 - (1 - x)^2), %g, %g), ", a, b);
	printf("delta = %g, eps = %g\n", delta, eps);
	Q = adapt_int (f4, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("my int  = %7.17f\n", Q);

	calls = 0;
	Q = gsl_int (f4, a, b, delta, eps);
	printf("calls = %10i, ", calls);
	printf("gsl int = %7.17f\n", Q);
	printf("M_PI from math.h            = %7.17f\n", M_PI);

	return 0;
}
