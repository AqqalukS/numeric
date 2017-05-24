#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

double adapt (double f2, double f3, double func(double), double a, double b, 
		double delta, double eps, int nrec) {
	assert(nrec < 100000000);
	double f1 = func( a + 1*(b - a)/6);
	double f4 = func( a + 5*(b - a)/6);
	double Q = (2*f1 + f2 + f3 + 2*f4) * (b - a)/6;
	double q = (f1 + f2 + f3 + f4) * (b - a)/4;
	double tolerance = delta + eps * fabs(Q);
	double error = fabs(Q-q);
	if (error < tolerance) {
		return Q;
	}
	else {
		double Q1 = adapt (f1, f2, func, a, (a+b)/2, delta/sqrt(2.), eps, nrec + 1);
		double Q2 = adapt (f3, f4, func, (a+b)/2, b, delta/sqrt(2.), eps, nrec + 1);
		return Q1 + Q2;
	}
}

double adapt_int (double func(double), double a, double b, double delta, double eps) {
	int nrec = 0;
	if (isinf(a) && isinf(b)) {
		if (a == -INFINITY && b == INFINITY) {
			a = -1; b = 1;
		}
		else {
			a = 1; b = -1;
		}
		double f_inf_inf (double t) {
			return func(t/(1-t*t)) * (1 + t*t)/pow(1-t*t,2);
		}

		double f2 = f_inf_inf (a + 2*(b - a)/6);
		double f3 = f_inf_inf (a + 4*(b - a)/6);
		return adapt (f2, f3, f_inf_inf, a, b, delta, eps, nrec);
	}
	else if ((!isinf(a) && b == INFINITY) ||
			(!isinf(b) && a == INFINITY)) {
		double aa; 
		if (!isinf(a) && b == INFINITY) {
			aa = a;
			a = 0; b = 1;
		}
		else {
			aa = b;
			a = 1; b = 0;
		}
		double f_a_inf (double t) {
			return func(aa + t/(1-t)) * 1/pow(1-t,2);
		}
		double f2 = f_a_inf (a + 2*(b - a)/6);
		double f3 = f_a_inf (a + 4*(b - a)/6);
		return adapt (f2, f3, f_a_inf, a, b, delta, eps, nrec);
	}
	else if ((a == -INFINITY && !isinf(b)) ||
			(!isinf(a) && b == -INFINITY)) {
		double aa;
		if (b != -INFINITY) {
			aa = b;
			a = -1; b = 0;
		}
		else {
			aa = -a;
			a = 0; b = -1;
		}
		double f_inf_a (double t) {
			return func(aa - t/(1+t)) * 1/pow(1+t,2);
		}
		double f2 = f_inf_a (a + 2*(b - a)/6);
		double f3 = f_inf_a (a + 4*(b - a)/6);
		return adapt (f2, f3, f_inf_a, a, b, delta, eps, nrec);
	}
	else {
		double f2 = func( a + 2*(b - a) / 6);
		double f3 = func( a + 4*(b - a) / 6);
		return adapt (f2, f3, func, a, b, delta, eps, nrec);
	}
}
