#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"
#include "funcs.h"

int quasi_newton (double f(vector *), vector *x, double dx, double eps);
int newton_optim (double f(vector *), void gradient(vector *, vector *),
		void hessian(vector *, matrix *), vector *x, double eps);

int main(int argc, const char *argv[])
{
	// Exercise B
	printf("##############\n");
	printf("# Exercise B #\n");
	printf("##############\n");

	// Allocations
	int n = 2, ds;
	vector *x = vector_alloc (n);
	vector_set (x, 0, 2);
	vector_set (x, 1, 2);

	double fx, eps;
	eps = 1e-8;

	double dx = 1e-5;
	vector_set (x, 0, -3);
	vector_set (x, 1, -3);
	ds = quasi_newton (himmel, x, dx, eps);
	printf("Number of steps = %d\n", ds);
	vector_print (x);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);

	double F (vector *v) {
		double A = vector_get (v, 0);
		double T = vector_get (v, 1);
		double B = vector_get (v, 2);
		double s = 0;
		for (int i = 0; i < N; i++) {
			s += pow((A * exp(-t[i] / T) + B) - y[i], 2) / pow(e[i], 2);
		}
		return s;
	}

	int M = 3;
	vector *u = vector_alloc (M);
	vector_set (u, 0, 1);
	vector_set (u, 1, 1);
	vector_set (u, 2, 1);

	ds = quasi_newton (F, u, dx, eps);

	for (int i = 0; i < N; i++) {
		fprintf(stderr,"%g %g %g %g\n", t[i], y[i], 
			vector_get (u, 0) * exp(-t[i] / vector_get (u, 1)) 
					+ vector_get (u, 2), e[i]);
	}

	vector_free (x);
	vector_free (u);
	return 0;
}
