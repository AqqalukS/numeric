#include <stdio.h>
#include <math.h>
#include "matrix.h"

int newton_root (void f(vector *, vector *), vector *x, double dx, double eps);

int main(int argc, const char *argv[])
{
	vector *x = vector_alloc (2);
	vector *fx = vector_alloc (2);
	int ds, calls;
	double dx, eps;

	// first system of equation
	calls = 0;
	dx = 1e-3, eps = 1e-5;
	double A = 1000;

	void f1 (vector *v, vector *df) {
		calls++;
		double x = vector_get (v, 0);
		double y = vector_get (v, 1);
		vector_set (df, 0, A*x*y -1);
		vector_set (df, 1, exp(-x) + exp(-y) - 1 - 1./A);
	}
	
	vector_set (x, 0, 1);
	vector_set (x, 1, 1);
	
	printf("Initial:\n"); printf("x = ( ");
	for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);

	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	ds = newton_root (f1, x, dx, eps);

	printf("number of steps = %d\n", ds);

	printf("x = ( "); for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);
	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	// Rosenbrock's valley function
	calls = 0;
	dx = 1e-9, eps = 1e-9;

	void f2 (vector *v, vector *df) {
		calls++;
		double x = vector_get (v, 0);
		double y = vector_get (v, 1);
		vector_set (df, 0, -2*(1-x) - 400*x*(y - x*x));
		vector_set (df, 1, 200*(y-x*x));
	}
	
	vector_set (x, 0, 5);
	vector_set (x, 1, 5);
	
	printf("Initial:\n"); printf("x = ( ");
	for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f2(x, fx);

	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	ds = newton_root (f2, x, dx, eps);

	printf("number of steps = %d\n", ds);

	printf("x = ( "); for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f2(x, fx);
	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	// Himmelblau's function
	calls = 0;
	dx = 1e-9, eps = 1e-9;

	void f3 (vector *v, vector *df) {
		calls++;
		double x = vector_get (v, 0);
		double y = vector_get (v, 1);
		vector_set (df, 0, 4*x*(x*x + y - 11) + 2*(x + y*y -7));
		vector_set (df, 1, 4*y*(y*y + x - 7) + 2*(y + x*x - 11));
	}
	
	vector_set (x, 0, 1);
	vector_set (x, 1, 1);
	
	printf("Initial:\n"); printf("x = ( ");
	for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f3(x, fx);

	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	ds = newton_root (f3, x, dx, eps);

	printf("number of steps = %d\n", ds);

	printf("x = ( "); for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f3(x, fx);
	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");
	vector_free (x);
	vector_free (fx);
	return 0;
}
