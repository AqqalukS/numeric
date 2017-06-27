#include <stdio.h>
#include <math.h>
#include "root.h"
#include "matrix.h"
#include "funs.h"

int main (int argc, const char *argv[]) {
	vector *x = vector_alloc (2);
	vector *fx = vector_alloc (2);

	int ds;
	double eps;
	
	// first system of equation
	eps = 1e-9;

	vector_set (x, 0, 2);
	vector_set (x, 1, 1);
	
	printf("Initial:\n"); printf("x = ( ");
	for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);

	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	// Newton root finding with analytic jacobian
	ds = newton_jacobian (f1, J_f1, x, eps);

	printf("number of function calls = %d\n", ds);

	printf("x = ( "); for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);
	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	return 0;
}
