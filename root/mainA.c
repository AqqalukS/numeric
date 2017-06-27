#include <stdio.h>
#include <math.h>
#include "root.h"
#include "matrix.h"
#include "funs.h"

int main (int argc, const char *argv[]) {
	vector *x = vector_alloc (2);
	vector *fx = vector_alloc (2);
	int ds;
	double dx, eps;
	dx = 1e-9, eps = 1e-9;
	
	// first system of equation
	printf("#######################\n");
	printf("# Newton root finding #\n");
	printf("#######################\n");
	vector_set (x, 0, 2);
	vector_set (x, 1, 1);
	
	printf("Initial:\n"); printf("x = ( ");
	for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);

	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	// Newton root finding
	ds = newton_root (f1, x, dx, eps);

	printf("number of function calls = %d\n", ds);

	printf("x = ( "); for (int i = 0; i < x->size; i++) printf("%g ", vector_get (x, i));
	printf(")\n");

	f1(x, fx);
	printf("f(x) = ( ");for (int i = 0; i < x->size; i++) printf("%g ", vector_get (fx, i));
	printf(")\n");

	return 0;
}
