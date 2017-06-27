#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lineq.h"
#include "funs.h"

void least_square_fit (vector *x, vector *y, vector *dy, vector *c, 
		matrix *S, double f(int, double));


int main(int argc, const char *argv[])
{
	// Exercise A
	double x1[] = {0.100, 0.145, 0.211, 0.307, 0.447,
			0.649, 0.944, 1.372, 1.995, 2.900};

	double y1[] = {12.644, 9.235, 7.377, 6.460, 5.555,
			5.896, 5.673, 6.964, 8.896, 11.355};

	double dy1[] = {0.858, 0.359, 0.505, 0.403, 0.683,
			0.605, 0.856, 0.351, 1.083, 1.002};

	// Allocate stuff
	int m = 3, n = 10, nx = 1000;
	vector *x = vector_alloc (n);
	vector *y = vector_alloc (n);
	vector *dy = vector_alloc (n);
	vector *c = vector_alloc (m);
	matrix *S = matrix_alloc (m, m);

	for (int i = 0; i < n; i++) {
		vector_set (x, i, x1[i]);
		vector_set (y, i, y1[i]);
		vector_set (dy, i, dy1[i]);
		printf("%g %g %g\n", x1[i], y1[i], dy1[i]);
	}
	printf("\n\n");

	least_square_fit (x, y, dy, c, S, funs);

	double dx =  (vector_get (x, n-1) - vector_get (x, 0)) / (nx-1);

	double xi = vector_get (x, 0) - vector_get (x, 0) / 5;
	while (xi < vector_get(x, n-1) + dx) {
		printf("%g %g %g %g\n", xi,
			fit(xi, funs, c),
			fit(xi, funs, c) - df(xi, funs, S),
			fit(xi, funs, c) + df(xi, funs, S));
		xi += dx;
	}
	printf("\n\n");

	vector_free (x);
	vector_free (y);
	vector_free (c);
	vector_free (dy);
	matrix_free (S);

	return 0;
}
