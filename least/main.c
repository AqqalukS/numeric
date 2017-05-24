#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lineq.h"

void least_square_fit (vector *x, vector *y, vector *dy, vector *c, 
		matrix *S, double f(int, double), int m);

// Exercise A
double funs (int i, double x) {
	switch (i) {
		case 0: return 1.0/x; 	break;
		case 1: return 1.0; 	break;
		case 2: return x; 	break;
		default:{fprintf(stderr, "function: wrong i: %d", i);
			return NAN;
		}
	}
}

double fit (double x, double f(int, double), vector *c) {
	double s = 0;
	for (int i = 0; i < c->size; i++) {
		s += vector_get (c, i) * f(i, x);
	}
	return s;
}

double df (double x, double f(int, double), matrix *S) {
	double s = 0;
	for (int i = 0; i < S->size1; i++) {
		for (int j = 0; j < S->size1; j++) {
			s += f(i, x) * matrix_get (S, i, j) * f(j, x);
		}
	}
	return sqrt(s);
}

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
	int m = 3, n = 10, nx = 100;
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

	least_square_fit (x, y, dy, c, S, funs, m);

	double dx = vector_get (x, 0) + (vector_get (x, n-1) -
				vector_get (x, 0)) / (nx-1);

	double xi = vector_get (x, 0) - dx/10;
	while (xi < vector_get(x, n-1) + dx) {
		printf("%g %g %g %g\n", xi,
			fit(xi, funs, c),
			fit(xi, funs, c) - df(xi, funs, S),
			fit(xi, funs, c) + df(xi, funs, S));
		xi += dx;
	}

	vector_free (x);
	vector_free (y);
	vector_free (c);
	vector_free (dy);
	matrix_free (S);

	return 0;
}
