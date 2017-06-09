#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "funs.h"
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
