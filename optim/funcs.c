#include <math.h>
#include "funcs.h"
#include "vector.h"
#include "matrix.h"

// Functions
double rosen (vector *v) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	return pow(1 - x, 2) + 100 * pow(y - x*x, 2);
}

double himmel (vector *v) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	return pow(x*x + y - 11, 2) + pow(x + y*y - 7, 2);
}

// gradients
void rosen_df (vector *v, vector *df) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	vector_set (df, 0, -2*(1-x) - 400*x*(y - x*x));
	vector_set (df, 1, 200*(y-x*x));
}

void himmel_df (vector *v, vector *df) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	vector_set (df, 0, 4*x*(x*x + y - 11) + 2*(x + y*y -7));
	vector_set (df, 1, 4*y*(y*y + x - 7) + 2*(y + x*x - 11));
}

// Hessian matrices
void rosen_H (vector *v, matrix *H) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	matrix_set (H, 0, 0, 1200*x*x - 400*y + 2);
	matrix_set (H, 1, 1, 200);
	matrix_set (H, 0, 1, -400*x);
	matrix_set (H, 1, 0, -400*x);
}

void himmel_H (vector *v, matrix *H) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	matrix_set (H, 0, 0, 2 * (6 *x*x - 21) + 4*y);
	matrix_set (H, 1, 1, 2 * (6 *y*y - 13) + 4*x);
	matrix_set (H, 0, 1, 4*x + 4*y);
	matrix_set (H, 1, 0, 4*x + 4*y);
}
