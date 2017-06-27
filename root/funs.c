#include <math.h>
#include "vector.h"
#include "matrix.h"
void f1 (vector *v, vector *df) {
	double A = 10000;
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	vector_set (df, 0, A*x*y -1);
	vector_set (df, 1, exp(-x) + exp(-y) - 1 - 1./A);
}
void J_f1 (vector *v, matrix *J) {
	double A = 10000;
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	matrix_set (J, 0, 0, A * y);
	matrix_set (J, 0, 1, A * x);
	matrix_set (J, 1, 0, -exp(-x));
	matrix_set (J, 1, 1, -exp(-y));
}
void rosen (vector *v, vector *df) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	vector_set (df, 0, -2*(1-x) - 400*x*(y - x*x));
	vector_set (df, 1, 200*(y-x*x));
}
void J_rosen (vector *v, matrix *J) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	matrix_set (J, 0, 0, 1200*x*x - 400*y + 2);
	matrix_set (J, 0, 1, -400*x);
	matrix_set (J, 1, 0, -400*x);
	matrix_set (J, 1, 1, 200);
}
void himmel (vector *v, vector *df) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	vector_set (df, 0, 4*x*(x*x + y - 11) + 2*(x + y*y -7));
	vector_set (df, 1, 4*y*(y*y + x - 7) + 2*(y + x*x - 11));
}
void J_himmel (vector *v, matrix *J) {
	double x = vector_get (v, 0);
	double y = vector_get (v, 1);
	matrix_set (J, 0, 0, 12*x*x + 4*y - 42);
	matrix_set (J, 0, 1, 4*x + 4*y);
	matrix_set (J, 1, 0, 4*x + 4*y);
	matrix_set (J, 1, 1, 12*y*y + 4*x - 26);
}
