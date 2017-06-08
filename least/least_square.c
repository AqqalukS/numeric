#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lineq.h"
#include "jacobi.h"


void least_square_fit (vector *x, vector *y, vector *dy, vector *c, 
		matrix *S, double f(int j, double x)) {
	int n = x->size, m = c->size;
	matrix *A = matrix_alloc (n, m);
	matrix *R = matrix_alloc (m, m);
	matrix *Rinv = matrix_alloc (m, m);
	vector *b = vector_alloc (n);
	vector *r = vector_alloc (m);

	for (int i = 0; i < n; i++) {
		vector_set (b, i, vector_get (y, i)/
				vector_get (dy, i));
		for (int j = 0; j < m; j++) {
			matrix_set (A, i, j, f(j, vector_get (x,i)) /
						vector_get (dy, i));
		}
	}

	qr_gs_decomp (A, R);
	qr_gs_solve (A, R, b, c);

	matrix_set_identity (Rinv);
	for (int i = 0; i < m; i++) {
		matrix_get_column (Rinv, i, r);
		backsub (R, r);
		matrix_set_column (Rinv, i, r);
	}

	matrix_x_matrixT (Rinv, Rinv, S);
	matrix_free (A);
	matrix_free (R);
	vector_free (b);
	vector_free (r);
}

void least_square_singular (vector *x, vector *y, vector *dy, vector *c, matrix *S,
		double, f(int, double)) {
	int n = x->size, m = c->size;
	matrix *A = matrix_alloc (n, m);
	matrix *D = matrix_alloc (m, m);
	matrix *S = matrix_alloc (m, m);
	matrix *V = matrix_alloc (m, m);
	matrix *U = matrix_alloc (m, m);
	vector *b = vector_alloc (n);
	vector *e = vector_alloc (m);

	for (int i = 0; i < n; i++) {
		vector_set (b, i, vector_get (y, i)/
				vector_get (dy, i));
		for (int j = 0; j < m; j++) {
			matrix_set (A, i, j, f(j, vector_get (x,i)) /
						vector_get (dy, i));
		}
	}
	matrixT_x_matrix (A, A, D);
	int *jacobi = jacobi_cyclic (D, e, V);
	matrix_set_identity (D); 
	for (int i = 0; i < m; i++) {
		matrix_set (D, i, i, vector_get (e, i));
	}

	free (jacobi);
	matrix_free (A);
	matrix_free (D);
	matrix_free (S);
	matrix_free (V);
	matrix_free (U);
	vector_free (b);
	vector_free (e);
}
