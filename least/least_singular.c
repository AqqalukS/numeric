#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "lineq.h"
#include "jacobi.h"

void least_square_singular (vector *x, vector *y, vector *dy, vector *c, matrix *S,
		double f(int, double)) {
	int n = x->size;
	int m = S->size1;
	matrix *A = matrix_alloc (n, m);
	matrix *A2 = matrix_alloc (n, m);
	matrix *R = matrix_alloc (m, m);
	matrix *D = matrix_alloc (m, m);
	vector *b = vector_alloc (n);

	for (int i = 0; i < n; i++) {
		vector_set (b, i, vector_get (y, i) / vector_get (dy, i));
		for (int j = 0; j < m; j++) {
			matrix_set (A, i, j, f(j, vector_get (x,i)) / vector_get (dy, i));
			matrix_set (A2, i, j, matrix_get (A, i, j));
		}
	}
fprintf(stderr, "1\n");

	vector *e = vector_alloc (m);
fprintf(stderr, "2\n");
	fprintf(stderr, "e(n) = e(%i)\n", e->size);
	matrix *B = matrix_alloc (m, m);
	matrix *V = matrix_alloc (m, m);

	fprintf(stderr, "D(n, m) = D(%i, %i)\n", D->size1, D->size2);
	fprintf(stderr, "e(n) = e(%i)\n", e->size);
	matrixT_x_matrix (A, A2, D);
	fprintf(stderr, "D(n, m) = D(%i, %i)\n", D->size1, D->size2);
	fprintf(stderr, "e(n) = e(%i)\n", e->size);
	int jacobi = jacobi_cyclic (D, e, V);
	matrix_set_identity (B); 
	for (int i = 0; i < m; i++) {
		matrix_set (B, i, i, pow(vector_get (e, i), -0.5));
		fprintf(stderr, "%g\n", vector_get (e, i));
	}
	for (int i = 0; i < m; i++) {
		matrix_set (R, i, i, pow(vector_get (e, i), 0.5));
	}

fprintf(stderr, "6\n");
	vector *u = vector_alloc (n);
	matrix *AV = matrix_alloc (n, m);
	matrix_x_matrix (A, V, AV);
fprintf(stderr, "3\n");
	matrix *U = matrix_alloc (n, m);
	matrix_x_matrix (AV, B, U);
fprintf(stderr, "4\n");

	qr_gs_solve (R, U, b, u);
	matrix_x_vector (V, u, c);
	for (int i = 0; i < m; i++) {
		matrix_set (R, i, i, pow(matrix_get (R, i, i), -2));
	}
	matrix_x_matrix (V, R, B);
	matrix_x_matrixT (B, V, S);

	matrix_free (A);
	matrix_free (A2);
	matrix_free (R);
	matrix_free (D);
	matrix_free (B);
	matrix_free (V);
	matrix_free (U);
	matrix_free (AV);
	vector_free (b);
	vector_free (e);
	vector_free (u);
}
