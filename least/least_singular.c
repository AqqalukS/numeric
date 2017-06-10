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
	matrix *B = matrix_alloc (m, m);
	matrix *V = matrix_alloc (m, m);
	matrix *AV = matrix_alloc (n, m);
	matrix *U = matrix_alloc (n, m);
	vector *u = vector_alloc (m);
	vector *e = vector_alloc (m);
	vector *b = vector_alloc (n);

	for (int i = 0; i < n; i++) {
		vector_set (b, i, vector_get (y, i) / vector_get (dy, i));
		for (int j = 0; j < m; j++) {
			matrix_set (A, i, j, f(j, vector_get (x,i)) / vector_get (dy, i));
			matrix_set (A2, i, j, matrix_get (A, i, j));
		}
	}


	matrixT_x_matrix (A, A2, D);
/*	double s;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			s = 0;
			for (int k = 0; k < n; k++) {
				s += matrix_get (A, k, i) * matrix_get (A2, k, j);
			}
			matrix_set (D, i, j, s);
		}
	}*/
	int jacobi = jacobi_cyclic (D, e, V);
	matrix_set_identity (D); 
	for (int i = 0; i < m; i++) {
		matrix_set (D, i, i, pow(vector_get (e, i), -0.5));
		fprintf(stderr, "%g\n", vector_get (e, i));
	}
	for (int i = 0; i < m; i++) {
		matrix_set (R, i, i, pow(vector_get (e, i), 0.5));
	}

	matrix_x_matrix (A, V, AV);
	matrix_x_matrix (AV, D, U);

	matrixT_x_vector (U, b, u);
	backsub (R, u);

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
