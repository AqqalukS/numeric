#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#define RND ((double)rand() / RAND_MAX)

void lanczos (matrix *A, matrix *u, matrix *V, matrix *T);

int main (int argc, const char *argv[]) {
	int n = atof(argv[1]);
	matrix *A = matrix_alloc (n, n);
	matrix *V = matrix_alloc (n, n);
	matrix *T = matrix_alloc (n, n);
	matrix *R = matrix_alloc (n, n);
	vector *v = vector_alloc (n);

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, RND);
			matrix_set (A, j, i, matrix_get (A, i, j));
		}
		vector_set (v, i, RND);
	}

	printf("########################################\n");
	printf("# Lanczos tridiagonalization algorithm #\n");
	printf("########################################\n");

	printf("Matrix A\n");
	matrix_print (A);

	lanczos (A, v, V, T);

	printf("Matrix V\n");
	matrix_print (V);

	printf("matrix T\n");
	matrix_print (T);

	printf("V A V*\n");
	matrix_x_matrix (V, T, R);
	matrix_x_matrixT (R, V, T);
	matrix_print (T);

	printf("V V*\n");
	matrix_x_matrixT (V, V, R);
	matrix_print (R);


	matrix_free (A);
	matrix_free (T);
	matrix_free (R);
	matrix_free (V);
	return 0;
}
