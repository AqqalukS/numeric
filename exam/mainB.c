#include <stdio.h>
#include "matrix.h"
#include "vector.h"
#define RND ((double)rand() / RAND_MAX)

void arnoldi (matrix *A, vector *q, matrix *Q, matrix *H);
void lanczos (matrix *A, vector *u, matrix *V, matrix *T);

int main (int argc, const char *argv[]) {
	int n = atof(argv[1]);
	matrix *A = matrix_alloc (n, n);
	matrix *Q = matrix_alloc (n, n);
	matrix *H = matrix_alloc (n, n);
	matrix *R = matrix_alloc (n, n);
	vector *q = vector_alloc (n);

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, RND);
			matrix_set (A, j, i, matrix_get (A, i, j));
		}
		vector_set (q, i, RND);
	}

	arnoldi (A, q, Q, H);
	printf("#####################\n");
	printf("# Arnoldi iteration #\n");
	printf("#####################\n");

	printf("Matrix A\n");
	matrix_print (A);

	printf("# Arnoldi iteration #\n");
	printf("#####################\n");
	printf("Matrix H =\n");
	matrix_print (H);

	printf("Matrix Q =\n");
	matrix_print (Q);
	
	lanczos (A, q, Q, H);

	printf("# Lanczos algorithm #\n");
	printf("#####################\n");
	printf("Matrix T =\n");
	matrix_print (H);

	printf("Matrix V =\n");
	matrix_print (Q);

	arnoldi (A, q, Q, H);

	printf("Q H Q* = A\n");
	matrix_x_matrix (Q ,H, R);
	matrix_x_matrixT (R, Q, H);
	matrix_print (H);

	matrix_free (A);
	matrix_free (Q);
	matrix_free (H);
	matrix_free (R);
	return 0;
}
