#include "lineq.h"
#define RND ((double)rand()/RAND_MAX)

int main() {
	printf("############################################################\n");
	printf("# Exercise A - QR-decomp by Gram-Schmidt orthogonalization #\n");
	printf("############################################################\n");
	printf("# Part one - Orthogonalization of a nxm (n>m) Matrix A #\n");
	printf("########################################################\n");
	int n = 6, m = 3;
	matrix *A = matrix_alloc (n,m);
	matrix *Q = matrix_alloc (n,m);
	matrix *C = matrix_alloc (n,m);
	matrix *R = matrix_alloc (m,m);
	matrix *I = matrix_alloc (m,m);

	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			matrix_set (A, i, j, RND);
		}
	}
	printf("# Matrix A\n");
	matrix_print (A);

	qr_gs_decomp (A, R);

	// copy A -> Q
	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			matrix_set (Q, i, j, matrix_get (A, i, j));
		}
	}
	matrixT_x_matrix (A, Q, I);
	printf("# Q^T * Q = \n");
	matrix_print (I);

	printf("# Matrix R\n");
	matrix_print (R);

	matrix_x_matrix (A, R, C);

	printf("# QR = A \n");
	matrix_print (C);

	matrix_free (A);
	matrix_free (I);
	matrix_free (C);
	matrix_free (R);

	// first task second exercise
	printf("############################\n");
	printf("# Part two - solve QRx = b #\n");
	printf("############################\n");

	matrix *D = matrix_alloc (n,n);
	matrix *E = matrix_alloc (n,n);
	matrix *RR = matrix_alloc (n,n);
	vector *b = vector_alloc (n);
	vector *x = vector_alloc (n);
	vector *x2 = vector_alloc (n);

	for (int i = 0; i < D->size1; i++) {
		for (int j = 0; j < D->size2; j++) {
			matrix_set (D, i, j, RND);
			matrix_set (E, i, j, matrix_get (D, i, j));
		}
		vector_set (b, i, RND);
	}

	printf("# matrix A(n,n)\n");
	matrix_print (D);

	printf("# Vector b(n)\n");
	vector_print (b);

	qr_gs_decomp (D, RR);
	qr_gs_solve (D, RR, b, x);
	printf("# Solving for x in QRx = b\n");
	vector_print (x);

	matrix_x_vector (E, x, x2);
	printf("# A x = b =\n");
	vector_print (x2);

	matrix_free (D);
	matrix_free (E);
	matrix_free (RR);
	vector_free (b);
	vector_free (x2);
	vector_free (x);
	return 0;
}
