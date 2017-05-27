#include "lineq.h"
#define RND ((double)rand()/RAND_MAX)

int main() {
	printf("#############################################\n");
	printf("# Exercise C - QR-decomp by Givens rotation #\n");
	printf("#############################################\n");
	printf("# Part one - Orthogonalization of a nxm (n>m) Matrix A #\n");
	printf("########################################################\n");
	int n = 6, m = 3;
	matrix *A = matrix_alloc (n,m);
	matrix *Q1 = matrix_alloc (n,m);
	matrix *Q2 = matrix_alloc (n,m);
	matrix *C = matrix_alloc (n,m);
	matrix *R = matrix_alloc (m,m);
	matrix *I = matrix_alloc (m,m);

	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			matrix_set (A, i, j, RND);
			matrix_set (Q2, i, j, matrix_get (A, i, j));
		}
	}
	printf("# Matrix A\n");
	matrix_print (A);

	givens (Q2);
	givens_unpack_Q (Q2, Q1);

	// copy A -> Q
	for (int i = 0; i < Q1->size1; i++) {
		for (int j = 0; j < Q1->size2; j++) {
			matrix_set (Q2, i, j, matrix_get (Q1, i, j));
		}
	}
	matrixT_x_matrix (Q1, Q2, I);
	printf("# Q^T * Q = \n");
	matrix_print (I);

	printf("# Matrix R isn't calculated as in Gram-Schmidt, so we get it by:\n");
	printf("# R = Q^T A since Q is an orthogonal matrix\n");
	matrixT_x_matrix (Q1, A, R);
	matrix_print (R);
	printf("# This shows that QR = A \n");

	matrix_x_matrix (A, R, C);

	matrix_print (C);

	matrix_free (A);
	matrix_free (Q1);
	matrix_free (Q2);
	matrix_free (I);
	matrix_free (C);
	matrix_free (R);

	
	// first task second exercise
	printf("############################\n");
	printf("# Part two - solve QRx = b #\n");
	printf("############################\n");

	matrix *D = matrix_alloc (n,n);
	matrix *E = matrix_alloc (n,n);
	matrix *invA = matrix_alloc (n,n);
	vector *b = vector_alloc (n);
	vector *x = vector_alloc (n);
	vector *b2 = vector_alloc (n);

	for (int i = 0; i < D->size1; i++) {
		for (int j = 0; j < D->size2; j++) {
			matrix_set (D, i, j, RND);
			matrix_set (E, i, j, matrix_get (D, i, j));
		}
		vector_set (b, i, RND);
		vector_set (x, i, vector_get (b, i));
	}

	printf("# matrix A(n,n) = \n");
	matrix_print (D);

	printf("# Vector b(n) = \n");
	vector_print (x);

	givens (D);
	givens_solve (D, x);
	printf("# Solving for x in QRx = b\n");
	vector_print (x);

	matrix_x_vector (E, x, b2);
	printf("# A x = b =\n");
	vector_print (b2);

	printf("######################################\n");
	printf("# Matrix inverse by Givens rotations #\n");
	printf("######################################\n");
	
	givens_inverse (D, invA);

	printf("# A =\n");
	matrix_print (E);

	printf("# A^-1 = \n");
	matrix_print (invA);

	matrix_x_matrix (E, invA, D);
	printf("# A * A^-1 = I = \n");
	matrix_print (D);

	matrix_free (D);
	matrix_free (E);
	matrix_free (invA);
	vector_free (b);
	vector_free (b2);
	vector_free (x);

	return 0;
}
