#include "lineq.h"

#define RND ((double)rand()/RAND_MAX)

int main() {
	// first task first exercise
	int n = 6, m = 3;
	matrix *A = matrix_alloc (n,m);
	matrix *C = matrix_alloc (n,m);
	matrix *R = matrix_alloc (m,m);

	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			matrix_set (A, i, j, RND);
		}
	}
	printf("# Matrix A\n");
	matrix_print (A);
	qr_gs_decomp (A, R);

	printf("# Matrix Q\n");
	matrix_print (A);

	printf("# Matrix R\n");
	matrix_print (R);

	matrix_x_matrix (A, R, C);

	printf("# QR = A \n");
	matrix_print (C);

	matrix_free (A);
	matrix_free (C);
	matrix_free (R);

	// first task second exercise
	
	matrix *D = matrix_alloc (n,n);
	matrix *E = matrix_alloc (n,n);
	matrix *RR = matrix_alloc (n,n);
	vector *b = vector_alloc (n);

	for (int i = 0; i < D->size1; i++) {
		for (int j = 0; j < D->size2; j++) {
			matrix_set (D, i, j, RND);
			matrix_set (E, i, j, matrix_get (D, i, j));
		}
		vector_set (b, i, RND);
	}

	printf("# matrix A(n,n)\n");
	matrix_print (D);
	matrix_print (E);

	printf("# Vector b(n)\n");
	vector_print (b);

	vector *x = vector_alloc (n);
	qr_gs_decomp (D, RR);
	qr_gs_solve (D, RR, b, x);
	vector_print (x);

	vector *x2 = vector_alloc (n);
	matrix_x_vector (E, x, x2);
	vector_print (x2);
	vector_print (b);

	// second task

	matrix *DT = matrix_alloc (n, n);	
	matrix *E2 = matrix_alloc (n, n);	
	qr_gs_inverse (D, RR, DT);

	matrix_x_matrix (DT, E, E2);
	printf("matrix times inverse matrix \n");
	matrix_print (E2);
	printf("inverse matrix \n");
	matrix_print (DT);

	matrix_free (D);
	matrix_free (E);
	matrix_free (E2);
	matrix_free (RR);
	matrix_free (DT);
	vector_free (b);
	vector_free (x2);
	vector_free (x);
	return 0;
}
