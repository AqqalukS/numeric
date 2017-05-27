#include "lineq.h"
#define RND ((double)rand()/RAND_MAX)

int main(int argc, const char *argv[]) {
	printf("################################################################\n");
	printf("# Exercise B - Matrix inverse by Gram-Schmidt QR factorization #\n");
	printf("################################################################\n");

	int n = 6;
	
	matrix *Q = matrix_alloc (n,n);
	matrix *A = matrix_alloc (n,n);
	matrix *invA = matrix_alloc (n,n);
	matrix *R = matrix_alloc (n,n);
	// vector *b = vector_alloc (n);

	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < A->size2; j++) {
			matrix_set (Q, i, j, RND);
			matrix_set (A, i, j, matrix_get (Q, i, j));
		}
		//vector_set (b, i, RND);
	}

	qr_gs_decomp (Q, R);
	qr_gs_inverse (Q, R, invA);

	printf("# A =\n");
	matrix_print (A);

	printf("# A^-1 \n");
	matrix_print (invA);

	matrix_x_matrix (A, invA, Q);
	printf("# A * A^-1 = I = \n");
	matrix_print (Q);

	matrix_free (Q);
	matrix_free (A);
	matrix_free (invA);
	matrix_free (R);
	return 0;
}
