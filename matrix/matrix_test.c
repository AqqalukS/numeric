#include <vector.h>
#include "matrix.h"
#define RND ((double)rand()/RAND_MAX)

int main()
{
	int n = 4, m = 4;
	matrix *A = matrix_alloc (n, m);
	matrix *B = matrix_alloc (n, n);
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			matrix_set (A, i, j, RND);
		}
		for (int j = 0; j < n; j++) {
			matrix_set (B, i, j, RND);
		}
	}
	printf("# matrix A\n");
	matrix_print (A);
	printf("# matrix B\n");
	matrix_print (B);

	matrix *C = matrix_alloc (n,m);
	matrix_x_matrix (A,B,C);
	
	printf("# A times B\n");
	matrix_print (C);

	matrix_set_identity (B);
	printf("# Identity \n");
	matrix_print (B);

	vector *vn = vector_alloc (n);
	vector *vm = vector_alloc (m);
	matrix_get_column (A, 0, vn);

	printf("# first column of A\n");
	vector_print (vn);

	for (int i = 0; i < n; i++) {
		vector_set (vn, i, 5);
	}
	for (int i = 0; i < m; i++) {
		vector_set (vm, i, RND);
	}

	matrix_set_column (A, 0, vn);

	vector *bn = vector_alloc (n);
	vector *bm = vector_alloc (m);

	matrix_x_vector (A, vm, bn);
	matrixT_x_vector (A, vn, bm);

	printf("# A times vector\n");
	vector_print (bn);

	printf("# A^T times vector\n");
	vector_print (bm);


	vector_free (vm);
	vector_free (vn);
	vector_free (bn);
	vector_free (bm);
	matrix_free (A);
	matrix_free (B);
	matrix_free (C);
	return 0;
}
