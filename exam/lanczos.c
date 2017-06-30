#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "vector.h"
#include "matrix.h"
#define RND ((double)rand() / RAND_MAX)

void orthonormalize (vector *v, matrix *V, int m);

void lanczos (matrix *A, vector *u, matrix *V, matrix *T) {
	/* Lanczos tridiagonalization algorithm:
	 * tridigonalizes A and output in to a orthonormal matrix V and trigiagonal real symmetric
	 * matrix T = V* AV size nxn in this case.
	 */
	assert(A->size1 == A->size2);
	assert(A->size1 == V->size1); assert(V->size1 == V->size2);
	assert(A->size1 == T->size1); assert(T->size1 == T->size2);	
	int n = A->size1;
	double beta, alpha;

	vector *w = vector_alloc (n);
	vector *v = vector_alloc (n);

	vector_cp (v, u);
	beta = vectorT_x_vector (v, v);
	vector_scale (v, 1/sqrt(beta));
	matrix_set_column (V, 0, v);

	matrix_x_vector (A, v, w); 			// w1 = Av1
	alpha =  vectorT_x_vector (w, v); 		// alpha1 = w1' * v1;
	matrix_set (T, 0, 0, alpha);
	vector_scale (v, alpha);
	vector_minus_vector (w, v);			// w1 = w1' - alpha * v1
	for (int i = 1; i < n; i++) {
		beta = vector_norm (w);			// beta_i = ||w_i-1||
		matrix_set (T, i, i-1, beta); matrix_set (T, i-1, i, beta);
		if (beta != 0) {
			vector_cp (v, w);
			vector_scale (v, 1 / beta);
		}
		else {
			for (int k = 0; k < n; k++) {
				vector_set (v, k, RND);
			}
			orthonormalize (v, V, i);
		}
		matrix_set_column (V, i, v);
		matrix_x_vector (A, v, w);
		alpha = vectorT_x_vector (w, v);
		matrix_set (T, i, i, alpha);
		
		vector_scale (v, alpha);
		vector_minus_vector (w, v);
		matrix_get_column (V, i-1, v);
		vector_scale (v, beta);
		vector_minus_vector (w, v);
	}

	vector_free (w);
	vector_free (v);
}
