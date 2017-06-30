#include <math.h>
#include "vector.h"
#include "matrix.h"

void arnoldi (matrix *A, vector *q, matrix *Q, matrix *H) {
	int n = Q->size1;
	vector *qi = vector_alloc (n);
	vector *qk = vector_alloc (n);
	double normq, QiQk;

	vector_cp (qi, q);
	normq = vector_norm (qi);
	vector_scale (qi, 1/normq);
	matrix_set_column (Q, 0, qi);

	for (int k = 1; k < n; k++) {
		matrix_get_column (Q, k - 1, qi);
		matrix_x_vector (A, qi, qk);
		for (int i = 0; i < k; i++) {
			matrix_get_column (Q, i, qi);
			QiQk = vectorT_x_vector (qi, qk);
			matrix_set (H, i, k - 1, QiQk);
			vector_scale (qi, QiQk);
			vector_minus_vector (qk, qi);
		}
		normq = vector_norm (qk);
		matrix_set (H, k, k - 1, normq);
		vector_scale (qk, 1/normq);
		matrix_set_column (Q, k, qk);
	}
	vector_cp (qi, qk);
	matrix_x_vector (A, qi, qk);
	for (int i = 0; i < n; i++) {
		matrix_get_column (Q, i, qi);
		QiQk = vectorT_x_vector (qi, qk);
		matrix_set (H, i, n - 1, QiQk);
		vector_scale (qi, QiQk);
		vector_minus_vector (qk, qi);
	}
	vector_free (qi);
	vector_free (qk);
}
