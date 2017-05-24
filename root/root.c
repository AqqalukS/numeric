#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "lineq.h"


int newton_root (void f(vector *, vector *), vector *x, double dx, double eps) {
	int n = x->size, ds = 0;
	matrix *J = matrix_alloc (n, n);
	matrix *R = matrix_alloc (n, n);
	vector *fx = vector_alloc (n);
	vector *fdx = vector_alloc (n);
	vector *xd = vector_alloc (n);
	vector *Dx = vector_alloc (n);
	double normfx, l;

	do {
		ds ++;
		f(x, fx);
		for (int j = 0; j < n; j++) {
			vector_set (x, j, vector_get (x, j) + dx);
			f(x, fdx); 
			//vector_add (fdx, -1, fx);
			for (int i = 0; i < n; i++) {
				matrix_set (J, i, j, (vector_get (fdx, i) -
							vector_get (fx, i)) / dx);
			}
			vector_set (x, j, vector_get (x, j) - dx); 
		}
		qr_gs_decomp (J, R);

		vector_scale (fx, -1);
		qr_gs_solve (J, R, fx, Dx);
		vector_scale (fx, -1);

		normfx = vector_norm (fx);
		l = 2;
		do {
			l /= 2;
			for (int i = 0; i < n; i++) {
				vector_set (xd, i, vector_get (x, i) + vector_get (Dx, i) * l);
			}
			f (xd, fdx);
		} while (vector_norm (fdx) > (1 - l/2) * normfx && l  > 0.02);
		for (int i = 0; i < n; i++) {
			vector_set (x, i, vector_get (xd, i));
			vector_set (fx, i, vector_get (fdx, i));
		}

	} while (vector_norm (Dx) > dx && normfx > eps);

	matrix_free (J);
	vector_free (fx);
	vector_free (fdx);
	vector_free (xd);
	vector_free (Dx);
	return ds;
}
