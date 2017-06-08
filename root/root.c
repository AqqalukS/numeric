#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "root.h"
#include "vector.h"
#include "matrix.h"
#include "lineq.h"


int newton_root (void f(vector *, vector *), vector *x, double dx, double eps) {
	int n = x->size, ds = 0;
	double normfx, l;
	matrix *J = matrix_alloc (n, n);
	vector *fx = vector_alloc (n);
	vector *fdx = vector_alloc (n);
	vector *xd = vector_alloc (n);
	vector *Dx = vector_alloc (n);

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
		givens_decomp (J);

		vector_scale (fx, -1);
		vector_cp (Dx, fx);
		vector_scale (fx, -1);

		givens_solve (J, Dx);

		normfx = vector_norm (fx);
		l = 2;
		do {
			l /= 2;
			vector_cp (xd, x);
			vector_scale (Dx, l);
			vector_plus_vector (xd, Dx);
			vector_scale (Dx, 1/l);
			for (int i = 0; i < n; i++) {
				vector_set (xd, i, vector_get (x, i) + vector_get (Dx, i) * l);
			}
			f (xd, fdx);
		} while (vector_norm (fdx) > (1 - l/2) * normfx && l  > 0.02);
		vector_cp (x, xd);
		vector_cp (fx, fdx);
		normfx = vector_norm (fx);
	} while (vector_norm (Dx) > dx && normfx > eps);

	matrix_free (J);
	vector_free (fx);
	vector_free (fdx);
	vector_free (xd);
	vector_free (Dx);
	return ds;
}

int newton_jacobian (void f(vector *, vector *), void jacobian(vector *, matrix *), vector *x,
				double eps) {
	int n = x->size, ds = 0;
	double normfx, l;
	matrix *J = matrix_alloc (n, n);
	vector *fx = vector_alloc (n);
	vector *fdx = vector_alloc (n);
	vector *xd = vector_alloc (n);
	vector *Dx = vector_alloc (n);

	do {
		ds ++;
		f(x, fx);
		jacobian(x, J);

		givens_decomp (J);

		vector_scale (fx, -1);
		vector_cp (Dx, fx);
		vector_scale (fx, -1);

		givens_solve (J, Dx);

		normfx = vector_norm (fx);
		l = 2;
		do {
			l /= 2;
			for (int i = 0; i < n; i++) {
				vector_set (xd, i, vector_get (x, i) + vector_get (Dx, i) * l);
			}
			f (xd, fdx);
		} while (vector_norm (fdx) > (1 - l/2) * normfx && l  > 0.02);
		vector_cp (x, xd);
		vector_cp (fx, fdx);
		normfx = vector_norm (fx);
	} while (normfx > eps);

	matrix_free (J);
	vector_free (fx);
	vector_free (fdx);
	vector_free (xd);
	vector_free (Dx);
	return ds;
}
	
