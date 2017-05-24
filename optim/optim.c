#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"
#include "lineq.h"

int quasi_newton (double f(vector *), vector *x, double dx, double eps) {
	int n = x->size, ds = 0;

	void gradient (vector *v, vector *dv) {
		double fv = f(v);
		for (int i = 0; i < v->size; i++) {
			vector_set (v, i, vector_get (v, i) + dx);
			vector_set (dv, i, (f(v) - fv) / dx);
			vector_set (v, i, vector_get (v, i) - dx);
		}
	}

	matrix *H = matrix_alloc (n, n);
	vector *df = vector_alloc (n);
	vector *dfdx = vector_alloc (n);
	vector *Dx = vector_alloc (n);
	vector *s = vector_alloc (n);
	vector *u = vector_alloc (n);
	vector *y = vector_alloc (n);
	vector *hd = vector_alloc (n);

	double fx = f(x), fdx, sdf;

	matrix_set_identity (H);
	gradient (x, df);
	do {
		ds++;
		vector_scale (df, -1);
		matrix_x_vector (H, df, s);
		vector_scale (df, -1);
		vector_scale (s, 2);
		do {
			vector_scale (s, 0.5);
			vector_cp (Dx, x);
			vector_plus_vector (Dx, s);
			fdx = f(Dx);
			sdf = vector_dot (s, df);
			if (fabs(fdx) < fabs(fx) + 0.1 * sdf) {
				break;
			}
			if (vector_norm (s) < dx) {
				matrix_set_identity (H);
				break;
			}
		} while (1);
		gradient (Dx, dfdx);
		vector_cp (y, dfdx);
		vector_minus_vector (y, df);
		matrix_x_vector (H, y, hd);
		vector_cp (u, s);
		vector_minus_vector (u, hd);

		for (int i = 0; i < H->size1; i++) {
			for (int j = 0; j < H->size2; j++) {
				matrix_set (H, i, j, matrix_get (H, i, j) 
					+ vector_get (u, i) * vector_get (u, j)
					* 1/ vector_dot (u, y));
			}
		}
		vector_cp (x, Dx);
		vector_cp (df, dfdx);
		fx = fdx;
	} while (vector_norm (Dx) > dx && vector_norm (df) > eps);

	matrix_free (H);
	vector_free (df);
	vector_free (dfdx);
	vector_free (Dx);
	vector_free (s);
	vector_free (u);
	vector_free (y);
	vector_free (hd);
	return ds;
}

int newton_optim (double f(vector *x), void gradient(vector *x, vector *df),
		void hessian(vector *x, matrix *H), vector *x, double eps) {
	int n = x->size, ds = 0;
	double fx, a = 1e-4, l, dxdf, normfx;
	vector *df = vector_alloc (n);
	vector *dx = vector_alloc (n);
	vector *dxx = vector_alloc (n);
	vector *xk = vector_alloc (n);
	matrix *H  = matrix_alloc (n, n);
	matrix *R  = matrix_alloc (n, n);

	gradient (x, df);

	do {
		ds++;
		fx = f(x);
		hessian (x, H);	
		
		qr_gs_decomp (H, R);
		vector_scale (df, -1);
		qr_gs_solve (H, R, df, dx);
		vector_scale (df, -1);

		l = 1;
		do {
			l /= 2;
			vector_cp (xk, x);
			vector_cp (dxx, dx);
			vector_scale (dxx, l);
			vector_plus_vector (xk, dxx);
			dxdf = vectorT_x_vector (dxx, df);
		} while (f (xk) < fx + a * l * dxdf && l > 1e-3);
		vector_cp (x, xk);
		gradient (x, df);
		normfx = vector_norm (df);
	} while (normfx > eps);

	vector_free (df);
	vector_free (dx);
	vector_free (dxx);
	vector_free (xk);
	matrix_free (H);
	matrix_free (R);
	return ds;
}
