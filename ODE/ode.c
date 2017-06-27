// e_i = ||dy|| < (epsilon ||y_i|| + tau) sqrt(h_i/(b-a))
// h_next = h (tau/e_i)^0.25 0.95
#include "vector.h"
#include "ode.h"
#define stepper rkstep23

void rkstep23 (double t, double h, vector *y, void f(double t, vector *y, vector *dydt), vector *yh, vector *err) {
	int n = y->size;
	vector *k1 = vector_alloc (n);
	vector *k2 = vector_alloc (n);
	vector *k3 = vector_alloc (n);
	vector *k4 = vector_alloc (n);
	vector *yt = vector_alloc (n);	

	f (t, y, k1);
	for (int i = 0; i < n; i++) {
		vector_set (yt, i, vector_get (y, i) + 1./2 * vector_get (k1, i) * h);
	}
	f (t + 1./2 * h, y, k2);
	for (int i = 0; i < n; i++) {
		vector_set (yt, i, vector_get (y, i) + 3./4 * vector_get (k2, i) * h);
	}
	f (t + 3./4 * h, y, k3);
	for (int i = 0; i < n; i++) {
		vector_set (yh, i, vector_get (y, i) + (
					2./9 * vector_get (k1, i) +
					1./3 * vector_get (k2, i) + 
					4./9 *vector_get (k3, i) 
				       )* h);
	}
	f (t + h, y, k4);
	for (int i = 0; i < n; i++) {
		vector_set (yt, i, vector_get (y, i) + (
					7./24 * vector_get (k1, i) +
					1./4 * vector_get (k2, i) + 
					1./3 *vector_get (k3, i) +
					1./8 * vector_get (k4, i) + 
				       ) * h);
		vector_set (err, i, vector_get (yh , i) - vector_get (yt, i));
	}

	vector_free (k1);
	vector_free (k2);
	vector_free (k3);
	vector_free (k4);
	vector_free (yt);
}

void driver (
	double *t, double b, double *h, vector *y, double acc, double eps, 
	void stepper (
		double t, double h, vector *y, 
		void f(double t, vector *y, vector *dydt), vector *yh, vector *err
		), 
	void f (double t, double *y, double *dydt)) {
	int n = y->size, k = 0;
	double tau, e, normy, s;
	vector *yh = vector_alloc (n);
	vector *err = vector_alloc (n);
	while (t < b) {
		if (t + h > b) {
			h = b - t;
		}
		stepper (t, h, y, f, yh, err);
		s = 0;
		for (int i = 0; i < n; i++) {
			s += vector_get (err, i)*vector_get (err, i);
		}
		e = sqrt (s);
		s = 0;
		for (int i = 0; i < n; i++) {
			s += vector_get (yh, i)*vector_get (yh, i);
		}
		normy = sqrt (s);

		tau = (normy * eps + acc) * sqrt (h/(b - a));
		if (err < tau) {
			k ++;
			t = t + h;
			vector_cp (y, yh);
		}
		if (err > 0) {
			h *= pow (tau/err, 0.25) * 0.95;
		else
			h *= 2;
		}
	}
}
