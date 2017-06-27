#include <math.h>
#include "vector.h"
#include "matrix.h"

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
