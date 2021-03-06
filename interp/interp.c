#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "interp.h"
#include "vector.h"

double linterp (vector *x, vector *y, double z) {
	assert(x->size == y->size);
	int n = x->size;
	assert(n>0 && z >= vector_get(x, 0) && z <= vector_get(x, n-1));
	int i = 0, j = n - 1;
	while (j - i > 1){
		int m = (i + j)/2;
		if (z > vector_get(x, m)) {
			i = m;
		}
		else j = m;
	}
	double dy = vector_get(y, i+1) - vector_get(y, i);
       	double dx = vector_get(x, i+1) - vector_get(x, i);
	double p = dy/dx;
	return vector_get(y, i) +p * (z - vector_get(x, i));
}

double linterp_integ (vector *x, vector *y, double z) {
	double s = 0, dx, dy;
	int j = x->size - 1, i = 0, m;
	while (j - i > 1) {
		m = (i + j) / 2;
		if (z > vector_get (x, m)) {
			i = m;
		}
		else {
			j = m;
		}
	}
	for (int j = 0; j < i; j++) {
		dx = vector_get (x, j+1) - vector_get (x, j);
		dy = vector_get (y, j+1) - vector_get (y, j);
		s += vector_get (y, j) * dx + 0.5 * dy/dx * dx * dx;
	}
	dx = vector_get (x, i+1) - vector_get (x, i);
	dy = vector_get (y, i+1) - vector_get (y, i);
	s = s + vector_get (y, i) * (z - vector_get (x, i)) 
			+ 0.5 * dy/dx * (z - vector_get (x, i)) * (z - vector_get (x, i));
	return s;
}


qspline * qspline_alloc (vector *vx, vector *vy) {
	assert(vx->size == vy->size);
	// Allocations
	int n = vx->size;
	qspline *s = (qspline*) malloc(sizeof(qspline));
	s->x = (double*) malloc(n * sizeof(double));
	s->y = (double*) malloc(n * sizeof(double));
	s->b = (double*) malloc((n-1) * sizeof(double));
	s->c = (double*) malloc((n-1) * sizeof(double));
	s->n = n;

	int i;
	double x[n], y[n];
	for (i = 0; i < n; i++) {
		s->x[i] = x[i] = vector_get(vx, i);
		s->y[i] = y[i] = vector_get(vy, i);
	}

	double p[n-1], dx[n-1], dy[n-1];
	
	// Calculate p
	for (i = 0; i < n-1; i++) {
		dx[i] = x[i+1] - x[i]; dy[i] = y[i+1] - y[i];
		p[i] = dy[i]/dx[i];
	}

	// forward recursive
	s->c[0] = 0;
	for (i = 0; i < n-2; i++) {
		s->c[i+1] = (p[i+1] - p[i] - s->c[i] * dx[i]) / dx[i+1];
	}
	// backwards recursive
	s->c[n-2] = 2;
	for (i = n-3; i >= 0; i--) {	
		s->c[i] = (p[i+1] - p[i] - s->c[i+1] * dx[i+1]) / dx[i];
	}
	// calculate b
	for (i = 0; i < n-1; i++) {
		s->b[i] = p[i] - s->c[i] * dx[i];
	}
	return s;
}

double qspline_evaluate (qspline *s, double z) {
	assert (z <= s->x[s->n - 1] && z >= s->x[0]);
	int i = 0, j = s->n - 1;
	
	while (j - i > 1){
		int m = (i + j)/2;
		if (z > s->x[m]) {
			i = m;
		}
		else j = m;
	}
	double dx = z - s->x[i];
	return s->y[i] + s->b[i] * dx + s->c[i] * dx*dx;
}

double qspline_integ (qspline *s, double z) {
	assert (z >= s->x[0]);
	assert (z <= s->x[s->n - 1]);
	double sum = 0, w;
	int i = 0, j = s->n -1, m;
	while (j - i > 1) {
		m = (i + j) / 2;
		if (z > s->x[m]) {
			i = m;
		}
		else {
			j = m;
		}
	}
	for (j = 0; j < i; j++) {
		w = s->x[j+1] - s->x[j];
		sum += s->y[j] * w 
			+ 0.5 * s->b[j] * pow (w, 2) 
			+ 1/3 * s->c[j] * pow (w, 3);
	}
	w = z - s->x[i];
	sum = sum + s->y[i] * w 
		+ 0.5 * s->b[i] * pow (w, 2) 
		+ 1/3 * s->c[i] * pow (w, 3);
	return sum;
}

double qspline_deriv (qspline *s, double z) {
	int i = 0, j = s->n-1;

	while (j - i > 1) {
		int m = (i + j) / 2;
		if (z > s->x[m]) {
			i = m;
		}
		else {
			j = m;
		}
	}
	double w = z - s->x[i];
	return s->b[i] + 2*s->c[i] * w;
}

void qspline_free (qspline *s) {
	free (s->x); free (s->y); free (s->b); free (s->c); free (s); 
}


cspline * cspline_alloc (vector *vx, vector *vy) {
	assert(vx->size == vy->size);
	// Allocations 
	int n  = vx->size;
	cspline *s = (cspline*) malloc(sizeof(cspline));
	s->x = (double*) malloc(n * sizeof(double));
	s->y = (double*) malloc(n * sizeof(double));
	s->b = (double*) malloc((n-1) * sizeof(double));
	s->c = (double*) malloc((n-1) * sizeof(double));
	s->d = (double*) malloc((n-1) * sizeof(double));
	s->n = n;

	int i = 0;
	double x[n], y[n];
	for (i = 0; i < n; i++) {
		s->x[i] = x[i] = vector_get (vx, i);
		s->y[i] = y[i] = vector_get (vy, i);
	}

	double p[n-1], dx[n-1], dy[n-1];
	
	// Calculate p
	for (i = 0; i < n-1; i++) {
		dx[i] = x[i+1] - x[i]; dy[i] = y[i+1] - y[i];
		p[i] = dy[i]/dx[i];
	}
	// the tridiagonal matrix
	double D[n], Q[n-1], B[n];

	D[0] = 2; D[n] = 2;
	Q[0] = 1;
	B[0] = 3 * p[1]; B[n] = 3 * p[n-1];	

	for (i = 0; i < n-2; i++) {
		D[i+1] = 2 + 2 * dx[i]/dx[i+1];
		Q[i+1] = dx[i]/dx[i+1];
		B[i+1] = 3 * (p[i] + p[i+1] * dx[i]/dx[i+1]);
	}

	// Gauss elimination

	for (i = 1; i < n; i++) {
		D[i] -= Q[i-1]/D[i-1];
		B[i] -= B[i-1]/D[i-1];
	}


	// Back-substitution
	s->b[n] = B[n]/D[n];
	for (i = n-2; i < 0; i--) {
		s->b[i] = (B[i] - Q[i] * s->b[i+1])/D[i];
	}

	// Calculate c and d
	
	for (i = 0; i < n-1; i++) {
		s->c[i] = (-2 * s->b[i] - s->b[i+1] + 3 * p[i])/dx[i];
		s->d[i] = (s->b[i] - s->b[i+1] - 2 * p[i])/dx[i]/dx[i];
	}
	return s;
}


double cspline_evaluate (cspline *s, double z) {
	assert(z <= s->x[s->n - 1] && z >= s->x[0]);
	int i = 0, j = s->n - 1;
	
	while (j - i > 1){
		int m = (i + j)/2;
		if (z > s->x[m]) {
			i = m;
		}
		else j = m;
	}
	double dx = z - s->x[i];
	return s->y[i] + s->b[i] * dx + s->c[i] * dx*dx + s->d[i] * dx*dx*dx;
}

double cspline_integ (cspline *s, double z) {
	assert (z >= s->x[0]);
	assert (z <= s->x[s->n - 1]);
	double sum = 0, w;
	int i = 0;
	while (i < s->n && s->x[i+1] < z) {
		w = s->x[i+1] - s->x[i];
		sum += s->y[i] * w 
			+ 0.5 * s->b[i] * pow (w, 2) 
			+ 1/3 * s->c[i] * pow (w, 3) 
			+ 1/4 * s->d[i] * pow (w, 4);
		i++;
	}
	w = z - s->x[i];
	return  sum + s->y[i] * w 
		+ 0.5 * s->b[i] * pow (w, 2) 
		+ 1/3 * s->c[i] * pow (w, 3) 
		+ 1/4 * s->d[i] * pow (w, 4);
}

double cspline_deriv (cspline *s, double z) {
	int i = 0, j = s->n-1;

	while (j - i > 1) {
		int m = (i + j) / 2;
		if (z > s->x[m]) {
			i = m;
		}
		else {
			j = m;
		}
	}
	double w = z - s->x[i];
	return s->b[i] + 2*s->c[i] * w + 3 * s->d[i] * w*w;
}

void cspline_free (cspline *s) {
	free (s->x); free (s->y); free (s->b);
	free (s->c); free (s->d); free (s); 
}

