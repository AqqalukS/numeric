#include "vector.h"
#include <math.h>

vector * vector_alloc (int n) {
	vector *v = (vector*) malloc(sizeof(vector));
	v->size = n;
	v->data = (double*) malloc(n *sizeof(double));
	assert(v!=NULL);

	return v;
}

void vector_free (vector *v) {
	free (v->data);
	free (v);
}

double vector_get (vector* v, int i) {
	assert(i>=0 && i < v->size);
	return v->data[i];
}

void vector_set (vector *v, int i, double x) {
	assert(i>=0 && i < v->size);
	v->data[i] = x;
}

void vector_print (vector *v) {
	for (int i = 0; i < v->size; i++) {
		printf("%7.5f\n", vector_get (v, i));
	}
	printf("\n");
}

void vector_add (vector *v, double x, vector *u) {
	assert(v->size == u->size);
	for (int i = 0; i < v->size; i++) {
		vector_set (u, i, vector_get (v, i) + x);
	}
}

double vector_dot (vector *v, vector *u) {
	assert(v->size == u->size);
	double s = 0;
	for (int i = 0; i < v->size; i++) {
		s += vector_get (v, i) * vector_get (u, i);
	}
	return s;
}

double vector_norm (vector *v) {
	return sqrt(vector_dot (v, v));
}

void vector_scale (vector *v, double x) {
	for (int i = 0; i < v->size; i++) {
		vector_set (v, i, vector_get (v, i) * x);
	}
}

void vector_plus_vector (vector *v, vector *u) {
	assert(v->size == u->size);
	for (int i = 0; i < v->size; i++) {
		vector_set (v, i, vector_get (v, i) + vector_get (u, i));
	}
}

void vector_minus_vector (vector *v, vector *u) {
	assert(v->size == u->size);
	for (int i = 0; i < v->size; i++) {
		vector_set (v, i, vector_get (v, i) - vector_get (u, i));
	}
}

void vector_cp (vector *v, vector *u) {
	assert(v->size == u->size);
	for (int i = 0; i < v->size; i++) {
		vector_set (v, i, vector_get (u, i));
	}
}

double vectorT_x_vector (vector *v, vector *u) {
	assert(v->size == u->size);
	double s = 0;
	for (int i = 0; i < v->size; i++) {
		s += vector_get (v, i) * vector_get (u, i);
	}
	return s;
}

