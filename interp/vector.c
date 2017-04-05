#include "vector.h"
#include <stdlib.h>
#include <assert.h>

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
