#ifndef VECTOR_H
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {int size; double *data;} vector;

vector * vector_alloc (int n);
void vector_free (vector *v);
double vector_get (vector *v, int i);
void vector_set (vector *v, int i, double x);
void vector_print (vector *v);
void vector_add (vector *v, double x, vector *u);
double vector_dot (vector *v, vector *u);
double vector_norm (vector *v);
void vector_scale (vector *v, double x);
void vector_plus_vector (vector *v, vector *u);
void vector_minus_vector (vector *v, vector *u);
void vector_cp (vector *v, vector *u);
double vectorT_x_vector (vector *v, vector *x);

#define VECTOR_H 
#endif
