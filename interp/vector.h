#ifndef VECTOR_H
typedef struct {int size; double *data;} vector;

vector * vector_alloc (int n);
void vector_free (vector *v);
double vector_get (vector *v, int i);
void vector_set (vector *v, int i, double x);
#define VECTOR_H 
#endif
