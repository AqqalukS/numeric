#ifndef HAVE_FUNCS_H
#include "vector.h"
#include "matrix.h"
double rosen (vector *v);
double himmel (vector *v);
void rosen_df (vector *v, vector *df);
void himmel_df (vector *v, vector *df);
void rosen_H (vector *v, matrix *H);
void himmel_H (vector *v, matrix *H);
#define HAVE_FUNCS_H
#endif
