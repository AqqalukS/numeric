#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"
#include "funcs.h"

int quasi_newton (double f(vector *), vector *x, double dx, double eps);
int newton_optim (double f(vector *), void gradient(vector *, vector *),
		void hessian(vector *, matrix *), vector *x, double eps);

int main(int argc, const char *argv[])
{
	// Exercise A
	printf("##############\n");
	printf("# Exercise A #\n");
	printf("##############\n");
	// Allocations
	int n = 2, ds;
	double fx, eps;
	vector *x = vector_alloc (n);
	vector_set (x, 0, 2);
	vector_set (x, 1, 2);

	eps = 1e-8;

	fx = rosen (x);
	printf("Initial x = \n");
	vector_print (x);
	printf("f(x) = %g\n", fx);

	ds = newton_optim (rosen, rosen_df, rosen_H, x, eps);
	fx = rosen (x);
	printf("number of steps = %d\n", ds);
	printf("after optimisation x = \n");
	vector_print (x);
	printf("f(x) = %g\n", fx);
	printf("\nabla f = \n");
	vector *df = vector_alloc (n);
	rosen_df (x, df);
	vector_print (df);

	vector_set (x, 0, 1);
	vector_set (x, 1, 1);
	fx = himmel (x);
	printf("Initial x = \n");
	vector_print (x);
	printf("f(x) = %g\n", fx);
	himmel_df (x, df);
	vector_print (df);

	ds = newton_optim (himmel, himmel_df, himmel_H, x, eps);
	fx = himmel (x);
	printf("number of steps = %d\n", ds);
	printf("after optimisation x = \n");
	vector_print (x);
	printf("f(x) = %g\n", fx);
	printf("\nabla f = \n");
	himmel_df (x, df);
	vector_print (df);

	vector_free (x);
	vector_free (df);
	return 0;
}
