#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

double gsl_int (double func(double), double a, double b, double epsabs, double epsrel) {
	double integrand (double x, void *params) {
		return func(x);
	}
	gsl_function f;
	f.function = integrand;
	f.params = NULL;
	int limit = 1000;
	double result = 0;

	if (!isinf(a) && !isinf(b)) {
		double res,error;
								    
		gsl_integration_workspace *w =
			gsl_integration_workspace_alloc (limit);

		int status = 
			gsl_integration_qags
				( &f, a, b, epsabs, epsrel, limit, w, &res, &error);
		gsl_integration_workspace_free (w);

		if(status!=GSL_SUCCESS) return NAN;
		else result = res;
	}
	else if (isinf(a) && isinf(b)) { 
		double resiu,resil,erriu, erril, a = 0, bb = 0;
		            
		gsl_integration_workspace 
			*w_iu =
				gsl_integration_workspace_alloc (limit),
			*w_il =
				gsl_integration_workspace_alloc (limit);

		int status_iu = 
			gsl_integration_qagiu
				( &f, a, epsabs, epsrel, limit, w_iu, &resiu, &erriu);
		int status_il = 
			gsl_integration_qagil
				( &f, bb, epsabs, epsrel, limit, w_il, &resil, &erril);

		gsl_integration_workspace_free (w_iu);
		gsl_integration_workspace_free (w_il);
		if ((status_iu != GSL_SUCCESS) | (status_il != GSL_SUCCESS)) {
			return NAN;
		}
		else {
			if (b == -INFINITY) {
				result = -(resiu+resil);
			}
			else {
				result = resiu+resil;
			}
		}
	}
	else if ((!isinf(a) && b == INFINITY) ||
			(a == INFINITY && !isinf(b))) {
		double aa;
		if (!isinf(b)) {
			aa = b;
		}
		else {
			aa = a;
		}
		double resiu, erriu;
		            
		gsl_integration_workspace 
			*w_iu =
				gsl_integration_workspace_alloc (limit);

		int status_iu = 
			gsl_integration_qagiu
				( &f, aa, epsabs, epsrel, limit, w_iu, &resiu, &erriu);

		gsl_integration_workspace_free (w_iu);
		if (status_iu != GSL_SUCCESS) return NAN;
		else {
			if (!isinf(b)) {
				result = -resiu;
			}
			else {
				result = resiu;
			}
		}
	}
	else if ((a == -INFINITY && !isinf(b)) ||
			(!isinf(a) && b == -INFINITY)) {
		double bb;
	       	if (a != -INFINITY) {
	       		bb = a;
	       	}
		else {
			bb = b;
		}

		double resil, erril;
		            
		gsl_integration_workspace 
			*w_il =
				gsl_integration_workspace_alloc (limit);

		int status_il = 
			gsl_integration_qagil
				( &f, bb, epsabs, epsrel, limit, w_il, &resil, &erril);

		gsl_integration_workspace_free (w_il);
		if (status_il != GSL_SUCCESS) return NAN;
		else {
			if (a != -INFINITY) {
				result = -resil;
			}
			else {
				result = resil;
			}
		}
	}
	return result;
}
