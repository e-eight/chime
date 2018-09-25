#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include "tdho.h"
#include "operator_r_sq.h"
//#include "cubature.h"

double integrand_l(double y, void *gsl_params) {
    integrand_l_params *params = (integrand_l_params *) gsl_params;
    return gsl_sf_laguerre_n(params->n, params->alpha, y)
	* gsl_sf_laguerre_n(params->np, params->alphap, y);
}

double integral_l1(wf_params *i_params, wf_params *f_params) {
    int n = i_params->n, l = i_params->l, np = f_params->n, lp = f_params->l;
    double GSL_A = 0, GSL_B = 1, GSL_ALPHA = (lp + l - 1) / 2;
    integrand_l_params params = { n, l+0.5, np, lp+0.5 };
    
    gsl_integration_fixed_workspace *workspace
    	= gsl_integration_workspace_alloc(1000);
    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;

    int num_nodes =  (n + np + 1) / 2 + 5;
    workspace = gsl_integration_fixed_alloc(T, num_nodes, GSL_A, GSL_B, GSL_ALPHA, 0.0);

    /* gsl_integration_workspace *workspace */
    /* 	= gsl_integration_workspace_alloc(1000); */
        
    gsl_function F;
    F.function = &integrand_l;
    F.params = &params;
    
    double result;
    gsl_integration_fixed(&F, &result, workspace);

    /* gsl_integration_qags (&F, 0, 1, 0, 1e-3, 1000, workspace, &result, &error); */

    gsl_integration_fixed_free(workspace);

    //  hcubature_v(1, &integrand_l, &integrand_params, 1, &tmin, &tmax,
//	      0, 0, 1e-4, ERROR_INDIVIDUAL, &result, &error);
    
    return result;
}
