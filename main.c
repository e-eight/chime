#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "tdho.h"
#include "utility.h"

double radial_1d(double p, void* gsl_params) {
    wf_params* params = (wf_params*) gsl_params;
    double b_factor = M_SQRT2 / gsl_pow_3(params->b);
    return b_factor * norm_nl(params) * gsl_sf_laguerre_n(params->n, 0.5, 2*p);
}

int main() {
    wf_params iparams = { 100, 0, 0, 1 };
    double N_nl;

    N_nl = norm_nl(&iparams);
    

    gsl_integration_fixed_workspace* w
	= gsl_integration_workspace_alloc(1000);
    const gsl_integration_fixed_type* T = gsl_integration_fixed_laguerre;

    int num_nodes = (iparams.n + 1) / 2 + 1;
    w = gsl_integration_fixed_alloc(T, num_nodes, 0.0, 1.0, 0.5, 0.0); 

    gsl_function F;
    F.function = &radial_1d;
    F.params = &iparams;

    double result;
    gsl_integration_fixed(&F, &result, w);

    printf("%f\n", result);

    gsl_integration_fixed_free(w);

/* #pragma omp parallel for */
/*     for (int i=0; i <= 4; ++i) */
/* 	for (int j=0; j <= 4; ++j) */
/* 	    if (i > j) */
/* 		printf("%d\t %d\t %d\n", i, j, omp_get_thread_num()); */
    
}
