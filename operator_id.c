#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include "tdho.h"
#include "operator_id.h"

double operator_id_1b(q_nums *i_nums, q_nums *f_nums) {
    return (!isequal_q_nums(i_nums, f_nums)) ? 0 : 2;
}

double radial_1d(double p, void *gsl_params) {
    wf_params *params = (wf_params *) gsl_params;
    double b_factor = M_SQRT2 / gsl_pow_3(params->b);
    return b_factor * norm_nl(params) * gsl_sf_laguerre_n(params->n, 0.5, 2*p);
}

double operator_id_2b(q_nums *i_nums, q_nums *f_nums, double *b) {
    double pi_factor = 2 / gsl_pow_2(2 * M_PI);

    if (!isequal_q_nums(i_nums, f_nums))
    	return 0;
    else if (i_nums->l != 0 || f_nums->l !=0)
    	return 0;

    wf_params i_params = { i_nums->n, 0, 0, *b };
    wf_params f_params = { f_nums->n, 0, 0, *b };

    gsl_integration_fixed_workspace *w_i
	= gsl_integration_workspace_alloc(1000);
    gsl_integration_fixed_workspace *w_f
	= gsl_integration_workspace_alloc(1000);
    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;

    int num_nodes_i = (i_params.n + 1) / 2 + 1;
    int	num_nodes_f = (f_params.n + 1) / 2 + 1;
    w_i = gsl_integration_fixed_alloc(T, num_nodes_i, 0.0, 1.0, 0.5, 0.0); 
    w_f = gsl_integration_fixed_alloc(T, num_nodes_f, 0.0, 1.0, 0.5, 0.0);
    
    gsl_function F_i, F_f;
    F_i.function = &radial_1d;
    F_i.params = &i_params;
    F_f.function = &radial_1d;
    F_f.params = &f_params;

    double result_i, result_f;
    gsl_integration_fixed(&F_i, &result_i, w_i);
    gsl_integration_fixed(&F_f, &result_f, w_f);

    gsl_integration_fixed_free(w_i);
    gsl_integration_fixed_free(w_f);

    return pi_factor * result_i * result_f;
}

double operator_id(q_nums *i_nums, q_nums *f_nums, double *b) {
    return operator_id_1b(i_nums, f_nums) + operator_id_2b(i_nums, f_nums, b);
}

	
