#include <gsl/gsl_math.h>
#include "cuba.h"
#include "tdho.c"
#include "utility.h"

double operator_id_1b(q_nums i_nums, q_nums f_nums) {
    return (!isequal_q_nums(i_nums, f_nums)) ? 0 : 2;
}

double operator_id_2b(q_nums i_nums, q_nums f_nums, double b) {
    double pi_factor = 2 / gsl_pow_2(2 * M_PI);

    if (!isequal_q_nums(i_nums, f_nums))
	return 0;
    else if (i_nums.l != 0 || f_nums.l !=0)
	return 0;

    wf_params i_params = { i_nums.n, 0, 0, b };
    wf_params f_params = { f_nums.n, 0, 0, b };
    
}

double operator_id(q_nums i_nums, q_nums f_nums) {
    
}

	
