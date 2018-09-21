#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include "utility.h"

double cg_coeff (int ja, int jb, int j, int mja, int mjb, int mj) {
    int triangle_cond, mj_cond;
    triangle_cond = abs(ja-jb) <= j && j <= (ja+jb);
    mj_cond = mja + mjb == mj;
    if (!triangle_cond || !mj_cond)
	return 0;
    else
	return (gsl_pow_int(-1, ja-jb+mj) * sqrt(2*j + 1)
		* gsl_sf_coupling_3j(2*ja, 2*jb, 2*j, 2*mja, 2*mjb, -2*mj));
}
